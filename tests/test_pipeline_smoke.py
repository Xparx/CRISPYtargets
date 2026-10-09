"""End-to-end smoke test for the crispytargets design pipeline.

Runs the documented invocation (``crispytargets design``, equivalently
``python -m crispytargets design``) against a synthetic 5kb chromosome with one gene model, through FlashFry
index/discover/score plus Doench 2016 rescoring, and checks the ranked
outputs are structurally sound. Skipped when java or the bundled FlashFry
JAR is missing. Per-guide flank geometry used in scoring is pinned
directly, without java, in test_flanks_pam.py.
"""

import os
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC = REPO_ROOT / "src"
TOOLS_DIR = REPO_ROOT / "tools"
JAR = TOOLS_DIR / "FlashFry-assembly-1.15.jar"

pytestmark = pytest.mark.skipif(
    shutil.which("java") is None or not JAR.exists(),
    reason="java or FlashFry JAR not available",
)

GFF3 = (
    "chr1\tcrispytest\tgene\t200\t900\t.\t+\t.\tID=gene1\n"
    "chr1\tcrispytest\tmRNA\t200\t900\t.\t+\t.\tID=tx1;Parent=gene1\n"
    "chr1\tcrispytest\texon\t200\t900\t.\t+\t.\tID=exon1;Parent=tx1\n"
)


def write_reference(path, length=5000, seed=1234):
    """Deterministic pseudo-random chromosome so runs are reproducible."""
    rng = random.Random(seed)
    seq = "".join(rng.choice("ACGT") for _ in range(length))
    with open(path, "w") as f:
        f.write(">chr1\n")
        for i in range(0, length, 60):
            f.write(seq[i:i + 60] + "\n")


@pytest.fixture(scope="module")
def pipeline_out(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("smoke")
    ref = tmp / "chr1.fa"
    gff = tmp / "test.gff3"
    out = tmp / "out"
    write_reference(ref)
    gff.write_text(GFF3)

    cmd = [
        sys.executable, "-m", "crispytargets", "design",
        "--gff", str(gff),
        "--ref", str(ref),
        "--outdir", str(out),
        "--flashfry-jar", str(JAR),
        "--tss-features", "mRNA",
        "--promoter-upstream", "100",
        "--promoter-downstream", "50",
        "--browser-bed",
    ]
    # Run against src/ without requiring an install first.
    env = {**os.environ, "PYTHONPATH": str(SRC)}
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=REPO_ROOT, env=env)
    assert proc.returncode == 0, f"pipeline failed:\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
    return out


def test_output_files_exist(pipeline_out):
    for name in ("final_sgrnas.tsv", "final_sgrnas.bed", "design.log",
                 "targets.output.scored", "targets.fasta"):
        assert (pipeline_out / name).exists(), f"missing {name}"


def test_ranked_table_structure(pipeline_out):
    df = pd.read_csv(pipeline_out / "final_sgrnas.tsv", sep="\t")
    assert len(df) > 0
    assert df.columns[0] == "gene_id"
    assert set(df["gene_id"]) == {"gene1"}
    for col in ("Doench2016OnTarget", "over_cutoff", "fwd", "rev", "target", "context"):
        assert col in df.columns, f"missing column {col}"

    scores = df["Doench2016OnTarget"]
    assert scores.notna().all()
    assert scores.between(0.0, 1.0).all()

    # T7 oligo convention from design_oligos: G-anchored core, atac/aaac tails.
    assert df["fwd"].str.startswith("agat").all()
    assert df["rev"].str.startswith("aaac").all()
    assert df["fwd"].str.len().eq(24).all()


def test_scores_show_plausible_spread(pipeline_out):
    """All-near-maximum scores indicate a feature-extraction bug, not luck."""
    df = pd.read_csv(pipeline_out / "final_sgrnas.tsv", sep="\t")
    scores = df["Doench2016OnTarget"]
    assert scores.max() < 1.0
    assert scores.max() - scores.min() > 0.05


def test_browser_bed(pipeline_out):
    lines = (pipeline_out / "final_sgrnas.bed").read_text().splitlines()
    assert lines[0].startswith("track name=CrispyGuides")
    rows = [l.split("\t") for l in lines[1:]]
    assert rows
    for row in rows:
        assert len(row) == 7
        assert row[6] in ("0,180,0", "255,140,0", "200,0,0")
    positions = [int(r[1]) for r in rows]
    assert positions == sorted(positions)
