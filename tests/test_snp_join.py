"""Guide-level SNP assignment and orientation tests.

SNP annotation must be computed inside each guide's own 23bp window and in the
guide's own orientation, not the target region's. These tests build a reference
with genes on both strands plus a dense SNP track, run the real CLI, and then
recompute every expectation directly from the raw coordinates (FASTA + BED),
independently of the pipeline code. They pin the defects this suite replaced:

1. SNPs were joined against the multi-kb REGION span, so any SNP landing in the
   region's first 23bp was painted onto every guide of that region, while real
   in-guide SNPs elsewhere went unannotated.
2. The guide-frame position and the ALT/IUPAC complement were taken from the
   REGION's strand, so guides facing the other way got mirrored positions and
   the wrong allele.
3. Flanks for RS2/RS3 were read on the region's strand, so every guide whose
   orientation disagreed with its region got context bases from the wrong side.

They also pin the cheap/expensive split: snp_pam_dist and snp_cut_dist arrive
with --snps alone, the variant score columns only with --rescore-snps.
"""

import os
import random
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pandas as pd
import pytest
from pyfaidx import Fasta

from crispytargets.design.doench2016 import score_guides

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC = REPO_ROOT / "src"
JAR = REPO_ROOT / "tools" / "FlashFry-assembly-1.15.jar"

pytestmark = pytest.mark.skipif(
    shutil.which("java") is None or not JAR.exists(),
    reason="java or FlashFry JAR not available",
)

CHROM_LEN = 2000
SNP_STEP = 40

GFF3 = (
    "chr1\tcrispytest\tgene\t300\t900\t.\t+\t.\tID=gene1\n"
    "chr1\tcrispytest\tmRNA\t300\t900\t.\t+\t.\tID=tx1;Parent=gene1\n"
    "chr1\tcrispytest\texon\t300\t900\t.\t+\t.\tID=exon1;Parent=tx1\n"
    "chr1\tcrispytest\tgene\t1100\t1700\t.\t-\t.\tID=gene2\n"
    "chr1\tcrispytest\tmRNA\t1100\t1700\t.\t-\t.\tID=tx2;Parent=gene2\n"
    "chr1\tcrispytest\texon\t1100\t1700\t.\t-\t.\tID=exon2;Parent=tx2\n"
)

_VARIANT_COLUMNS = [c for c in (
    "Doench2016OnTarget_avg1snp", "Doench2016OnTarget_worst1snp",
    "Doench2016OnTarget_allsnps", "Doench2016_delta_avg",
    "Doench2016_delta_worst", "Doench2016_delta_all",
)]


def revcomp(seq):
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


def write_reference(path, length=CHROM_LEN, seed=1234):
    """Deterministic pseudo-random chromosome so runs are reproducible."""
    rng = random.Random(seed)
    seq = "".join(rng.choice("ACGT") for _ in range(length))
    with open(path, "w") as f:
        f.write(">chr1\n")
        for i in range(0, length, 60):
            f.write(seq[i:i + 60] + "\n")


def snp_positions():
    """SNP positions (0-based, as BED/pyranges report them), every 40bp."""
    return list(range(200, CHROM_LEN - 199, SNP_STEP))


def write_snps(path):
    with open(path, "w") as f:
        for p in snp_positions():
            f.write(f"chr1\t{p}\t{p + 1}\tchr1:{p}:A/G\t1000\t.\n")


def expected_guide_hits(window_start, orientation, positions):
    """SNPs inside [g, g+23), as guide-frame positions.

    FWD guides read the window forward, RVS guides read it reverse
    complemented, so position p maps to g+22-p on the minus strand.
    """
    hits = [p for p in positions if window_start <= p < window_start + 23]
    if orientation == "RVS":
        return sorted(window_start + 22 - p for p in hits)
    return sorted(p - window_start for p in hits)


def expected_flanks(genome, window_start, orientation):
    """Guide-orientation flanks: 5bp 5' and 3bp 3' of the guide+PAM window."""
    g = window_start
    if orientation == "RVS":
        return revcomp(str(genome["chr1"][g + 23:g + 28].seq)), \
               revcomp(str(genome["chr1"][g - 3:g].seq))
    return str(genome["chr1"][g - 5:g].seq), str(genome["chr1"][g + 23:g + 26].seq)


def run_design(tmp_path, extra_args, out_name):
    ref = tmp_path / "chr1.fa"
    gff = tmp_path / "test.gff3"
    snps = tmp_path / "snps.bed"
    out = tmp_path / out_name
    if not ref.exists():
        write_reference(ref)
        gff.write_text(GFF3)
        write_snps(snps)
    cmd = [
        sys.executable, "-m", "crispytargets", "design",
        "--gff", str(gff),
        "--ref", str(ref),
        "--snps", str(snps),
        "--outdir", str(out),
        "--flashfry-jar", str(JAR),
        "--tss-features", "mRNA",
        "--promoter-upstream", "100",
        "--promoter-downstream", "50",
        "--browser-bed",
    ] + list(extra_args)
    env = {**os.environ, "PYTHONPATH": str(SRC)}
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=REPO_ROOT, env=env)
    assert proc.returncode == 0, f"pipeline failed:\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
    return out, ref


@pytest.fixture(scope="module")
def cheap_run(tmp_path_factory):
    """--snps only: annotation, entropy and distances, no variant scoring."""
    return run_design(tmp_path_factory.mktemp("snp_cheap"), [], "out")


@pytest.fixture(scope="module")
def rescore_run(tmp_path_factory):
    """--snps --rescore-snps: the variant score columns as well."""
    return run_design(tmp_path_factory.mktemp("snp_rescore"), ["--rescore-snps"], "out")


def load(out, ref):
    df = pd.read_csv(out / "final_sgrnas.tsv", sep="\t")
    return df, Fasta(str(ref))


def assert_matches_expectation(df, genome):
    """Every guide's annotation must equal the position-only ground truth."""
    positions = snp_positions()
    annotated_fwd = annotated_rvs = annotated_mismatched = 0
    for row in df.itertuples():
        g = int(row.Start) + int(row.start)
        exp_rels = expected_guide_hits(g, row.orientation, positions)
        if exp_rels:
            assert row.num_snps == len(exp_rels), (
                f"{row.gene_id}@{g}: num_snps {row.num_snps} != {len(exp_rels)} true hits")
            got_rels = sorted(int(x.split(":")[0]) for x in row.snp_overlaps.split("|"))
            assert got_rels == exp_rels, (
                f"{row.gene_id}@{g} ({row.orientation}): annotated {got_rels}, "
                f"true guide-frame positions {exp_rels}")
            assert row.entropy > 0
            assert row.snp_pam_dist == min(20 - r for r in exp_rels)
            assert row.snp_cut_dist == min(abs(r - 17) for r in exp_rels)
            if row.orientation == "FWD":
                annotated_fwd += 1
            else:
                annotated_rvs += 1
            region_matches = (row.Strand == "+") == (row.orientation == "FWD")
            if not region_matches:
                annotated_mismatched += 1
        else:
            assert row.num_snps == 0, f"{row.gene_id}@{g}: annotated without a true hit"
            assert row.snp_overlaps == "."
            assert row.entropy == 0
            assert pd.isna(row.snp_pam_dist) and pd.isna(row.snp_cut_dist)
    # The fixture must actually exercise both strands and the mismatch class,
    # otherwise these tests would pass against the region-based join.
    assert annotated_fwd > 0 and annotated_rvs > 0
    assert annotated_mismatched > 0, "no annotated guide disagrees with its region strand"


def test_annotation_is_guide_level(cheap_run):
    out, ref = cheap_run
    df, genome = load(out, ref)
    assert_matches_expectation(df, genome)


def test_distances_need_only_snps(cheap_run):
    out, ref = cheap_run
    df, _ = load(out, ref)
    assert "snp_pam_dist" in df.columns and "snp_cut_dist" in df.columns
    for col in _VARIANT_COLUMNS:
        assert col not in df.columns, f"{col} must require --rescore-snps"


def test_reference_scores_use_guide_orientation_flanks(cheap_run):
    """RS2 context bases must follow the guide's orientation, not the region's."""
    out, ref = cheap_run
    df, genome = load(out, ref)
    mismatch = ((df["Strand"] == "+") & (df["orientation"] == "RVS")) | \
               ((df["Strand"] == "-") & (df["orientation"] == "FWD"))
    assert mismatch.sum() > 0, "fixture exercises no orientation/region mismatch"
    for row in df.itertuples():
        g = int(row.Start) + int(row.start)
        f5, f3 = expected_flanks(genome, g, row.orientation)
        want = float(score_guides([row.target[:20]], flanks=[(f5, f3)],
                                  pams=[row.target[20:23]])[0])
        assert row.Doench2016OnTarget == pytest.approx(want, abs=1e-9), (
            f"{row.gene_id}@{g} ({row.Strand}/{row.orientation}): "
            f"RS2 {row.Doench2016OnTarget} != {want} with orientation-correct flanks")


def test_browser_bed_reports_guide_strand(cheap_run):
    out, ref = cheap_run
    df, _ = load(out, ref)
    want = Counter()
    for row in df.itertuples():
        strand = "-" if row.orientation == "RVS" else "+"
        want[(int(row.Start) + int(row.start), int(row.Start) + int(row.stop), strand)] += 1
    got = Counter()
    lines = (out / "final_sgrnas.bed").read_text().splitlines()[1:]
    for line in lines:
        chrom, start, end, _name, _score, strand, _rgb = line.split("\t")
        got[(int(start), int(end), strand)] += 1
    assert got == want


def test_variant_scores_match_reconstruction(rescore_run):
    """Single-SNP guides: avg/worst/all must equal the reconstructed variant."""
    out, ref = rescore_run
    df, genome = load(out, ref)
    positions = snp_positions()
    single = [r for r in df.itertuples() if r.num_snps == 1]
    assert single, "fixture produced no single-SNP guides"
    for row in single:
        g = int(row.Start) + int(row.start)
        rel = int(row.snp_overlaps.split(":")[0])
        alt = "G"  # ALT of the bed name chr1:<p>:A/G, in guide orientation below
        if row.orientation == "RVS":
            alt = revcomp(alt)
        guide, pam = row.target[:20], row.target[20:23]
        if rel < 20:
            guide = guide[:rel] + alt + guide[rel + 1:]
        else:  # SNP inside the PAM
            pam = pam[:rel - 20] + alt + pam[rel - 19:]
        f5, f3 = expected_flanks(genome, g, row.orientation)
        want = float(score_guides([guide], flanks=[(f5, f3)], pams=[pam])[0])
        assert row.Doench2016OnTarget_avg1snp == pytest.approx(want, abs=1e-9)
        assert row.Doench2016OnTarget_worst1snp == pytest.approx(want, abs=1e-9)
        assert row.Doench2016OnTarget_allsnps == pytest.approx(want, abs=1e-9)
        assert row.Doench2016_delta_avg == pytest.approx(
            want - row.Doench2016OnTarget, abs=1e-9)
        # sanity: the guide really is expected to carry that SNP
        assert expected_guide_hits(g, row.orientation, positions) == [rel]
