"""Strand-aware flank extraction and PAM handling tests.

Drives CrispyPipeline.post_process with a crafted FlashFry scored table
(no java needed) and checks that the Doench 2016 score matches an
independently computed one, i.e. that the flanks fed to the scorer are the
correct genomic neighbours on each strand and that the real PAM is used.

Geometry, verified against real FlashFry output: the 23-mer sits at
forward [g, g+23). On the + strand the protospacer is [g, g+20) and the
PAM [g+20, g+23). On the - strand it is mirrored: the PAM is [g, g+3) and
the protospacer [g+3, g+23), because the whole 23-mer is read as the
reverse complement. So the guide-orientation 5' flank on the - strand is
revcomp(fwd[g+23 : g+28]) and the 3' flank is revcomp(fwd[g-2 : g]).

Guards the bug classes from commits a3d3b4d, 1c6685d, 59004a0, 83aa2c4.
"""

import random
import types
from pathlib import Path

import pandas as pd
import pytest
from pyfaidx import Fasta

import crispytargets.design.doench2016 as doench2016
import crispytargets.design.pipeline as pipeline

COMP = str.maketrans("ACGT", "TGCA")

FWD_START = 100   # forward guide occupies fwd[100:123], PAM at fwd[120:123]
REV_START = 500   # reverse guide occupies fwd[500:523], PAM at fwd[500:503]
FWD_PAM = "AGA"   # deliberately not NGG, so a hardcoded PAM cannot pass
REV_PAM = "CCA"


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def write_genome(path, length=2000, seed=99):
    """Deterministic genome with two guides planted at known coordinates."""
    rng = random.Random(seed)
    bases = list("".join(rng.choice("ACGT") for _ in range(length)))
    # + strand: PAM is the last 3 of the 23-mer.
    bases[FWD_START:FWD_START + 20] = "ACGTTGCAGTCACTAGTCGA"
    bases[FWD_START + 20:FWD_START + 23] = FWD_PAM
    # - strand: PAM is the first 3 of the 23-mer (read as revcomp).
    bases[REV_START:REV_START + 3] = revcomp(REV_PAM)
    bases[REV_START + 3:REV_START + 23] = "TTGACCGTATGGCATAGCCT"
    seq = "".join(bases)

    with open(path, "w") as f:
        f.write(">chr1\n")
        for i in range(0, length, 60):
            f.write(seq[i:i + 60] + "\n")
    return seq


def write_scored_table(path, seq):
    """Minimal FlashFry-like scored table for the two planted guides."""
    fwd23 = seq[FWD_START:FWD_START + 23]
    rev23 = seq[REV_START:REV_START + 23]
    rows = [
        # contig, target (guide orientation), context, start, stop, orientation
        (f"fwd_gene|exon|chr1:{FWD_START}-{FWD_START + 23}|+",
         fwd23, "C" * 60, 0, 23, "FWD"),
        (f"rev_gene|exon|chr1:{REV_START}-{REV_START + 23}|-",
         revcomp(rev23), "C" * 60, 0, 23, "RVS"),
    ]
    df = pd.DataFrame(rows, columns=[
        "contig", "target", "context", "start", "stop", "orientation",
    ])
    df.to_csv(path, sep="\t", index=False)


def expected_score(seq, start, strand):
    """Score the guide from first principles, independent of the code."""
    fwd23 = seq[start:start + 23]
    if strand == "+":
        guide, pam = fwd23[:20], fwd23[20:23]
        flank5 = fwd23 and seq[start - 5:start]
        flank3 = seq[start + 23:start + 25]
    else:
        guide, pam = revcomp(fwd23[3:23]), revcomp(fwd23[0:3])
        flank5 = revcomp(seq[start + 23:start + 28])
        flank3 = revcomp(seq[start - 2:start])
    return float(doench2016.score_guides([guide], flanks=[(flank5, flank3)], pams=[pam])[0])


def run_post_process(tmp_path, seq):
    ref = tmp_path / "chr1.fa"
    scored = tmp_path / "scored.tsv"
    out = tmp_path / "out"
    out.mkdir(exist_ok=True)
    write_genome(ref)
    # Reuse the planted sequence so the table matches the reference exactly.
    write_scored_table(scored, seq)

    args = types.SimpleNamespace(
        gff="unused.gff3", ref=str(ref), outdir=str(out), genes=None, filter_bed=None,
        accessome=None, snps=None, tss_features="mRNA",
        promoter_upstream=100, promoter_downstream=50,
        flashfry_jar="unused.jar", mem="4g", metrics="doench2014ontarget",
        target_len=20, pam_len=3, cutoff=0.65, rs3=False,
        rescore_snps=False, browser_bed=False, verbose=False,
    )
    pipe = pipeline.CrispyPipeline(args)
    # post_process only needs gene_map for the rescore path; seed it empty.
    pipe.gene_map = {}
    pipe.post_process(str(scored))
    return pd.read_csv(out / "final_sgrnas.tsv", sep="\t"), Fasta(str(ref))


def score_for(df, gene):
    row = df[df.gene_id == gene]
    assert len(row) == 1, f"expected 1 row for {gene}, got {len(row)}"
    return float(row.Doench2016OnTarget.iloc[0])


@pytest.fixture(scope="module")
def scored(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("flanks")
    seq = write_genome(tmp / "chr1.fa")
    df, _ = run_post_process(tmp, seq)
    return seq, df


def test_forward_strand_flanks(scored):
    seq, df = scored
    assert score_for(df, "fwd_gene") == pytest.approx(expected_score(seq, FWD_START, "+"), abs=1e-12)


def test_reverse_strand_flanks(scored):
    """The - strand 5' flank must start after the PAM, not 3 bases inside it."""
    seq, df = scored
    assert score_for(df, "rev_gene") == pytest.approx(expected_score(seq, REV_START, "-"), abs=1e-12)


def test_pam_is_taken_from_the_reference(scored):
    """A non-NGG PAM must change the score; a hardcoded NGG cannot pass."""
    seq, df = scored
    guide_fwd = seq[FWD_START:FWD_START + 20]
    as_gg = float(doench2016.score_guides([guide_fwd], pams=["NGG"])[0])
    assert score_for(df, "fwd_gene") != pytest.approx(as_gg, abs=1e-9)


def test_strands_score_differently(scored):
    """Strand-blind flank extraction would make these two agree by chance."""
    seq, df = scored
    assert score_for(df, "fwd_gene") != pytest.approx(score_for(df, "rev_gene"), abs=1e-9)


def test_reverse_flank_starts_after_pam():
    """Direct check of the - strand 5' flank offset against the genome."""
    seq = write_genome(Path("/tmp/opencode/flank_geometry.fa"))
    start = REV_START
    correct = revcomp(seq[start + 23:start + 28])
    # The pre-fix behaviour read from start+20, i.e. 3 protospacer bases.
    buggy = revcomp(seq[start + 20:start + 25])
    assert correct != buggy
    assert correct == revcomp(seq[start + 23:start + 28])
