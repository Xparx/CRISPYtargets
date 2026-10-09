"""SNP liftover round-trip tests for crispytargets.mapping.lift_snps.

Builds a tiny old/new assembly pair that differs by a known 100bp
insertion, lifts a VCF through the real minimap2 + pysam path, and checks
the positions land where the insertion dictates. The block-mapping maths
is also pinned directly from a hand-written PAF, which needs no minimap2.

Skip conditions: the round trip needs minimap2; the PAF test needs only
pysam (imported at module load by lift_snps itself).
"""

import os
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

from crispytargets.mapping.lift_snps import AlignmentBlockMapper, blocks_from_cigar

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC = REPO_ROOT / "src"

pytestmark = pytest.mark.skipif(
    shutil.which("minimap2") is None, reason="minimap2 not available"
)

CHR_LEN = 5000
INSERT_AT = 2500     # old-coordinate boundary of the insertion
INSERT_LEN = 100
BEFORE_POS = 500     # 1-based old position before the insertion
AFTER_POS = 4000     # 1-based old position after the insertion


def make_pair(tmp_path, seed=7):
    """Write old.fa and new.fa, where new = old with a 100bp insertion."""
    rng = random.Random(seed)
    old = "".join(rng.choice("ACGT") for _ in range(CHR_LEN))
    insertion = "".join(rng.choice("ACGT") for _ in range(INSERT_LEN))
    new = old[:INSERT_AT] + insertion + old[INSERT_AT:]

    old_fa = tmp_path / "old.fa"
    new_fa = tmp_path / "new.fa"
    for path, seq in ((old_fa, old), (new_fa, new)):
        with open(path, "w") as f:
            f.write(">chr1\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")
    return old_fa, new_fa, old


def write_vcf(path, old):
    """VCF on old-assembly coordinates: one SNP each side of the insertion."""
    header = pysam.VariantHeader()
    header.add_line("##fileformat=VCFv4.2")
    header.add_line("##contig=<ID=chr1>")
    header.add_sample("S1")
    with pysam.VariantFile(str(path), "w", header=header) as vcf:
        for pos in (BEFORE_POS, AFTER_POS):
            ref = old[pos - 1]
            alt = "A" if ref != "A" else "C"
            rec = vcf.new_record(contig="chr1", start=pos - 1, stop=pos,
                                 alleles=(ref, alt), id=f"rs{pos}")
            vcf.write(rec)


def lift(tmp_path):
    """Run the documented CLI against the constructed assemblies."""
    old_fa, new_fa, old = make_pair(tmp_path)
    vcf = tmp_path / "snps.vcf"
    out = tmp_path / "lifted"
    write_vcf(vcf, old)

    proc = subprocess.run(
        [sys.executable, "-m", "crispytargets", "lift",
         "--old-genome", str(old_fa),
         "--new-genome", str(new_fa),
         "--vcf", str(vcf),
         "--output", str(out),
         "--bed", "--threads", "2"],
        capture_output=True, text=True,
        env={**os.environ, "PYTHONPATH": str(SRC)},
    )
    assert proc.returncode == 0, f"lift failed:\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"

    out_vcf = out / "snps_lifted.vcf.gz"
    assert out_vcf.exists(), f"missing {out_vcf}"
    with pysam.VariantFile(str(out_vcf)) as vf:
        recs = {r.id: r for r in vf}
    return recs, out


def test_positions_lift_across_known_insertion(tmp_path):
    recs, _ = lift(tmp_path)
    assert set(recs) == {f"rs{BEFORE_POS}", f"rs{AFTER_POS}"}
    # Before the insertion: unchanged.
    assert recs[f"rs{BEFORE_POS}"].pos == BEFORE_POS
    # After it: shifted by exactly the insertion length.
    assert recs[f"rs{AFTER_POS}"].pos == AFTER_POS + INSERT_LEN


def test_lifted_records_keep_alleles_and_contig(tmp_path):
    recs, _ = lift(tmp_path)
    for rid, rec in recs.items():
        assert rec.chrom == "chr1"
        assert rec.alts is not None and len(rec.alts) == 1


def test_bed_matches_lifted_positions(tmp_path):
    recs, out = lift(tmp_path)
    lines = (out / "snps_lifted.bed").read_text().splitlines()
    rows = [l.split("\t") for l in lines if not l.startswith("#")]
    assert len(rows) == 2
    # BED start is the lifted 1-based position minus one.
    starts = sorted(int(r[1]) for r in rows)
    assert starts == sorted(rec.pos - 1 for rec in recs.values())


def test_paf_block_mapping_maths(tmp_path):
    """Pin AlignmentBlockMapper coordinate math without minimap2.

    Hand-written PAF: query is the OLD genome, target the NEW genome
    (lift_snps aligns new as reference, old as query).
    """

    paf = tmp_path / "hand.paf"
    # qname qlen qstart qend strand tname tlen tstart tend ...
    paf.write_text(
        "chr1\t5000\t0\t2500\t+\tchr1\t5100\t0\t2500\t2500\t2500\t60\n"
        "chr1\t5000\t2500\t5000\t+\tchr1\t5100\t2600\t5100\t2500\t2500\t60\n"
    )
    mapper = AlignmentBlockMapper()
    mapper.load_from_paf(str(paf))

    assert mapper.map_position("chr1", 1) == ("chr1", 1)
    assert mapper.map_position("chr1", 2000) == ("chr1", 2000)
    # Old 2500 is the last base BEFORE the insertion, so it is unaffected.
    assert mapper.map_position("chr1", 2500) == ("chr1", 2500)
    # Old 2501 is the first base of the second block, shifted by the insertion.
    assert mapper.map_position("chr1", 2501) == ("chr1", 2601)
    assert mapper.map_position("chr1", 5000) == ("chr1", 5100)
    # A position on a chromosome with no alignment fails cleanly.
    assert mapper.map_position("chrZ", 100) == (None, None)


def test_cigar_splits_record_into_blocks(tmp_path):
    """A single record with an internal insertion must yield two blocks.

    This is the real minimap2 output for a 100bp insertion: one record
    carrying cg:Z:2500M100D2500M, not two records. Treating it as one
    collinear block is the bug that shifted every downstream position.
    """

    paf = tmp_path / "cigar.paf"
    paf.write_text(
        "chr1\t5000\t0\t5000\t+\tchr1\t5100\t0\t5100\t5000\t5081\t60"
        "\tcg:Z:2500M100D2500M\n"
    )
    mapper = AlignmentBlockMapper()
    mapper.load_from_paf(str(paf))

    assert mapper.records_without_cigar == 0
    # Before the insertion: unchanged.
    assert mapper.map_position("chr1", 1) == ("chr1", 1)
    assert mapper.map_position("chr1", 2500) == ("chr1", 2500)
    # After it: shifted by the 100bp insertion, which linear interpolation
    # over the whole record would have missed.
    assert mapper.map_position("chr1", 2501) == ("chr1", 2601)
    assert mapper.map_position("chr1", 5000) == ("chr1", 5100)


def test_cigar_handles_insertion_and_deletion():
    """blocks_from_cigar: I advances query, D advances target."""

    # 100 target-only bases (D) after 2500 aligned.
    blocks = list(blocks_from_cigar(0, "chr1", 0, "2500M100D2500M"))
    assert blocks == [
        (0, 2500, "chr1", 0, 2500),
        (2500, 5000, "chr1", 2600, 5100),
    ]

    # 50 query-only bases (I) after 1000 aligned.
    blocks = list(blocks_from_cigar(0, "chr1", 0, "1000M50I1000M"))
    assert blocks == [
        (0, 1000, "chr1", 0, 1000),
        (1050, 2050, "chr1", 1000, 2000),
    ]

    # Soft clip consumes query but emits no block.
    blocks = list(blocks_from_cigar(5, "chr1", 5, "10S100M"))
    assert blocks == [(5, 105, "chr1", 5, 105)]

    # Adjacent aligned runs (M/=/X) each emit a block; they are contiguous,
    # so lookups give the same answer as one merged block.
    blocks = list(blocks_from_cigar(0, "chr1", 0, "10M5=3X10M"))
    assert blocks == [
        (0, 10, "chr1", 0, 10),
        (10, 15, "chr1", 10, 15),
        (15, 18, "chr1", 15, 18),
        (18, 28, "chr1", 18, 28),
    ]


def test_record_without_cigar_warns_and_falls_back(tmp_path, caplog):
    """A record with no cg:Z: logs a warning and keeps the old behaviour."""
    import logging as _logging


    paf = tmp_path / "nocigar.paf"
    paf.write_text("chr1\t5000\t0\t5000\t+\tchr1\t5100\t0\t5100\t5000\t5081\t60\n")
    mapper = AlignmentBlockMapper()
    with caplog.at_level(_logging.WARNING):
        mapper.load_from_paf(str(paf))

    assert mapper.records_without_cigar == 1
    assert "minimap2 -c" in caplog.text
    assert mapper.map_position("chr1", 2500) == ("chr1", 2500)


def test_paf_mapping_rejects_unmapped_positions(tmp_path):

    paf = tmp_path / "gapped.paf"
    # One block covering only [0, 1000) of the old chromosome.
    paf.write_text("chr1\t5000\t0\t1000\t+\tchr1\t5100\t0\t1000\t1000\t1000\t60\n")
    mapper = AlignmentBlockMapper()
    mapper.load_from_paf(str(paf))

    assert mapper.map_position("chr1", 500) == ("chr1", 500)
    assert mapper.map_position("chr1", 1001) == (None, None)
