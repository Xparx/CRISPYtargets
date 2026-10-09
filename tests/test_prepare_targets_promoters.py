"""Promoter extraction strand-handling tests for CrispyPipeline.prepare_targets.

Regression for the crash where a GFF3 whose TSS features appear on only one
strand (or not at all) raised ``KeyError`` from pandas internals. PyRanges
returns a DataFrame with zero columns for an empty subset, so ``tpts.df['Strand']``
and ``neg_df['End']`` both failed. Any single-strand annotation hit this:
draft or region-filtered GFF3s, organellar genomes, and any GFF3 whose
feature labels do not match ``--tss-features``.

These tests call prepare_targets directly and need no java or FlashFry.
"""

import logging
import types

import pyranges as pr

import crispytargets.design.pipeline as pipeline

UPSTREAM = 100
DOWNSTREAM = 50
# prepare_targets extends every target by this before returning, so promoter
# coordinates as returned are shifted by the same amount.
FLANK = max(UPSTREAM, DOWNSTREAM, 100)


def write_gff3(path, rows):
    path.write_text("".join(
        f"chr1\ttest\t{feat}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n"
        for feat, start, end, strand, attrs in rows
    ))


def prepare(tmp_path, rows, tss_features="mRNA"):
    """Run prepare_targets against a GFF3 of the given feature rows."""
    gff = tmp_path / "test.gff3"
    write_gff3(gff, rows)
    out = tmp_path / "out"
    out.mkdir(exist_ok=True)
    args = types.SimpleNamespace(
        gff=str(gff), ref="unused.fa", outdir=str(out), genes=None, filter_bed=None,
        accessome=None, snps=None, tss_features=tss_features,
        promoter_upstream=UPSTREAM, promoter_downstream=DOWNSTREAM,
        flashfry_jar="unused.jar", mem="4g", metrics="doench2014ontarget",
        target_len=20, pam_len=3,
    )
    return pipeline.CrispyPipeline(args).prepare_targets()


def promoters(targets):
    df = targets.df
    return df[df.Feature == "promoter"] if len(df) else df


def expected_promoter(gff_path, strand):
    """Restate the promoter rule independently of the code under test.

    TSS is Start on the plus strand and End on the minus strand; upstream
    and downstream swap sides accordingly. Coordinates come from the same
    PyRanges reader so the GFF3 1-based to 0-based conversion is shared.
    """
    mrna = pr.read_gff3(str(gff_path)).df
    row = mrna[mrna.Feature == "mRNA"]
    row = row[row.Strand == strand].iloc[0]
    if strand == "+":
        start, end = row.Start - UPSTREAM, row.Start + DOWNSTREAM
    else:
        start, end = row.End - DOWNSTREAM, row.End + UPSTREAM
    return max(0, start - FLANK), end + FLANK


def test_plus_strand_only_builds_promoter(tmp_path):
    targets = prepare(tmp_path, [
        ("gene", 200, 900, "+", "ID=gene1"),
        ("mRNA", 200, 900, "+", "ID=tx1;Parent=gene1"),
        ("exon", 200, 900, "+", "ID=exon1;Parent=tx1"),
    ])
    prom = promoters(targets)
    assert len(prom) == 1
    assert prom.Strand.iloc[0] == "+"
    expected = expected_promoter(tmp_path / "test.gff3", "+")
    assert (int(prom.Start.iloc[0]), int(prom.End.iloc[0])) == expected


def test_minus_strand_only_builds_promoter(tmp_path):
    targets = prepare(tmp_path, [
        ("gene", 200, 900, "-", "ID=gene1"),
        ("mRNA", 200, 900, "-", "ID=tx1;Parent=gene1"),
        ("exon", 200, 900, "-", "ID=exon1;Parent=tx1"),
    ])
    prom = promoters(targets)
    assert len(prom) == 1
    assert prom.Strand.iloc[0] == "-"
    expected = expected_promoter(tmp_path / "test.gff3", "-")
    assert (int(prom.Start.iloc[0]), int(prom.End.iloc[0])) == expected


def test_minus_promoter_is_not_mirrored_from_plus(tmp_path):
    """Minus promoter must sit at the transcript's End, not its Start.

    This is the swap that strand-blind code would produce; without it the
    two strand rules above could both pass if a shared mistake were made.
    """
    targets = prepare(tmp_path, [
        ("gene", 200, 900, "-", "ID=gene1"),
        ("mRNA", 200, 900, "-", "ID=tx1;Parent=gene1"),
        ("exon", 200, 900, "-", "ID=exon1;Parent=tx1"),
    ])
    prom = promoters(targets)
    start, end = int(prom.Start.iloc[0]), int(prom.End.iloc[0])
    # Anchored at End=900 (0-based 900): promoter straddles 900 and extends
    # further upstream than downstream.
    assert start < 900 < end
    assert 900 - start == DOWNSTREAM + FLANK
    assert end - 900 == UPSTREAM + FLANK


def test_both_strands_build_promoters_independently(tmp_path):
    rows = [
        ("gene", 200, 900, "+", "ID=gene1"),
        ("mRNA", 200, 900, "+", "ID=tx1;Parent=gene1"),
        ("exon", 200, 900, "+", "ID=exon1;Parent=tx1"),
        ("gene", 2000, 2700, "-", "ID=gene2"),
        ("mRNA", 2000, 2700, "-", "ID=tx2;Parent=gene2"),
        ("exon", 2000, 2700, "-", "ID=exon2;Parent=tx2"),
    ]
    targets = prepare(tmp_path, rows)
    prom = promoters(targets)
    assert len(prom) == 2
    assert set(prom.Strand) == {"+", "-"}

    gff = tmp_path / "test.gff3"
    for strand in ("+", "-"):
        got = prom[prom.Strand == strand].iloc[0]
        assert (int(got.Start), int(got.End)) == expected_promoter(gff, strand)


def test_unmatched_tss_features_degrade_to_exons_with_warning(tmp_path, caplog):
    """A --tss-features/GFF3 label mismatch must warn, not raise KeyError."""
    rows = [
        ("gene", 200, 900, "+", "ID=gene1"),
        ("transcript", 200, 900, "+", "ID=tx1;Parent=gene1"),
        ("exon", 200, 900, "+", "ID=exon1;Parent=tx1"),
    ]
    with caplog.at_level(logging.WARNING):
        targets = prepare(tmp_path, rows, tss_features="mRNA")
    assert len(targets) > 0
    assert "promoter" not in set(targets.df.Feature)
    assert "No TSS features matched" in caplog.text
    assert "--tss-features" in caplog.text


def test_exons_only_gff3_still_yields_targets(tmp_path):
    """No TSS match at all must still return the exon-derived targets."""
    targets = prepare(tmp_path, [
        ("gene", 200, 900, "+", "ID=gene1"),
        ("CDS", 300, 700, "+", "ID=cds1;Parent=gene1"),
    ])
    assert len(targets) == 1
    assert set(targets.df.Feature) == {"CDS"}
