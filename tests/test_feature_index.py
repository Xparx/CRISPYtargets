"""feature_index tests: transcription-order block ordinals.

assign_feature_index numbers each distinct target region within its
(gene_id, feature) group so guides can be sorted first-to-last along a
gene: ascending coordinates on the + strand, descending on the - strand.
Promoters number separately from CDS/exon blocks, genes number
independently of each other, and every guide inside one block shares the
block's index.

Pure DataFrame logic: no java, FlashFry, or reference files needed.
"""

import pandas as pd

from crispytargets.design.utils import assign_feature_index


def _row(gene, feature, start, strand, end=None, chrom="Chr1"):
    return {
        'gene_id': gene, 'feature': feature, 'Chromosome': chrom,
        'Start': start, 'End': start + 50 if end is None else end,
        'Strand': strand, 'target': 'ACGTACGTACGTACGTACGTNGG',
    }


def _index_map(df):
    """(gene, feature, Start) -> feature_index for one row per region."""
    out = assign_feature_index(df)
    regions = out.drop_duplicates(
        subset=['gene_id', 'feature', 'Start'])[
        ['gene_id', 'feature', 'Start', 'feature_index']]
    return {(r.gene_id, r.feature, r.Start): int(r.feature_index)
            for r in regions.itertuples()}


def test_plus_strand_orders_by_ascending_start():
    df = pd.DataFrame([
        _row('G1', 'CDS', 300, '+'),
        _row('G1', 'CDS', 100, '+'),
        _row('G1', 'CDS', 200, '+'),
    ])
    out = assign_feature_index(df)
    got = dict(zip(out['Start'], out['feature_index']))
    assert got == {300: 3, 100: 1, 200: 2}


def test_minus_strand_orders_by_descending_start():
    df = pd.DataFrame([
        _row('G1', 'CDS', 100, '-'),
        _row('G1', 'CDS', 300, '-'),
        _row('G1', 'CDS', 200, '-'),
    ])
    out = assign_feature_index(df)
    got = dict(zip(out['Start'], out['feature_index']))
    # 5' end of a - strand gene is the HIGH coordinate end.
    assert got == {100: 3, 300: 1, 200: 2}


def test_guides_in_one_block_share_the_index():
    df = pd.DataFrame([
        _row('G1', 'CDS', 100, '+'),  # two guides from the same block
        _row('G1', 'CDS', 100, '+'),
        _row('G1', 'CDS', 200, '+'),
    ])
    out = assign_feature_index(df)
    assert list(out['feature_index']) == [1, 1, 2]
    assert len(out) == len(df)


def test_promoters_number_separately_from_cds():
    df = pd.DataFrame([
        _row('G1', 'CDS', 100, '+'),
        _row('G1', 'CDS', 200, '+'),
        _row('G1', 'promoter', 50, '+'),
        _row('G1', 'promoter', 40, '+'),
    ])
    idx = _index_map(df)
    assert idx == {
        ('G1', 'CDS', 100): 1, ('G1', 'CDS', 200): 2,
        ('G1', 'promoter', 50): 2, ('G1', 'promoter', 40): 1,
    }


def test_genes_number_independently():
    df = pd.DataFrame([
        _row('G1', 'CDS', 100, '+'),
        _row('G2', 'CDS', 100, '+'),
    ])
    out = assign_feature_index(df)
    assert list(out['feature_index']) == [1, 1]


def test_non_plus_non_minus_strand_ascending():
    df = pd.DataFrame([
        _row('G1', 'CDS', 200, '.'),
        _row('G1', 'CDS', 100, '.'),
    ])
    out = assign_feature_index(df)
    got = dict(zip(out['Start'], out['feature_index']))
    assert got == {200: 2, 100: 1}
