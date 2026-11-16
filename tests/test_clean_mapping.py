"""Tests for isoform-to-locus cleaning and secondary-column hygiene.

clean-mapping collapses isoform IDs to locus IDs. That collapse can make a
secondary list repeat the primary locus, which reads as a self-ortholog in
every downstream join, so the primary is dropped from the secondary.
"""

import pandas as pd

from crispytargets.mapping.clean_mapping import clean_locus_id, drop_primary


def test_clean_locus_id_collapses_isoform_to_locus():
    assert clean_locus_id("KY.Chr1.1.v1.nonSL13-1") == "KY.Chr1.1"
    assert clean_locus_id("KH.C2.208.v2.A.SL4-1") == "KH.C2.208"


def test_clean_locus_id_keeps_short_ids():
    assert clean_locus_id("KH.L28.8") == "KH.L28.8"
    assert clean_locus_id("KY.Chr1.1") == "KY.Chr1.1"


def test_clean_locus_id_deduplicates_a_pipe_list():
    value = "KY.Chr1.1.v1-1|KY.Chr1.1.v2-1|KY.Chr1.2.v1-1"
    assert clean_locus_id(value) == "KY.Chr1.1|KY.Chr1.2"


def test_drop_primary_removes_the_primary_locus():
    assert drop_primary("KH.C1.1|KH.C1.7", "KH.C1.1") == "KH.C1.7"


def test_drop_primary_leaves_empty_when_nothing_remains():
    assert drop_primary("KH.C1.1", "KH.C1.1") == ""
    assert drop_primary("KH.C1.1|KH.C1.1", "KH.C1.1") == ""


def test_drop_primary_keeps_missing_values_missing():
    assert pd.isna(drop_primary(pd.NA, "KH.C1.1"))
    # A missing primary cannot match anything, so the secondary survives.
    assert drop_primary("KH.C1.1", pd.NA) == "KH.C1.1"
