"""Tests for the species mapping QA gate.

The gate decides whether a cluster run may be accepted, so its arithmetic
and its failure behaviour are checked here: a preserved backbone passes, a
changed backbone fails, and a species hit that lands in the secondary column
still counts as agreement with the reference. One case runs a non-KH species
label to keep the column naming generic.
"""

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

TOOL = Path(__file__).resolve().parents[1] / "tools" / "qa_species_mapping.py"


def load_tool():
    spec = importlib.util.spec_from_file_location("qa_species_mapping", TOOL)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_inputs(tmp_path, species="KH2013", primary="KH.C1.1", secondary=None,
                 base_change=False):
    reference_hit = {"KH2013": "KH.C1.1",
                     "ENSEMBL": "ENSCING00000004107"}.get(species)
    base = pd.DataFrame({
        "gene_id": ["KY21.Chr1.1", "KY21.Chr1.2", "KY21.Chr1.3"],
        "KY19_primary": ["KY.Chr1.1", "KY.Chr1.2", "KY.Chr1.3"],
        "Human_primary": ["FGF13_HUMAN", None, "TAF9_HUMAN"],
    })
    extended = base.copy()
    if base_change:
        extended.loc[0, "Human_primary"] = "OTHER_HUMAN"
    extended[f"{species}_primary"] = [primary, None, None]
    extended[f"{species}_secondary"] = [secondary, None, None]

    base_path = tmp_path / "base.tsv"
    ext_path = tmp_path / "extended.tsv"
    reference_path = tmp_path / "reference.tsv"
    report_path = tmp_path / "report.tsv"
    projection_path = tmp_path / "projection.tsv"

    base.to_csv(base_path, sep="\t", index=False)
    extended.to_csv(ext_path, sep="\t", index=False)
    pd.DataFrame({
        "gene_id": ["KY21.Chr1.1", "KY21.Chr1.2"],
        species: [reference_hit, "KH.C1.9"],
    }).to_csv(reference_path, sep="\t", index=False)
    return {
        "base": base_path, "extended": ext_path, "reference": reference_path,
        "report": report_path, "projection": projection_path,
        "species": species,
    }


def run_tool(module, paths, monkeypatch):
    argv = ["qa_species_mapping.py",
            "--base", str(paths["base"]),
            "--extended", str(paths["extended"]),
            "--species", paths["species"],
            "--reference", str(paths["reference"]),
            "--projection", str(paths["projection"]),
            "--report", str(paths["report"])]
    monkeypatch.setattr(sys, "argv", argv)
    return module.main()


def test_preserved_backbone_passes_and_writes_report(tmp_path, monkeypatch):
    module = load_tool()
    paths = write_inputs(tmp_path)
    assert run_tool(module, paths, monkeypatch) == 0

    report = pd.read_csv(paths["report"], sep="\t").set_index("metric")["value"]
    assert int(report["rows_extended"]) == 3
    assert int(report["base_column_mismatches"]) == 0
    assert int(report["kh2013_primary_rows"]) == 1
    assert float(report["kh2013_primary_row_coverage"]) == pytest.approx(1 / 3, abs=1e-4)
    assert int(report["kh2013_distinct_primary"]) == 1
    assert float(report["reference_pair_recall"]) == pytest.approx(0.5, abs=1e-4)

    projection = pd.read_csv(paths["projection"], sep="\t")
    assert list(projection.columns) == [
        "gene_id", "KH2013_primary", "KH2013_secondary"]
    assert len(projection) == 3


def test_changed_base_column_fails_the_gate(tmp_path, monkeypatch):
    module = load_tool()
    paths = write_inputs(tmp_path, base_change=True)
    with pytest.raises(SystemExit):
        run_tool(module, paths, monkeypatch)

    report = pd.read_csv(paths["report"], sep="\t").set_index("metric")["value"]
    assert int(report["base_column_mismatches"]) == 1


def test_secondary_hit_counts_as_reference_agreement(tmp_path, monkeypatch):
    module = load_tool()
    paths = write_inputs(tmp_path, primary="KH.C1.7", secondary="KH.C1.1|KH.C1.4")
    assert run_tool(module, paths, monkeypatch) == 0

    report = pd.read_csv(paths["report"], sep="\t").set_index("metric")["value"]
    assert float(report["reference_pair_recall"]) == pytest.approx(0.5, abs=1e-4)
    assert float(report["reference_primary_precision"]) == pytest.approx(0.0, abs=1e-4)


def test_non_kh_species_uses_its_own_columns(tmp_path, monkeypatch):
    module = load_tool()
    paths = write_inputs(tmp_path, species="ENSEMBL",
                         primary="ENSCING00000004107")
    assert run_tool(module, paths, monkeypatch) == 0

    report = pd.read_csv(paths["report"], sep="\t").set_index("metric")["value"]
    assert int(report["ensembl_primary_rows"]) == 1
    assert float(report["reference_pair_recall"]) == pytest.approx(0.5, abs=1e-4)

    projection = pd.read_csv(paths["projection"], sep="\t")
    assert list(projection.columns) == [
        "gene_id", "ENSEMBL_primary", "ENSEMBL_secondary"]
