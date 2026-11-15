"""Tests for the KH2013 dictionary QA gate.

The gate decides whether a cluster run may be accepted, so its arithmetic
and its failure behaviour are checked here: a preserved backbone passes, a
changed backbone fails, and a KH2013 hit that lands in the secondary column
still counts as agreement with GHOST.
"""

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

TOOL = Path(__file__).resolve().parents[1] / "tools" / "qa_kh2013_mapping.py"


def load_tool():
    spec = importlib.util.spec_from_file_location("qa_kh2013_mapping", TOOL)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_inputs(tmp_path, primary="KH.C1.1", secondary=None, base_change=False):
    base = pd.DataFrame({
        "gene_id": ["KY21.Chr1.1", "KY21.Chr1.2", "KY21.Chr1.3"],
        "KY19_primary": ["KY.Chr1.1", "KY.Chr1.2", "KY.Chr1.3"],
        "Human_primary": ["FGF13_HUMAN", None, "TAF9_HUMAN"],
    })
    extended = base.copy()
    if base_change:
        extended.loc[0, "Human_primary"] = "OTHER_HUMAN"
    extended["KH2013_primary"] = [primary, None, None]
    extended["KH2013_secondary"] = [secondary, None, None]

    base_path = tmp_path / "base.tsv"
    ext_path = tmp_path / "extended.tsv"
    ghost_path = tmp_path / "ghost.txt"
    report_path = tmp_path / "report.tsv"
    projection_path = tmp_path / "projection.tsv"

    base.to_csv(base_path, sep="\t", index=False)
    extended.to_csv(ext_path, sep="\t", index=False)
    ghost_path.write_text(
        "Mapping of KH to KY21 [Eval>=1e-5,Identity>=97%]\n"
        "KH.C1.1\tKY21.Chr1.1\n"
        "KH.C1.9\tKY21.Chr1.2\n")
    return {
        "base": base_path, "extended": ext_path, "ghost": ghost_path,
        "report": report_path, "projection": projection_path,
    }


def run_tool(module, paths, monkeypatch):
    argv = ["qa_kh2013_mapping.py",
            "--base", str(paths["base"]),
            "--extended", str(paths["extended"]),
            "--ghost", str(paths["ghost"]),
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
    assert int(report["kh2013_primary_nonnull"]) == 1
    assert float(report["kh2013_primary_coverage"]) == pytest.approx(1 / 3, abs=1e-4)
    assert float(report["ghost_pair_recall"]) == pytest.approx(0.5, abs=1e-4)

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


def test_secondary_hit_counts_as_ghost_agreement(tmp_path, monkeypatch):
    module = load_tool()
    paths = write_inputs(tmp_path, primary="KH.C1.7", secondary="KH.C1.1|KH.C1.4")
    assert run_tool(module, paths, monkeypatch) == 0

    report = pd.read_csv(paths["report"], sep="\t").set_index("metric")["value"]
    assert float(report["ghost_pair_recall"]) == pytest.approx(0.5, abs=1e-4)
    assert float(report["ghost_primary_precision"]) == pytest.approx(0.0, abs=1e-4)
