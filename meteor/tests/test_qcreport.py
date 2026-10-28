# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Test the mapped reads report of meteor merge"""

import json
from pathlib import Path

import pandas as pd
import pytest

from ..qcreport import QCReport


def _census(total, counted, msp_count, layout="2.0.22") -> dict:
    """census_stage_2 content: counts are in different sections across versions"""
    if layout == "2.0.22":
        return {
            "mapping": {"total_read_count": total, "mapped_read_count": total},
            "counting": {"counted_reads": counted},
            "profiling_stats": {"msp_count": msp_count, "msp_signal": 0.5},
        }
    return {
        "mapping": {"total_read_count": str(total), "counted_reads": str(counted)},
        "profiling": {"msp_count": str(msp_count), "msp_signal": "0.5"},
    }


CENSUS = {
    "s1": _census(100, 80, 10),
    "s2": _census(200, 50, 20, layout="2.0.23"),
    "s3": _census(0, 0, 0),
    "s4": {"mapping": {}},
}


def _mqc(tmp_path: Path, name: str) -> dict:
    return json.loads((tmp_path / f"run_{name}_mqc.json").read_text())


def test_sample_qc(tmp_path: Path) -> None:
    qc = QCReport(tmp_path, "run").sample_qc(CENSUS)
    assert qc.columns.tolist() == QCReport.QC_COLUMNS
    assert qc["sample"].tolist() == ["s1", "s2", "s3", "s4"]
    assert qc["total_read_count"].tolist()[:3] == [100, 200, 0]
    assert qc["mapped_reads"].tolist()[:3] == [80, 50, 0]
    assert qc["mapped_read_rate"].tolist()[:2] == [80.0, 25.0]
    assert qc.loc[2:, "mapped_read_rate"].isna().all()  # no read, no census
    assert qc["msp_richness"].tolist()[:3] == [10, 20, 0]
    assert qc["msp_signal_rate"].tolist()[:3] == [50.0, 50.0, 50.0]


MSP = pd.DataFrame(
    {
        "msp_name": ["msp_1", "msp_2", "msp_3", "msp_4"],
        "s1": [1.0, 1.0, 2.0, None],
        "s2": [0.0, 3.0, None, 1.0],
        "s3": [0.0, 0.0, 0.0, 0.0],
    }
)
TAXONOMY = pd.DataFrame(
    {
        "msp_name": ["msp_1", "msp_2", "msp_3"],
        "phylum": ["Bacillota", "Bacteroidota", "unclassified "],
    }
)


def test_execute(tmp_path: Path, caplog) -> None:
    caplog.set_level("INFO")
    QCReport(tmp_path, "run").execute(CENSUS)
    assert "no phylum composition" in caplog.text
    assert not (tmp_path / "run_phylum.tsv").exists()
    qc = pd.read_csv(tmp_path / "run_qc.tsv", sep="\t")
    assert qc.columns.tolist() == QCReport.QC_COLUMNS
    assert len(qc) == 4
    general = _mqc(tmp_path, "general_stats")
    assert general["plot_type"] == "generalstats"
    assert [next(iter(header)) for header in general["headers"]] == [
        "total_read_count", "mapped_reads", "mapped_read_rate", "msp_richness", "msp_signal_rate",
    ]
    assert general["data"]["s1"] == {
        "total_read_count": 100, "mapped_reads": 80, "mapped_read_rate": 80.0,
        "msp_richness": 10, "msp_signal_rate": 50.0,
    }
    assert general["data"]["s3"]["mapped_read_rate"] is None
    reads = _mqc(tmp_path, "reads")
    assert reads["plot_type"] == "bargraph"
    assert reads["parent_id"] == "meteor"
    assert reads["data"] == {
        "s1": {"Mapped": 80, "Not mapped": 20},
        "s2": {"Mapped": 50, "Not mapped": 150},
        "s3": {"Mapped": 0, "Not mapped": 0},
    }
    richness = _mqc(tmp_path, "richness")
    assert richness["plot_type"] == "scatter"
    # samples without reads or without census are not plotted
    assert richness["data"] == {"s1": {"x": 100, "y": 10}, "s2": {"x": 200, "y": 20}}


def test_phylum_composition(tmp_path: Path) -> None:
    composition = QCReport(tmp_path, "run").phylum_composition(MSP, TAXONOMY)
    # s3 has no species; msp_4 is not in the taxonomy; phyla by mean abundance
    assert composition.index.tolist() == ["s1", "s2"]
    assert composition.columns.tolist() == ["Bacteroidota", "unclassified", "Bacillota"]
    assert composition.loc["s1"].to_dict() == {"Bacteroidota": 25.0, "unclassified": 50.0, "Bacillota": 25.0}
    assert composition.loc["s2"].to_dict() == {"Bacteroidota": 75.0, "unclassified": 25.0, "Bacillota": 0.0}


def test_write_composition(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(QCReport, "TOP_PHYLA", 2)
    QCReport(tmp_path, "run").execute(CENSUS, MSP, TAXONOMY)
    table = pd.read_csv(tmp_path / "run_phylum.tsv", sep="\t", index_col=0)
    assert table.columns.tolist() == ["Bacteroidota", "unclassified", "Bacillota"]
    phylum = _mqc(tmp_path, "phylum")
    assert phylum["plot_type"] == "bargraph"
    assert phylum["data"]["s1"] == {"Bacteroidota": 25.0, "unclassified": 50.0, "Other phyla": 25.0}
    for values in phylum["data"].values():
        assert sum(values.values()) == pytest.approx(100.0)

