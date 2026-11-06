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

"""Test meteor resistance_merge"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ..resistancemerging import ResistanceMerging
from ..session import Component

PARAMETERS = {
    "annotations": ["resfinder", "mustard", "kegg"],
    "min_depth": 10,
    "min_alt_count": 3,
    "min_frequency": 0.1,
    "min_gene_coverage": 0.8,
    "base_quality": 20,
    "core_size": 100,
}


def _sample(root: Path, name: str, rows: list[dict], reference: str = "cat", **changed) -> None:
    """Output of meteor resistance for one sample"""
    sample_dir = root / name
    sample_dir.mkdir(parents=True)
    columns = {
        "msp_name": "M1",
        "resfinder": "tet(W)",
        "resfinderfg": pd.NA,
        "mustard": pd.NA,
        "mustard_prediction": pd.NA,
        "target": pd.NA,
        "gene_length": 900,
        "alignments": 100,
        "interpretable": True,
        "major_per_kb": 1.0,
        "minor_per_kb": 0.0,
        "control_genes": 2,
        "control_major_per_kb": 2.0,
        "control_minor_per_kb": 0.5,
    }
    table = pd.DataFrame([{**columns, **row} for row in rows])
    table = table[["gene_id", *columns]]
    table.to_csv(sample_dir / f"{name}_resistance.tsv", sep="\t", index=False, na_rep="NA")
    (sample_dir / f"{name}_census_stage_resistance.json").write_text(
        json.dumps(
            {
                "meteor_version": "2.0.23",
                "sample_info": {"sample_name": name},
                "mapping": {"reference_name": reference},
                "resistance": {
                    **PARAMETERS,
                    **changed,
                    "resistance_file": f"{name}_resistance.tsv",
                },
            }
        )
    )


def _merging(tmp_path: Path, sample_file: Path | None = None) -> ResistanceMerging:
    meteor = Component
    meteor.resistance_dir = tmp_path / "resistance"
    meteor.merging_dir = tmp_path / "merged"
    return ResistanceMerging(meteor, "out", sample_file)


@pytest.fixture
def three_samples(tmp_path: Path) -> Path:
    root = tmp_path / "resistance"
    # gene 1: interpretable in S1 and S2, more variants than its controls in S1 only
    # gene 2 (gyrA, no MSP): reads in S1 (not interpretable) and S3
    _sample(root, "S1", [
        {"gene_id": 1, "major_per_kb": 4.0, "minor_per_kb": 1.0},
        {"gene_id": 2, "msp_name": pd.NA, "resfinder": pd.NA, "target": "gyrA",
         "interpretable": False, "control_genes": 0,
         "control_major_per_kb": np.nan, "control_minor_per_kb": np.nan},
    ])
    _sample(root, "S2", [{"gene_id": 1, "major_per_kb": 1.0, "minor_per_kb": 0.0}])
    _sample(root, "S3", [
        {"gene_id": 2, "msp_name": pd.NA, "resfinder": pd.NA, "target": "gyrA",
         "major_per_kb": 6.0, "control_genes": 0,
         "control_major_per_kb": np.nan, "control_minor_per_kb": np.nan},
    ])
    return root


def test_execute(three_samples: Path, tmp_path: Path) -> None:
    _merging(tmp_path).execute()
    out = tmp_path / "merged"
    long = pd.read_table(out / "out_resistance.tsv")
    assert list(long.columns[:2]) == ["sample", "gene_id"]
    assert list(zip(long["sample"], long["gene_id"])) == [("S1", 1), ("S1", 2), ("S2", 1), ("S3", 2)]
    major = pd.read_table(out / "out_resistance_major_per_kb.tsv", index_col=0)
    assert list(major.columns) == ["S1", "S2", "S3"]
    assert major.loc[1, "S1"] == 4.0 and major.loc[1, "S2"] == 1.0
    assert pd.isna(major.loc[1, "S3"])  # no reads
    assert pd.isna(major.loc[2, "S1"]) and major.loc[2, "S3"] == 6.0  # not interpretable in S1
    interpretable = pd.read_table(out / "out_resistance_interpretable.tsv", index_col=0)
    assert interpretable.to_dict("index") == {1: {"S1": 1, "S2": 1, "S3": 0}, 2: {"S1": 0, "S2": 0, "S3": 1}}
    genes = pd.read_table(out / "out_resistance_genes.tsv").set_index("gene_id")
    assert genes.loc[1, "resfinder"] == "tet(W)" and genes.loc[2, "target"] == "gyrA"
    assert genes.loc[1, "samples_with_reads"] == 2 and genes.loc[1, "samples_interpretable"] == 2
    assert genes.loc[2, "samples_with_reads"] == 2 and genes.loc[2, "samples_interpretable"] == 1
    assert genes.loc[1, "interpretable_fraction"] == pytest.approx(2 / 3)
    assert genes.loc[1, "median_major_per_kb"] == pytest.approx(2.5)
    assert genes.loc[1, "median_control_major_per_kb"] == pytest.approx(2.0)
    # S1: 4 > 2, S2: 1 < 2
    assert genes.loc[1, "samples_with_controls"] == 2
    assert genes.loc[1, "major_above_controls_fraction"] == pytest.approx(0.5)
    assert genes.loc[1, "minor_above_controls_fraction"] == pytest.approx(0.5)
    # gene 2 has no controls (no MSP)
    assert genes.loc[2, "samples_with_controls"] == 0
    assert pd.isna(genes.loc[2, "major_above_controls_fraction"])
    assert genes.loc[2, "median_major_per_kb"] == pytest.approx(6.0)
    census = json.loads((out / "out_resistance_census.json").read_text())["resistance_merging"]
    assert census["samples"] == 3 and census["genes"] == 2 and census["interpretable_genes"] == 2
    assert census["reference_name"] == "cat" and census["min_depth"] == 10


def test_execute_sample_list(three_samples: Path, tmp_path: Path) -> None:
    sample_file = tmp_path / "samples.txt"
    sample_file.write_text("# selection\nS3\nS1\nS9\n")
    _merging(tmp_path, sample_file).execute()
    major = pd.read_table(tmp_path / "merged" / "out_resistance_major_per_kb.tsv", index_col=0)
    assert list(major.columns) == ["S3", "S1"]
    long = pd.read_table(tmp_path / "merged" / "out_resistance.tsv")
    assert long["sample"].unique().tolist() == ["S3", "S1"]


@pytest.mark.parametrize(
    "changed",
    [{"reference": "other"}, {"min_depth": 5}, {"annotations": ["resfinder"]}],
    ids=["catalogue", "min_depth", "annotations"],
)
def test_execute_incompatible_samples(tmp_path: Path, changed: dict) -> None:
    root = tmp_path / "resistance"
    _sample(root, "S1", [{"gene_id": 1}])
    _sample(root, "S2", [{"gene_id": 1}], **changed)
    with pytest.raises(SystemExit):
        _merging(tmp_path).execute()


def test_execute_errors(tmp_path: Path) -> None:
    (tmp_path / "resistance").mkdir()
    with pytest.raises(SystemExit):  # no sample
        _merging(tmp_path).execute()
    _sample(tmp_path / "resistance", "S1", [{"gene_id": 1}])
    _sample(tmp_path / "resistance" / "again", "S1", [{"gene_id": 1}])
    with pytest.raises(SystemExit):  # the same sample twice
        _merging(tmp_path).execute()


def test_execute_missing_table(tmp_path: Path) -> None:
    _sample(tmp_path / "resistance", "S1", [{"gene_id": 1}])
    (tmp_path / "resistance" / "S1" / "S1_resistance.tsv").unlink()
    with pytest.raises(SystemExit):
        _merging(tmp_path).execute()
