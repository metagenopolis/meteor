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

"""Test strain sharing"""

import logging
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ..strainsharing import StrainSharing

COLUMNS = ["sample1", "sample2", "overlap_noN_info_pc", "distance"]


def _comparison(path: Path, rows: list[tuple]) -> Path:
    pd.DataFrame(rows, columns=COLUMNS).to_csv(path, sep="\t", index=False)
    return path


def _read(output_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    return tuple(
        pd.read_csv(output_dir / name, sep="\t")
        for name in (StrainSharing.EVENTS_FILE, StrainSharing.SAMPLES_FILE, StrainSharing.MSP_FILE)
    )


@pytest.fixture
def cohort(tmp_path: Path) -> tuple[dict[str, Path], Path]:
    """12 samples of 6 individuals (s0-s1, s2-s3, ...: two time points each).

    msp_a: pairs of the same individual at 0, the others at 0.01 + 0.0001 * k
    except s0-s2 and s0-s3 at 0 (transmission between two individuals)
    msp_b: 3 pairs only, one below the threshold, one with a low overlap
    """
    samples = [f"s{i}" for i in range(12)]
    rows = []
    for k, (first, second) in enumerate(combinations(samples, 2)):
        same = int(first[1:]) // 2 == int(second[1:]) // 2
        transmission = (first, second) in {("s0", "s2"), ("s0", "s3")}
        rows.append((first, second, 80.0, 0.0 if same or transmission else 0.01 + 0.0001 * k))
    files = {
        "msp_a": _comparison(tmp_path / "msp_a_comparison.tsv", rows),
        "msp_b": _comparison(
            tmp_path / "msp_b_comparison.tsv",
            [("s0", "s1", 50.0, 0.00005), ("s0", "s2", 5.0, 0.0), ("s1", "s2", 50.0, 0.02)],
        ),
        "msp_c": _comparison(tmp_path / "msp_c_comparison.tsv", []),
    }
    individuals = tmp_path / "individuals.tsv"
    individuals.write_text(
        "sample\tindividual\n# comment\n\n" + "".join(f"s{i}\tind{i // 2}\r\n" for i in range(12))
    )
    return files, individuals


def test_strain_sharing_individuals(cohort, tmp_path: Path) -> None:
    files, individuals = cohort
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    StrainSharing(output_dir, individuals).execute(files)
    events, samples, msp = _read(output_dir)
    msp = msp.set_index("msp_name")
    assert (msp["threshold"] == StrainSharing.THRESHOLD).all()
    # msp_a: 60 pairs of different individuals, 2 of them share (s0-s2, s0-s3);
    # the 6 pairs of the same individual share but are no transmission
    assert msp.loc["msp_a", ["n_sample_pairs", "n_shared_pairs"]].tolist() == [60, 2]
    assert msp.loc["msp_a", "transmissibility"] == pytest.approx(2 / 60)
    assert msp.loc["msp_a", ["n_same_individual_pairs", "n_same_individual_shared"]].tolist() == [6, 6]
    # msp_b: low-overlap pair not compared; s0-s1 is the same individual; s1-s2
    # differ: transmissibility computed on 1 pair, however few
    assert msp.loc["msp_b", ["n_sample_pairs", "n_shared_pairs"]].tolist() == [1, 0]
    assert msp.loc["msp_b", "transmissibility"] == 0
    assert msp.loc["msp_b", ["n_same_individual_pairs", "n_same_individual_shared"]].tolist() == [1, 1]
    # empty comparison table
    assert msp.loc["msp_c", "n_sample_pairs"] == 0 and np.isnan(msp.loc["msp_c", "transmissibility"])
    # events = every shared pair, with its relation
    assert len(events) == 6 + 2 + 1
    assert set(events.loc[events["relation"] == "different_individuals", "sample2"]) == {"s2", "s3"}
    assert (events.loc[events["msp_name"] == "msp_b", "relation"] == "same_individual").all()
    # sample pairs: one row per compared pair, counts summed over MSPs
    pairs = samples.set_index(["sample1", "sample2"])
    assert len(samples) == 66
    assert list(samples.columns) == ["sample1", "sample2", "relation", "sharing_rate", "n_msp", "n_msp_shared"]
    assert pairs.loc[("s0", "s1"), ["relation", "n_msp", "n_msp_shared"]].tolist() == ["same_individual", 2, 2]
    assert pairs.loc[("s1", "s2"), ["relation", "n_msp", "n_msp_shared"]].tolist() == ["different_individuals", 2, 0]
    assert pairs.loc[("s0", "s2"), ["n_msp", "n_msp_shared", "sharing_rate"]].tolist() == [1, 1, 1.0]
    assert samples["sharing_rate"].tolist() == pytest.approx((samples["n_msp_shared"] / samples["n_msp"]).tolist())
    assert samples["n_msp"].sum() == (msp["n_sample_pairs"] + msp["n_same_individual_pairs"]).sum()


def test_strain_sharing_mother_infant(tmp_path: Path) -> None:
    """mother and infant are different individuals: their shared strains are
    transmissions; samples given in either order"""
    rows = [(f"infant{i}", f"mother{i}", 50.0, 0.0 if i < 8 else 0.5) for i in range(12)]
    rows += [(f"mother{i}", f"mother{j}", 50.0, 0.1 + i * 0.01 + j * 0.001) for i, j in combinations(range(12), 2)]
    files = {"msp_a": _comparison(tmp_path / "msp_a_comparison.tsv", rows)}
    individuals = tmp_path / "individuals.tsv"
    individuals.write_text("".join(f"mother{i}\tM{i}\ninfant{i}\tI{i}\n" for i in range(12)))
    StrainSharing(tmp_path, individuals).execute(files)
    _, samples, msp = _read(tmp_path)
    assert msp.loc[0, "n_sample_pairs"] == 12 + 66
    assert msp.loc[0, "n_shared_pairs"] == 8
    assert msp.loc[0, "transmissibility"] == pytest.approx(8 / 78)
    assert msp.loc[0, "n_same_individual_pairs"] == 0
    assert ("infant0", "mother0") in set(zip(samples["sample1"], samples["sample2"]))
    assert (samples["sample1"] < samples["sample2"]).all()


def test_strain_sharing_without_individuals(cohort, tmp_path: Path, caplog) -> None:
    caplog.set_level(logging.INFO)
    files, _ = cohort
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    StrainSharing(output_dir, None, min_overlap=1.0).execute(files)
    assert "No individual file" in caplog.text
    _, samples, msp = _read(output_dir)
    msp = msp.set_index("msp_name")
    assert (samples["relation"] == "different_individuals").all()
    assert msp["n_same_individual_pairs"].sum() == 0
    # min_overlap 1: the 5 % pair of msp_b is compared and shared
    assert msp.loc["msp_b", ["n_sample_pairs", "n_shared_pairs"]].tolist() == [3, 2]
    assert msp.loc["msp_b", "transmissibility"] == pytest.approx(2 / 3)
    assert msp.loc["msp_a", ["n_sample_pairs", "n_shared_pairs"]].tolist() == [66, 8]


def test_strain_sharing_reduce(cohort, tmp_path: Path, monkeypatch) -> None:
    """Partial pair counts summed in several steps give the same table"""
    files, individuals = cohort
    for name, reduce_pairs in (("once", 10**9), ("often", 1)):
        output_dir = tmp_path / name
        output_dir.mkdir()
        monkeypatch.setattr(StrainSharing, "REDUCE_PAIRS", reduce_pairs)
        StrainSharing(output_dir, individuals).execute(files)
    once, often = (
        (tmp_path / name / StrainSharing.SAMPLES_FILE).read_text() for name in ("once", "often")
    )
    assert once == often


def test_strain_sharing_no_comparison(tmp_path: Path) -> None:
    StrainSharing(tmp_path).execute({})
    events, samples, msp = _read(tmp_path)
    assert events.empty and samples.empty and msp.empty
    assert list(samples.columns) == ["sample1", "sample2", "relation", "sharing_rate", "n_msp", "n_msp_shared"]


@pytest.mark.parametrize(
    "content",
    ["s1\n", "s1\tind1\ns1\tind2\n", "# only a comment\n"],
)
def test_read_individuals_errors(tmp_path: Path, content: str) -> None:
    individuals = tmp_path / "individuals.tsv"
    individuals.write_text(content)
    with pytest.raises(SystemExit):
        StrainSharing.read_individuals(individuals)
