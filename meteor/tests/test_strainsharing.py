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
    """12 samples in 6 pairs of related samples (s0-s1, s2-s3, ...).

    msp_a: related pairs at 0, unrelated pairs at 0.01 + 0.0001 * k
    msp_b: 3 pairs only, one below the fixed threshold, one with a low overlap
    """
    samples = [f"s{i}" for i in range(12)]
    rows = []
    for k, (first, second) in enumerate(combinations(samples, 2)):
        related = int(first[1:]) // 2 == int(second[1:]) // 2
        rows.append((first, second, 80.0, 0.0 if related else 0.01 + 0.0001 * k))
    files = {
        "msp_a": _comparison(tmp_path / "msp_a_comparison.tsv", rows),
        "msp_b": _comparison(
            tmp_path / "msp_b_comparison.tsv",
            [("s0", "s1", 50.0, 0.00005), ("s0", "s2", 5.0, 0.0), ("s1", "s2", 50.0, 0.02)],
        ),
        "msp_c": _comparison(tmp_path / "msp_c_comparison.tsv", []),
    }
    groups = tmp_path / "groups.tsv"
    groups.write_text("sample\tgroup\n# comment\n\n" + "".join(f"s{i}\tfam{i // 2}\r\n" for i in range(12)))
    return files, groups


def test_strain_sharing_groups(cohort, tmp_path: Path) -> None:
    files, groups = cohort
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    StrainSharing(output_dir, groups).execute(files)
    events, samples, msp = _read(output_dir)
    msp = msp.set_index("msp_name")
    # msp_a: 6 related + 60 unrelated pairs; threshold = 3rd percentile of the unrelated
    unrelated = [0.01 + 0.0001 * k for k, (a, b) in enumerate(combinations(range(12), 2)) if a // 2 != b // 2]
    assert msp.loc["msp_a", "threshold"] == pytest.approx(np.percentile(unrelated, 3))
    assert msp.loc["msp_a", "threshold_source"] == "unrelated_p3"
    assert msp.loc["msp_a", ["pairs_compared", "related_compared", "unrelated_compared"]].tolist() == [66, 6, 60]
    assert msp.loc["msp_a", "related_shared"] == 6
    assert np.isnan(msp.loc["msp_a", "transmissibility"])  # 6 related pairs < min_pairs
    assert msp.loc["msp_a", "unrelated_shared"] == int(np.sum(np.array(unrelated) <= np.percentile(unrelated, 3)))
    assert msp.loc["msp_a", "unrelated_sharing_rate"] == pytest.approx(msp.loc["msp_a", "unrelated_shared"] / 60)
    # msp_b: too few unrelated pairs for a percentile, low-overlap pair not compared
    assert msp.loc["msp_b", "threshold"] == StrainSharing.FIXED_THRESHOLD
    assert msp.loc["msp_b", "threshold_source"] == "fixed"
    assert msp.loc["msp_b", ["pairs_compared", "pairs_shared"]].tolist() == [2, 1]
    assert np.isnan(msp.loc["msp_b", "sharing_rate"])
    # empty comparison table
    assert msp.loc["msp_c", "pairs_compared"] == 0
    # events = shared pairs, with their relation
    assert len(events) == msp["pairs_shared"].sum()
    assert set(events.loc[events["msp_name"] == "msp_b", "sample2"]) == {"s1"}
    assert set(events.loc[events["relation"] == "related", "sample1"]) == {"s0", "s2", "s4", "s6", "s8", "s10"}
    # sample pairs: one row per compared pair, counts summed over MSPs
    pairs = samples.set_index(["sample1", "sample2"])
    assert len(samples) == 66
    assert pairs.loc[("s0", "s1"), ["relation", "msp_compared", "msp_shared"]].tolist() == ["related", 2, 2]
    assert pairs.loc[("s1", "s2"), ["relation", "msp_compared", "msp_shared"]].tolist() == ["unrelated", 2, 0]
    assert pairs.loc[("s0", "s2"), "msp_compared"] == 1
    assert samples["sharing_rate"].tolist() == pytest.approx((samples["msp_shared"] / samples["msp_compared"]).tolist())
    assert samples["msp_compared"].sum() == msp["pairs_compared"].sum()


def test_strain_sharing_transmissibility(tmp_path: Path) -> None:
    """related pairs >= min_pairs: transmissibility; samples given in either order"""
    rows = [(f"b{i}", f"a{i}", 50.0, 0.0 if i < 8 else 0.5) for i in range(12)]
    rows += [(f"a{i}", f"a{j}", 50.0, 0.1 + i * 0.01 + j * 0.001) for i, j in combinations(range(12), 2)]
    files = {"msp_a": _comparison(tmp_path / "msp_a_comparison.tsv", rows)}
    groups = tmp_path / "groups.tsv"
    groups.write_text("".join(f"a{i}\tfam{i}\nb{i}\tfam{i}\n" for i in range(12)))
    StrainSharing(tmp_path, groups, min_pairs=10).execute(files)
    _, samples, msp = _read(tmp_path)
    assert msp.loc[0, "related_compared"] == 12
    assert msp.loc[0, "transmissibility"] == pytest.approx(8 / 12)
    assert msp.loc[0, "unrelated_shared"] == 0 or msp.loc[0, "unrelated_sharing_rate"] < 0.05
    assert ("a0", "b0") in set(zip(samples["sample1"], samples["sample2"]))
    assert (samples["sample1"] < samples["sample2"]).all()


def test_strain_sharing_without_groups(cohort, tmp_path: Path, caplog) -> None:
    files, _ = cohort
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    StrainSharing(output_dir, None, min_overlap=1.0, min_pairs=2).execute(files)
    assert "No group file" in caplog.text
    _, samples, msp = _read(output_dir)
    msp = msp.set_index("msp_name")
    assert (msp["threshold_source"] == "fixed").all()
    assert (samples["relation"] == "unknown").all()
    assert msp["related_compared"].sum() == 0
    assert msp["transmissibility"].isna().all()
    # min_overlap 1: the 5 % pair of msp_b is compared and shared
    assert msp.loc["msp_b", ["pairs_compared", "pairs_shared"]].tolist() == [3, 2]
    assert msp.loc["msp_b", "sharing_rate"] == pytest.approx(2 / 3)
    assert msp.loc["msp_a", "pairs_shared"] == 6


def test_strain_sharing_reduce(cohort, tmp_path: Path, monkeypatch) -> None:
    """Partial pair counts summed in several steps give the same table"""
    files, groups = cohort
    for name, reduce_pairs in (("once", 10**9), ("often", 1)):
        output_dir = tmp_path / name
        output_dir.mkdir()
        monkeypatch.setattr(StrainSharing, "REDUCE_PAIRS", reduce_pairs)
        StrainSharing(output_dir, groups).execute(files)
    once, often = (
        (tmp_path / name / StrainSharing.SAMPLES_FILE).read_text() for name in ("once", "often")
    )
    assert once == often


def test_strain_sharing_no_comparison(tmp_path: Path) -> None:
    StrainSharing(tmp_path).execute({})
    events, samples, msp = _read(tmp_path)
    assert events.empty and samples.empty and msp.empty
    assert list(samples.columns) == ["sample1", "sample2", "relation", "msp_compared", "msp_shared", "sharing_rate"]


@pytest.mark.parametrize(
    "content",
    ["s1\n", "s1\tfam1\ns1\tfam2\n", "# only a comment\n"],
)
def test_read_groups_errors(tmp_path: Path, content: str) -> None:
    groups = tmp_path / "groups.tsv"
    groups.write_text(content)
    with pytest.raises(SystemExit):
        StrainSharing.read_groups(groups)
