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

"""Test phylogeny"""

from ..session import Component
from ..phylogeny import Phylogeny
import pytest
from pathlib import Path
import pickle


@pytest.fixture
def phylogeny_builder(datadir: Path, tmp_path: Path) -> Phylogeny:
    meteor = Component
    meteor.tmp_dir = tmp_path
    meteor.tree_dir = tmp_path / "tree"
    meteor.tree_dir.mkdir()
    meteor.threads = 1
    meteor.DEFAULT_GAP_CHAR = "?"
    return Phylogeny(meteor, [Path(datadir / "msp_0864.fasta")], 0.5, 4)


def test_compute_site_info(phylogeny_builder: Phylogeny):
    """Test compute_site_info"""
    assert phylogeny_builder.compute_site_info(["?GCT", "A?CT", "A?CT", "A?CT"]) == [
        0.25,
        0.75,
        0.0,
        0.0,
    ]
    assert phylogeny_builder.compute_site_info(["?G", "A?", "A?", "A?"]) == [0.25, 0.75]
    assert phylogeny_builder.compute_site_info(["?CG", "A?T", "A?G", "AC?"]) == [
        0.25,
        0.5,
        0.25,
    ]
    assert phylogeny_builder.compute_site_info(["?G", "?A", "?T", "GC"]) == [0.75, 0]


def test_clean_sites(phylogeny_builder: Phylogeny, datadir: Path, tmpdir: Path):
    msp = tmpdir / "msp_0864_clean.fasta"
    msp_expected_file = datadir / "msp_0864_dict.pck"
    print(msp_expected_file)
    with open(msp_expected_file, "rb") as msp_file:
        msp_expected = pickle.load(msp_file)
        with msp.open("w") as f:
            result_dict, _ = phylogeny_builder.clean_sites(
                phylogeny_builder.msp_file_list[0], f
            )
            for key in msp_expected:
                assert key in result_dict
                assert len(result_dict[key]) == len(msp_expected[key])
                assert result_dict[key] == msp_expected[key]


def test_remove_edge_labels(phylogeny_builder: Phylogeny):
    assert phylogeny_builder.remove_edge_labels("edge.0:1") == ":1"
    assert (
        phylogeny_builder.remove_edge_labels("(A:0.1, B:0.2)edge.1:0.3")
        == "(A:0.1, B:0.2):0.3"
    )


def test_execute(phylogeny_builder: Phylogeny):
    phylogeny_builder.execute()
    result = phylogeny_builder.meteor.tree_dir / "msp_0864.tree"
    assert result.exists()
