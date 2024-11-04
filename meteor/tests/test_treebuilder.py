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

"""Test treebuilder analysis"""

from ..session import Component
from ..treebuilder import TreeBuilder
import pytest
from pathlib import Path
from hashlib import md5
from ete3 import Tree
import pandas as pd


@pytest.fixture
def treebuilder_builder(datadir: Path, tmp_path: Path) -> TreeBuilder:
    meteor = Component
    meteor.tmp_path = tmp_path
    meteor.tree_dir = tmp_path / "tree"
    meteor.tree_dir.mkdir()
    meteor.threads = 1
    meteor.strain_dir = datadir / "strain"
    return TreeBuilder(meteor, 0.5, 4, 500, 500, None)


def test_concatenate(treebuilder_builder: TreeBuilder, datadir: Path):
    msp_file_dict = {}
    msp_file_dict["D15B1"] = [
        Path(datadir / "strain" / "D15B1" / "msp_0864.fasta.xz"),
        Path(datadir / "strain" / "D15B2" / "msp_0864.fasta.xz"),
        Path(datadir / "strain" / "D15B3" / "msp_0864.fasta.xz"),
    ]
    msp_list = treebuilder_builder.concatenate(msp_file_dict)
    with msp_list[0].open("rb") as msp:
        assert md5(msp.read()).hexdigest() == "719ba21ae065b2700a56e307d578c12e"


def test_get_msp_distance(treebuilder_builder: TreeBuilder, datadir: Path):
    tree_file = datadir / "tree" / "msp_0864.raxml.bestTree"
    tree_expected_file = datadir / "expected_output" / "msp_0864_dist.tsv"
    distance_matrix = treebuilder_builder.get_msp_distance(
        Tree(str(tree_file.resolve()))
    )
    pd.read_csv(
        tree_expected_file,
        sep="\t",
        index_col=0,
        header=0,
    ).equals(distance_matrix)


def test_execute(treebuilder_builder: TreeBuilder, tmp_path: Path):
    """Test execute"""
    matrix = treebuilder_builder.meteor.tree_dir / "msp_0864.tsv"
    treebuilder_builder.execute()
    assert matrix.exists()
