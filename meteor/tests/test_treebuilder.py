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
import json
import logging
import shutil


@pytest.fixture
def treebuilder_builder(datadir: Path, tmp_path: Path) -> TreeBuilder:
    meteor = Component
    meteor.tmp_path = tmp_path
    meteor.tree_dir = tmp_path / "tree"
    meteor.tree_dir.mkdir()
    meteor.threads = 1
    meteor.strain_dir = datadir / "strain"
    return TreeBuilder(meteor, 0.5, 4, False, 500, 500, None)


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


# def test_get_msp_distance(treebuilder_builder: TreeBuilder, datadir: Path):
#     tree_file = datadir / "tree" / "msp_0864.raxml.bestTree"
#     tree_expected_file = datadir / "expected_output" / "msp_0864_dist.tsv"
#     distance_matrix = treebuilder_builder.get_msp_distance(
#         Tree(str(tree_file.resolve()))
#     )
#     pd.read_csv(
#         tree_expected_file,
#         sep="\t",
#         index_col=0,
#         header=0,
#     ).equals(distance_matrix)


# def test_execute(treebuilder_builder: TreeBuilder, tmp_path: Path):
#     """Test execute"""
#     matrix = treebuilder_builder.meteor.tree_dir / "msp_0864.tsv"
#     treebuilder_builder.execute()
#     assert matrix.exists()


def _treebuilder(strain_dir: Path, tree_dir: Path, tmp_path: Path, sample_file: Path | None) -> TreeBuilder:
    meteor = Component
    meteor.tmp_path = tmp_path
    meteor.tree_dir = tree_dir
    meteor.threads = 1
    meteor.strain_dir = strain_dir
    return TreeBuilder(meteor, 0.5, 4, False, 500, 500, None, sample_file=sample_file)


def _samples(fasta: Path) -> list[str]:
    return sorted(line[1:].strip() for line in fasta.read_text().splitlines() if line.startswith(">"))


def test_execute_sample_list(datadir: Path, tmp_path: Path, caplog):
    """Trees of the samples of the list; a previous tree of other samples is recomputed"""
    caplog.set_level(logging.INFO)
    strain_dir = tmp_path / "strain"
    shutil.copytree(datadir / "strain", strain_dir)
    tree_dir = tmp_path / "tree"
    sample_file = tmp_path / "samples.txt"
    sample_file.write_text("D15B1\nD15B3\n")
    _treebuilder(strain_dir, tree_dir, tmp_path, sample_file).execute()
    assert _samples(tree_dir / "msp_0864.fasta") == ["D15B1", "D15B3"]
    assert "D15B2" not in (tree_dir / "msp_0864.tree").read_text()
    assert "D15B1" in (tree_dir / "msp_0864.tree").read_text()
    # every sample: the concatenated fasta and the tree outputs are recomputed
    (tree_dir / "msp_0864_qcflags.tsv").write_text("stale\n")
    _treebuilder(strain_dir, tree_dir, tmp_path, None).execute()
    assert "msp_0864: samples changed, the tree is recomputed." in caplog.text
    assert _samples(tree_dir / "msp_0864.fasta") == ["D15B1", "D15B2", "D15B3"]
    assert "D15B2" in (tree_dir / "msp_0864.tree").read_text()
    assert not (tree_dir / "msp_0864_qcflags.tsv").exists()
    # same samples: the fasta and the tree are reused
    caplog.clear()
    outputs = [tree_dir / name for name in ("msp_0864.fasta", "msp_0864.tree", "msp_0864.tsv")]
    written = [path.stat().st_mtime_ns for path in outputs]
    _treebuilder(strain_dir, tree_dir, tmp_path, None).execute()
    assert "samples changed" not in caplog.text
    assert [path.stat().st_mtime_ns for path in outputs] == written


def test_execute_sample_list_no_sample(datadir: Path, tmp_path: Path):
    sample_file = tmp_path / "samples.txt"
    sample_file.write_text("unknown\n")
    with pytest.raises(SystemExit):
        _treebuilder(datadir / "strain", tmp_path / "tree", tmp_path, sample_file).execute()


def test_execute_strain_sharing_and_stale_comparisons(datadir: Path, tmp_path: Path, caplog):
    """meteor tree writes the strain sharing tables; a GTR tree directory of
    meteor <= 2.0.22 gets its comparison tables recomputed (issue #133)"""
    caplog.set_level(logging.INFO)
    tree_dir = tmp_path / "tree"
    groups = tmp_path / "groups.tsv"
    groups.write_text("D15B1\tfam1\nD15B2\tfam1\nD15B3\tfam2\n")
    builder = _treebuilder(datadir / "strain", tree_dir, tmp_path, None)
    builder.groups_file = groups
    builder.execute()
    for name in ("strain_sharing_events.tsv", "strain_sharing_samples.tsv", "strain_sharing_msp.tsv"):
        assert (tree_dir / name).exists()
    assert (tree_dir / "strain_sharing_msp.tsv").read_text().splitlines()[1].startswith("msp_0864\t")
    comparison = tree_dir / "msp_0864_comparison.tsv"
    written = comparison.stat().st_mtime_ns
    # same meteor: tables reused
    _treebuilder(datadir / "strain", tree_dir, tmp_path, None).execute()
    assert comparison.stat().st_mtime_ns == written
    # GTR tables of meteor 2.0.22: removed, the MSP is recomputed
    census = tree_dir / "census_stage_4.json"
    config = json.loads(census.read_text())
    config["phylogeny"].update(tree_mode="GTR", phylogeny_meteor_version="2.0.22")
    census.write_text(json.dumps(config))
    _treebuilder(datadir / "strain", tree_dir, tmp_path, None).execute()
    assert "computed in GTR by meteor 2.0.22" in caplog.text
    assert comparison.exists() and comparison.stat().st_mtime_ns != written


@pytest.mark.parametrize(
    "phylogeny",
    [
        {"tree_mode": "TN93", "phylogeny_meteor_version": "2.0.22"},
        {"tree_mode": "GTR", "phylogeny_meteor_version": "2.0.23"},
        {"tree_mode": "GTR"},
        {"tree_mode": "GTR", "phylogeny_meteor_version": "not a version"},
    ],
)
def test_remove_stale_comparisons_kept(treebuilder_builder: TreeBuilder, phylogeny: dict):
    tree_dir = treebuilder_builder.meteor.tree_dir
    (tree_dir / "census_stage_4.json").write_text(json.dumps({"phylogeny": phylogeny}))
    (tree_dir / "msp_1_comparison.tsv").write_text("x\n")
    treebuilder_builder.remove_stale_comparisons()
    assert (tree_dir / "msp_1_comparison.tsv").exists()

