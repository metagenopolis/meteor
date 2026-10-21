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

"""Test the TreeShrink QC of strain trees"""

import random
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from .. import treeshrink_qc
from ..phylogeny import Phylogeny
from ..session import Component
from ..treebuilder import TreeBuilder
from ..treeshrink_qc import TreeShrinkQC, apply, check_available, per_gene_flag

try:
    treeshrink_qc._load_reference()
    HAS_TREESHRINK = True
except ImportError:
    HAS_TREESHRINK = False
needs_treeshrink = pytest.mark.skipif(
    not HAS_TREESHRINK, reason="TreeShrink >= 1.4.0 is not installed"
)


def random_tree(seed: int, n: int, long: int | None = None) -> str:
    """Random binary tree of n samples S0..Sn-1, sample `long` on a long branch"""
    rng = random.Random(seed)
    nodes = [
        f"S{i}:{1.0 if i == long else rng.uniform(0.005, 0.02):.4f}" for i in range(n)
    ]
    while len(nodes) > 1:
        rng.shuffle(nodes)
        left, right = nodes.pop(), nodes.pop()
        nodes.append(f"({left},{right}):{rng.uniform(0.001, 0.01):.4f}")
    return nodes[0].rsplit(":", 1)[0] + ";"


def distance_matrix(n: int) -> pd.DataFrame:
    names = [f"S{i}" for i in range(n)]
    return pd.DataFrame(
        [[abs(i - j) / n for j in range(n)] for i in range(n)], index=names, columns=names
    )


@needs_treeshrink
def test_per_gene_flag() -> None:
    res = per_gene_flag(random_tree(100, 100, long=7))
    assert res["n"] == 100
    assert res["k"] == 20  # min(100 / 5, 2 * sqrt(100))
    assert set(res["flagged"]) == {"S7"}
    assert res["flagged"]["S7"] == max(res["features"].values())
    assert res["flagged"]["S7"] > res["threshold"]
    # no long branch: nothing flagged
    assert per_gene_flag(random_tree(50, 50))["flagged"] == {}


@needs_treeshrink
def test_per_gene_flag_small_tree() -> None:
    assert per_gene_flag("(A:1,B:2,C:3);") == {
        "features": {}, "threshold": None, "flagged": {}, "k": 0, "n": 3
    }
    # diameter 0 once E is removed: no ratio, nothing flagged
    res = per_gene_flag("(A:0,B:0,C:0,D:0,E:1);")
    assert (res["n"], res["k"], res["features"], res["flagged"]) == (5, 1, {}, {})


@needs_treeshrink
def test_per_gene_flag_degenerate_threshold(monkeypatch) -> None:
    """The reference threshold fails on degenerate features: nothing flagged"""
    def fail(*_, **__):
        raise ValueError("degenerate")
    real = treeshrink_qc._load_reference
    monkeypatch.setattr(
        treeshrink_qc, "_load_reference", lambda: (*real()[:2], fail)
    )
    res = per_gene_flag(random_tree(100, 100, long=7))
    assert res["threshold"] is None and res["flagged"] == {}
    assert res["features"]


@pytest.mark.skipif(
    not HAS_TREESHRINK or shutil.which("run_treeshrink.py") is None,
    reason="TreeShrink >= 1.4.0 command line is not installed",
)
@pytest.mark.parametrize(
    "seed,n,long", [(100, 100, 7), (50, 50, 3), (60, 60, None), (200, 200, 150)]
)
def test_per_gene_flag_matches_treeshrink_cli(
    tmp_path: Path, seed: int, n: int, long: int | None
) -> None:
    """Same flags as run_treeshrink.py -m per-gene"""
    newick = random_tree(seed, n, long)
    (tmp_path / "input.tree").write_text(newick + "\n")
    subprocess.run(
        [sys.executable, "-W", "ignore", str(shutil.which("run_treeshrink.py")),
         "-t", str(tmp_path / "input.tree"), "-m", "per-gene",
         "-o", str(tmp_path / "out"), "-O", "out", "-f"],
        check=True, capture_output=True,
    )
    official = (tmp_path / "out" / "out.txt").read_text().split()
    assert sorted(per_gene_flag(newick)["flagged"]) == sorted(official)
    assert (long is None) == (official == [])


@needs_treeshrink
def test_apply(tmp_path: Path) -> None:
    tree_file = tmp_path / "msp_0001.tree"
    tree_file.write_text(random_tree(100, 100, long=7))
    res = apply(tree_file, distance_matrix(100))
    flags = pd.read_csv(tmp_path / "msp_0001_qcflags.tsv", sep="\t")
    assert list(flags.columns) == ["sample", "feature", "threshold", "flagged"]
    assert len(flags) == len(res["features"])
    assert flags.loc[flags["flagged"] == 1, "sample"].tolist() == ["S7"]
    shrink = pd.read_csv(tmp_path / "msp_0001_shrink.tsv", sep="\t", index_col=0)
    assert shrink.shape == (99, 99)
    assert "S7" not in shrink.index and "S7" not in shrink.columns
    # the tree is left untouched
    assert tree_file.read_text() == random_tree(100, 100, long=7)
    # a new run without flags removes the matrix of the previous run
    tree_file.write_text(random_tree(50, 50))
    assert apply(tree_file, distance_matrix(50), TreeShrinkQC())["flagged"] == {}
    assert not (tmp_path / "msp_0001_shrink.tsv").exists()
    flags = pd.read_csv(tmp_path / "msp_0001_qcflags.tsv", sep="\t")
    assert flags["flagged"].sum() == 0


def test_check_available_missing(monkeypatch) -> None:
    def missing():
        raise ImportError("No module named 'treeshrink'")
    monkeypatch.setattr(treeshrink_qc, "_load_reference", missing)
    with pytest.raises(SystemExit) as error:
        check_available()
    assert error.value.code == 1


@needs_treeshrink
def test_check_available() -> None:
    check_available()


@pytest.fixture
def phylogeny_tshrink(tmp_path: Path) -> Phylogeny:
    meteor = Component
    meteor.tmp_dir = tmp_path
    meteor.tree_dir = tmp_path / "tree"
    meteor.tree_dir.mkdir()
    meteor.threads = 1
    meteor.DEFAULT_GAP_CHAR = "?"
    msp = Path(__file__).parent / "test_phylogeny" / "msp_0864.fasta"
    return Phylogeny(meteor, [msp], 0.5, 4, False, tshrink=True)


@needs_treeshrink
def test_phylogeny_with_treeshrink(phylogeny_tshrink: Phylogeny) -> None:
    meteor = phylogeny_tshrink.meteor
    tree_dir = meteor.tree_dir
    msp = phylogeny_tshrink.msp_file_list[0]
    # in this process (execute runs it in worker processes)
    assert phylogeny_tshrink.process_msp_file(msp, 1, 1, tree_dir, meteor.tmp_dir)[1]
    flags_file = tree_dir / "msp_0864_qcflags.tsv"
    first = flags_file.read_text()
    assert first.startswith("sample\tfeature\tthreshold\tflagged\n")
    tree = (tree_dir / "msp_0864.tree").read_text()
    matrix = (tree_dir / "msp_0864.tsv").read_text()
    # existing tree: not rebuilt, QC run from the tree and its distance matrix
    flags_file.unlink()
    assert phylogeny_tshrink.process_msp_file(msp, 1, 1, tree_dir, meteor.tmp_dir)[1]
    assert flags_file.read_text() == first
    flags_file.unlink()
    phylogeny_tshrink.execute()
    assert flags_file.read_text() == first
    assert (tree_dir / "msp_0864.tree").read_text() == tree
    assert (tree_dir / "msp_0864.tsv").read_text() == matrix


def test_treeshrink_qc_failure_is_not_fatal(
    phylogeny_tshrink: Phylogeny, tmp_path: Path, monkeypatch, caplog
) -> None:
    def fail(*_):
        raise RuntimeError("boom")
    monkeypatch.setattr("meteor.phylogeny.treeshrink_apply", fail)
    tree_file = tmp_path / "msp_0001.tree"
    tree_file.write_text("(A:1,B:2,C:3,D:4);")
    distance_matrix(4).to_csv(tree_file.with_suffix(".tsv"), sep="\t")
    phylogeny_tshrink.treeshrink_qc(tree_file)
    assert "TreeShrink QC failed for msp_0001: boom" in caplog.text


def test_treebuilder_requires_treeshrink(tmp_path: Path, monkeypatch) -> None:
    def missing():
        raise ImportError("No module named 'treeshrink'")
    monkeypatch.setattr(treeshrink_qc, "_load_reference", missing)
    meteor = Component
    meteor.tmp_path = tmp_path
    meteor.tree_dir = tmp_path / "tree"
    with pytest.raises(SystemExit):
        TreeBuilder(meteor, 0.5, 4, False, 500, 500, None, tshrink=True)
    assert not meteor.tree_dir.exists()
