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

import io
import math
import random
import shutil
import subprocess
import sys
import warnings
from contextlib import redirect_stdout
from pathlib import Path

import pandas as pd
import pytest

from ..phylogeny import Phylogeny
from ..session import Component
from ..treeshrink_qc import (
    KShrink,
    TreeShrinkQC,
    apply,
    diameter_ratios,
    loglnorm_threshold,
    per_gene_test,
)

# TreeShrink 1.4.0 (reference implementation, test dependency only)
try:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from treeshrink._vendor.dendropy import Tree  # type: ignore
        from treeshrink.optimal_filter_lib import TreeFilter  # type: ignore
        from treeshrink.threshold_lib import threshold_loglnorm  # type: ignore
    HAS_TREESHRINK = True
except ImportError:
    HAS_TREESHRINK = False
needs_treeshrink = pytest.mark.skipif(
    not HAS_TREESHRINK, reason="TreeShrink 1.4.0 is not installed"
)


def random_tree(
    seed: int, n: int, long: int | None = None, digits: int = 4, zero: bool = False
) -> str:
    """Random binary tree of n samples S0..Sn-1

    `long`: sample on a long branch; `digits`: rounding of the branch lengths
    (few digits = ties); `zero`: some zero branch lengths.
    """
    rng = random.Random(seed)

    def length(low: float, high: float) -> float:
        if zero and rng.random() < 0.3:
            return 0.0
        return round(rng.uniform(low, high), digits)

    nodes = [f"S{i}:{1.0 if i == long else length(0.005, 0.02)}" for i in range(n)]
    while len(nodes) > 1:
        rng.shuffle(nodes)
        left, right = nodes.pop(), nodes.pop()
        nodes.append(f"({left},{right}):{length(0.0, 0.01)}")
    return nodes[0].rsplit(":", 1)[0] + ";"


def distance_matrix(n: int) -> pd.DataFrame:
    names = [f"S{i}" for i in range(n)]
    return pd.DataFrame(
        [[abs(i - j) / n for j in range(n)] for i in range(n)], index=names, columns=names
    )


def test_k_shrink() -> None:
    """Hand-checked star tree: the longest path is D-C, then C-B, then B-A"""
    tree = KShrink("(A:1,B:2,C:3,D:10,E:0.5,F:0.5,G:0.5,H:0.5,I:0.5,J:0.5);")
    diameters, removals = tree.run(2)
    assert diameters == [13.0, 5.0, 3.0]
    assert [sorted(r) for r in removals] == [[], ["D"], ["C", "D"]]
    assert diameter_ratios(diameters, removals) == {"D": 13 / 5, "C": 5 / 3}
    # the search stops when a single sample is left
    assert KShrink("(A:1,B:2,C:3);").run(5) == ([5.0, 3.0], [[], ["C"]])


def test_k_shrink_switches_side() -> None:
    """Best removal set taking both ends of the longest path"""
    # A and B are both far from the rest: remove A then B, not A then C
    tree = KShrink("((A:10,B:9):0.1,(C:1,D:1,E:1,F:1,G:1):0.1);")
    diameters, removals = tree.run(2)
    assert diameters == [19.0, 10.2, 2.0]
    assert [sorted(r) for r in removals] == [[], ["A"], ["A", "B"]]


def test_diameter_ratios_stop_at_zero_diameter() -> None:
    assert diameter_ratios([4.0, 2.0, 0.0, 0.0], [[], ["A"], ["A", "B"], ["A", "B", "C"]]) == {
        "A": 2.0
    }


def test_loglnorm_threshold() -> None:
    assert loglnorm_threshold([3.0]) is None
    assert loglnorm_threshold([1.0, 1.0, 3.0]) is None  # ratios of 1 are ignored
    assert loglnorm_threshold([3.0, 3.0]) == 3.0  # no spread
    # mean 0 and sd 1 of log(log(r)): threshold exp(exp(z(0.95)))
    values = [math.exp(math.exp(x)) for x in (-1.0, 1.0)]
    expected = round(math.exp(math.exp(1.6448536269514722 * math.sqrt(2))), 6)
    assert loglnorm_threshold(values) == pytest.approx(expected, rel=1e-12)
    assert loglnorm_threshold([1.0001, 1e300]) == math.inf


def test_per_gene_test() -> None:
    res = per_gene_test(random_tree(100, 100, long=7))
    assert (res.n, res.k) == (100, 20)  # k = min(100 / 5, 2 * sqrt(100))
    assert len(res.diameters) == len(res.removals) == 21
    assert set(res.flagged) == {"S7"}
    assert res.flagged["S7"] == max(res.ratios.values()) > res.threshold
    assert per_gene_test(random_tree(50, 50)).flagged == {}


def test_per_gene_test_small_or_flat_trees() -> None:
    assert per_gene_test("(A:1,B:2,C:3);").ratios == {}
    # diameter 0 once E is removed: no ratio, nothing flagged
    res = per_gene_test("(A:0,B:0,C:0,D:0,E:1);")
    assert (res.k, res.diameters, res.ratios, res.flagged) == (1, [1.0, 0.0], {}, {})
    # no branch length at all
    res = per_gene_test("(A:0,B:0,C:0,D:0,E:0);")
    assert (res.diameters, res.flagged) == ([], {})
    res = per_gene_test("(A,B,(C,D),E);")
    assert (res.diameters, res.flagged) == ([], {})


def reference(newick: str) -> tuple:
    """Diameters, removal sets, ratios, threshold and flags of TreeShrink 1.4.0"""
    tree = Tree.get(data=newick, schema="newick", preserve_underscores=True)
    k = TreeShrinkQC().k(sum(1 for _ in tree.leaf_node_iter()))
    tree_filter = TreeFilter(ddpTree=tree, centroid_reroot=False, scaling=(5, 2))
    with redirect_stdout(io.StringIO()):
        tree_filter.optFilter(d=k)
    levels = [i for i, e in enumerate(tree_filter.best_entries) if e is not None]
    removals = [sorted(tree_filter.list_removals(d=i)) for i in levels]
    diameters = [tree_filter.min_diams[i] for i in levels]
    ratios: dict[str, float] = {}
    for i in range(1, len(diameters)):
        if diameters[i] == 0:
            break
        for sample in tree_filter.list_removals(d=i):
            ratios[sample] = max(ratios.get(sample, 0.0), diameters[i - 1] / diameters[i])
    threshold, flagged = None, set()
    if len(ratios) > 1:
        try:
            threshold = threshold_loglnorm(list(ratios.values()), e=0.05)
            flagged = {s for s, r in ratios.items() if r > threshold}
        except ValueError:  # fewer than two ratios above 1
            pass
    return diameters, removals, ratios, threshold, flagged


@needs_treeshrink
@pytest.mark.parametrize("seed", range(40))
def test_same_as_treeshrink(seed: int) -> None:
    """Same diameters, removal sets, ratios, threshold and flags as TreeShrink"""
    rng = random.Random(seed)
    n = rng.choice([5, 12, 40, 100, 250])
    newick = random_tree(
        seed,
        n,
        long=rng.choice([None, rng.randrange(n)]),
        digits=rng.choice([1, 2, 4, 8]),  # 1-2 digits: many ties
        zero=rng.random() < 0.3,
    )
    diameters, removals, ratios, threshold, flagged = reference(newick)
    res = per_gene_test(newick)
    assert res.diameters == diameters
    assert [sorted(r) for r in res.removals] == removals
    assert res.ratios == ratios
    if threshold is None or math.isinf(threshold):
        assert res.threshold == threshold
    else:  # same up to the last bits of the mean and sd
        assert res.threshold == pytest.approx(threshold, rel=1e-9)
    assert set(res.flagged) == flagged


@needs_treeshrink
def test_same_as_treeshrink_on_strain_tree(phylogeny_tshrink: Phylogeny) -> None:
    meteor = phylogeny_tshrink.meteor
    msp = phylogeny_tshrink.msp_file_list[0]
    phylogeny_tshrink.process_msp_file(msp, 1, 1, meteor.tree_dir, meteor.tmp_dir)
    newick = (meteor.tree_dir / "msp_0864.tree").read_text()
    diameters, removals, ratios, _, flagged = reference(newick)
    res = per_gene_test(newick)
    assert (res.diameters, [sorted(r) for r in res.removals]) == (diameters, removals)
    assert (res.ratios, set(res.flagged)) == (ratios, flagged)


@pytest.mark.skipif(
    not HAS_TREESHRINK or shutil.which("run_treeshrink.py") is None,
    reason="TreeShrink 1.4.0 command line is not installed",
)
@pytest.mark.parametrize(
    "seed,n,long", [(100, 100, 7), (50, 50, 3), (60, 60, None), (200, 200, 150)]
)
def test_same_flags_as_treeshrink_cli(
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
    assert sorted(per_gene_test(newick).flagged) == sorted(official)
    assert (long is None) == (official == [])


def test_apply(tmp_path: Path) -> None:
    tree_file = tmp_path / "msp_0001.tree"
    tree_file.write_text(random_tree(100, 100, long=7))
    res = apply(tree_file, distance_matrix(100))
    flags = pd.read_csv(tmp_path / "msp_0001_qcflags.tsv", sep="\t")
    assert list(flags.columns) == ["sample", "diameter_ratio", "threshold", "flagged"]
    assert len(flags) == len(res.ratios)
    assert flags.loc[flags["flagged"] == 1, "sample"].tolist() == ["S7"]
    shrink = pd.read_csv(tmp_path / "msp_0001_shrink.tsv", sep="\t", index_col=0)
    assert shrink.shape == (99, 99)
    assert "S7" not in shrink.index and "S7" not in shrink.columns
    # the tree is left untouched
    assert tree_file.read_text() == random_tree(100, 100, long=7)
    # a new run without flags removes the matrix of the previous run
    tree_file.write_text(random_tree(50, 50))
    assert apply(tree_file, distance_matrix(50), TreeShrinkQC()).flagged == {}
    assert not (tmp_path / "msp_0001_shrink.tsv").exists()
    flags = pd.read_csv(tmp_path / "msp_0001_qcflags.tsv", sep="\t")
    assert flags["flagged"].sum() == 0
    # small tree: header only, threshold not available
    tree_file.write_text("(A:1,B:2,C:3);")
    assert apply(tree_file, None).threshold is None
    assert (tmp_path / "msp_0001_qcflags.tsv").read_text() == (
        "sample\tdiameter_ratio\tthreshold\tflagged\n"
    )


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


def test_phylogeny_with_treeshrink(phylogeny_tshrink: Phylogeny) -> None:
    meteor = phylogeny_tshrink.meteor
    tree_dir = meteor.tree_dir
    msp = phylogeny_tshrink.msp_file_list[0]
    # in this process (execute runs it in worker processes)
    assert phylogeny_tshrink.process_msp_file(msp, 1, 1, tree_dir, meteor.tmp_dir)[1]
    flags_file = tree_dir / "msp_0864_qcflags.tsv"
    first = flags_file.read_text()
    assert first.startswith("sample\tdiameter_ratio\tthreshold\tflagged\n")
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
