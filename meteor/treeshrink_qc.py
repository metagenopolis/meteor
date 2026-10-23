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

"""Long-branch (rogue taxon) QC of the strain trees.

Per-gene test of TreeShrink (Mai and Mirarab, BMC Genomics 2018,
doi:10.1186/s12864-018-4620-2), written for meteor: the tree is held in
integer-indexed lists and neither TreeShrink, dendropy nor scipy is needed.
The flags are those of TreeShrink 1.4.0 (``run_treeshrink.py -m per-gene``),
tie-breaking included.

For a tree of n samples:

1. k-shrink. For each j <= k = min(n/5, 2*sqrt(n)), find j samples whose
   removal makes the diameter of the tree (longest sample-to-sample path) as
   small as possible. Removing a sample that is not an end of the longest
   path leaves the diameter unchanged, so each step removes an end of the
   current longest path. A candidate removal set keeps one end and removes
   the other one at each step, switching at most once to the end it kept:
   there are j + 1 candidates with j samples removed (one per switching step,
   plus no switch), and each step only updates the ancestors of the removed
   sample.
2. Diameter ratios. r_j = D(j-1) / D(j), where D(j) is the smallest diameter
   with j samples removed. The ratio of a sample is the largest r_j over the
   levels j whose best removal set contains it.
3. Threshold. log(log(r)) of the ratios above 1 is fitted by a normal
   distribution, and the samples whose ratio is above its 1 - alpha quantile
   are flagged.

The test is additive: ``<msp>.tree`` and ``<msp>.tsv`` are not modified.

Outputs (per MSP, in the tree dir):
  <msp>_qcflags.tsv  sample, diameter_ratio, threshold, flagged (0/1) for the
                     samples of the best removal sets
  <msp>_shrink.tsv   tip-to-tip distance matrix without the flagged samples
                     (only written when at least one sample is flagged)
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from pathlib import Path
from statistics import NormalDist, fmean, stdev
from typing import Callable

import pandas as pd
from cogent3 import make_tree

ALPHA_DEFAULT = 0.05  # false-positive tolerance of the threshold
KSCALING_DEFAULT = (5, 2)  # k = min(n / a, b * sqrt(n))
MIN_TIPS = 4  # smaller trees are not tested

# Record of a node for the live (not removed) samples below it:
# (deepest sample, deepest sample below another child, height, diameter of
# the paths bending at the node). None when no live sample is below the node;
# the second sample and the diameter are None when a single child has some.
Record = tuple[int, int | None, float, float | None] | None


@dataclass(frozen=True)
class TreeShrinkQC:
    """Parameters of the per-gene test"""

    alpha: float = ALPHA_DEFAULT
    kscaling: tuple[int, int] = KSCALING_DEFAULT

    def k(self, n: int) -> int:
        """Largest number of samples removed from a tree of n samples"""
        a, b = self.kscaling
        return min(n // a, int(b * math.sqrt(n)))


@dataclass
class ShrinkResult:
    """Outcome of the per-gene test on one tree"""

    n: int
    k: int
    # smallest diameter after removing 0, 1, ... samples
    diameters: list[float] = field(default_factory=list)
    # best removal set for each number of removed samples (removals[0] = [])
    removals: list[list[str]] = field(default_factory=list)
    ratios: dict[str, float] = field(default_factory=dict)
    threshold: float | None = None
    flagged: dict[str, float] = field(default_factory=dict)


class _Candidate:
    """A candidate removal set: the records it changed and how it was built"""

    __slots__ = ("bend", "records", "removed", "kept", "previous")

    def __init__(
        self,
        bend: int | None,
        records: dict[int, Record],
        removed: int | None = None,
        kept: int | None = None,
        previous: _Candidate | None = None,
    ) -> None:
        self.bend = bend  # node where the longest path bends
        self.records = records  # records that differ from the full tree
        self.removed = removed  # sample removed at the last step
        self.kept = kept  # other end of the path it was removed from
        self.previous = previous

    def removed_samples(self) -> list[int]:
        """Samples removed, last first"""
        samples = []
        candidate: _Candidate | None = self
        while candidate is not None and candidate.removed is not None:
            samples.append(candidate.removed)
            candidate = candidate.previous
        return samples


class KShrink:
    """k-shrink search on a tree stored as integer-indexed lists"""

    def __init__(self, newick: str) -> None:
        tree = make_tree(treestring=newick)
        nodes = list(tree.postorder(include_self=True))
        index = {id(node): i for i, node in enumerate(nodes)}
        self.names: list[str | None] = [
            str(node.name) if node.is_tip() else None for node in nodes
        ]
        self.length = [float(node.length or 0.0) for node in nodes]
        self.parent = [-1] * len(nodes)
        self.children: list[list[int]] = [[] for _ in nodes]
        for i, node in enumerate(nodes):
            for child in node.children:
                self.children[i].append(index[id(child)])
                self.parent[index[id(child)]] = i
        self.n = sum(name is not None for name in self.names)
        # records of the full tree, nodes in postorder (children first)
        self.base: list[Record] = []
        bend, diameter = None, -1.0
        for node in range(len(nodes)):
            record = self._record(node, self.base.__getitem__)
            self.base.append(record)
            # first node in postorder with the largest non-zero diameter
            if record is not None and record[3] and record[3] > diameter:
                bend, diameter = node, record[3]
        self.full = _Candidate(bend, {})

    def _record(self, node: int, lookup: Callable[[int], Record]) -> Record:
        """Record of a node from the records of its children"""
        if self.names[node] is not None:
            return (node, node, 0.0, 0.0)
        best = second = -1.0
        deepest = other = None
        for child in self.children[node]:
            record = lookup(child)
            if record is None:
                continue
            height = record[2] + self.length[child]
            if height > best:
                second, other = best, deepest
                best, deepest = height, record[0]
            elif height > second:
                second, other = height, record[0]
        if deepest is None:
            return None
        return (deepest, other, best, best + second if other is not None else None)

    def _lookup(self, records: dict[int, Record]) -> Callable[[int], Record]:
        base = self.base
        return lambda node: records[node] if node in records else base[node]

    def _bend_record(self, candidate: _Candidate) -> tuple[int, int, float]:
        assert candidate.bend is not None
        record = self._lookup(candidate.records)(candidate.bend)
        assert record is not None and record[1] is not None and record[3] is not None
        return record[0], record[1], record[3]

    def _remove(
        self, candidate: _Candidate, sample: int, records: dict[int, Record]
    ) -> _Candidate:
        """Candidate without `sample`, an end of the longest path

        `records` (the records of `candidate` or a copy) is updated in place.
        """
        lookup = self._lookup(records)
        first, second, _ = self._bend_record(candidate)
        kept = second if sample == first else first
        records[sample] = None
        bend, diameter = None, -math.inf

        def consider(node: int) -> None:
            nonlocal bend, diameter
            record = lookup(node)
            if record is not None and record[3] is not None and record[3] > diameter:
                bend, diameter = node, record[3]

        # only the ancestors of the removed sample change
        node = self.parent[sample]
        while node != -1:
            records[node] = self._record(node, lookup)
            consider(node)
            node = self.parent[node]
        # the kept end is an end of a longest path of the new tree: the path
        # bends on an ancestor of the removed sample (above) or between the
        # kept end and the previous bending node
        node = self.parent[kept]
        while node != candidate.bend:
            consider(node)
            node = self.parent[node]
        consider(node)
        return _Candidate(bend, records, sample, kept, candidate)

    def run(self, k: int) -> tuple[list[float], list[list[str]]]:
        """Smallest diameter and best removal set for 0..k removed samples"""
        diameters: list[float] = []
        best_sets: list[list[str]] = []
        level = [self.full] if self.full.bend is not None else []
        for removed in range(k + 1):
            following: list[_Candidate] = []
            best: tuple[float, _Candidate] | None = None
            for i, candidate in enumerate(level):
                if candidate.bend is None:  # fewer than two live samples
                    level = []
                    break
                first, second, diameter = self._bend_record(candidate)
                if best is None or diameter < best[0]:
                    best = (diameter, candidate)
                if removed == k:
                    continue
                far, kept = (second, first) if first == candidate.kept else (first, second)
                if i == 0:
                    # the first candidate may also switch to the end it kept
                    following.append(self._remove(candidate, far, dict(candidate.records)))
                    following.append(self._remove(candidate, kept, candidate.records))
                else:
                    following.append(self._remove(candidate, far, candidate.records))
            if best is None:
                break
            diameters.append(best[0])
            best_sets.append([str(self.names[s]) for s in best[1].removed_samples()])
            level = following
        return diameters, best_sets


def diameter_ratios(diameters: list[float], removals: list[list[str]]) -> dict[str, float]:
    """Largest diameter ratio D(j-1) / D(j) of each removed sample"""
    ratios: dict[str, float] = {}
    for j in range(1, len(diameters)):
        if diameters[j] == 0:  # only zero branch lengths left
            break
        ratio = diameters[j - 1] / diameters[j]
        for sample in removals[j]:
            ratios[sample] = max(ratios.get(sample, 0.0), ratio)
    return ratios


def loglnorm_threshold(ratios: list[float], alpha: float = ALPHA_DEFAULT) -> float | None:
    """1 - alpha quantile of a log-lognormal fit of the ratios above 1

    Rounded to 6 decimals, as TreeShrink does. None when fewer than two
    ratios are above 1.
    """
    loglog = [math.log(math.log(r)) for r in ratios if 1 < r < math.inf]
    if len(loglog) < 2:
        return None
    mean, sd = fmean(loglog), stdev(loglog)
    # NormalDist needs sd > 0 on Python < 3.12
    quantile = NormalDist(mean, sd).inv_cdf(1 - alpha) if sd > 0 else mean
    try:
        return round(math.exp(math.exp(quantile)), 6)
    except OverflowError:  # widely spread ratios: nothing can be above
        return math.inf


def per_gene_test(newick: str, qc: TreeShrinkQC | None = None) -> ShrinkResult:
    """Per-gene TreeShrink test of one tree"""
    qc = qc or TreeShrinkQC()
    tree = KShrink(newick)
    if tree.n < MIN_TIPS:
        return ShrinkResult(n=tree.n, k=0)
    k = qc.k(tree.n)
    diameters, removals = tree.run(k)
    ratios = diameter_ratios(diameters, removals)
    threshold = loglnorm_threshold(list(ratios.values()), qc.alpha)
    flagged = {}
    if threshold is not None:
        flagged = {sample: r for sample, r in ratios.items() if r > threshold}
    return ShrinkResult(tree.n, k, diameters, removals, ratios, threshold, flagged)


def apply(
    tree_file: Path, dist_df: pd.DataFrame | None, qc: TreeShrinkQC | None = None
) -> ShrinkResult:
    """Test one tree and write the QC files next to it

    :param tree_file: (Path) newick tree of the MSP
    :param dist_df: tip-to-tip distance matrix (index = columns = samples),
                    used to write the matrix without the flagged samples
    """
    res = per_gene_test(tree_file.read_text(), qc)
    stem = tree_file.stem
    tree_dir = tree_file.parent
    threshold = res.threshold if res.threshold is not None else math.nan
    with (tree_dir / f"{stem}_qcflags.tsv").open("w") as out:
        out.write("sample\tdiameter_ratio\tthreshold\tflagged\n")
        for sample in sorted(res.ratios):
            out.write(
                f"{sample}\t{res.ratios[sample]:.6f}\t{threshold:.6f}\t"
                f"{int(sample in res.flagged)}\n"
            )
    shrink_file = tree_dir / f"{stem}_shrink.tsv"
    # a matrix left by an earlier run would no longer match the flags
    shrink_file.unlink(missing_ok=True)
    if res.flagged and dist_df is not None:
        keep = [sample for sample in dist_df.index if sample not in res.flagged]
        if len(keep) >= 2:
            dist_df.loc[keep, keep].to_csv(shrink_file, sep="\t")
    logging.info(
        "TreeShrink QC %s: %d/%d samples flagged (k=%d, threshold=%s)",
        stem,
        len(res.flagged),
        res.n,
        res.k,
        f"{res.threshold:.4f}" if res.threshold is not None else "n/a",
    )
    return res
