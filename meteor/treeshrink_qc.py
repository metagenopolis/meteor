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

"""Per-MSP long-branch (rogue taxon) QC for meteor strain trees.

Applies the TreeShrink algorithm (Mai & Mirarab, BMC Genomics 2018,
doi:10.1186/s12864-018-4620-2; reference implementation
https://github.com/uym2/TreeShrink, GPL-3.0) to each per-MSP strain tree
built by :mod:`meteor.phylogeny`.

For each MSP tree the *per-gene* outlier test is run:

  1. k-shrink: find up to ``k = min(n/5, 2*sqrt(n))`` leaves whose removal
     maximally reduces the tree diameter (``TreeFilter``).
  2. feature: for removal step ``i`` the ratio ``r_i = diam(i-1)/diam(i)``;
     each leaf's feature is the max ``r_i`` over the steps where it is in the
     optimal removal set.
  3. test: fit a log-lognormal on the features and flag the leaves above the
     ``(1 - alpha)`` quantile (``threshold_loglnorm``).

The step is additive: ``<msp>.tree`` and ``<msp>.tsv`` are left untouched.

Outputs (per MSP, in the tree dir):
  <msp>_qcflags.tsv  sample, feature (max diameter-reduction ratio),
                     threshold, flagged (0/1)
  <msp>_shrink.tsv   tip-to-tip distance matrix without the flagged samples
                     (only written when at least one sample is flagged)

Requires TreeShrink >= 1.4.0 (``conda install -c bioconda treeshrink``, or
``pip install git+https://github.com/uym2/TreeShrink.git@v1.4.0``); it is not
on PyPI.
"""

from __future__ import annotations

import io
import logging
import math
import sys
import warnings
from contextlib import redirect_stdout
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

# TreeShrink >= 1.3 defaults
KSCALING_DEFAULT = (5, 2)  # k = min(n/a, b*sqrt(n))
ALPHA_DEFAULT = 0.05  # false-positive tolerance
MIN_TIPS = 4  # trees with fewer tips are too small to QC


@dataclass
class TreeShrinkQC:
    """Config for the per-MSP TreeShrink QC step."""

    alpha: float = ALPHA_DEFAULT
    kscaling: tuple[int, int] = KSCALING_DEFAULT
    centroid: bool = False  # reroot at unit-weight centroid (deep trees)


def _load_reference():
    """Import the reference TreeShrink components (TreeShrink >= 1.4.0)."""
    with warnings.catch_warnings():
        # the dendropy copy vendored by TreeShrink has invalid escape sequences
        warnings.simplefilter("ignore", SyntaxWarning)
        from treeshrink._vendor.dendropy import Tree  # type: ignore
        from treeshrink.optimal_filter_lib import TreeFilter  # type: ignore
        from treeshrink.threshold_lib import threshold_loglnorm  # type: ignore
    return Tree, TreeFilter, threshold_loglnorm


def check_available() -> None:
    """Stop with an error when TreeShrink >= 1.4.0 cannot be imported."""
    try:
        _load_reference()
    except ImportError as error:
        logging.error(
            "--treeshrink requires TreeShrink >= 1.4.0 (conda install -c bioconda "
            "treeshrink): %s",
            error,
        )
        sys.exit(1)


def per_gene_flag(newick: str, qc: TreeShrinkQC | None = None) -> dict:
    """Run the per-gene TreeShrink test on a single strain tree.

    Returns dict with keys:
      features : {sample: feature}   feature = max diameter-reduction ratio
      threshold: float | None        log-lognormal (1-alpha) threshold
      flagged  : {sample: feature}   samples above threshold (may be empty)
      k        : int                 removals solved
      n        : int                 number of tips
    """
    qc = qc or TreeShrinkQC()
    Tree, TreeFilter, threshold_loglnorm = _load_reference()
    a_tree = Tree.get(data=newick, schema="newick", preserve_underscores=True)
    n = sum(1 for _ in a_tree.leaf_node_iter())
    if n < MIN_TIPS:
        return {"features": {}, "threshold": None, "flagged": {}, "k": 0, "n": n}

    a, b = qc.kscaling
    k = min(n // a, int(b * math.sqrt(n)))
    tree_filter = TreeFilter(
        ddpTree=a_tree, centroid_reroot=qc.centroid, scaling=qc.kscaling
    )
    with redirect_stdout(io.StringIO()):  # "Solving k-shrink with k = ..."
        tree_filter.optFilter(d=k)

    # per-gene feature of a leaf: max ratio over the optimal removal steps
    features: dict[str, float] = {}
    for i in range(1, len(tree_filter.min_diams)):
        if tree_filter.min_diams[i] == 0:
            break
        ratio = tree_filter.min_diams[i - 1] / tree_filter.min_diams[i]
        for sample in tree_filter.list_removals(d=i):
            features[sample] = max(features.get(sample, 0.0), ratio)

    threshold = None
    flagged: dict[str, float] = {}
    if len(features) > 1:
        try:
            threshold = float(threshold_loglnorm(list(features.values()), e=qc.alpha))
            flagged = {s: r for s, r in features.items() if r > threshold}
        except Exception:  # pylint: disable=broad-except
            # degenerate distribution (zero spread...): the reference raises,
            # TreeShrink then removes nothing
            threshold, flagged = None, {}
    return {
        "features": features,
        "threshold": threshold,
        "flagged": flagged,
        "k": k,
        "n": n,
    }


def apply(
    tree_file: Path, dist_df: pd.DataFrame | None, qc: TreeShrinkQC | None = None
) -> dict:
    """QC one tree (additive outputs) and return the result of per_gene_flag.

    ``dist_df`` is the tip-to-tip distance matrix (index = columns = sample
    names) used to write the shrunk matrix when samples are flagged.
    """
    qc = qc or TreeShrinkQC()
    stem = tree_file.stem
    tree_dir = tree_file.parent
    res = per_gene_flag(tree_file.read_text(), qc)

    with (tree_dir / f"{stem}_qcflags.tsv").open("w") as out:
        out.write("sample\tfeature\tthreshold\tflagged\n")
        thr = res["threshold"] if res["threshold"] is not None else float("nan")
        for sample in sorted(res["features"]):
            out.write(
                f"{sample}\t{res['features'][sample]:.6f}\t{thr:.6f}\t"
                f"{int(sample in res['flagged'])}\n"
            )

    shrink_file = tree_dir / f"{stem}_shrink.tsv"
    # a matrix left by an earlier run would no longer match the flags
    shrink_file.unlink(missing_ok=True)
    if res["flagged"] and dist_df is not None:
        keep = [s for s in dist_df.index if s not in res["flagged"]]
        if len(keep) >= 2:
            dist_df.loc[keep, keep].to_csv(shrink_file, sep="\t")

    logging.info(
        "TreeShrink QC %s: %d/%d samples flagged (k=%d, threshold=%s)",
        stem,
        len(res["flagged"]),
        res["n"],
        res["k"],
        f"{res['threshold']:.4f}" if res["threshold"] is not None else "n/a",
    )
    return res
