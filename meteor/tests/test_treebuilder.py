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
import pandas as pd
import pytest
from pathlib import Path
from hashlib import md5
import json
import logging
import lzma
import random
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


def _treebuilder(
    strain_dir: Path, tree_dir: Path, tmp_path: Path, sample_file: Path | None,
    ref_dir: Path | None = None,
) -> TreeBuilder:
    meteor = Component
    meteor.tmp_path = tmp_path
    meteor.tree_dir = tree_dir
    meteor.threads = 1
    meteor.strain_dir = strain_dir
    meteor.ref_dir = ref_dir
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
    raw = tree_dir / "msp_0864_raw_distances.tsv"
    assert sorted(raw.read_text().splitlines()[0].split("\t")[1:]) == ["D15B1", "D15B3"]
    # every sample: the concatenated fasta and the tree outputs are recomputed
    (tree_dir / "msp_0864_qcflags.tsv").write_text("stale\n")
    _treebuilder(strain_dir, tree_dir, tmp_path, None).execute()
    assert "msp_0864: samples changed, the tree is recomputed." in caplog.text
    assert _samples(tree_dir / "msp_0864.fasta") == ["D15B1", "D15B2", "D15B3"]
    assert "D15B2" in (tree_dir / "msp_0864.tree").read_text()
    assert not (tree_dir / "msp_0864_qcflags.tsv").exists()
    assert sorted(raw.read_text().splitlines()[0].split("\t")[1:]) == ["D15B1", "D15B2", "D15B3"]
    # same samples: the fasta and the tree are reused
    caplog.clear()
    outputs = [tree_dir / name for name in ("msp_0864.fasta", "msp_0864.tree", "msp_0864.tsv")]
    written = [path.stat().st_mtime_ns for path in outputs]
    _treebuilder(strain_dir, tree_dir, tmp_path, None).execute()
    assert "samples changed" not in caplog.text
    assert [path.stat().st_mtime_ns for path in outputs] == written


def test_execute_sample_list_no_sample(datadir: Path, tmp_path: Path, caplog):
    """No sample of the list in the strain directory: a warning, no tree"""
    sample_file = tmp_path / "samples.txt"
    sample_file.write_text("unknown\n")
    builder = _treebuilder(datadir / "strain", tmp_path / "tree", tmp_path, sample_file)
    builder.execute()
    assert f"No sample of {sample_file} was found." in caplog.text
    assert "No tree built" in caplog.text
    assert not list((tmp_path / "tree").iterdir())
    assert not builder.meteor.tmp_dir.exists()


def test_execute_strain_sharing_and_stale_comparisons(datadir: Path, tmp_path: Path, caplog):
    """meteor tree writes the strain sharing tables; a GTR tree directory of
    meteor <= 2.0.22 gets its comparison tables recomputed (issue #133)"""
    caplog.set_level(logging.INFO)
    tree_dir = tmp_path / "tree"
    individuals = tmp_path / "individuals.tsv"
    individuals.write_text("D15B1\tind1\nD15B2\tind1\nD15B3\tind2\n")
    builder = _treebuilder(datadir / "strain", tree_dir, tmp_path, None)
    builder.individuals_file = individuals
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
    builder = _treebuilder(datadir / "strain", tree_dir, tmp_path, None)
    builder.gtr = True
    builder.execute()
    assert "computed in GTR by meteor 2.0.22" in caplog.text
    assert comparison.exists() and comparison.stat().st_mtime_ns != written


def test_execute_model_change(tmp_path: Path, caplog):
    """A tree directory of the other model (GTR, the former default, run again in
    TN93, or the reverse) gets its trees and distances recomputed"""
    strain_dir = _strain_dir(tmp_path / "strain")
    tree_dir = tmp_path / "tree"
    outputs = [tree_dir / name for name in ("msp_0001.tree", "msp_0001.tsv",
               "msp_0001_raw_distances.tsv", "msp_0001_comparison.tsv")]

    def run(gtr: bool) -> list[int]:
        builder = _treebuilder(strain_dir, tree_dir, tmp_path, None)
        builder.gtr = gtr
        builder.execute()
        return [path.stat().st_mtime_ns for path in outputs]

    def mode() -> str:
        return json.loads((tree_dir / "census_stage_4.json").read_text())["phylogeny"]["tree_mode"]

    written = run(True)
    assert mode() == "GTR"
    gtr_distances = (tree_dir / "msp_0001_raw_distances.tsv").read_text()
    # same model: reused
    assert run(True) == written
    # TN93: recomputed
    caplog.clear()
    tn93 = run(False)
    assert "were computed in GTR, not in TN93" in caplog.text
    assert all(after != before for after, before in zip(tn93, written))
    assert mode() == "TN93"
    assert (tree_dir / "msp_0001_raw_distances.tsv").read_text() != gtr_distances
    # TN93 again: reused, no warning
    caplog.clear()
    assert run(False) == tn93
    assert "not in" not in caplog.text


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



def _catalogue(ref_dir: Path, thresholds: dict[str, float]) -> Path:
    """Minimal catalogue with a strain_thresholds annotation"""
    database = ref_dir / "database"
    database.mkdir(parents=True)
    pd.DataFrame(
        {"msp_name": list(thresholds), "threshold": list(thresholds.values()),
         "exact": [True] * len(thresholds)}
    ).to_feather(database / "cat_strain_thresholds.feather")
    (ref_dir / "cat_reference.json").write_text(
        json.dumps(
            {
                "reference_file": {"database_dir": "database"},
                "annotation": {"strain_thresholds": {"filename": "cat_strain_thresholds.feather"}},
            }
        )
    )
    return ref_dir


def _strain_dir(path: Path) -> Path:
    """Strain directory of 4 samples, one MSP (msp_0001) of 2000 bases: s1, s2
    and s3 differ from s0 at 1, 4 and 40 sites"""
    rng = random.Random(1)
    base = [rng.choice("ACGT") for _ in range(2000)]
    for sample, n_diff in {"s0": 0, "s1": 1, "s2": 4, "s3": 40}.items():
        seq = list(base)
        for position in range(0, n_diff * 50, 50):
            seq[position] = "A" if seq[position] != "A" else "C"
        (path / sample).mkdir(parents=True)
        (path / sample / f"{sample}_census_stage_3.json").write_text(
            json.dumps({"sample_info": {"sample_name": sample}})
        )
        with lzma.open(path / sample / "msp_0001.fasta.xz", "wt") as fasta:
            fasta.write(f">{sample}\n{''.join(seq)}\n")
    return path


def test_execute_catalogue_strain_thresholds(tmp_path: Path, caplog):
    """meteor tree -r: the MSPs of the strain_thresholds annotation of the catalogue
    use their own threshold, in the comparison tables (also when they are reused
    from a run without -r) and in the strain sharing"""
    strain_dir = _strain_dir(tmp_path / "strain")
    tree_dir = tmp_path / "tree"
    comparison = tree_dir / "msp_0001_comparison.tsv"

    def same_strain() -> set[tuple[str, str]]:
        table = pd.read_csv(comparison, sep="\t")
        return {
            tuple(sorted(pair))
            for pair in table.loc[table["distance_category"] == "same_strain", ["sample1", "sample2"]]
            .itertuples(index=False)
        }

    def shared() -> tuple[float, int]:
        msp = pd.read_csv(tree_dir / "strain_sharing_msp.tsv", sep="\t").set_index("msp_name")
        return msp.loc["msp_0001", "threshold"], msp.loc["msp_0001", "n_shared_pairs"]

    # default threshold (5.9e-4): s0-s1 only (1 site / 2000)
    _treebuilder(strain_dir, tree_dir, tmp_path, None).execute()
    assert same_strain() == {("s0", "s1")}
    assert shared() == (pytest.approx(0.000594809910301849), 1)
    distances = pd.read_csv(comparison, sep="\t")["distance"].tolist()
    # threshold of the catalogue (3e-3): s0-s2, s1-s2 too; the tree and the
    # distances are reused, the categories updated
    ref_dir = _catalogue(tmp_path / "catalogue", {"msp_0001": 0.003, "msp_other": 1e-9})
    caplog.clear()
    _treebuilder(strain_dir, tree_dir, tmp_path, None, ref_dir).execute()
    assert same_strain() == {("s0", "s1"), ("s0", "s2"), ("s1", "s2")}
    assert shared() == (0.003, 3)
    # TN93 (the default, as _treebuilder): the distances of the thresholds, no warning
    assert "were computed on TN93 distances" not in caplog.text
    assert pd.read_csv(comparison, sep="\t")["distance"].tolist() == distances
    census = json.loads((tree_dir / "census_stage_4.json").read_text())
    assert census["phylogeny"]["per_msp_strain_thresholds"] == 2
    # without -r again: back to the default threshold
    _treebuilder(strain_dir, tree_dir, tmp_path, None).execute()
    assert same_strain() == {("s0", "s1")}
    # a new tree directory with -r: the threshold of the catalogue from the start;
    # GTR distances with thresholds computed on TN93: a warning
    tree_dir = tmp_path / "tree_ref"
    comparison = tree_dir / "msp_0001_comparison.tsv"
    caplog.clear()
    builder = _treebuilder(strain_dir, tree_dir, tmp_path, None, ref_dir)
    builder.gtr = True
    builder.execute()
    assert same_strain() == {("s0", "s1"), ("s0", "s2"), ("s1", "s2")}
    assert "were computed on TN93 distances" in caplog.text
