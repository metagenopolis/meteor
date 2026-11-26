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

"""Test meteor resistance_summary"""

import json
import lzma
import random
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ..resistancesummary import ResistanceSummary, residue_classes_map
from ..session import Component

PARAMETERS = {
    "annotations": ["resfinder", "mustard", "kegg"],
    "min_depth": 10,
    "min_alt_count": 3,
    "min_frequency": 0.1,
    "min_gene_coverage": 0.8,
    "base_quality": 20,
    "core_size": 100,
}


def _sample(
    root: Path,
    name: str,
    rows: list[dict],
    reference: str = "cat",
    consensus: dict[int, str] | None = None,
    positions: list[dict] | None = None,
    **changed,
) -> None:
    """Output of meteor resistance for one sample (consensus: {gene_id: sequence},
    positions: rows of its table of the known resistance positions)"""
    sample_dir = root / name
    sample_dir.mkdir(parents=True)
    columns = {
        "msp_name": "M1",
        "resfinder": "tet(W)",
        "resfinderfg": pd.NA,
        "mustard": pd.NA,
        "mustard_prediction": pd.NA,
        "target": pd.NA,
        "gene_length": 900,
        "alignments": 100,
        "interpretable": True,
        "major_per_kb": 1.0,
        "minor_per_kb": 0.0,
        "control_genes": 2,
        "control_major_per_kb": 2.0,
        "control_minor_per_kb": 0.5,
        "reads_resistant": np.nan,
        "reads_wt": np.nan,
        "reads_other": np.nan,
        "reads_linked": np.nan,
    }
    table = pd.DataFrame([{**columns, **row} for row in rows])
    table = table[["gene_id", *columns]]
    table.to_csv(sample_dir / f"{name}_resistance.tsv", sep="\t", index=False, na_rep="NA")
    (sample_dir / f"{name}_census_stage_resistance.json").write_text(
        json.dumps(
            {
                "meteor_version": "2.0.23",
                "sample_info": {"sample_name": name},
                "mapping": {"reference_name": reference},
                "resistance": {
                    **PARAMETERS,
                    **changed,
                    "resistance_file": f"{name}_resistance.tsv",
                    **({"consensus_file": f"{name}_resistance_consensus.fasta.xz"}
                       if consensus is not None else {}),
                    **({"positions_file": f"{name}_resistance_positions.tsv"}
                       if positions is not None else {}),
                },
            }
        )
    )
    if positions is not None:
        pd.DataFrame([{**POSITION, **row} for row in positions]).to_csv(
            sample_dir / f"{name}_resistance_positions.tsv", sep="\t", index=False, na_rep="NA"
        )
    if consensus is not None:
        with lzma.open(sample_dir / f"{name}_resistance_consensus.fasta.xz", "wt") as fasta:
            for gene, sequence in consensus.items():
                fasta.write(f">{gene}\n{sequence}\n")


# a read gyrA 83 position with its catalogue residue
POSITION = {
    "gene_id": 10, "msp_name": "M1", "family": "gyrA", "codon_pos": 83, "ref_residue": "S",
    "resistance_residues": "L,W", "gene_aa_pos": 80, "gene_nt_pos": 238, "gene_residue": "S",
    "window_identity": 0.9, "depth": 20, "readable": True, "major_residue": "S",
    "major_frequency": 1.0, "major_is_resistance": False, "resistance_reads": 0,
    "resistance_frequency": 0.0, "known_resistance": False, "residue_counts": "S:20",
}


def _merging(tmp_path: Path, sample_file: Path | None = None, **options) -> ResistanceSummary:
    meteor = Component
    meteor.resistance_dir = tmp_path / "resistance"
    meteor.merging_dir = tmp_path / "merged"
    meteor.tmp_path = tmp_path
    meteor.threads = 1
    return ResistanceSummary(meteor, "out", sample_file, **options)


@pytest.fixture
def three_samples(tmp_path: Path) -> Path:
    root = tmp_path / "resistance"
    # gene 1: interpretable in S1 and S2, more variants than its controls in S1 only
    # gene 2 (gyrA, no MSP): reads in S1 (not interpretable) and S3
    _sample(root, "S1", [
        {"gene_id": 1, "major_per_kb": 4.0, "minor_per_kb": 1.0},
        {"gene_id": 2, "msp_name": pd.NA, "resfinder": pd.NA, "target": "gyrA",
         "interpretable": False, "control_genes": 0,
         "control_major_per_kb": np.nan, "control_minor_per_kb": np.nan},
    ])
    _sample(root, "S2", [{"gene_id": 1, "major_per_kb": 1.0, "minor_per_kb": 0.0}])
    _sample(root, "S3", [
        {"gene_id": 2, "msp_name": pd.NA, "resfinder": pd.NA, "target": "gyrA",
         "major_per_kb": 6.0, "control_genes": 0,
         "control_major_per_kb": np.nan, "control_minor_per_kb": np.nan},
    ])
    return root


def test_execute(three_samples: Path, tmp_path: Path) -> None:
    _merging(tmp_path).execute()
    out = tmp_path / "merged"
    long = pd.read_table(out / "out_resistance.tsv")
    assert list(long.columns[:2]) == ["sample", "gene_id"]
    assert list(zip(long["sample"], long["gene_id"])) == [("S1", 1), ("S1", 2), ("S2", 1), ("S3", 2)]
    major = pd.read_table(out / "out_resistance_major_per_kb.tsv", index_col=0)
    assert list(major.columns) == ["S1", "S2", "S3"]
    assert major.loc[1, "S1"] == 4.0 and major.loc[1, "S2"] == 1.0
    assert pd.isna(major.loc[1, "S3"])  # no reads
    assert pd.isna(major.loc[2, "S1"]) and major.loc[2, "S3"] == 6.0  # not interpretable in S1
    interpretable = pd.read_table(out / "out_resistance_interpretable.tsv", index_col=0)
    assert interpretable.to_dict("index") == {1: {"S1": 1, "S2": 1, "S3": 0}, 2: {"S1": 0, "S2": 0, "S3": 1}}
    genes = pd.read_table(out / "out_resistance_genes.tsv").set_index("gene_id")
    assert genes.loc[1, "resfinder"] == "tet(W)" and genes.loc[2, "target"] == "gyrA"
    assert genes.loc[1, "samples_with_reads"] == 2 and genes.loc[1, "samples_interpretable"] == 2
    assert genes.loc[2, "samples_with_reads"] == 2 and genes.loc[2, "samples_interpretable"] == 1
    assert genes.loc[1, "interpretable_fraction"] == pytest.approx(2 / 3)
    assert genes.loc[1, "median_major_per_kb"] == pytest.approx(2.5)
    assert genes.loc[1, "median_control_major_per_kb"] == pytest.approx(2.0)
    # S1: 4 > 2, S2: 1 < 2
    assert genes.loc[1, "samples_with_controls"] == 2
    assert genes.loc[1, "major_above_controls_fraction"] == pytest.approx(0.5)
    assert genes.loc[1, "minor_above_controls_fraction"] == pytest.approx(0.5)
    # gene 2 has no controls (no MSP)
    assert genes.loc[2, "samples_with_controls"] == 0
    assert pd.isna(genes.loc[2, "major_above_controls_fraction"])
    assert genes.loc[2, "median_major_per_kb"] == pytest.approx(6.0)
    census = json.loads((out / "out_resistance_census.json").read_text())["resistance_summary"]
    assert census["samples"] == 3 and census["genes"] == 2 and census["interpretable_genes"] == 2
    assert census["reference_name"] == "cat" and census["min_depth"] == 10
    assert census["trees"] is None
    assert not (out / "out_resistance_trees.tsv").exists()


def test_execute_sample_list(three_samples: Path, tmp_path: Path) -> None:
    sample_file = tmp_path / "samples.txt"
    sample_file.write_text("# selection\nS3\nS1\nS9\n")
    _merging(tmp_path, sample_file).execute()
    major = pd.read_table(tmp_path / "merged" / "out_resistance_major_per_kb.tsv", index_col=0)
    assert list(major.columns) == ["S3", "S1"]
    long = pd.read_table(tmp_path / "merged" / "out_resistance.tsv")
    assert long["sample"].unique().tolist() == ["S3", "S1"]


@pytest.mark.parametrize(
    "changed",
    [
        {"reference": "other"},
        {"min_depth": 5},
        {"annotations": ["resfinder"]},
        {"target_kos": ["K02469"]},
    ],
    ids=["catalogue", "min_depth", "annotations", "target_kos"],
)
def test_execute_incompatible_samples(tmp_path: Path, changed: dict) -> None:
    root = tmp_path / "resistance"
    _sample(root, "S1", [{"gene_id": 1}])
    _sample(root, "S2", [{"gene_id": 1}], **changed)
    with pytest.raises(SystemExit):
        _merging(tmp_path).execute()


def test_execute_errors(tmp_path: Path) -> None:
    (tmp_path / "resistance").mkdir()
    with pytest.raises(SystemExit):  # no sample
        _merging(tmp_path).execute()
    _sample(tmp_path / "resistance", "S1", [{"gene_id": 1}])
    _sample(tmp_path / "resistance" / "again", "S1", [{"gene_id": 1}])
    with pytest.raises(SystemExit):  # the same sample twice
        _merging(tmp_path).execute()


def test_execute_missing_table(tmp_path: Path) -> None:
    _sample(tmp_path / "resistance", "S1", [{"gene_id": 1}])
    (tmp_path / "resistance" / "S1" / "S1_resistance.tsv").unlink()
    with pytest.raises(SystemExit):
        _merging(tmp_path).execute()


def test_execute_positions(tmp_path: Path) -> None:
    """Cohort residue of each position and, per sample and family, the genes
    read with a resistance residue: a residue carried by most samples is the
    wild type of the species, not an acquired mutation"""
    root = tmp_path / "resistance"
    leucine = {"major_residue": "L", "major_is_resistance": True, "resistance_reads": 20,
               "resistance_frequency": 1.0, "known_resistance": True, "residue_counts": "L:20"}
    # rpsL gene 11: K43R everywhere (wild type of its species)
    rpsl = {"gene_id": 11, "family": "rpsL", "codon_pos": 43, "resistance_residues": "R",
            "major_residue": "R", "major_is_resistance": True, "resistance_reads": 20,
            "resistance_frequency": 1.0, "known_resistance": True, "residue_counts": "R:20"}
    _sample(root, "S1", [{
        "gene_id": 10, "target": "gyrA", "reads_resistant": 4, "reads_wt": 16,
        "reads_other": 0, "reads_linked": 0,
    }], positions=[
        # a minor resistant population, S in the majority
        {"resistance_reads": 4, "resistance_frequency": 0.2, "known_resistance": True,
         "major_frequency": 0.8, "residue_counts": "S:16,L:4"},
        rpsl,
    ])
    # gene 12: another catalogue allele of gyrA in M1, carrying L83 at a
    # lower depth: the sample votes S (20 reads) over L (10 reads)
    _sample(root, "S2", [
        {"gene_id": 10, "target": "gyrA", "reads_resistant": 0, "reads_wt": 20,
         "reads_other": 0, "reads_linked": 0},
        {"gene_id": 12, "target": "gyrA", "reads_resistant": 10, "reads_wt": 0,
         "reads_other": 0, "reads_linked": 0},
    ], positions=[
        {}, rpsl, {**leucine, "gene_id": 12, "depth": 10, "resistance_reads": 10,
                   "residue_counts": "L:10"}
    ])
    _sample(root, "S4", [{
        "gene_id": 10, "target": "gyrA", "reads_resistant": 0, "reads_wt": 20,
        "reads_other": 0, "reads_linked": 0,
    }], positions=[{}])
    _sample(root, "S3", [{
        "gene_id": 10, "target": "gyrA", "reads_resistant": 20, "reads_wt": 0,
        "reads_other": 0, "reads_linked": 0,
    }], positions=[
        leucine,
        {**rpsl, "depth": 5, "readable": False, "known_resistance": False},  # not read: ignored
        {"codon_pos": 87, "depth": 0, "readable": False, "residue_counts": ""},
    ])
    _merging(tmp_path).execute()
    out = tmp_path / "merged"
    long = pd.read_table(out / "out_resistance_positions.tsv")
    assert list(long.columns[:2]) == ["sample", "gene_id"]
    gyra = long[(long["gene_id"] == 10) & (long["codon_pos"] == 83)].set_index("sample")
    # one vote per individual (issue 156): S for S1, S2 and S4, L for S3
    assert set(gyra["cohort_residue"]) == {
        "S"
    } and set(gyra["cohort_samples"]) == {4}
    assert (gyra["cohort_share"] == 0.75).all()
    assert gyra["differs_from_cohort"].to_dict() == {
        "S1": False, "S2": False, "S3": True, "S4": False
    }
    # the allele of gene 12 carries L, but the sample votes S: acquired
    gene_12 = long[long["gene_id"] == 12].iloc[0]
    assert gene_12["cohort_residue"] == "S" and gene_12["differs_from_cohort"]
    per_sample = pd.read_table(out / "out_resistance_positions_samples.tsv")
    per_sample = per_sample.set_index(["sample", "family"])
    assert per_sample.loc[("S1", "gyrA"), "genes_readable"] == 1
    assert per_sample.loc[("S1", "gyrA"), "fraction_resistance"] == 1  # minor population
    assert per_sample.loc[("S1", "gyrA"), "fraction_major_resistance"] == 0
    assert per_sample.loc[("S3", "gyrA"), "fraction_acquired"] == 1
    assert per_sample.loc[("S2", "gyrA"), "genes_readable"] == 2
    assert per_sample.loc[("S2", "gyrA"), "fraction_acquired"] == 0.5
    # read totals of the family genes, summed over the genes (issue 157)
    assert per_sample.loc[("S1", "gyrA"), "reads_resistant"] == 4
    assert per_sample.loc[("S1", "gyrA"), "resistant_reads_fraction"] == 4 / 20
    assert per_sample.loc[("S3", "gyrA"), "resistant_reads_fraction"] == 1
    assert per_sample.loc[("S4", "gyrA"), "resistant_reads_fraction"] == 0
    assert per_sample.loc[("S2", "gyrA"), "reads_resistant"] == 10  # gene 12 only
    assert per_sample.loc[("S2", "gyrA"), "reads_wt"] == 20  # gene 10 only
    assert per_sample.loc[("S2", "gyrA"), "resistant_reads_fraction"] == pytest.approx(10 / 30)
    # rpsL: resistance residue as the cohort residue, never acquired
    assert per_sample.loc[("S1", "rpsL"), "fraction_major_resistance"] == 1
    assert per_sample.loc[("S1", "rpsL"), "genes_acquired"] == 0
    assert ("S3", "rpsL") not in per_sample.index
    census = json.loads((out / "out_resistance_census.json").read_text())["resistance_summary"]
    assert census["positions_read"] == 7 and census["positions_with_resistance"] == 5


def test_execute_positions_several_families(tmp_path: Path) -> None:
    """One row per sample and family when several target families are read,
    with the reads of the genes of that family only; a gene of two KOs
    (target "gyrB,parE") counts in the family of its positions"""
    root = tmp_path / "resistance"
    gyrb = {"gene_id": 13, "family": "gyrB", "codon_pos": 426, "ref_residue": "D",
            "resistance_residues": "N", "major_residue": "D"}
    rpsl = {"gene_id": 11, "family": "rpsL", "codon_pos": 43, "ref_residue": "K",
            "resistance_residues": "R", "major_residue": "K"}
    for name in ("S1", "S2"):
        _sample(root, name, [
            {"gene_id": 10, "target": "gyrA", "reads_resistant": 4, "reads_wt": 16,
             "reads_other": 0, "reads_linked": 0},
            {"gene_id": 11, "target": "rpsL", "reads_resistant": 0, "reads_wt": 30,
             "reads_other": 1, "reads_linked": 0},
            {"gene_id": 13, "target": "gyrB,parE", "reads_resistant": 5, "reads_wt": 5,
             "reads_other": 0, "reads_linked": 0},
        ], positions=[{}, rpsl, gyrb])
    _merging(tmp_path).execute()
    per_sample = pd.read_table(tmp_path / "merged" / "out_resistance_positions_samples.tsv")
    assert "target" not in per_sample.columns
    assert len(per_sample) == 6 and not per_sample.duplicated(["sample", "family"]).any()
    per_sample = per_sample.set_index(["sample", "family"])
    assert per_sample.loc[("S1", "gyrA"), "reads_resistant"] == 4
    assert per_sample.loc[("S1", "gyrA"), "reads_wt"] == 16
    assert per_sample.loc[("S1", "rpsL"), "reads_wt"] == 30
    assert per_sample.loc[("S1", "rpsL"), "resistant_reads_fraction"] == 0
    assert per_sample.loc[("S1", "gyrB"), "resistant_reads_fraction"] == 0.5


def test_execute_positions_missing_file(tmp_path: Path) -> None:
    _sample(tmp_path / "resistance", "S1", [{"gene_id": 10}], positions=[{}])
    (tmp_path / "resistance" / "S1" / "S1_resistance_positions.tsv").unlink()
    with pytest.raises(SystemExit):
        _merging(tmp_path).execute()


def test_execute_positions_old_tables(tmp_path: Path) -> None:
    """Sample tables of an older meteor resistance, without the read totals:
    the positions still merge, the read columns are left out"""
    root = tmp_path / "resistance"
    for name in ("S1", "S2"):
        _sample(root, name, [{"gene_id": 10, "target": "gyrA"}], positions=[{}])
        path = root / name / f"{name}_resistance.tsv"
        table = pd.read_table(path)
        table = table.drop(
            columns=["reads_resistant", "reads_wt", "reads_other", "reads_linked"]
        )
        table.to_csv(path, sep="\t", index=False, na_rep="NA")
    _merging(tmp_path).execute()
    per_sample = pd.read_table(tmp_path / "merged" / "out_resistance_positions_samples.tsv")
    assert "reads_resistant" not in per_sample.columns
    assert int(per_sample["genes_readable"].sum()) == 2


def test_residue_classes_map() -> None:
    """The antibiotic classes of the resistance residues (issue 155)"""
    assert residue_classes_map("S:QUINOLONE;L:QUINOLONE/AMINOCOUMARIN") == {
        "S": frozenset({"QUINOLONE"}),
        "L": frozenset({"QUINOLONE", "AMINOCOUMARIN"}),
    }
    assert residue_classes_map("S:QUINOLONE") == {"S": frozenset({"QUINOLONE"})}
    assert residue_classes_map("") == {}
    assert residue_classes_map(float("nan")) == {}
    assert residue_classes_map("S:") == {}


def test_execute_positions_residue_classes(tmp_path: Path) -> None:
    """Per antibiotic class counts: the cohort residue is S, one sample carries
    the QUINOLONE residue L as its major: it is the only one with a quinolone
    major resistance, and it is acquired"""
    root = tmp_path / "resistance"
    quinolone = "L:QUINOLONE;W:QUINOLONE"
    _sample(root, "S1", [{"gene_id": 10, "target": "gyrA"}], positions=[
        {"major_residue": "L", "major_is_resistance": True, "resistance_reads": 20,
         "resistance_frequency": 1.0, "known_resistance": True, "residue_counts": "L:20",
         "residue_classes": quinolone},
    ])
    _sample(root, "S2", [{"gene_id": 10, "target": "gyrA"}], positions=[{"residue_classes": quinolone}])
    _sample(root, "S3", [{"gene_id": 10, "target": "gyrA"}], positions=[{"residue_classes": quinolone}])
    _merging(tmp_path).execute()
    out = tmp_path / "merged"
    long = pd.read_table(out / "out_resistance_positions.tsv")
    assert long["residue_classes"].tolist() == [quinolone] * 3
    assert set(long["cohort_residue"]) == {"S"}
    per_sample = pd.read_table(out / "out_resistance_positions_samples.tsv")
    per_sample = per_sample.set_index(["sample", "family"])
    assert per_sample.loc[("S1", "gyrA"), "genes_QUINOLONE_readable"] == 1
    assert per_sample.loc[("S1", "gyrA"), "genes_QUINOLONE_major_resistance"] == 1
    assert per_sample.loc[("S1", "gyrA"), "genes_QUINOLONE_acquired"] == 1
    assert per_sample.loc[("S1", "gyrA"), "fraction_QUINOLONE_major_resistance"] == 1
    assert per_sample.loc[("S2", "gyrA"), "genes_QUINOLONE_readable"] == 1
    assert per_sample.loc[("S2", "gyrA"), "genes_QUINOLONE_major_resistance"] == 0
    assert per_sample.loc[("S2", "gyrA"), "fraction_QUINOLONE_major_resistance"] == 0


def _consensus_cohort(root: Path, n: int = 6) -> None:
    """n samples: gene 1 (tet(W)) interpretable in all of them, in two allele
    groups; gene 2 (gyrA) in 3 samples only (no tree with -m 4)"""
    rng = random.Random(1)
    base = "".join(rng.choice("ACGT") for _ in range(600))
    other = "".join(
        ("A" if b != "A" else "C") if i % 20 == 0 else b for i, b in enumerate(base)
    )
    for k in range(n):
        allele = base if k < n // 2 else other
        # a few sample-specific sites and an uncovered end
        sequence = list(allele)
        sequence[100 + k] = "G" if sequence[100 + k] != "G" else "T"
        sequence[-30:] = "N" * 30
        consensus = {1: "".join(sequence)}
        rows = [{"gene_id": 1, "gene_length": 600}]
        if k < 3:
            consensus[2] = base
            rows.append({"gene_id": 2, "msp_name": "M2", "resfinder": pd.NA, "target": "gyrA"})
        _sample(root, f"S{k}", rows, consensus=consensus)


def test_execute_trees(tmp_path: Path) -> None:
    _consensus_cohort(tmp_path / "resistance")
    _merging(tmp_path, trees=True, min_samples=4, gtr=False).execute()
    out = tmp_path / "merged"
    trees = pd.read_table(out / "out_resistance_trees.tsv")
    # gene 2 is interpretable in 3 samples only
    assert trees["gene_id"].tolist() == [1]
    row = trees.iloc[0]
    assert row["family"] == "tet(W)" and row["msp_name"] == "M1"
    assert row["samples"] == 6 and row["tree"] == "gene_1.tree"
    tree_dir = out / "out_resistance_trees"
    newick = (tree_dir / "gene_1.tree").read_text()
    assert all(f"S{k}" in newick for k in range(6))
    assert (tree_dir / "gene_1_raw_distances.tsv").exists()
    distances = pd.read_table(tree_dir / "gene_1.tsv", index_col=0)
    # the two allele groups are further apart than the samples of a group
    assert distances.loc["S0", "S1"] < distances.loc["S0", "S5"]
    census = json.loads((out / "out_resistance_census.json").read_text())["resistance_summary"]
    assert census["trees"] == 1 and census["tree_mode"] == "TN93"
    assert census["tree_min_samples"] == 4
    # same sequences: the tree is reused; a sample less: recomputed
    before = (tree_dir / "gene_1.tree").stat().st_mtime_ns
    _merging(tmp_path, trees=True, min_samples=4, gtr=False).execute()
    assert (tree_dir / "gene_1.tree").stat().st_mtime_ns == before
    sample_file = tmp_path / "samples.txt"
    sample_file.write_text("\n".join(f"S{k}" for k in range(5)))
    _merging(tmp_path, sample_file, trees=True, min_samples=4, gtr=False).execute()
    assert "S5" not in (tree_dir / "gene_1.tree").read_text()


def test_execute_trees_family_priority(tmp_path: Path) -> None:
    root = tmp_path / "resistance"
    for k in range(4):
        _sample(
            root, f"S{k}",
            [{"gene_id": 7, "resfinder": pd.NA, "resfinderfg": pd.NA, "mustard": "tetM"}],
            consensus={7: "ACGT" * 50},
        )
    # fewer sites than -s: no tree
    _merging(tmp_path, trees=True, min_samples=4, min_info_sites=1000, gtr=False).execute()
    trees = pd.read_table(tmp_path / "merged" / "out_resistance_trees.tsv")
    assert trees.loc[0, "family"] == "tetM"
    assert pd.isna(trees.loc[0, "tree"])


def test_execute_trees_without_consensus(three_samples: Path, tmp_path: Path) -> None:
    """Samples of an older meteor resistance (no consensus file)"""
    with pytest.raises(SystemExit):
        _merging(tmp_path, trees=True).execute()
