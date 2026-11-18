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

"""Tables and trees of the variants of the antibiotic resistance genes of
several samples (meteor resistance_summary)"""

import fnmatch
import logging
import lzma
import os
import sys
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from shutil import rmtree
from tempfile import mkdtemp
from time import perf_counter
from typing import ClassVar

import numpy as np
import pandas as pd

from meteor.session import Component, Session

CENSUS_PATTERN = "*_census_stage_resistance.json"
# parameters that must be the same in every sample: rates are not comparable otherwise
SHARED_PARAMETERS = (
    "annotations",
    "target_kos",
    "min_depth",
    "min_alt_count",
    "min_frequency",
    "min_gene_coverage",
    "base_quality",
    "core_size",
)
# per-gene columns, the same in every sample
GENE_COLUMNS = (
    "msp_name",
    "resfinder",
    "resfinderfg",
    "mustard",
    "mustard_prediction",
    "target",
    "gene_length",
)
# annotations giving the family of a gene, by priority
FAMILY_COLUMNS = ("resfinder", "resfinderfg", "mustard", "target")


def residue_classes_map(value) -> dict[str, frozenset[str]]:
    """`residue_classes` of the resistance_positions annotation, e.g.
    "S:QUINOLONE;L:QUINOLONE/AMINOCOUMARIN" ->
    {"S": frozenset({"QUINOLONE"}), "L": frozenset({"QUINOLONE", "AMINOCOUMARIN"})}"""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return {}
    classes: dict[str, set[str]] = {}
    for part in str(value).split(";"):
        if ":" not in part:
            continue
        residue, names = part.split(":", 1)
        if residue:
            classes.setdefault(residue, set()).update(
                name for name in (n.strip() for n in names.split("/")) if name
            )
    return {residue: frozenset(names) for residue, names in classes.items() if names}


# outputs of a gene tree, removed when its sequences change
TREE_OUTPUTS = (
    ".tree", ".tsv", "_raw_distances.tsv", "_comparison.tsv", "_qcflags.tsv", "_shrink.tsv",
)


@dataclass
class ResistanceSummary(Session):
    """Merge the meteor resistance tables of several samples and, with trees,
    build one tree per resistance gene from the consensus sequences of the samples"""

    DEFAULT_PREFIX: ClassVar[str] = "output"
    DEFAULT_MIN_SAMPLES: ClassVar[int] = 4
    DEFAULT_MAX_GAP: ClassVar[float] = 1.0
    DEFAULT_MIN_INFO_SITES: ClassVar[int] = 10
    DEFAULT_NUM_THREADS: ClassVar[int] = 1
    MATRICES: ClassVar[tuple[str, ...]] = ("major_per_kb", "minor_per_kb")

    meteor: type[Component]
    prefix: str = DEFAULT_PREFIX
    sample_file: Path | None = None
    trees: bool = False
    min_samples: int = DEFAULT_MIN_SAMPLES
    max_gap: float = DEFAULT_MAX_GAP
    min_info_sites: int = DEFAULT_MIN_INFO_SITES
    gtr: bool = False

    def find_samples(self) -> dict[str, tuple[Path, dict]]:
        """{sample name: (resistance table, census)} of the samples of the
        input directory (and of the sample list)"""
        census_files = [
            Path(root, filename)
            for root, _, files in os.walk(self.meteor.resistance_dir, followlinks=True)
            for filename in fnmatch.filter(files, CENSUS_PATTERN)
        ]
        if not census_files:
            logging.error("No %s file found in %s", CENSUS_PATTERN, self.meteor.resistance_dir)
            sys.exit(1)
        samples: dict[str, tuple[Path, dict]] = {}
        for census_file in sorted(census_files):
            census = self.read_json(census_file)
            name = census["sample_info"]["sample_name"]
            if name in samples:
                logging.error(
                    "Sample %s is found twice: %s and %s",
                    name,
                    samples[name][0].parent,
                    census_file.parent,
                )
                sys.exit(1)
            table = census_file.parent / census["resistance"]["resistance_file"]
            if not table.exists():
                logging.error("%s is missing", table)
                sys.exit(1)
            samples[name] = (table, census)
        logging.info("%d samples found in %s", len(samples), self.meteor.resistance_dir)
        return self.select_samples(samples, self.sample_file)

    @staticmethod
    def check_samples(samples: dict[str, tuple[Path, dict]]) -> dict:
        """Same catalogue and same parameters in every sample"""
        first_name, (_, first) = next(iter(samples.items()))
        reference = first["mapping"]["reference_name"]
        parameters = {key: first["resistance"].get(key) for key in SHARED_PARAMETERS}
        for name, (_, census) in samples.items():
            if census["mapping"]["reference_name"] != reference:
                logging.error(
                    "%s was analysed with the catalogue %s, %s with %s",
                    name,
                    census["mapping"]["reference_name"],
                    first_name,
                    reference,
                )
                sys.exit(1)
            different = [
                key
                for key in SHARED_PARAMETERS
                if census["resistance"].get(key) != parameters[key]
            ]
            if different:
                logging.error(
                    "%s and %s were analysed with different parameters (%s): "
                    "their rates are not comparable",
                    name,
                    first_name,
                    ", ".join(different),
                )
                sys.exit(1)
        return {"reference_name": reference, **parameters}

    @staticmethod
    def gene_summary(table: pd.DataFrame, n_samples: int) -> pd.DataFrame:
        """One row per gene: its annotations, the samples where it has reads and
        where it is interpretable, and its median rates and those of its controls
        over the samples where it is interpretable"""
        annotations = table.groupby("gene_id")[
            [c for c in GENE_COLUMNS if c in table]
        ].first()
        summary = annotations.assign(
            samples_with_reads=table.groupby("gene_id").size(),
            samples_interpretable=table.groupby("gene_id")["interpretable"].sum(),
        )
        summary["interpretable_fraction"] = summary["samples_interpretable"] / n_samples
        interpretable = table[table["interpretable"]]
        medians = interpretable.groupby("gene_id")[
            ["major_per_kb", "minor_per_kb", "control_major_per_kb", "control_minor_per_kb"]
        ].median()
        summary = summary.join(medians.add_prefix("median_"))
        # samples where the gene has more variants than its controls
        compared = interpretable[interpretable["control_major_per_kb"].notna()]
        summary["samples_with_controls"] = compared.groupby("gene_id").size()
        for kind in ("major", "minor"):
            above = compared[f"{kind}_per_kb"] > compared[f"control_{kind}_per_kb"]
            summary[f"{kind}_above_controls_fraction"] = (
                above.groupby(compared["gene_id"]).mean()
            )
        summary["samples_with_controls"] = (
            summary["samples_with_controls"].fillna(0).astype(int)
        )
        summary["samples_interpretable"] = summary["samples_interpretable"].astype(int)
        return summary.reset_index().sort_values("gene_id")

    @staticmethod
    def merge_positions(samples: dict[str, tuple[Path, dict]]) -> pd.DataFrame | None:
        """Residues of every sample at the known resistance positions, with the
        cohort residue of each position: the individual residue of most of the
        individuals reading the position. The individual residue of a sample is
        the residue of its reads at the position, over the genes of the family
        in the MSP, weighted by depth (several catalogue genes can be alleles
        of the gene of one species, and one of them carrying a resistance
        residue must not outvote the sample; a gene of no MSP stands alone). A
        resistance residue that is also the cohort residue is the wild type of
        the species carrying the gene (the positions come from E. coli and
        K. pneumoniae), not an acquired mutation.

        :return: long table, None when the samples have no positions table
        """
        files = {}
        for name, (table, census) in samples.items():
            positions_name = census["resistance"].get("positions_file")
            if positions_name is None:
                continue
            path = table.parent / positions_name
            if not path.exists():
                logging.error("%s is missing", path)
                sys.exit(1)
            files[name] = path
        if not files:
            return None
        residue_columns = ("ref_residue", "resistance_residues", "gene_residue", "major_residue")
        tables = [
            pd.read_table(path, dtype=dict.fromkeys(residue_columns, str)).assign(sample=name)
            for name, path in files.items()
        ]
        # samples without any position read: left out of the concatenation
        long = pd.concat([t for t in tables if not t.empty] or tables[:1], ignore_index=True)
        for column in ("readable", "major_is_resistance", "known_resistance"):
            long[column] = long[column].astype(bool)
        key = ["cohort_group", "family", "codon_pos"]
        long["cohort_group"] = long["msp_name"].where(
            long["msp_name"].notna(), "gene_" + long["gene_id"].astype(str)
        )
        readable = long[long["readable"]]
        # one vote per individual: the residue of its reads at the position,
        # over the genes of the family in the MSP, weighted by depth (most
        # reads first, then alphabetical: a stable individual residue)
        individual = (
            readable.groupby(["sample", *key, "major_residue"], sort=False)["depth"].sum()
            .reset_index()
            .rename(columns={"depth": "individual_reads"})
            .sort_values(
                ["sample", *key, "individual_reads", "major_residue"],
                ascending=[True] * 4 + [False, True],
            )
            .groupby(["sample", *key], sort=False)["major_residue"]
            .first()
            .rename("individual_residue")
            .reset_index()
        )
        # the cohort residue: the individual residue of most individuals
        votes = (
            individual.assign(major_residue=individual["individual_residue"])
            .groupby([*key, "major_residue"], sort=False)
            .size()
            .rename("individuals")
            .reset_index()
            .sort_values(
                [*key, "individuals", "major_residue"], ascending=[True] * 3 + [False, True]
            )
        )
        cohort = (
            votes.groupby(key, sort=False)["major_residue"]
            .first()
            .rename("cohort_residue")
            .to_frame()
        )
        cohort["cohort_samples"] = individual.groupby(key, sort=False)["sample"].nunique()
        cohort["cohort_share"] = (
            votes.groupby(key, sort=False)["individuals"].max() / cohort["cohort_samples"]
        )
        long = long.join(cohort, on=key).drop(columns="cohort_group")
        long["cohort_samples"] = long["cohort_samples"].fillna(0).astype(int)
        long["cohort_share"] = long["cohort_share"].fillna(0.0)
        long["differs_from_cohort"] = long["readable"] & (
            long["major_residue"] != long["cohort_residue"]
        )
        return long[["sample", *[c for c in long.columns if c != "sample"]]]

    @staticmethod
    def position_samples(
        long: pd.DataFrame, table: pd.DataFrame | None = None
    ) -> pd.DataFrame:
        """Per sample and target family: genes read at >= 1 known resistance
        position, and how many carry a resistance residue (minor population
        included), as their major residue, and as a major residue that differs
        from the cohort residue (candidate acquired mutation). With `table`
        (the merged per-gene table): the reads of the family genes at the
        known positions, resistant (a resistance residue on at least one),
        wt (a codon read, no resistance residue), other (no codon read),
        linked (a resistance residue on at least two positions, e.g. gyrA
        S83L + D87N) and the fraction of the resistant reads over the reads
        reading a position"""
        readable = long[long["readable"]].assign(
            acquired=lambda t: t["major_is_resistance"] & t["differs_from_cohort"]
        )
        genes = readable.groupby(["sample", "family", "gene_id"])[
            ["known_resistance", "major_is_resistance", "acquired"]
        ].any()
        per_sample = genes.groupby(["sample", "family"]).agg(
            genes_readable=("known_resistance", "size"),
            genes_resistance=("known_resistance", "sum"),
            genes_major_resistance=("major_is_resistance", "sum"),
            genes_acquired=("acquired", "sum"),
        )
        for column in ("resistance", "major_resistance", "acquired"):
            per_sample[f"fraction_{column}"] = (
                per_sample[f"genes_{column}"] / per_sample["genes_readable"]
            )
        read_totals = ("reads_resistant", "reads_wt", "reads_other", "reads_linked")
        if table is not None and "target" in table and all(c in table for c in read_totals):
            reads = (
                table[table["target"].notna()]
                .groupby(["sample", "target"])[list(read_totals)]
                .sum(min_count=1)
                .rename(columns={"target": "family"})
            )
            per_sample = per_sample.join(reads, how="left")
            per_sample["resistant_reads_fraction"] = np.where(
                per_sample["reads_resistant"].notna()
                & (per_sample["reads_resistant"] + per_sample["reads_wt"] > 0),
                per_sample["reads_resistant"] / (per_sample["reads_resistant"] + per_sample["reads_wt"]),
                np.nan,
            )
        if "residue_classes" in long.columns:
            # per antibiotic class of the known residues: the genes read at one
            # position of the class at least, and among them the genes whose
            # major residue is a resistance residue of the class (and differs
            # from the cohort residue); without the column, the family is its
            # only class (the current columns)
            entries: dict[tuple, dict[str, set]] = {}
            for row in readable.itertuples(index=False):
                position_classes = residue_classes_map(row.residue_classes)
                if not position_classes:
                    continue
                entry = entries.setdefault(
                    (row.sample, row.family, row.gene_id),
                    {"positions": set(), "major": set(), "acquired": set()},
                )
                entry["positions"].update(set().union(*position_classes.values()))
                entry["major"].update(position_classes.get(row.major_residue, ()))
                if row.differs_from_cohort:
                    entry["acquired"].update(position_classes.get(row.major_residue, ()))
            all_classes = (
                sorted(set().union(*[e["positions"] for e in entries.values()])) if entries else []
            )
            classes_table = pd.DataFrame(
                [
                    {
                        "sample": sample,
                        "family": family,
                        **{c: int(c in e["positions"]) for c in all_classes},
                        **{f"major_{c}": int(c in e["major"]) for c in all_classes},
                        **{f"acquired_{c}": int(c in e["acquired"]) for c in all_classes},
                    }
                    for (sample, family, _), e in entries.items()
                ]
            )
            per_class = classes_table.groupby(["sample", "family"]).sum()
            per_sample = per_sample.join(
                per_class.rename(
                    columns={
                        **{c: f"genes_{c}_readable" for c in all_classes},
                        **{f"major_{c}": f"genes_{c}_major_resistance" for c in all_classes},
                        **{f"acquired_{c}": f"genes_{c}_acquired" for c in all_classes},
                    }
                ),
                how="left",
            )
            class_columns = [
                f"genes_{c}_{suffix}"
                for c in all_classes
                for suffix in ("readable", "major_resistance", "acquired")
            ]
            per_sample[class_columns] = per_sample[class_columns].fillna(0).astype(int)
            for c in all_classes:
                for prefix in ("major_resistance", "acquired"):
                    per_sample[f"fraction_{c}_{prefix}"] = (
                        per_sample[f"genes_{c}_{prefix}"] / per_sample[f"genes_{c}_readable"]
                    )
        return per_sample.reset_index()

    @staticmethod
    def read_consensus(path: Path) -> dict[int, str]:
        """{gene_id: consensus sequence} of one sample"""
        sequences: dict[int, str] = {}
        with lzma.open(path, "rt", encoding="UTF-8") as fasta:
            gene = None
            for line in fasta:
                line = line.strip()
                if line.startswith(">"):
                    gene = int(line[1:])
                elif gene is not None and line:
                    sequences[gene] = sequences.get(gene, "") + line
        return sequences

    def build_trees(
        self, samples: dict[str, tuple[Path, dict]], summary: pd.DataFrame
    ) -> pd.DataFrame:
        """One tree per resistance gene interpretable in at least min_samples
        samples, from the consensus sequences of these samples (the sequences
        are on the coordinates of the catalogue gene: no alignment needed).

        :return: one row per gene with a FASTA: family, annotations, samples, tree
        """
        from meteor.phylogeny import Phylogeny

        tree_dir = self.meteor.merging_dir / f"{self.prefix}_resistance_trees"
        tree_dir.mkdir(exist_ok=True, parents=True)
        per_gene: dict[int, dict[str, str]] = defaultdict(dict)
        for name, (table, census) in samples.items():
            consensus_name = census["resistance"].get("consensus_file")
            consensus_file = table.parent / consensus_name if consensus_name else None
            if consensus_file is None or not consensus_file.exists():
                logging.error(
                    "%s has no consensus sequences: run meteor resistance again "
                    "to build the trees",
                    name,
                )
                sys.exit(1)
            for gene, sequence in self.read_consensus(consensus_file).items():
                per_gene[gene][name] = sequence
        fasta_files = []
        index = []
        for gene in sorted(per_gene):
            sequences = per_gene[gene]
            if len(sequences) < self.min_samples:
                continue
            fasta_file = tree_dir / f"gene_{gene}.fasta"
            content = "".join(f">{name}\n{seq}\n" for name, seq in sequences.items())
            # sequences changed (samples, parameters): the tree is recomputed
            if not fasta_file.exists() or fasta_file.read_text(encoding="UTF-8") != content:
                for suffix in TREE_OUTPUTS:
                    (tree_dir / f"{fasta_file.stem}{suffix}").unlink(missing_ok=True)
                fasta_file.write_text(content, encoding="UTF-8")
            fasta_files.append(fasta_file)
            index.append({"gene_id": gene, "samples": len(sequences)})
        logging.info(
            "%d resistance genes interpretable in at least %d samples: one tree each",
            len(fasta_files),
            self.min_samples,
        )
        self.meteor.tree_dir = tree_dir
        self.meteor.tmp_dir = Path(mkdtemp(dir=self.meteor.tmp_path))
        try:
            phylogeny = Phylogeny(
                self.meteor, fasta_files, self.max_gap, self.min_info_sites, self.gtr
            )
            if fasta_files:
                phylogeny.execute()
        finally:
            rmtree(self.meteor.tmp_dir, ignore_errors=True)
        built = {tree.name for tree in phylogeny.tree_files}
        annotations = summary.set_index("gene_id")
        columns = ["msp_name", *[c for c in FAMILY_COLUMNS if c in annotations]]
        trees = pd.DataFrame(index, columns=["gene_id", "samples"])
        trees = trees.join(annotations[columns], on="gene_id")
        family = pd.Series(pd.NA, index=trees.index, dtype=object)
        for column in reversed([c for c in FAMILY_COLUMNS if c in trees]):
            family = trees[column].where(trees[column].notna(), family)
        trees.insert(1, "family", family)
        trees["tree"] = [
            f"gene_{gene}.tree" if f"gene_{gene}.tree" in built else pd.NA
            for gene in trees["gene_id"]
        ]
        return trees

    def execute(self) -> None:
        """Merge the resistance tables of the samples"""
        start = perf_counter()
        samples = self.find_samples()
        parameters = self.check_samples(samples)
        names = list(samples)
        table = pd.concat(
            [
                pd.read_table(path).assign(sample=name)
                for name, (path, _) in samples.items()
            ],
            ignore_index=True,
        )
        table["interpretable"] = table["interpretable"].astype(bool)
        table = table[["sample", *[c for c in table.columns if c != "sample"]]]
        self.meteor.merging_dir.mkdir(exist_ok=True, parents=True)
        out = self.meteor.merging_dir
        long_file = out / f"{self.prefix}_resistance.tsv"
        table.to_csv(long_file, sep="\t", index=False, na_rep="NA")
        # gene x sample matrices: rates where the gene is interpretable, NA otherwise
        for value in self.MATRICES:
            matrix = (
                table[table["interpretable"]]
                .pivot(index="gene_id", columns="sample", values=value)
                .reindex(columns=names)
                .sort_index()
            )
            matrix.to_csv(
                out / f"{self.prefix}_resistance_{value}.tsv", sep="\t", na_rep="NA"
            )
        interpretable = (
            table.assign(interpretable=table["interpretable"].astype(np.int8))
            .pivot(index="gene_id", columns="sample", values="interpretable")
            .reindex(columns=names)
            .fillna(0)
            .astype(np.int8)
            .sort_index()
        )
        interpretable.to_csv(out / f"{self.prefix}_resistance_interpretable.tsv", sep="\t")
        summary = self.gene_summary(table, len(names))
        summary.to_csv(
            out / f"{self.prefix}_resistance_genes.tsv", sep="\t", index=False, na_rep="NA"
        )
        positions = self.merge_positions(samples)
        positions_census = {}
        if positions is not None:
            positions.to_csv(
                out / f"{self.prefix}_resistance_positions.tsv", sep="\t", index=False, na_rep="NA"
            )
            per_sample = self.position_samples(positions, table)
            per_sample.to_csv(
                out / f"{self.prefix}_resistance_positions_samples.tsv",
                sep="\t",
                index=False,
                na_rep="NA",
            )
            positions_census = {
                "positions_read": int(positions["readable"].sum()),
                "positions_with_resistance": int(positions["known_resistance"].sum()),
            }
        tree_count = None
        if self.trees:
            trees = self.build_trees(samples, summary)
            trees.to_csv(
                out / f"{self.prefix}_resistance_trees.tsv", sep="\t", index=False, na_rep="NA"
            )
            tree_count = int(trees["tree"].notna().sum())
        self.save_config(
            {
                "meteor_version": self.meteor.version,
                "resistance_summary": {
                    "summary_date": datetime.now().strftime("%Y-%m-%d"),
                    "samples": len(names),
                    "genes": int(table["gene_id"].nunique()),
                    "interpretable_genes": int(
                        table.loc[table["interpretable"], "gene_id"].nunique()
                    ),
                    **parameters,
                    **positions_census,
                    "trees": tree_count,
                    **(
                        {
                            "tree_min_samples": self.min_samples,
                            "tree_max_gap": self.max_gap,
                            "tree_min_info_sites": self.min_info_sites,
                            "tree_mode": "GTR" if self.gtr else "TN93",
                        }
                        if self.trees
                        else {}
                    ),
                },
            },
            out / f"{self.prefix}_resistance_census.json",
        )
        logging.info(
            "%d samples, %d resistance genes with reads (%d interpretable in at least "
            "one sample) merged in %f seconds: %s",
            len(names),
            table["gene_id"].nunique(),
            table.loc[table["interpretable"], "gene_id"].nunique(),
            perf_counter() - start,
            long_file,
        )
