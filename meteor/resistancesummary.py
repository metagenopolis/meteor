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
