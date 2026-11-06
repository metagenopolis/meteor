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

"""Tables of the variants of the antibiotic resistance genes of several samples
(meteor resistance_merge)"""

import fnmatch
import logging
import os
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
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


@dataclass
class ResistanceMerging(Session):
    """Merge the meteor resistance tables of several samples"""

    DEFAULT_PREFIX: ClassVar[str] = "output"
    MATRICES: ClassVar[tuple[str, ...]] = ("major_per_kb", "minor_per_kb")

    meteor: type[Component]
    prefix: str = DEFAULT_PREFIX
    sample_file: Path | None = None

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
        self.save_config(
            {
                "meteor_version": self.meteor.version,
                "resistance_merging": {
                    "merging_date": datetime.now().strftime("%Y-%m-%d"),
                    "samples": len(names),
                    "genes": int(table["gene_id"].nunique()),
                    "interpretable_genes": int(
                        table.loc[table["interpretable"], "gene_id"].nunique()
                    ),
                    **parameters,
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
