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

"""Mapped reads, MSP richness and phylum composition of each sample, as
tables and MultiQC custom content files"""

import json
import logging
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, ClassVar

import pandas as pd


@dataclass
class QCReport:
    """<prefix>_qc.tsv (one row per sample), <prefix>_phylum.tsv and MultiQC
    custom content files <prefix>_*_mqc.json, read by `multiqc <merging dir>`
    (MultiQC is not a dependency of meteor)."""

    QC_COLUMNS: ClassVar[list[str]] = [
        "sample", "total_read_count", "mapped_reads", "mapped_read_rate",
        "msp_richness", "msp_signal_rate",
    ]
    # phyla shown in the composition plot, the others are summed
    TOP_PHYLA: ClassVar[int] = 10

    merging_dir: Path
    prefix: str

    @staticmethod
    def census_value(config: dict, name: str) -> Any:
        """Value of a field in any section of a census (sections differ
        between meteor versions); None when absent"""
        for section in config.values():
            if isinstance(section, dict) and name in section:
                return section[name]
        return None

    def sample_qc(self, census: dict[str, dict]) -> pd.DataFrame:
        """One row per sample from its census. Mapped reads are the reads
        counted on the catalogue genes (after the identity filter); the MSP
        richness is the number of MSPs detected in the sample; the MSP signal
        rate is the % of the gene signal carried by the MSPs.

        :param census: {sample name: census_stage_2 content}, in output order
        """
        rows = []
        for sample, config in census.items():
            rows.append(
                {
                    "sample": sample,
                    "total_read_count": self.census_value(config, "total_read_count"),
                    "mapped_reads": self.census_value(config, "counted_reads"),
                    "msp_richness": self.census_value(config, "msp_count"),
                    "msp_signal_rate": self.census_value(config, "msp_signal"),
                }
            )
        qc = pd.DataFrame(rows, columns=[c for c in self.QC_COLUMNS if c != "mapped_read_rate"])
        for column in qc.columns.drop("sample"):
            qc[column] = pd.to_numeric(qc[column], errors="coerce")
        total = qc["total_read_count"].where(qc["total_read_count"] > 0)
        qc["mapped_read_rate"] = (100 * qc["mapped_reads"] / total).round(2)
        qc["msp_signal_rate"] = (100 * qc["msp_signal_rate"]).round(2)
        return qc[self.QC_COLUMNS]

    def phylum_composition(self, msp_table: pd.DataFrame, taxonomy: pd.DataFrame) -> pd.DataFrame:
        """Relative abundance (%) of each phylum in each sample (samples as rows),
        from the species table before filters; samples without species are left out

        :param msp_table: merged species table (msp_name, then one column per sample)
        :param taxonomy: catalogue taxonomy (msp_name, phylum)
        """
        abundance = msp_table.set_index("msp_name").astype(float).fillna(0.0)
        phylum = (
            taxonomy.drop_duplicates("msp_name").set_index("msp_name")["phylum"]
            .reindex(abundance.index).fillna("").astype(str).str.strip()
            .replace("", "unclassified")
        )
        by_phylum = abundance.groupby(phylum.to_numpy()).sum()
        totals = by_phylum.sum(axis=0)
        by_phylum = 100 * by_phylum.loc[:, totals > 0] / totals[totals > 0]
        return by_phylum.T.reindex(
            columns=by_phylum.mean(axis=1).sort_values(ascending=False).index
        )

    @staticmethod
    def _json_value(value: Any) -> Any:
        """JSON value: NaN as null, numpy numbers as Python numbers"""
        if value is None or (isinstance(value, float) and math.isnan(value)):
            return None
        return value.item() if hasattr(value, "item") else value

    def _write_mqc(self, name: str, content: dict) -> None:
        content = {"parent_id": "meteor", "parent_name": "Meteor", **content}
        with (self.merging_dir / f"{self.prefix}_{name}_mqc.json").open("w", encoding="utf-8") as out:
            json.dump(content, out, indent=1)

    def write(self, qc: pd.DataFrame) -> None:
        """<prefix>_qc.tsv and the MultiQC files"""
        qc_file = self.merging_dir / f"{self.prefix}_qc.tsv"
        qc.to_csv(qc_file, sep="\t", index=False)
        records = {
            row["sample"]: {key: self._json_value(value) for key, value in row.items() if key != "sample"}
            for row in qc.to_dict(orient="records")
        }
        # General statistics table
        headers = {
            "total_read_count": {"title": "Reads", "description": "Reads of the sample", "format": "{:,.0f}", "scale": "Blues"},
            "mapped_reads": {"title": "Mapped reads", "description": "Reads mapped and counted on the catalogue genes", "format": "{:,.0f}", "scale": "Purples"},
            "mapped_read_rate": {"title": "% Mapped", "description": "Reads mapped and counted on the catalogue genes", "suffix": "%", "min": 0, "max": 100, "scale": "RdYlGn"},
            "msp_richness": {"title": "MSP richness", "description": "Species (MSPs) detected", "format": "{:,.0f}", "scale": "Greens"},
            "msp_signal_rate": {"title": "% MSP signal", "description": "Gene signal carried by the MSPs", "suffix": "%", "min": 0, "max": 100, "scale": "YlGn"},
        }
        self._write_mqc(
            "general_stats",
            {
                "id": "meteor_general_stats",
                "plot_type": "generalstats",
                "headers": [{column: header} for column, header in headers.items()],
                "data": {sample: {column: values[column] for column in headers} for sample, values in records.items()},
            },
        )
        # Mapped and unmapped reads per sample
        reads = {}
        for sample, values in records.items():
            total, mapped = values["total_read_count"], values["mapped_reads"]
            if total is not None and mapped is not None:
                reads[sample] = {"Mapped": mapped, "Not mapped": max(total - mapped, 0)}
        self._write_mqc(
            "reads",
            {
                "id": "meteor_reads",
                "section_name": "Meteor: mapped reads",
                "description": "Reads mapped and counted on the genes of the catalogue (after the identity filter) among the reads of each sample.",
                "plot_type": "bargraph",
                "pconfig": {"id": "meteor_reads_plot", "title": "Meteor: mapped reads", "ylab": "Reads"},
                "categories": {"Mapped": {"color": "#437bb1"}, "Not mapped": {"color": "#b1084c"}},
                "data": reads,
            },
        )
        # MSP richness against the number of reads
        self._write_mqc(
            "richness",
            {
                "id": "meteor_richness",
                "section_name": "Meteor: MSP richness",
                "description": "Species (MSPs) detected in each sample against its number of reads: "
                "samples far below the others at a similar depth, or at the low end of the reads, "
                "may be too shallow.",
                "plot_type": "scatter",
                "pconfig": {"id": "meteor_richness_plot", "title": "Meteor: MSP richness", "xlab": "Reads", "ylab": "MSP richness", "xlog": True},
                "data": {
                    sample: {"x": values["total_read_count"], "y": values["msp_richness"]}
                    for sample, values in records.items()
                    if values["total_read_count"] and values["msp_richness"] is not None
                },
            },
        )
        logging.info("Sample QC saved as %s (MultiQC files: %s_*_mqc.json)", qc_file, self.prefix)

    def write_composition(self, composition: pd.DataFrame) -> None:
        """<prefix>_phylum.tsv (every phylum) and the MultiQC composition plot
        (the TOP_PHYLA most abundant phyla on average, the others summed)"""
        composition.round(4).to_csv(
            self.merging_dir / f"{self.prefix}_phylum.tsv", sep="\t", index_label="sample"
        )
        shown = composition.iloc[:, : self.TOP_PHYLA].copy()
        if composition.shape[1] > self.TOP_PHYLA:
            shown["Other phyla"] = composition.iloc[:, self.TOP_PHYLA :].sum(axis=1)
        self._write_mqc(
            "phylum",
            {
                "id": "meteor_phylum",
                "section_name": "Meteor: phylum composition",
                "description": f"Relative abundance of the phyla in each sample (species table before filters; "
                f"the {self.TOP_PHYLA} most abundant phyla on average, the others summed). "
                f"Every phylum is in {self.prefix}_phylum.tsv.",
                "plot_type": "bargraph",
                "pconfig": {"id": "meteor_phylum_plot", "title": "Meteor: phylum composition", "ylab": "Relative abundance (%)", "ymax": 100, "cpswitch": False},
                "data": {
                    sample: {phylum: round(float(value), 4) for phylum, value in row.items()}
                    for sample, row in shown.iterrows()
                },
            },
        )

    def execute(
        self,
        census: dict[str, dict],
        msp_table: pd.DataFrame | None = None,
        taxonomy: pd.DataFrame | None = None,
    ) -> None:
        """Write the report

        :param census: {sample name: census_stage_2 content}, in output order
        :param msp_table: merged species table before filters (None when the
            species table is not merged: no composition)
        :param taxonomy: catalogue taxonomy (msp_name, phylum)
        """
        self.write(self.sample_qc(census))
        if msp_table is None or taxonomy is None:
            logging.info("The species table is not merged: no phylum composition in the report.")
            return
        self.write_composition(self.phylum_composition(msp_table, taxonomy))
