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

"""Strain sharing between samples, from the comparison tables of meteor tree"""

import logging
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import ClassVar

import numpy as np
import pandas as pd


@dataclass
class StrainSharing:
    """Strain-sharing events, sharing rate of each sample pair and
    transmissibility of each MSP.

    Two samples are compared on an MSP when the sites informative in both
    cover at least min_overlap % of the alignment. They share its strain when
    their distance is at most THRESHOLD. The transmissibility of an MSP is the
    fraction of the compared pairs of samples of different individuals that
    share its strain: two samples of the same individual (time points) are not
    a transmission.
    """

    DEFAULT_MIN_OVERLAP: ClassVar[float] = 10.0
    # strain identity threshold of every MSP (same_strain of the comparison tables)
    THRESHOLD: ClassVar[float] = 0.0001
    # partial counts are summed every REDUCE_PAIRS pairs (bounded memory)
    REDUCE_PAIRS: ClassVar[int] = 5_000_000
    EVENTS_FILE: ClassVar[str] = "strain_sharing_events.tsv"
    SAMPLES_FILE: ClassVar[str] = "strain_sharing_samples.tsv"
    MSP_FILE: ClassVar[str] = "strain_sharing_msp.tsv"
    SAME: ClassVar[str] = "same_individual"
    DIFFERENT: ClassVar[str] = "different_individuals"

    output_dir: Path
    individuals_file: Path | None = None
    min_overlap: float = DEFAULT_MIN_OVERLAP
    individuals: dict[str, str] = field(default_factory=dict, init=False)

    def __post_init__(self) -> None:
        if self.individuals_file is not None:
            self.individuals = self.read_individuals(self.individuals_file)

    @staticmethod
    def read_individuals(individuals_file: Path) -> dict[str, str]:
        """{sample: individual} from a two-column tab-separated file (sample, individual).

        Empty lines and lines starting with # are ignored; a header line is
        harmless (its "sample" matches no sample).
        """
        individuals: dict[str, str] = {}
        with individuals_file.open(encoding="utf-8") as lines:
            for number, line in enumerate(lines, start=1):
                fields = [value.strip() for value in line.rstrip("\r\n").split("\t")]
                if not fields[0] or fields[0].startswith("#"):
                    continue
                if len(fields) < 2 or not fields[1]:
                    logging.error(
                        "Line %d of %s: expected a sample and its individual separated by a tab.",
                        number,
                        individuals_file,
                    )
                    sys.exit(1)
                sample, individual = fields[0], fields[1]
                if individuals.get(sample, individual) != individual:
                    logging.error(
                        "Sample %s has two individuals in %s: %s and %s.",
                        sample,
                        individuals_file,
                        individuals[sample],
                        individual,
                    )
                    sys.exit(1)
                individuals[sample] = individual
        if not individuals:
            logging.error("No sample and individual found in %s.", individuals_file)
            sys.exit(1)
        return individuals

    def relation(self, sample1: pd.Series, sample2: pd.Series) -> np.ndarray:
        """same_individual or different_individuals; a sample missing from the
        individual file is an individual of its own"""
        individual1 = sample1.map(self.individuals).fillna(sample1).to_numpy()
        individual2 = sample2.map(self.individuals).fillna(sample2).to_numpy()
        return np.where(individual1 == individual2, self.SAME, self.DIFFERENT).astype(object)

    @staticmethod
    def rate(shared: int, compared: int) -> float:
        """shared / compared, NaN when nothing is compared"""
        return shared / compared if compared else np.nan

    def execute(self, comparison_files: dict[str, Path]) -> None:
        """Write the events, sample pair and MSP tables.

        :param comparison_files: {msp name: <msp>_comparison.tsv}
        """
        if self.individuals_file is None:
            logging.info(
                "No individual file (--ind): every sample is a distinct individual."
            )
        codes: dict[str, int] = {}
        pair_counts: pd.DataFrame | None = None
        pending: list[pd.DataFrame] = []
        pending_rows = 0
        msp_rows = []
        events_file = self.output_dir / self.EVENTS_FILE
        with events_file.open("w", encoding="utf-8") as events:
            events.write(
                "msp_name\tsample1\tsample2\trelation\tdistance\toverlap_noN_info_pc\tthreshold\n"
            )
            for msp_name, comparison_file in sorted(comparison_files.items()):
                table = pd.read_csv(
                    comparison_file,
                    sep="\t",
                    usecols=["sample1", "sample2", "overlap_noN_info_pc", "distance"],
                    dtype={"sample1": str, "sample2": str},
                )
                table = table.loc[
                    (table["overlap_noN_info_pc"] >= self.min_overlap)
                    & table["distance"].notna()
                ]
                relation = self.relation(table["sample1"], table["sample2"])
                distance = table["distance"].to_numpy(dtype=float)
                shared = distance <= self.THRESHOLD
                different = relation == self.DIFFERENT
                msp_rows.append(
                    {
                        "msp_name": msp_name,
                        "threshold": self.THRESHOLD,
                        "n_sample_pairs": int(different.sum()),
                        "n_shared_pairs": int((shared & different).sum()),
                        "n_same_individual_pairs": int((~different).sum()),
                        "n_same_individual_shared": int((shared & ~different).sum()),
                    }
                )
                if not len(table):
                    continue
                # events
                event_table = table.loc[shared].assign(
                    relation=relation[shared], threshold=self.THRESHOLD
                )
                event_table.insert(0, "msp_name", msp_name)
                event_table[
                    ["msp_name", "sample1", "sample2", "relation", "distance",
                     "overlap_noN_info_pc", "threshold"]
                ].to_csv(events, sep="\t", index=False, header=False)
                # sample pairs, as integer keys (smaller code first)
                for name in pd.unique(pd.concat([table["sample1"], table["sample2"]])):
                    codes.setdefault(name, len(codes))
                code1 = table["sample1"].map(codes).to_numpy(dtype=np.int64)
                code2 = table["sample2"].map(codes).to_numpy(dtype=np.int64)
                keys = np.minimum(code1, code2) << 32 | np.maximum(code1, code2)
                pending.append(pd.DataFrame({"key": keys, "compared": 1, "shared": shared.astype(np.int64)}))
                pending_rows += len(keys)
                if pending_rows >= self.REDUCE_PAIRS:
                    pair_counts = self._reduce(pair_counts, pending)
                    pending, pending_rows = [], 0
        pair_counts = self._reduce(pair_counts, pending)
        self._write_samples(pair_counts, codes)
        self._write_msp(msp_rows)
        logging.info(
            "Strain sharing of %d MSPs saved in %s, %s and %s",
            len(msp_rows),
            self.EVENTS_FILE,
            self.SAMPLES_FILE,
            self.MSP_FILE,
        )

    @staticmethod
    def _reduce(pair_counts: pd.DataFrame | None, pending: list[pd.DataFrame]) -> pd.DataFrame | None:
        """Sum the pending pair counts into pair_counts (indexed by pair key)"""
        if not pending:
            return pair_counts
        counts = pd.concat(pending).groupby("key")[["compared", "shared"]].sum()
        if pair_counts is None:
            return counts
        return pair_counts.add(counts, fill_value=0).astype(np.int64)

    def _write_samples(self, pair_counts: pd.DataFrame | None, codes: dict[str, int]) -> None:
        columns = ["sample1", "sample2", "relation", "sharing_rate", "n_msp", "n_msp_shared"]
        if pair_counts is None:
            pd.DataFrame(columns=columns).to_csv(
                self.output_dir / self.SAMPLES_FILE, sep="\t", index=False
            )
            return
        names = np.empty(len(codes), dtype=object)
        names[list(codes.values())] = list(codes.keys())
        keys = pair_counts.index.to_numpy(dtype=np.int64)
        first, second = names[keys >> 32], names[keys & 0xFFFFFFFF]
        # sample1 before sample2 in alphabetical order
        swap = first > second
        sample1 = pd.Series(np.where(swap, second, first))
        sample2 = pd.Series(np.where(swap, first, second))
        samples = pd.DataFrame(
            {
                "sample1": sample1,
                "sample2": sample2,
                "relation": self.relation(sample1, sample2),
                "n_msp": pair_counts["compared"].to_numpy(),
                "n_msp_shared": pair_counts["shared"].to_numpy(),
            }
        )
        samples["sharing_rate"] = samples["n_msp_shared"] / samples["n_msp"]
        samples.sort_values(["sample1", "sample2"]).to_csv(
            self.output_dir / self.SAMPLES_FILE, sep="\t", index=False, columns=columns
        )

    def _write_msp(self, msp_rows: list[dict]) -> None:
        msp = pd.DataFrame(
            msp_rows,
            columns=[
                "msp_name", "threshold", "n_sample_pairs", "n_shared_pairs",
                "n_same_individual_pairs", "n_same_individual_shared",
            ],
        )
        msp["transmissibility"] = [
            self.rate(s, c) for s, c in zip(msp["n_shared_pairs"], msp["n_sample_pairs"])
        ]
        msp.to_csv(
            self.output_dir / self.MSP_FILE,
            sep="\t",
            index=False,
            columns=[
                "msp_name", "threshold", "transmissibility", "n_sample_pairs", "n_shared_pairs",
                "n_same_individual_pairs", "n_same_individual_shared",
            ],
        )
