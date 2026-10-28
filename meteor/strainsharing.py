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
    """Strain-sharing events, sharing rate of each sample pair and of each MSP.

    Two samples are compared on an MSP when the sites informative in both
    cover at least min_overlap % of the alignment. They share its strain when
    their distance is at most the threshold of the MSP: the
    THRESHOLD_PERCENTILE-th percentile of the distances between samples of
    different groups (unrelated samples) when a group file gives enough of
    them, as StrainPhlAn's strain identity thresholds (Valles-Colomer et al.,
    Nature 2023), otherwise FIXED_THRESHOLD.
    """

    DEFAULT_MIN_OVERLAP: ClassVar[float] = 10.0
    DEFAULT_MIN_PAIRS: ClassVar[int] = 10
    THRESHOLD_PERCENTILE: ClassVar[float] = 3.0
    # unrelated pairs needed to set the threshold of an MSP from their distances
    MIN_CALIBRATION_PAIRS: ClassVar[int] = 50
    # same_strain limit of the comparison tables
    FIXED_THRESHOLD: ClassVar[float] = 0.0001
    # partial counts are summed every REDUCE_PAIRS pairs (bounded memory)
    REDUCE_PAIRS: ClassVar[int] = 5_000_000
    EVENTS_FILE: ClassVar[str] = "strain_sharing_events.tsv"
    SAMPLES_FILE: ClassVar[str] = "strain_sharing_samples.tsv"
    MSP_FILE: ClassVar[str] = "strain_sharing_msp.tsv"

    output_dir: Path
    groups_file: Path | None = None
    min_overlap: float = DEFAULT_MIN_OVERLAP
    min_pairs: int = DEFAULT_MIN_PAIRS
    groups: dict[str, str] = field(default_factory=dict, init=False)

    def __post_init__(self) -> None:
        if self.groups_file is not None:
            self.groups = self.read_groups(self.groups_file)

    @staticmethod
    def read_groups(groups_file: Path) -> dict[str, str]:
        """{sample: group} from a two-column tab-separated file (sample, group).

        Empty lines and lines starting with # are ignored; a header line is
        harmless (its "sample" matches no sample).
        """
        groups: dict[str, str] = {}
        with groups_file.open(encoding="utf-8") as lines:
            for number, line in enumerate(lines, start=1):
                fields = [value.strip() for value in line.rstrip("\r\n").split("\t")]
                if not fields[0] or fields[0].startswith("#"):
                    continue
                if len(fields) < 2 or not fields[1]:
                    logging.error(
                        "Line %d of %s: expected a sample and its group separated by a tab.",
                        number,
                        groups_file,
                    )
                    sys.exit(1)
                sample, group = fields[0], fields[1]
                if groups.get(sample, group) != group:
                    logging.error(
                        "Sample %s has two groups in %s: %s and %s.",
                        sample,
                        groups_file,
                        groups[sample],
                        group,
                    )
                    sys.exit(1)
                groups[sample] = group
        if not groups:
            logging.error("No sample and group found in %s.", groups_file)
            sys.exit(1)
        return groups

    def relation(self, sample1: pd.Series, sample2: pd.Series) -> np.ndarray:
        """related (same group), unrelated (different groups) or unknown"""
        if not self.groups:
            return np.full(len(sample1), "unknown", dtype=object)
        group1 = sample1.map(self.groups)
        group2 = sample2.map(self.groups)
        known = group1.notna().to_numpy() & group2.notna().to_numpy()
        same = (group1.to_numpy() == group2.to_numpy()) & known
        return np.where(same, "related", np.where(known, "unrelated", "unknown")).astype(object)

    def threshold(self, distance: np.ndarray, relation: np.ndarray) -> tuple[float, str]:
        """Strain identity threshold of an MSP and how it was set"""
        unrelated = distance[relation == "unrelated"]
        if len(unrelated) >= self.MIN_CALIBRATION_PAIRS:
            return (
                float(np.percentile(unrelated, self.THRESHOLD_PERCENTILE)),
                f"unrelated_p{self.THRESHOLD_PERCENTILE:g}",
            )
        return self.FIXED_THRESHOLD, "fixed"

    def rate(self, shared: int, compared: int) -> float:
        """Sharing rate of an MSP, NaN under min_pairs compared pairs"""
        return shared / compared if compared >= self.min_pairs else np.nan

    def execute(self, comparison_files: dict[str, Path]) -> None:
        """Write the events, sample pair and MSP tables.

        :param comparison_files: {msp name: <msp>_comparison.tsv}
        """
        if self.groups_file is None:
            logging.warning(
                "No group file (--groups): strain sharing uses the fixed threshold "
                "%g and no transmissibility is computed.",
                self.FIXED_THRESHOLD,
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
                threshold, source = self.threshold(distance, relation)
                shared = distance <= threshold
                row: dict[str, object] = {
                    "msp_name": msp_name,
                    "threshold": threshold,
                    "threshold_source": source,
                    "pairs_compared": len(table),
                    "pairs_shared": int(shared.sum()),
                }
                for name in ("related", "unrelated"):
                    in_relation = relation == name
                    row[f"{name}_compared"] = int(in_relation.sum())
                    row[f"{name}_shared"] = int((shared & in_relation).sum())
                msp_rows.append(row)
                if not len(table):
                    continue
                # events
                event_table = table.loc[shared].assign(relation=relation[shared], threshold=threshold)
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
        columns = ["sample1", "sample2", "relation", "msp_compared", "msp_shared", "sharing_rate"]
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
                "msp_compared": pair_counts["compared"].to_numpy(),
                "msp_shared": pair_counts["shared"].to_numpy(),
            }
        )
        samples["sharing_rate"] = samples["msp_shared"] / samples["msp_compared"]
        samples.sort_values(["sample1", "sample2"]).to_csv(
            self.output_dir / self.SAMPLES_FILE, sep="\t", index=False, columns=columns
        )

    def _write_msp(self, msp_rows: list[dict]) -> None:
        msp = pd.DataFrame(
            msp_rows,
            columns=[
                "msp_name", "threshold", "threshold_source", "pairs_compared", "pairs_shared",
                "related_compared", "related_shared", "unrelated_compared", "unrelated_shared",
            ],
        )
        msp["sharing_rate"] = [
            self.rate(s, c) for s, c in zip(msp["pairs_shared"], msp["pairs_compared"])
        ]
        msp["transmissibility"] = [
            self.rate(s, c) for s, c in zip(msp["related_shared"], msp["related_compared"])
        ]
        msp["unrelated_sharing_rate"] = [
            self.rate(s, c) for s, c in zip(msp["unrelated_shared"], msp["unrelated_compared"])
        ]
        msp.to_csv(
            self.output_dir / self.MSP_FILE,
            sep="\t",
            index=False,
            columns=[
                "msp_name", "threshold", "threshold_source",
                "pairs_compared", "pairs_shared", "sharing_rate",
                "related_compared", "related_shared", "transmissibility",
                "unrelated_compared", "unrelated_shared", "unrelated_sharing_rate",
            ],
        )
