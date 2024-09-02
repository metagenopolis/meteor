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

"""Run variant calling and performs counting"""

import sys
import logging
import lzma
import pandas as pd
from dataclasses import dataclass, field
from meteor.variantcalling import VariantCalling
from meteor.session import Session, Component
from tempfile import mkdtemp
from pathlib import Path
from time import perf_counter
from io import StringIO
import pysam
from typing import ClassVar
from shutil import rmtree


@dataclass
class Strain(Session):
    """Counter session map and count"""

    DEFAULT_MAX_DEPTH: ClassVar[int] = 100
    MIN_MIN_SNP_DEPTH: ClassVar[int] = 1
    MAX_MIN_SNP_DEPTH: ClassVar[int] = 10000
    DEFAULT_MIN_SNP_DEPTH: ClassVar[int] = 3
    DEFAULT_MIN_FREQUENCY: ClassVar[float] = 0.01
    DEFAULT_PLOIDY: ClassVar[int] = 1
    MIN_MIN_MSP_COVERAGE: ClassVar[int] = 1
    MAX_MIN_MSP_COVERAGE: ClassVar[int] = 100
    DEFAULT_MIN_MSP_COVERAGE: ClassVar[int] = 50
    DEFAULT_MIN_GENE_COVERAGE: ClassVar[float] = 0.5
    DEFAULT_NUM_THREADS: ClassVar[int] = 1
    DEFAULT_MIN_DEPTH: ClassVar[int] = 3
    MIN_DEPTH: ClassVar[int] = 1

    meteor: type[Component]
    max_depth: int
    min_depth: int
    min_snp_depth: int
    min_frequency: float
    ploidy: int
    min_msp_coverage: int
    min_gene_coverage: float
    core_size: int
    keep_consensus: bool
    json_data: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.meteor.tmp_dir = Path(mkdtemp(dir=self.meteor.tmp_path))
        self.meteor.strain_dir.mkdir(exist_ok=True, parents=True)

    def filter_coverage(
        self, cram_file: Path, bed_file: Path, reference_file: Path
    ) -> pd.DataFrame:
        """Filter gene coverage
        :param cram_file: (Path) Path to the cram file
        :param bed_file: (Path) Path to the bed file
        :return: (pd.DataFrame) Return a matrix with gene having a coverage above threshold
        """
        # Read the abundance file
        gene_interest = pd.read_csv(
            bed_file,
            sep="\t",
            names=[
                "gene_id",
                "startpos",
                "endpos",
            ],
            # header=1,
        )
        cov_df = pd.read_csv(
            StringIO(
                pysam.coverage(
                    "--reference",
                    str(reference_file.resolve()),
                    "-d",
                    str(self.max_depth),
                    str(cram_file.resolve()),
                )  # type: ignore[attr-defined]
            ),
            sep="\t",
            header=0,
            names=[
                "gene_id",
                "startpos",
                "endpos",
                "numreads",
                "covbases",
                "coverage",
                "meandepth",
                "meanbaseq",
                "meanmapq",
            ],
        )
        merged_df = cov_df.merge(gene_interest, on="gene_id", how="inner")
        filtered_df = merged_df[
            merged_df["coverage"] >= round(self.min_gene_coverage * 100.0)
        ]
        return filtered_df

        # for gene, startpos, endpos in gene_interest.itertuples(index=False):
        #     a = cramdesc.count_coverage(str(gene), start=startpos, stop=endpos)
        #     print(a)

    def is_only_question_marks(self, s: str) -> bool:
        """Check if a string is made of question marks
        :param s: (str) A consensus string of MSPs
        :return: (bool) Return True if the string is made of question marks only
        """
        return set(s) == {self.meteor.DEFAULT_GAP_CHAR}

    def get_msp_variant(
        self,
        consensus_file: Path,
        # count_file: Path,
        msp_file: Path,
        cram_file: Path,
        bed_file: Path,
        reference_file: Path,
    ) -> None:
        """Produce meaning full msp variants
        :param consensus_file: (Path) A path to consensus file
        :param count_file: (Path) A path to count file
        :param msp_file: (Path) A path to msp file
        :param cram_file: (Path) A path to cram file
        :param bed_file: (Path) A path to bed file
        """
        # Read the abundance file
        # count = pd.read_csv(
        #     count_file, sep="\t", names=["gene_id", "gene_length", "value"], header=1
        # )
        # Read the msp content description file
        msp_content = pd.read_csv(
            msp_file,
            sep="\t",
            # names=["msp_name", "gene_id", "gene_name", "gene_category"],
            header=0,
        )
        msp_content = msp_content.loc[msp_content["gene_category"] == "core"]
        msp_content = msp_content.groupby("msp_name").head(self.core_size).reset_index()
        # Filter for gene with a minimum count
        # filtered_count = count[count["value"] >= self.min_gene_count]
        # Filter for gene with a minimum count
        if self.min_gene_coverage:
            filtered_coverage = self.filter_coverage(
                cram_file, bed_file, reference_file
            )
            joined_df = msp_content.merge(filtered_coverage, on="gene_id")
            # filtered_coverage = filtered_count.groupby("gene_id")
            # Join the two DataFrames based on gene_id
            # joined_df = msp_content.merge(filtered_count, on="gene_id")
            # Count the number of overlapping gene IDs for each MSP
            overlapping_gene_counts = (
                joined_df.groupby("msp_name")["gene_id"]
                .nunique()
                .to_frame(name="overlapping_gene_count")
                .reset_index()
            )
        else:
            overlapping_gene_counts = (
                msp_content.groupby("msp_name")["gene_id"]
                .nunique()
                .to_frame(name="overlapping_gene_count")
                .reset_index()
            )
        # Filter the DataFrame to only include MSPs with at least min_msp_coverage overlapping gene IDs
        msp_with_overlapping_genes = overlapping_gene_counts[
            overlapping_gene_counts["overlapping_gene_count"] >= self.min_msp_coverage
        ]
        # msp_covered = joined_df.merge(msp_with_overlapping_genes, on="msp_name")
        # clear first ?
        if not consensus_file.exists():
            logging.error(
                """Consensus file %s is not available to perform a new strain extraction.
                Please consider to use strain with --kc option.""",
                consensus_file,
            )
            sys.exit(1)

        gene_dict = dict(
            self.get_sequences(consensus_file, use_lzma=True, id_as_int=True)
        )
        logging.info(
            "%s MSPs have sufficient signal for SNP analysis ",
            len(msp_with_overlapping_genes["msp_name"].values),
        )
        for msp_name in msp_with_overlapping_genes["msp_name"].values:
            msp_file = self.json_data["directory"] / Path(msp_name + ".fasta.xz")
            msp_seq = ""
            for gene_id in msp_content[msp_content["msp_name"] == msp_name][
                "gene_id"
            ].values:
                # if gene_id in msp_covered["gene_id"].values:
                msp_seq += gene_dict[gene_id]
                # else:
                #     msp_seq += "?" * len(gene_dict[gene_id])

            if not self.is_only_question_marks(msp_seq):
                with lzma.open(msp_file, "wt", preset=0) as msp:
                    print(
                        f">{self.json_data['census']['sample_info']['sample_name']}\n{msp_seq}\n",
                        file=msp,
                    )
            else:
                logging.info(
                    "%s contains only question mark characters according to depth thresholds. Skipping...",
                    msp_name,
                )
        # Remove consensus file
        if not self.keep_consensus:
            logging.info(
                "Consensus catalogue is not kept (--kc). Re-counting strain will need to be performed from scratch."
            )
            consensus_file.unlink(missing_ok=True)
        # for gene_id, seq in self.get_gene_consensus(consensus_file):
        #     if gene_id in msp_covered["gene_id"].values:
        #         msp_name = msp_covered[msp_covered["gene_id"] == gene_id][
        #             "msp_name"
        #         ].values[0]
        #         msp_file = self.json_data["directory"] / Path(msp_name + ".fasta")
        #         with msp_file.open("at") as msp:
        #             print(f">{gene_id}\n{fill(seq, width=80)}\n", file=msp)
        #     print(gene_id)
        #     print(seq)

    def execute(self) -> None:
        """Compute the strain analysis"""
        logging.info("Launch strain analysis")
        try:
            ref_json = self.get_reference_info(self.meteor.ref_dir)
            self.meteor.ref_name = ref_json["reference_info"]["reference_name"]
        except AssertionError:
            logging.error(
                "No *_reference.json file found in %s. "
                "One *_reference.json is expected",
                self.meteor.ref_dir,
            )
            sys.exit(1)
        try:
            census_json = self.get_census_stage(self.meteor.mapped_sample_dir, 1)
            sample_info = census_json["sample_info"]
            stage3_dir = self.meteor.strain_dir / sample_info["sample_name"]
            stage3_dir.mkdir(exist_ok=True, parents=True)
            self.json_data["mapped_sample_dir"] = self.meteor.mapped_sample_dir
            self.json_data["directory"] = stage3_dir
            self.json_data["census"] = census_json
            self.json_data["reference"] = ref_json
            self.json_data["Stage3FileName"] = (
                stage3_dir / f"{sample_info['sample_name']}_census_stage_3.json"
            )

            # Get the cram
            # Get the variant calling
            # Variant calling this library on the reference
            variant_calling_process = VariantCalling(
                self.meteor,
                self.json_data,
                self.max_depth,
                self.min_depth,
                self.min_snp_depth,
                self.min_frequency,
                self.ploidy,
                self.core_size,
            )
            if self.json_data["Stage3FileName"].exists():
                logging.info(
                    "Variant calling already done for sample: %s",
                    sample_info["sample_name"],
                )
                logging.info("Skipped !")
            else:
                variant_calling_process.execute()
                self.update_json(
                    ref_json,
                    "variant_calling",
                    {
                        "min_msp_coverage": str(self.min_msp_coverage),
                        "min_gene_coverage": str(self.min_gene_coverage),
                        # "min_gene_comptage": str(self.min_gene_count),
                    },
                )
            consensus_file = (
                self.json_data["directory"]
                / f"{sample_info['sample_name']}_consensus.fasta.xz"
            )
            # count_file = (
            #     self.json_data["mapped_sample_dir"] / f"{sample_info['sample_name']}.tsv"
            # )
            msp_file = (
                self.meteor.ref_dir
                / self.json_data["reference"]["reference_file"]["database_dir"]
                / self.json_data["reference"]["annotation"]["msp"]["filename"]
            )
            cram_file = (
                self.json_data["mapped_sample_dir"]
                / f"{sample_info['sample_name']}.cram"
            )
            reference_file = (
                self.meteor.ref_dir
                / self.json_data["reference"]["reference_file"]["fasta_dir"]
                / self.json_data["reference"]["reference_file"]["fasta_filename"]
            )
            bed_file = (
                self.meteor.ref_dir
                / self.json_data["reference"]["reference_file"]["database_dir"]
                / self.json_data["reference"]["annotation"]["bed"]["filename"]
            )
            start = perf_counter()
            # count_file,
            self.get_msp_variant(
                consensus_file, msp_file, cram_file, bed_file, reference_file
            )
            logging.info(
                "Completed strain analysis in %f seconds", perf_counter() - start
            )
        except AssertionError:
            logging.error(
                "No *_census_stage_0.json file found in %s",
                self.meteor.fastq_dir,
            )
            sys.exit(1)
        else:
            rmtree(self.meteor.tmp_dir, ignore_errors=True)
