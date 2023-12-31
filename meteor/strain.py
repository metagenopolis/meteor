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
import pandas as pd
from dataclasses import dataclass, field
from configparser import ConfigParser
from meteor.variantcalling import VariantCalling
from meteor.session import Session, Component
from typing import Type
from tempfile import mkdtemp
from pathlib import Path
from time import perf_counter
from io import StringIO
from pysam import coverage


@dataclass
class Strain(Session):
    """Counter session map and count"""

    meteor: Type[Component]
    max_depth: int
    # min_gene_count: int
    min_snp_depth: int
    min_frequency_non_reference: float
    min_msp_coverage: int
    min_gene_coverage: float
    keep_consensus: bool
    ini_data: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.meteor.tmp_dir = Path(mkdtemp(dir=self.meteor.tmp_path))
        self.meteor.strain_dir.mkdir(exist_ok=True, parents=True)

    def filter_coverage(self, bam_file: Path, bed_file: Path) -> pd.DataFrame:
        """Filter gene coverage
        :param bam_file: (Path) Path to the bam file
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
            header=1,
        )
        cov_df = pd.read_csv(
            StringIO(coverage("-d", str(self.max_depth), str(bam_file.resolve()))),
            sep="\t",
            header=1,
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
        #     a = bamdesc.count_coverage(str(gene), start=startpos, stop=endpos)
        #     print(a)

    def get_msp_variant(
        self,
        consensus_file: Path,
        # count_file: Path,
        msp_file: Path,
        bam_file: Path,
        bed_file: Path,
    ) -> None:
        """Produce meaning full msp variants
        :param consensus_file: (Path) A path to consensus file
        :param count_file: (Path) A path to count file
        :param msp_file: (Path) A path to msp file
        :param bam_file: (Path) A path to bam file
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
            names=["msp_name", "gene_id", "gene_name", "gene_category"],
            header=1,
        )
        msp_content = msp_content.loc[msp_content["gene_category"] == "core"]
        msp_content = msp_content.groupby("msp_name").head(100).reset_index()
        # Filter for gene with a minimum count
        # filtered_count = count[count["value"] >= self.min_gene_count]
        # Filter for gene with a minimum count
        if self.min_gene_coverage:
            filtered_coverage = self.filter_coverage(bam_file, bed_file)
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
        msp_covered = joined_df.merge(msp_with_overlapping_genes, on="msp_name")
        # clear first ?
        if not consensus_file.exists():
            logging.error(
                "Consensus file %s is not available to perform a new strain extraction. Please consider to use strain with --kc option.",
                consensus_file,
            )
            sys.exit()
        gene_dict = {
            gene_id: seq for gene_id, seq in self.get_sequences(consensus_file)
        }
        logging.info(
            "%s MSPs have sufficient signal for SNP analysis ",
            len(msp_with_overlapping_genes["msp_name"].values),
        )
        for msp_name in msp_with_overlapping_genes["msp_name"].values:
            msp_file = self.ini_data["directory"] / Path(msp_name + ".fasta")
            msp_seq = ""
            for gene_id in msp_content[msp_content["msp_name"] == msp_name][
                "gene_id"
            ].values:
                if gene_id in msp_covered["gene_id"].values:
                    msp_seq += gene_dict[gene_id]
                else:
                    msp_seq += "-" * len(gene_dict[gene_id])
            with msp_file.open("wt", encoding="UTF-8") as msp:
                print(
                    f">{self.ini_data['census']['sample_info']['sample_name']}\n{msp_seq}\n",
                    file=msp,
                )
        # Remove consensus file
        if not self.keep_consensus:
            logging.info(
                "Consensus catalogue is not kept. Re-counting strain will need to be performed from scratch."
            )
            consensus_file.unlink(missing_ok=True)
        # for gene_id, seq in self.get_gene_consensus(consensus_file):
        #     if gene_id in msp_covered["gene_id"].values:
        #         msp_name = msp_covered[msp_covered["gene_id"] == gene_id][
        #             "msp_name"
        #         ].values[0]
        #         msp_file = self.ini_data["directory"] / Path(msp_name + ".fasta")
        #         with msp_file.open("at") as msp:
        #             print(f">{gene_id}\n{fill(seq, width=80)}\n", file=msp)
        #     print(gene_id)
        #     print(seq)

    def execute(self) -> bool:
        """Compute the strain analysis"""
        logging.info("Launch strain analysis")
        try:
            # Get the ini ref
            ref_ini_file_list = list(self.meteor.ref_dir.glob("**/*_reference.ini"))
            assert len(ref_ini_file_list) == 1
            ref_ini_file = ref_ini_file_list[0]
            ref_ini = ConfigParser()
            with open(ref_ini_file, "rt", encoding="UTF-8") as ref:
                ref_ini.read_file(ref)
            self.meteor.ref_name = ref_ini["reference_info"]["reference_name"]
        except AssertionError:
            logging.error(
                "Error, no *_reference.ini file found in %s. "
                "One *_reference.ini is expected",
                self.meteor.ref_dir,
            )
            sys.exit()
        try:
            census_ini_file_list = list(
                self.meteor.mapped_sample_dir.glob("*_census_stage_1.ini")
            )
            assert len(census_ini_file_list) == 1
            census_ini_file = census_ini_file_list[0]
            census_ini = ConfigParser()
            with open(census_ini_file, "rt", encoding="UTF-8") as cens:
                census_ini.read_file(cens)
                sample_info = census_ini["sample_info"]
                stage3_dir = self.meteor.strain_dir / sample_info["sample_name"]
                stage3_dir.mkdir(exist_ok=True, parents=True)
                self.ini_data["mapped_sample_dir"] = self.meteor.mapped_sample_dir
                self.ini_data["directory"] = stage3_dir
                self.ini_data["census"] = census_ini
                self.ini_data["reference"] = ref_ini
                self.ini_data["reference"] = ref_ini
                self.ini_data["Stage3FileName"] = (
                    stage3_dir / f"{sample_info['sample_name']}_census_stage_3.ini"
                )

                # Get the bam
                # Get the variant calling
                # Variant calling this library on the reference
                variant_calling_process = VariantCalling(
                    self.meteor,
                    self.ini_data,
                    self.max_depth,
                    self.min_snp_depth,
                    self.min_frequency_non_reference,
                )
            if self.ini_data["Stage3FileName"].exists():
                logging.info(
                    "Variant calling already done for sample: %s",
                    sample_info["sample_name"],
                )
                logging.info("Skipped !")
            else:
                variant_calling_process.execute()
                ref_ini = ConfigParser()
                with open(ref_ini_file, "rt", encoding="UTF-8") as ref:
                    ref_ini.read_file(ref)
                self.update_ini(
                    ref_ini,
                    "variant_calling",
                    {
                        "min_msp_coverage": str(self.min_msp_coverage),
                        "min_gene_coverage": str(self.min_gene_coverage),
                        # "min_gene_comptage": str(self.min_gene_count),
                    },
                )
            consensus_file = (
                self.ini_data["directory"]
                / f"{sample_info['sample_name']}_consensus.fasta"
            )
            # count_file = (
            #     self.ini_data["mapped_sample_dir"] / f"{sample_info['sample_name']}.tsv"
            # )
            msp_file = (
                self.meteor.ref_dir
                / self.ini_data["reference"]["reference_file"]["database_dir"]
                / self.ini_data["reference"]["annotation"]["msp"]
            )
            bam_file = (
                self.ini_data["mapped_sample_dir"] / f"{sample_info['sample_name']}.bam"
            )
            bed_file = (
                self.meteor.ref_dir
                / self.ini_data["reference"]["reference_file"]["database_dir"]
                / self.ini_data["reference"]["annotation"]["bed"]
            )
            start = perf_counter()
            # count_file,
            self.get_msp_variant(consensus_file, msp_file, bam_file, bed_file)
            logging.info(
                "Completed strain analysis in %f seconds", perf_counter() - start
            )
        except AssertionError:
            logging.error(
                "Error, no *_census_stage_0.ini file found in %s", self.meteor.fastq_dir
            )
            sys.exit()
        return True
