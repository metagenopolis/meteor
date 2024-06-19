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

"""Effective variant calling"""

import logging
import sys
import tempfile
import lzma
from subprocess import check_call, run, Popen, PIPE
from dataclasses import dataclass
from pathlib import Path
from datetime import datetime
from meteor.session import Session, Component
from time import perf_counter
from tempfile import NamedTemporaryFile
from packaging.version import parse

# import memory_profiler
# import psutil
from pysam import AlignmentFile, FastaFile
from collections import defaultdict
import pandas as pd
from io import StringIO


@dataclass
class VariantCalling(Session):
    """Run bcftools"""

    meteor: type[Component]
    census: dict
    max_depth: int
    min_depth: int
    min_snp_depth: int
    min_frequency_non_reference: float

    def set_variantcalling_config(
        self,
        cram_file: Path,
        vcf_file: Path,
        consensus_file: Path,
        bcftool_version: str,
    ) -> dict:  # pragma: no cover
        """Define the census 1 configuration

        :param cmd: A string of the specific parameters
        :param cram_file: A path to the sam file
        :return: (Dict) A dict object with the census 1 config
        """
        config = {
            "meteor_version": self.meteor.version,
            "sample_info": self.census["census"]["sample_info"],
            "sample_file": self.census["census"]["sample_file"],
            "mapping": {
                "reference_name": self.census["reference"]["reference_info"][
                    "reference_name"
                ],
                "cram_name": cram_file.name,
            },
            "variant_calling": {
                "variant_calling_tool": "bcftools",
                "variant_calling_version": bcftool_version,
                "variant_calling_date": datetime.now().strftime("%Y-%m-%d"),
                "vcf_name": vcf_file.name,
                "consensus_name": consensus_file.name,
                "min_snp_depth": str(self.min_snp_depth),
                "min_frequency_non_reference": str(self.min_frequency_non_reference),
                "max_depth": str(self.max_depth),
            },
        }
        return config

    def group_consecutive_positions(
        self, position_count_dict: dict, gene_name: str, gene_length: int
    ):
        # Initialize the result list
        result = []

        # Initialize variables for tracking ranges
        start = 0
        current_count = position_count_dict.get(0, 0)

        # Create a sorted list of all positions from 0 to gene_length
        for pos in range(1, gene_length + 1):
            # Get the count for the current position, defaulting to 0 if not present in the dictionary
            count = position_count_dict.get(pos, 0)
            if count != current_count:
                # Append the current range to the result
                result.append((start, pos, current_count))
                # Start a new range
                start = pos
                current_count = count

        # Append the last range
        if start != pos:
            result.append((start, pos, current_count))
        return pd.DataFrame(
            [
                (gene_name, start, end, count)
                for start, end, count in result
                if count < self.min_depth
            ],
            columns=["gene_id", "startpos", "endpos", "coverage"],
        )

    def count_reads_in_gene(
        self,
        cram: AlignmentFile,
        gene_name: str,
        gene_length: int,
        reference_file: Path,
    ):
        """
        Counts the number of reads at each position in a specified gene from a BAM file,
        ignoring reads with gaps or deletions at that position.

        Parameters:
        cram (pysam.AlignmentFile): A pysam AlignmentFile object representing the BAM file
        gene_name (str): The name of the gene for which reads are to be counted

        Returns:
        dict: A dictionary with positions (int) as keys and the number of reads (int) at each position as values.
        """
        reads_dict: dict[int, int] = defaultdict(int)
        # gene_length = bam.lengths[bam.references.index(gene_name)]
        # Extraction of positions having at least one read
        with FastaFile(filename=str(reference_file.resolve())) as Fasta:
            for pileupcolumn in cram.pileup(
                contig=gene_name,
                start=0,
                end=gene_length,
                stepper="all",
                # stepper="samtools",
                max_depth=self.max_depth,
                fastafile=Fasta,
                multiple_iterators=False,
            ):
                read_count = sum(
                    1
                    for pileupread in pileupcolumn.pileups
                    if not pileupread.is_del and not pileupread.is_refskip
                )
                reads_dict[pileupcolumn.reference_pos] = read_count

        return reads_dict

    # @memory_profiler.profile
    def filter_low_cov_sites_python(
        self,
        cram_file: Path,
        reference_file: Path,
        temp_low_cov_sites: tempfile._TemporaryFileWrapper,
    ) -> None:
        """Create a bed file reporting a list of positions below coverage threshold

        :param cram_file:   Path to the input cram file
        :param temp_low_cov_sites: File handle to write low coverage sites
        """
        # Get list of genes
        gene_interest = pd.read_csv(
            self.matrix_file,
            sep="\t",
            names=[
                "gene_id",
                "gene_length",
                "coverage",
            ],
            header=0,
            compression="xz",
        )
        # Round the coverage column to 0 decimal places
        gene_interest["coverage"] = gene_interest["coverage"].round(0)

        # Optionally, if you want to convert the rounded values to integers
        # gene_interest["coverage"] = gene_interest["coverage"].astype(int)
        # Convert the 'coverage' column to a sparse column using SparseDtype
        gene_interest["coverage"] = gene_interest["coverage"].astype(
            pd.SparseDtype(int, 0)
        )
        # We need to do more work for these genes
        # their total count is above the min_depth
        # but we need to find on which positions they do.
        gene_tofilter = gene_interest[gene_interest["coverage"] >= self.min_depth][
            ["gene_id", "gene_length"]
        ]
        dfs = []
        with AlignmentFile(str(cram_file.resolve()), "rc") as cram:
            # For genes with a count
            for _, row in gene_tofilter.iterrows():
                reads_dict = self.count_reads_in_gene(
                    cram, str(row["gene_id"]), row["gene_length"], reference_file
                )
                df = self.group_consecutive_positions(
                    reads_dict, str(row["gene_id"]), row["gene_length"]
                )
                if len(df) > 0:
                    dfs.append(df)
        # All these genes are going to be replaced by gaps
        # Their count is below the threshold level
        gene_ignore = gene_interest[gene_interest["coverage"] < self.min_depth]
        dfs.append(
            pd.DataFrame(
                {
                    "gene_id": gene_ignore["gene_id"],
                    "startpos": 0,
                    "endpos": gene_ignore["gene_length"],
                    "coverage": gene_ignore["coverage"],
                }
            )
        )
        sum_cov_bed = pd.concat(dfs, ignore_index=True)
        # .query(f"coverage < {self.min_depth}")
        sum_cov_bed.to_csv(temp_low_cov_sites, sep="\t", header=False, index=False)

    def filter_low_cov_sites(
        self, output: str, temp_low_cov_sites: tempfile._TemporaryFileWrapper
    ) -> None:
        """Write sites with a low coverage"""
        sum_cov_bed = pd.read_csv(
            StringIO(output),
            sep="\t",
            names=[
                "gene_id",
                "startpos",
                "endpos",
                "coverage",
            ],
        )
        sum_cov_bed.query(f"coverage < {self.min_depth}").to_csv(
            temp_low_cov_sites, sep="\t", header=False, index=False
        )

    def execute(self) -> None:
        """Call variants reads"""
        # Start mapping
        cram_file = (
            self.census["mapped_sample_dir"]
            / f"{self.census['census']['sample_info']['sample_name']}.cram"
        )
        self.matrix_file = (
            self.census["mapped_sample_dir"]
            / f"{self.census['census']['sample_info']['sample_name']}.tsv.xz"
        )
        vcf_file = (
            self.census["directory"]
            / f"{self.census['census']['sample_info']['sample_name']}.vcf.gz"
        )
        consensus_file = (
            self.census["directory"]
            / f"{self.census['census']['sample_info']['sample_name']}_consensus.fasta.xz"
        )
        reference_file = (
            self.meteor.ref_dir
            / self.census["reference"]["reference_file"]["fasta_dir"]
            / self.census["reference"]["reference_file"]["fasta_filename"]
        )
        bed_file = (
            self.meteor.ref_dir
            / self.census["reference"]["reference_file"]["database_dir"]
            / self.census["reference"]["annotation"]["bed"]["filename"]
        )
        bcftools_exec = run(["bcftools", "--version"], check=False, capture_output=True)
        if bcftools_exec.returncode != 0:
            logging.error(
                "Checking bcftools failed:\n%s", bcftools_exec.stderr.decode("utf-8")
            )
            sys.exit(1)
        bcftools_version = (
            bcftools_exec.stdout.decode("utf-8").split("\n")[0].split(" ")[1]
        )
        if parse(bcftools_version) < self.meteor.MIN_BCFTOOLS_VERSION:
            logging.error(
                "The bcftools version %s is outdated for meteor. Please update bcftools to >= %s.",
                bcftools_version,
                self.meteor.MIN_BCFTOOLS_VERSION,
            )
            sys.exit(1)
        bedtools_exec = run(["bedtools", "--version"], check=False, capture_output=True)
        if bedtools_exec.returncode != 0:
            logging.error(
                "Checking bedtools failed:\n%s", bedtools_exec.stderr.decode("utf-8")
            )
            sys.exit(1)
        bedtools_version = bedtools_exec.stdout.decode("utf-8").split(" ")[1][1:]
        if parse(bedtools_version) < self.meteor.MIN_BEDTOOLS_VERSION:
            logging.error(
                "The bedtools version %s is outdated for meteor. Please update bedtools to >= %s.",
                bedtools_version,
                self.meteor.MIN_BEDTOOLS_VERSION,
            )
            sys.exit(1)
        start = perf_counter()
        with NamedTemporaryFile(
            mode="wt",
            dir=self.meteor.tmp_dir,  # delete=False
        ) as temp_vcf_pileup:
            with NamedTemporaryFile(
                mode="wt",
                dir=self.meteor.tmp_dir,  # delete=False
            ) as temp_vcf:
                startpileup = perf_counter()
                check_call(
                    [
                        "bcftools",
                        "mpileup",
                        "-d",
                        str(self.max_depth),
                        "-Ob",
                        "-R",
                        str(bed_file.resolve()),
                        "-f",
                        str(reference_file.resolve()),
                        str(cram_file.resolve()),
                        "--threads",
                        str(self.meteor.threads),
                        "-o",
                        temp_vcf_pileup.name,
                    ]
                )
                logging.info(
                    "Completed pileup step in %f seconds", perf_counter() - startpileup
                )
                startcall = perf_counter()
                check_call(
                    [
                        "bcftools",
                        "call",
                        "-v",
                        "-c",
                        "--ploidy",
                        str(1),
                        "-V",
                        "indels",
                        "--threads",
                        str(self.meteor.threads),
                        "-Oz",
                        "-o",
                        temp_vcf.name,
                        temp_vcf_pileup.name,
                    ]
                )
                logging.info(
                    "Completed calling step in %f seconds", perf_counter() - startcall
                )
                check_call(["bcftools", "index", temp_vcf.name])
                # Only SNP from 100 core genes
                logging.info(
                    "Completed index step in %f seconds", perf_counter() - startpileup
                )
                startview = perf_counter()
                check_call(
                    [
                        "bcftools",
                        "view",
                        "-R",
                        str(bed_file.resolve()),
                        "-f",
                        f"MIN(FMT/DP >= {str(self.min_snp_depth)})",
                        "--threads",
                        str(self.meteor.threads),
                        "-q",
                        str(self.min_frequency_non_reference),
                        "-Oz",
                        "-o",
                        str(vcf_file.resolve()),
                        temp_vcf.name,
                    ]
                )
                logging.info(
                    "Completed SNP filtering step in %f seconds",
                    perf_counter() - startview,
                )
                check_call(["bcftools", "index", str(vcf_file.resolve())])
                # The columns of the tab-delimited BED file are also CHROM, POS
                # and END (trailing columns are ignored), but coordinates are
                # 0-based, half-open. To indicate that a file be treated as BED
                # rather than the 1-based tab-delimited file, the file must have
                # the ".bed" or ".bed.gz" suffix (case-insensitive).
                # startlowcovpython = perf_counter()
                # with NamedTemporaryFile(
                #     mode="wt", dir=self.meteor.tmp_dir, delete=False, suffix=".bed"
                # ) as temp_low_cov_sites:
                #     self.filter_low_cov_sites_python(
                #         cram_file, reference_file, temp_low_cov_sites
                #     )
                # logging.info(
                #     "Completed low coverage python regions filtering step in %f seconds",
                #     perf_counter() - startlowcovpython,
                # )
                startlowcovbed = perf_counter()
                with NamedTemporaryFile(
                    mode="wt", dir=self.meteor.tmp_dir, delete=False, suffix=".bed"
                ) as temp_low_cov_sites:
                    # process = psutil.Process()
                    # mem_before = (
                    #     process.memory_info().rss / 1024 / 1024
                    # )  # Convert to MB
                    output = run(
                        [
                            "bedtools",
                            "genomecov",
                            "-bga",
                            "-ibam",
                            str(cram_file.resolve()),
                        ],
                        check=False,
                        capture_output=True,
                    ).stdout.decode("utf-8")
                    # mem_after = process.memory_info().rss / 1024 / 1024  # Convert to MB
                    # print(f"Memory usage before: {mem_before} MB")
                    # print(f"Memory usage after: {mem_after} MB")
                    # print(f"Memory usage during: {mem_after - mem_before} MB")

                    self.filter_low_cov_sites(output, temp_low_cov_sites)
                    logging.info(
                        "Completed low coverage bedtools regions filtering step in %f seconds",
                        perf_counter() - startlowcovbed,
                    )
                    startlowcov = perf_counter()
                    with Popen(
                        [
                            "bcftools",
                            "consensus",
                            "--mask",
                            temp_low_cov_sites.name,
                            "--mask-with",
                            self.meteor.DEFAULT_GAP_CHAR,
                            "-f",
                            str(reference_file.resolve()),
                            str(vcf_file.resolve()),
                        ],
                        stdout=PIPE,
                    ) as bcftools_process:
                        # capture output of bcftools_process
                        bcftools_output = bcftools_process.communicate()[0]
                        # compress output using lzma
                        compressed_output = lzma.compress(bcftools_output)

                        # write compressed output to file
                        with open(str(consensus_file.resolve()), "wb") as f:
                            f.write(compressed_output)
                logging.info(
                    "Completed low coverage regions filtering step in %f seconds",
                    perf_counter() - startlowcov,
                )
        logging.info("Completed SNP calling in %f seconds", perf_counter() - start)
        config = self.set_variantcalling_config(
            cram_file, vcf_file, consensus_file, bcftools_version
        )
        self.save_config(config, self.census["Stage3FileName"])
