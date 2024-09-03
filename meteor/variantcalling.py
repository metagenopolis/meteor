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
import lzma
import bgzip
from subprocess import CalledProcessError, run, Popen, PIPE
from dataclasses import dataclass
from pathlib import Path
from datetime import datetime
from meteor.session import Session, Component
from time import perf_counter
from tempfile import NamedTemporaryFile
from packaging.version import parse

# import psutil
from pysam import AlignmentFile, FastaFile, VariantFile, faidx, tabix_index
from collections import defaultdict
import pandas as pd
from typing import ClassVar
import numpy as np
import os


@dataclass
class VariantCalling(Session):
    """Run bcftools"""

    # from https://www.bioinformatics.org/sms/iupac.html
    IUPAC: ClassVar[dict] = {
        ("A",): "A",
        ("T",): "T",
        ("G",): "G",
        ("C",): "C",
        ("A", "G"): "R",
        ("C", "T"): "Y",
        ("C", "G"): "S",
        ("A", "T"): "W",
        ("G", "T"): "K",
        ("A", "C"): "M",
        ("C", "G", "T"): "B",
        ("A", "G", "T"): "D",
        ("A", "C", "T"): "H",
        ("A", "C", "G"): "V",
        ("A", "C", "G", "T"): "N",
    }

    meteor: type[Component]
    census: dict
    max_depth: int
    min_depth: int
    min_snp_depth: int
    min_frequency: float
    ploidy: int
    core_size: int

    def set_variantcalling_config(
        self,
        cram_file: Path,
        vcf_file: Path,
        consensus_file: Path,
        freebayes_version: str,
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
                "variant_calling_tool": "freebayes",
                "variant_calling_version": freebayes_version,
                "variant_calling_date": datetime.now().strftime("%Y-%m-%d"),
                "vcf_name": vcf_file.name,
                "consensus_name": consensus_file.name,
                "min_snp_depth": str(self.min_snp_depth),
                "min_frequency": str(self.min_frequency),
                "max_depth": str(self.max_depth),
            },
        }
        return config

    def group_consecutive_positions(
        self, position_count_dict: dict, gene_name: str, gene_length: int
    ):
        # Initialize the result list
        result = []
        # result = {gene_name: []}

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
                # result[gene_name] += [[start, pos, current_count]]
                # Start a new range
                start = pos
                current_count = count

        # Append the last range
        if start != pos:
            result.append((start, pos, current_count))
            # result[gene_name] += [[start, pos, current_count]]
        return pd.DataFrame(
            [
                (gene_name, start, end, count)
                for start, end, count in result
                if count < self.min_depth
            ],
            columns=["gene_id", "startpos", "endpos", "coverage"],
        )
        # return result

    def count_reads_in_gene(
        self,
        cram: AlignmentFile,
        gene_name: str,
        gene_length: int,
        Fasta: FastaFile,
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
    def filter_low_cov_sites(
        self,
        cram_file: Path,
        reference_file: Path,
        # temp_low_cov_sites: tempfile._TemporaryFileWrapper,
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
        # all_genes_dict = {}
        with AlignmentFile(
            str(cram_file.resolve()),
            "rc",
            reference_filename=str(reference_file.resolve()),
            threads=self.meteor.threads,
        ) as cram:
            with FastaFile(filename=str(reference_file.resolve())) as Fasta:
                # For genes with a count
                for _, row in gene_tofilter.iterrows():
                    reads_dict = self.count_reads_in_gene(
                        cram, str(row["gene_id"]), row["gene_length"], Fasta
                    )
                    df = self.group_consecutive_positions(
                        reads_dict, str(row["gene_id"]), row["gene_length"]
                    )
                    # genes_dict = self.group_consecutive_positions(reads_dict, str(row["gene_id"]), row["gene_length"])
                    # if len(genes_dict[str(row["gene_id"])]) > 0:
                    #     all_genes_dict.update(genes_dict)
                    if len(df) > 0:
                        dfs.append(df)
        # All these genes are going to be replaced by gaps
        # Their count is below the threshold level
        gene_ignore = gene_interest[
            gene_interest["coverage"] < self.min_depth
        ].set_index("gene_id")
        # gene_ignore = {row["gene_id"]: row["gene_length"] for _, row in gene_interest[gene_interest["coverage"] < self.min_depth].iterrows()}
        # dfs.append(
        #     pd.DataFrame(
        #         {
        #             "gene_id": gene_ignore["gene_id"],
        #             "startpos": 0,
        #             "endpos": gene_ignore["gene_length"],
        #             "coverage": gene_ignore["coverage"],
        #         }
        #     )
        # )
        # for _, row in gene_interest[gene_interest["coverage"] < self.min_depth].iterrows():
        #     gene_id = row["gene_id"]
        #     # If the gene_id is not already a key in the dict, create a new list
        #     if gene_id not in all_genes_dict:
        #         all_genes_dict[gene_id] = []
        #     # Append the new information to the list corresponding to the gene_id
        #     all_genes_dict[gene_id].append({
        #         "startpos": 0,
        #         "endpos": row["gene_length"],
        #         "coverage": row["coverage"]
        #     })
        sum_cov_bed = pd.concat(dfs, ignore_index=True).set_index("gene_id")
        # .query(f"coverage < {self.min_depth}")
        # sum_cov_bed.to_csv(temp_low_cov_sites, sep="\t", header=False, index=False)
        return sum_cov_bed, gene_ignore
        # return all_genes_dict

    # @memory_profiler.profile
    def create_consensus(
        self,
        reference_file,
        consensus_file,
        low_cov_sites,
        gene_ignore,
        vcf_file,
        bed_file,
    ):
        """Generate a consensus sequence by applying VCF variants to the provided reference genome."""
        # Read the CSV file using pandas
        bed_set = sorted(
            set(
                pd.read_csv(bed_file, usecols=[0], sep="\t", header=None)
                .iloc[:, 0]
                .astype(int)
            )
        )
        # low_cov_sites_dict = low_cov_sites.groupby(low_cov_sites.index).apply(lambda x: x.to_dict(orient='records')).to_dict()
        with VariantFile(vcf_file) as vcf:
            with FastaFile(filename=str(reference_file.resolve())) as Fasta:
                with lzma.open(consensus_file, "wt", preset=0) as consensus_f:
                    # Iterate over all reference sequences in the fasta file
                    len_bed_set = len(bed_set)
                    for i, gene_id in enumerate(bed_set):
                        if i % 1000 == 0:
                            logging.info(f"{i}/{len_bed_set}")
                        ref = str(gene_id)
                        if gene_id in gene_ignore.index:
                            consensus = [
                                self.meteor.DEFAULT_GAP_CHAR
                            ] * gene_ignore.loc[gene_id]["gene_length"]
                        else:
                            consensus = np.array(list(Fasta.fetch(ref)), dtype="<U1")
                            # Apply variants from VCF
                            # startvcf = perf_counter()
                            for record in vcf.fetch(ref):
                                # print(record.info.keys())
                                # print(record.info["AF"])
                                # print(record.alleles)
                                # print(record.alts)
                                ##INFO=<ID=RO,Number=1,Type=Integer,Description="Count of full observations of the reference haplotype.">
                                ##INFO=<ID=AO,Number=A,Type=Integer,Description="Count of full observations of this alternate haplotype.">
                                reference_frequency = record.info["RO"] / (
                                    record.info["RO"] + np.sum(record.info["AO"])
                                )
                                if reference_frequency >= self.min_frequency:
                                    keep_alts = tuple(sorted(list(record.alleles)))
                                else:
                                    keep_alts = tuple(sorted(list(record.alts)))
                                max_len = max(map(len, keep_alts))
                                # MNV vase
                                if max_len > 1:
                                    for i in range(max_len):
                                        mnv = tuple(
                                            sorted(
                                                set(
                                                    keep_alts[k][i]
                                                    for k in range(len(keep_alts))
                                                )
                                            )
                                        )
                                        consensus[record.pos + i - 1] = self.IUPAC[mnv]
                                else:
                                    consensus[record.pos - 1] = self.IUPAC[keep_alts]
                            # Update consensus array for each matching range
                            if gene_id in low_cov_sites.index:
                                selection = low_cov_sites.loc[gene_id]
                                if isinstance(selection, pd.Series):
                                    consensus[
                                        selection["startpos"] : selection["endpos"]
                                    ] = self.meteor.DEFAULT_GAP_CHAR
                                else:
                                    for _, row in selection.iterrows():
                                        # Mark as uncertain
                                        consensus[row["startpos"] : row["endpos"]] = (
                                            self.meteor.DEFAULT_GAP_CHAR
                                        )
                        consensus_f.write(f">{gene_id}\n")
                        consensus_f.write("".join(consensus) + "\n")
                        del consensus

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
        ).resolve()
        msp_file = (
            self.meteor.ref_dir
            / self.census["reference"]["reference_file"]["database_dir"]
            / self.census["reference"]["annotation"]["msp"]["filename"]
        )
        # print(self.census)
        annotation_file = (
            self.meteor.ref_dir
            / self.census["reference"]["reference_file"]["database_dir"]
            / self.census["reference"]["annotation"]["gene_id"]["filename"]
        )

        msp_content = pd.read_csv(
            msp_file,
            sep="\t",
            # names=["msp_name", "gene_id", "gene_name", "gene_category"],
            header=0,
        )
        gene_details = pd.read_csv(
            annotation_file,
            sep="\t",
            # names=["gene_id", "gene_name", "gene_length"],
            header=0,
        )
        if self.census["reference"]["reference_info"]["database_type"] == "complete":
            msp_content = msp_content.loc[msp_content["gene_category"] == "core"]
            msp_content = (
                msp_content.groupby("msp_name")
                .head(self.core_size)
                .reset_index(drop=True)
            )
            # print(msp_content)
            # print(gene_details)
            merged_df = pd.merge(msp_content, gene_details, on="gene_id")
            # Add a constant column with value 0
            merged_df["startpos"] = 0
            # Extract the required columns
            result_df = merged_df[["gene_id", "startpos", "gene_length"]]
            with NamedTemporaryFile(suffix=".bed", delete=False) as temp_bed_file:
                result_df.to_csv(
                    temp_bed_file.name, sep="\t", index=False, header=False
                )
            # print("Why is empty ?")
            # print(temp_bed_file.name)
            bed_file = temp_bed_file.name
        freebayes_exec = run(
            ["freebayes", "--version"], check=False, capture_output=True
        )
        if freebayes_exec.returncode != 0:
            logging.error(
                "Checking freebayes failed:\n%s", freebayes_exec.stderr.decode("utf-8")
            )
            sys.exit(1)
        freebayes_version = (
            freebayes_exec.stdout.decode("utf-8").strip().split("  ")[1][1:]
        )
        if parse(freebayes_version) < self.meteor.MIN_FREEBAYES_VERSION:
            logging.error(
                "The freebayes version %s is outdated for meteor. Please update freebayes to >= %s.",
                freebayes_version,
                self.meteor.MIN_FREEBAYES_VERSION,
            )
            sys.exit(1)

        start = perf_counter()
        logging.info("Run freebayes")
        startfreebayes = perf_counter()
        with reference_file.open("rb") as ref_fh:
            with bgzip.BGZipReader(ref_fh, num_threads=self.meteor.threads) as reader:
                decompressed_reference = reader.read()
        with NamedTemporaryFile(suffix=".fasta", delete=False) as temp_ref_file:
            temp_ref_file.write(decompressed_reference)
            temp_ref_file_path = temp_ref_file.name
        # index on the fly
        faidx(temp_ref_file.name)
        try:
            with Popen(
                [
                    "freebayes",
                    "-i",  # no indel
                    # "-X",
                    "-u",  # no complex observation that may include ins
                    "--min-alternate-count",
                    str(self.min_snp_depth),
                    "--min-alternate-fraction",
                    str(self.min_frequency),
                    "-t",
                    str(bed_file),
                    "-p",
                    str(self.ploidy),
                    "-f",
                    temp_ref_file_path,
                    "-b",
                    str(cram_file.resolve()),
                ],
                stdin=PIPE,
                stdout=PIPE,
            ) as freebayes_process:
                # capture output of bcftools_process
                freebayes_output = freebayes_process.communicate(
                    input=decompressed_reference
                )[0]
                # print(freebayes_output)
                # compress output using bgzip
                with open(str(vcf_file.resolve()), "wb") as raw:
                    with bgzip.BGZipWriter(raw) as fh:
                        fh.write(freebayes_output)
        except CalledProcessError as e:
            logging.error("Freebayes failed with return code %d", e.returncode)
            logging.error("Output: %s", e.output)
            sys.exit()
        finally:
            # Ensure the temporary file is removed after use
            if os.path.isfile(temp_ref_file_path):
                os.remove(temp_ref_file_path)
            if os.path.isfile(temp_ref_file_path + ".fai"):
                os.remove(temp_ref_file_path + ".fai")
        logging.info(
            "Completed freebayes step in %f seconds", perf_counter() - startfreebayes
        )
        # Index the vcf file
        logging.info("Indexing")
        startindexing = perf_counter()
        tabix_index(str(vcf_file.resolve()), preset="vcf", force=True)
        logging.info(
            "Completed indexing step in %f seconds", perf_counter() - startindexing
        )
        # The columns of the tab-delimited BED file are also CHROM, POS
        # and END (trailing columns are ignored), but coordinates are
        # 0-based, half-open. To indicate that a file be treated as BED
        # rather than the 1-based tab-delimited file, the file must have
        # the ".bed" or ".bed.gz" suffix (case-insensitive).
        startlowcovpython = perf_counter()
        logging.info("Detecting low coverage regions")
        low_cov_sites, gene_ignore = self.filter_low_cov_sites(
            cram_file, reference_file
        )
        logging.info(
            "Completed low coverage regions filtering step in %f seconds",
            perf_counter() - startlowcovpython,
        )
        logging.info("Consensus creation")
        startconsensuspython = perf_counter()
        self.create_consensus(
            reference_file,
            consensus_file,
            low_cov_sites,
            gene_ignore,
            vcf_file,
            bed_file,
        )
        logging.info(
            "Completed consensus step with python in %f seconds",
            perf_counter() - startconsensuspython,
        )
        logging.info("Completed SNP calling in %f seconds", perf_counter() - start)
        config = self.set_variantcalling_config(
            cram_file, vcf_file, consensus_file, freebayes_version
        )
        self.save_config(config, self.census["Stage3FileName"])
        if os.path.isfile(temp_bed_file.name):
            os.remove(temp_bed_file.name)
