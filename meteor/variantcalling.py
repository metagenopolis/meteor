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

from subprocess import check_call, run
from dataclasses import dataclass
from pathlib import Path
from configparser import ConfigParser
from datetime import datetime
from typing import Type
from meteor.session import Session, Component
from time import perf_counter
from tempfile import NamedTemporaryFile
from packaging.version import Version, parse
import logging
import sys


@dataclass
class VariantCalling(Session):
    """Run bcftools"""

    meteor: Type[Component]
    census: dict
    depth: int
    min_snp_depth: int
    min_frequency_non_reference: float

    def set_variantcalling_config(
        self, bam_file: Path, vcf_file: Path, consensus_file: Path, bcftool_version: str
    ) -> ConfigParser:  # pragma: no cover
        """Define the census 1 configuration

        :param cmd: A string of the specific parameters
        :param bam_file: A path to the sam file
        :return: (ConfigParser) A configparser object with the census 1 config
        """
        config = ConfigParser()
        config["sample_info"] = self.census["census"]["sample_info"]
        config["sample_file"] = self.census["census"]["sample_file"]
        config["mapping"] = {
            "reference_name": self.census["reference"]["reference_info"][
                "reference_name"
            ],
            "bam_name": bam_file.name,
        }
        config["variant_calling"] = {
            "variant_calling_tool": "bcftools",
            "variant_calling_version": bcftool_version,
            "variant_calling_date": datetime.now().strftime("%Y-%m-%d"),
            "vcf_name": vcf_file.name,
            "consensus_name": consensus_file.name,
        }
        return config

    def execute(self) -> bool:
        """Call variants reads"""
        # Start mapping
        bam_file = (
            self.census["mapped_sample_dir"]
            / f"{self.census['census']['sample_info']['sample_name']}.bam"
        )
        vcf_file = (
            self.census["directory"]
            / f"{self.census['census']['sample_info']['sample_name']}.vcf.gz"
        )
        consensus_file = (
            self.census["directory"]
            / f"{self.census['census']['sample_info']['sample_name']}_consensus.fasta"
        )
        reference_file = (
            self.meteor.ref_dir
            / self.census["reference"]["reference_file"]["fasta_dir"]
            / self.census["reference"]["reference_file"]["fasta_filename"]
        )
        bed_file = (
            self.meteor.ref_dir
            / self.census["reference"]["reference_file"]["database_dir"]
            / self.census["reference"]["annotation"]["bed"]
        )
        bcftools_version = (
            str(run(["bcftools", "--version"], capture_output=True).stdout)
            .split("\\n")[0]
            .split(" ")[1]
        )
        if parse(bcftools_version) <= Version("0.1.19"):
            logging.error(
                "Error, the bcftools version %s is outdated for meteor. Please update bcftools.",
                bcftools_version,
            )
            sys.exit()
        start = perf_counter()
        # "--skip-indels", ?
        with NamedTemporaryFile(
            mode="wt", dir=self.meteor.tmp_dir
        ) as temp_vcf_pileup_file:
            with NamedTemporaryFile(
                mode="wt", dir=self.meteor.tmp_dir
            ) as temp_vcf_file:
                check_call(
                    [
                        "bcftools",
                        "mpileup",
                        "-d",
                        str(self.depth),
                        "-Ob",
                        "-f",
                        str(reference_file.resolve()),
                        str(bam_file.resolve()),
                        "--threads",
                        str(self.meteor.threads),
                        "-o",
                        temp_vcf_pileup_file.name,
                    ]
                )
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
                        temp_vcf_file.name,
                        temp_vcf_pileup_file.name,
                    ]
                )
                check_call(["bcftools", "index", temp_vcf_file.name])
                # Only SNP from 100 core genes
                check_call(
                    [
                        "bcftools",
                        "view",
                        "-R",
                        str(bed_file.resolve()),
                        "-f",
                        f"MIN(FMT/DP>{str(self.min_snp_depth)})",
                        "--threads",
                        str(self.meteor.threads),
                        "-q",
                        str(self.min_frequency_non_reference),
                        "-Oz",
                        "-o",
                        str(vcf_file.resolve()),
                        temp_vcf_file.name,
                    ]
                )
                check_call(["bcftools", "index", str(vcf_file.resolve())])
                check_call(
                    [
                        "bcftools",
                        "consensus",
                        "-f",
                        str(reference_file.resolve()),
                        "-o",
                        str(consensus_file.resolve()),
                        str(vcf_file.resolve()),
                    ]
                )
        # gatk AddOrReplaceReadGroups -I Zymo_6300.bam -O test_grp.bam -RGID 4 -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM 20
        # gatk SplitNCigarReads -R ../../catalogue/mock/fasta/mock.fasta -I test_grp.bam -O test.bam
        # gatk --java-options "-Xmx8g" HaplotypeCaller -R ../../catalogue/mock/fasta/mock.fasta -I test.bam -O test.vcf.gz
        logging.info("Completed SNP calling in %f seconds", perf_counter() - start)
        config = self.set_variantcalling_config(
            bam_file, vcf_file, consensus_file, bcftools_version
        )
        self.save_config(config, self.census["Stage3FileName"])
        return True
