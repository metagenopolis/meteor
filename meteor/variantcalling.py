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

from subprocess import check_call
from dataclasses import dataclass
from pathlib import Path
from configparser import ConfigParser
from datetime import datetime
from typing import Type
from meteor.session import Session, Component
from time import perf_counter
from tempfile import NamedTemporaryFile
import logging


@dataclass
class VariantCalling(Session):
    """Run bcftools"""

    meteor: Type[Component]
    census: dict
    depth: int

    def set_variantcalling_config(
        self, cmd: str, sam_file: Path
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
            "variant_calling_tool": "bcftools",
            "variant_calling_version": "NA",
            "mapping_date": datetime.now().strftime("%Y-%m-%d"),
            "reference_name": self.census["reference"]["reference_info"][
                "reference_name"
            ],
            # "processed_read_count": self.census["census"]["sample_info"]["sequenced_read_count"]
        }
        config["mapping_file"] = {
            "mapping_file_count": "1",
            "bowtie_file_1": sam_file.name,
            "mapping_file_format": "sam",
            # "mapping_file_format": "bam",
        }
        return config

    def execute(self) -> bool:
        """Call variants reads"""
        # Start mapping
        bam_file = (
            self.census["directory"]
            / f"{self.census['census']['sample_info']['sample_name']}.bam"
        )
        vcf_file = (
            self.meteor.strain_dir
            / f"{self.census['census']['sample_info']['sample_name']}.vcf.gz"
        )
        reference = (
            self.meteor.ref_dir
            / self.census["reference"]["reference_file"]["fasta_dir"]
            / self.census["reference"]["reference_file"]["fasta_filename_1"]
        )
        start = perf_counter()
        # "--skip-indels", ?
        with NamedTemporaryFile(mode="wt", dir=self.meteor.tmp_dir) as temp_vcf_file:
            check_call(
                [
                    "bcftools",
                    "mpileup",
                    "-d",
                    str(self.depth),
                    "-Ob",
                    "-f",
                    str(reference.resolve()),
                    str(bam_file.resolve()),
                    "--threads",
                    str(self.meteor.threads),
                    "-o",
                    temp_vcf_file.name,
                ]
            )
            check_call(
                [
                    "bcftools",
                    "call",
                    "-v",
                    "-c",
                    "--threads",
                    str(self.meteor.threads),
                    "-Oz",
                    "-o",
                    str(vcf_file.resolve()),
                    temp_vcf_file.name,
                ]
            )
        # gatk AddOrReplaceReadGroups -I Zymo_6300.bam -O test_grp.bam -RGID 4 -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM 20
        # gatk SplitNCigarReads -R ../../catalogue/mock/fasta/mock.fasta -I test_grp.bam -O test.bam
        # gatk --java-options "-Xmx8g" HaplotypeCaller -R ../../catalogue/mock/fasta/mock.fasta -I test.bam -O test.vcf.gz
        logging.info("Completed mapping creation in %f seconds", perf_counter() - start)
        # config = self.set_variantcalling_config(parameters, sam_file)
        # self.save_config(config, self.census["Stage1FileName"])
        return True
