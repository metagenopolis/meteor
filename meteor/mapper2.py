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

"""Effective mapping"""

from subprocess import check_call
from dataclasses import dataclass
from pathlib import Path
from configparser import ConfigParser
from datetime import datetime
from typing import Type
from distutils.version import LooseVersion
from meteor.session import Session, Component
# from pysam import view, sort, index  # type: ignore[attr-defined]
# from tempfile import NamedTemporaryFile
# from memory_profiler import profile
from time import perf_counter
import logging
# import os
# fp=open('memory_profiler.log','w+')


@dataclass
class Mapper2(Session):
    """Run the bowtie"""
    meteor: Type[Component]
    census: dict
    fastq_list: list[str]
    mapping_type: str
    trim: int
    alignment_number: int
    counting_type: str

    def set_mapping_config(self, cmd: str, sam_file: Path) -> ConfigParser:  # pragma: no cover
        """Define the census 1 configuration

        :param cmd: A string of the specific parameters
        :param bam_file: A path to the sam file
        :return: (ConfigParser) A configparser object with the census 1 config
        """
        config = ConfigParser()
        config["sample_info"] = self.census["census"]["sample_info"]
        config["sample_file"] = self.census["census"]["sample_file"]
        config["mapping"] = {
            "mapping_tool": "bowtie2",
            "mapping_tool_version": "NA",
            "mapping_date":  datetime.now().strftime("%Y-%m-%d"),
            "reference_name": self.census["reference"]["reference_info"]["reference_name"],
            "mapping_cmdline": cmd,
            "parameters": "l-1-m5",
            "mapped_read_length": "-1",
            "mapped_read_length_type": "overall",
            "mismatches": "5",
            "is_mismatches_percentage": "1",
            "matches": "10000",
            "is_local_mapping": str(int(self.mapping_type == "local")),
            "mapping_software": "Meteor",
            "mapping_software_version": "3.3",
            # "processed_read_count": self.census["census"]["sample_info"]["sequenced_read_count"]
        }
        config["mapping_file"] = {
            "mapping_file_count": "1",
            "bowtie_file_1": sam_file.name,
            "mapping_file_format": "sam",
            # "mapping_file_format": "bam",
        }
        return config

    # # @profile(stream=fp)
    # def create_bam(self, samfile: str, bam_file: str) -> None:
    #     """Function that create a BAM file from a SAM file.

    #     :param samfile: (str) SAM filename
    #     """
    #     # convert sam to bam using pysam
    #     view("-@", str(self.meteor.threads), "-1", "-u", "-S", "-b", "-o",
    #          bam_file, samfile, catch_stdout=False)

    # # @profile(stream=fp)
    # def sort_bam(self, bamfile: str, sorted_bamfile: str) -> None:
    #     """Function that sort a BAM file using pysam

    #     :param bamfile: (str) BAM filename
    #     :return: bamfile (str) BAM filename
    #     """
    #     # sort the bam file
    #     sort("-o", sorted_bamfile, "-@", str(self.meteor.threads),
    #          "-O", "bam", bamfile, catch_stdout=False)
    #     # index the bam file
    #     index(sorted_bamfile)

    # @profile(stream=fp)
    def execute(self) -> bool:
        """Map reads"""
        # bam_file = self.census["directory"] / f"{self.census['census']['sample_info']['sample_name']}.bam"
        # test
        sam_file = self.census["directory"] / f"{self.census['census']['sample_info']['sample_name']}.sam"
        bowtie_index = (self.meteor.ref_dir /
                        self.census["reference"]["reference_file"]["fasta_dir"] /
                        self.census["reference"]["bowtie2_index"]["dna_space_bowtie_index_prefix_name_1"])
        # bowtie2 parameters
        if self.mapping_type == "local":
            parameters = f"-p {self.meteor.threads} --local --sensitive-local "
        else:
            parameters = f"-p {self.meteor.threads} --end-to-end --sensitive "
        if self.trim > 0:
            parameters += f"--trim-to {self.trim} "
        if self.alignment_number > 1:
            # and self.counting_type != "best"
            parameters += f"-k {self.alignment_number} "
        # Start mapping
        start = perf_counter()
        # if LooseVersion(biom.__version__) < LooseVersion("2.0.0"):
        check_call(["bowtie2", parameters, "--mm --no-unal",
                    "-x", str(bowtie_index.resolve()), "-U", ",".join(self.fastq_list),
                    "-S", str(sam_file.resolve())])
        logging.info("Completed mapping creation in %f seconds", perf_counter() - start)
        # with NamedTemporaryFile(mode="wt", dir=self.meteor.tmp_dir) as temp_sam_file:
        #     # execute command
        #     start = perf_counter()
        #     print(" ".join(["bowtie2", parameters, "--mm --no-unal",
        #                 "-x", str(bowtie_index.resolve()), "-U", ",".join(self.fastq_list),
        #                 "-S", temp_sam_file.name]))
        #     check_call(["bowtie2", parameters, "--mm --no-unal",
        #                 "-x", bowtie_index, "-U", ",".join(self.fastq_list),
        #                 "-S", temp_sam_file.name])
        #     logging.info("Completed mapping in %f seconds", perf_counter() - start)
        #     start = perf_counter()
        #     with NamedTemporaryFile(mode="wt", dir=self.meteor.tmp_dir) as temp_bam_file:
        #         self.create_bam(temp_sam_file.name, temp_bam_file.name)
        #         self.sort_bam(temp_bam_file.name, str(bam_file.resolve()))
        #     if not bam_file.exists():
        #         raise ValueError("Failed to create the bam")
        #     logging.info("Completed bam creation in %f seconds", perf_counter() - start)
        # config = self.set_mapping_config(parameters, bam_file)
        config = self.set_mapping_config(parameters, sam_file)
        self.save_config(config, self.census["Stage1FileName"])
        return True
