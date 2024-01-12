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

from subprocess import run
from dataclasses import dataclass
from pathlib import Path
from configparser import ConfigParser
from datetime import datetime
from typing import Type, List
from packaging.version import Version, parse
from re import findall
from meteor.session import Session, Component

# from memory_profiler import profile
from time import perf_counter
import logging
import sys

# import os
# fp=open('memory_profiler.log','w+')


@dataclass
class Mapper(Session):
    """Run the bowtie"""

    meteor: Type[Component]
    census: dict
    fastq_list: list[str]
    mapping_type: str
    trim: int
    alignment_number: int
    counting_type: str

    def set_mapping_config(
        self,
        cmd: str,
        sam_file: Path,
        bowtie_version: str,
        mapping_data: List[int],
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
            "mapping_tool": "bowtie2",
            "mapping_tool_version": bowtie_version,
            "mapping_date": datetime.now().strftime("%Y-%m-%d"),
            "reference_name": self.census["reference"]["reference_info"][
                "reference_name"
            ],
            "mapping_options": cmd,
            "trim": str(self.trim),
            "alignment_number": str(self.alignment_number),
            "mapping_type": self.mapping_type,
            "total_read_count": str(mapping_data[0]),
            "mapped_read_count": str(mapping_data[2] + mapping_data[3]),
            "overall_alignment_rate": str(
                round((mapping_data[2] + mapping_data[3]) / mapping_data[0] * 100, 2)
            ),
            "fastq_files": ",".join(self.fastq_list),
        }
        config["mapping_file"] = {
            "bowtie_file": sam_file.name,
        }
        return config

    # @profile(stream=fp)
    def execute(self) -> bool:
        """Map reads"""
        sam_file = (
            self.census["directory"]
            / f"{self.census['census']['sample_info']['sample_name']}.sam"
        )
        bowtie_index = (
            self.meteor.ref_dir
            / self.census["reference"]["reference_file"]["fasta_dir"]
            / self.census["reference"]["reference_info"]["reference_name"]
        )
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
        # Check the bowtie2 version
        bowtie_version = (
            str(run(["bowtie2", "--version"], capture_output=True).stdout)
            .split("\\n")[0]
            .split(" ")[2]
        )
        if bowtie_version.returncode != 0:
            logging.error(
                "bowtie2 failed with error code: {}".format(bowtie_version.returncode)
            )
            sys.exit()
        elif parse(bowtie_version) < Version("2.3.5"):
            logging.error(
                "Error, the bowtie2 version %s is outdated for meteor. Please update bowtie2 to >=2.3.5.",
                bowtie_version,
            )
            sys.exit()
        # Start mapping
        start = perf_counter()
        mapping_result = run(
            [
                "bowtie2",
                parameters,
                "--mm --no-unal",
                "-x",
                str(bowtie_index.resolve()),
                "-U",
                ",".join(self.fastq_list),
                "-S",
                str(sam_file.resolve()),
            ],
            capture_output=True,
        ).stderr.decode("utf-8")
        if mapping_result.returncode != 0:
            logging.error(
                "bowtie2 failed with error code: {}".format(mapping_result.returncode)
            )
            sys.exit()
        try:
            mapping_log = findall(r"([0-9]+)\s+\(", mapping_result)
            assert len(mapping_log) == 4
            mapping_data = [int(i) for i in mapping_log]
        except AssertionError:
            logging.error("Error, could not access the mapping result from bowtie2")
            sys.exit()
        logging.info("Completed mapping creation in %f seconds", perf_counter() - start)
        config = self.set_mapping_config(
            parameters, sam_file, bowtie_version, mapping_data
        )
        self.save_config(config, self.census["Stage1FileName"])
        return True
