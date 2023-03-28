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
from meteor.session import Session, Component
from time import perf_counter
import logging


@dataclass
class Mapper(Session):
    """Run the bowtie"""
    meteor: Type[Component]
    census: dict
    fastq_reindex: Path
    mapping_type: str
    trim: int
    alignment_number: int

    def set_mapping_config(self, cmd: str, sam_file: Path) -> ConfigParser:
        """Define the census 1 configuration

        :param cmd: A string of the specific parameters
        :param sam_file: A path to the sam file
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
            "processed_read_count": self.census["census"]["sample_info"]["sequenced_read_count"]
        }
        config["mapping_file"] = {
            "mapping_file_count": "1",
            "bowtie_file_1": sam_file.name,
            "mapping_file_format": "sam"
        }
        return config

    def execute(self) -> bool:
        bowtie_index = (self.meteor.ref_dir / self.census["reference"]["reference_info"]["reference_name"] /
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
            parameters += f"-k {self.alignment_number} "
        sam_file = self.census["directory"] / f"{self.census['census']['sample_info']['full_sample_name']}_1.sam"
        # execute command
        start = perf_counter()
        check_call(["bowtie2", parameters, "--mm --no-head --no-sq --no-unal --omit-sec-seq",
                    "-x", bowtie_index, "-U", self.fastq_reindex, "-S", sam_file])
        logging.info("Completed mapping creation in %f seconds", perf_counter() - start)
        config = self.set_mapping_config(parameters, sam_file)
        self.save_config(config, self.census["Stage1FileName"])
        return True
