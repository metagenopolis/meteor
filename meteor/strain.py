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
from dataclasses import dataclass, field
from pysam import VariantFile
from configparser import ConfigParser
from meteor.variantcalling import VariantCalling
from meteor.session import Session, Component
from typing import Type, Dict, Generator, List, Tuple
from tempfile import mkdtemp
from pathlib import Path

@dataclass
class Strain(Session):
    """Counter session map and count
    """
    meteor: Type[Component]
    ini_data: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.meteor.tmp_dir = Path(mkdtemp(dir=self.meteor.tmp_path))
        self.meteor.strain_dir.mkdir(exist_ok=True)

    # def find_strain(self) -> None:
    #     """Compute strain ???"""


    def execute(self) -> bool:
        """Compute the mapping"""
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
            logging.error("Error, no *_reference.ini file found in %s. "
                          "One *_reference.ini is expected", self.meteor.ref_dir)
            sys.exit()
        try:
            census_ini_file_list = list(self.meteor.mapped_sample_dir.glob("*_census_stage_1.ini"))
            assert len(census_ini_file_list) == 1
            census_ini_file = census_ini_file_list[0]
            census_ini = ConfigParser()
            with open(census_ini_file, "rt", encoding="UTF-8") as cens:
                census_ini.read_file(cens)
                sample_info = census_ini["sample_info"]
                self.ini_data["directory"] = self.meteor.mapped_sample_dir
                self.ini_data["census"] = census_ini
                self.ini_data["reference"] = ref_ini
                # Get the bam
                # Get the variant calling
                # Variant calling this library on the reference
                variant_calling_process = VariantCalling(self.meteor, self.ini_data)
                vcf_file = self.ini_data["directory"] / f"{sample_info['sample_name']}.vcf.gz"
                # fasta_file = self.get_strain(bcf_file)
            if not variant_calling_process.execute():
                raise ValueError(f"Error, TaskMainMapping failed: {census_ini_file}")
            # analyse the bcf file
        except AssertionError:
            logging.error("Error, no *_census_stage_0.ini file found in %s", self.meteor.fastq_dir)
            sys.exit()