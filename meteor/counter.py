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

"""Run mapping and performs counting"""

import gzip
import bz2
import lzma
import logging
import sys
from dataclasses import dataclass, field
from tempfile import mkdtemp, NamedTemporaryFile
import tempfile
from configparser import ConfigParser
from pathlib import Path
from subprocess import check_call
from meteor.mapper import Mapper
from meteor.session import Session, Component
from typing import Type


@dataclass
class Counter(Session):
    """Counter session map and count
    """
    meteor: Type[Component]
    counting_type: str
    mapping_type: str
    counting_only: bool
    mapping_only: bool
    ini_data: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.meteor.tmp_dir = Path(mkdtemp(dir=self.meteor.tmp_path))
        self.meteor.mapping_dir.mkdir(exist_ok=True)

    def count_index_fastq(self, fastq_file: Path, output_desc: tempfile._TemporaryFileWrapper) -> tuple:
        """Count the number of bases

        :param fastq_file: Path object of the fastq_file
        :param output_desc: Output tempfile._TemporaryFileWrapper descriptor
        :return: (Tuple) A tuple giving the read count and base count
        """
        read_count = 0
        base_count = 0
        if fastq_file.suffix == ".gz":
            in_fq = gzip.open(fastq_file, "rt")
        elif fastq_file.suffix == ".bz2":
            in_fq = bz2.open(fastq_file, "rt")
        elif fastq_file.suffix == ".xz":
            in_fq = lzma.open(fastq_file, "rt")
        else:
            in_fq = open(fastq_file, "rt", encoding="UTF-8")
        # read input fastq line by line
        with in_fq:
            for read_count, line in enumerate(in_fq, start=1):  # pylint: disable=unused-variable
                output_desc.write(f"@{read_count}\n")
                # read the sequence
                read = next(in_fq)
                output_desc.write(f"{read}")
                base_count += len(read.strip())
                # pass the plus
                output_desc.write(next(in_fq))
                # pass the quality
                output_desc.write(next(in_fq))
            output_desc.flush()
        return read_count, base_count

    def set_workflow_config(self, ref_ini) -> ConfigParser:
        """Write configuration file for reference genome

        :return: (ConfigParser) A configparser object
        """
        config = ConfigParser()
        config["worksession"] = {
            "meteor.reference.dir": self.meteor.ref_dir.name,
            "meteor.db.type": "binary",
            "meteor.mapping.program": "bowtie2",
            "meteor.mapping.file.format": "sam",
            "meteor.is.cpu.percentage": "0",
            "meteor.cpu.count": str(self.meteor.threads),
            "meteor.excluded.reference.count": "0"
        }
        config["main_reference"] = {
            "meteor.reference.name": self.meteor.ref_name,
            "meteor.matches": "10000",
            "meteor.mismatches": "5",
            "meteor.is.perc.mismatches": "1",
            "meteor.bestalignment": "1",
            "meteor.mapping.prefix.name": f"mapping_vs_{ref_ini['reference_info']['reference_name']}",
            "meteor.counting.prefix.name": f"vs_{ref_ini['reference_info']['reference_name']}"
        }
        return config

    def launch_mapping(self) -> None:
        """Create temporary indexed files and map aga
        """
        logging.info("Launch mapping")
        # loop on each library
        for library, dict_data in self.ini_data.items():
            census = dict_data["census"]
            sample_info = census["sample_info"]  # reference
            sample_file = census["sample_file"]
            logging.info("Meteor Mapping task description")
            logging.info("Sample name = %s",  sample_info["sample_name"])
            logging.info("Library name = %s", sample_info["full_sample_name"])
            logging.info("Project name = %s", sample_info["project_name"])
            logging.info("Sequencing device = %s", sample_info["sequencing_device"])
            logging.info("Workflow = %s", library.name)

            # reindexing this library reads and fill FLibraryIndexerReport
            fastq_path = self.meteor.fastq_dir / sample_file["fastq_file"]
            try:
                with NamedTemporaryFile(mode="wt", dir=self.meteor.tmp_dir) as output_desc:
                    read_count, base_count = self.count_index_fastq(fastq_path, output_desc)
                    census.set("sample_info", "census_status",  str(1))
                    census.set("sample_info", "indexed_read_length", str(1))
                    census.set("sample_info", "sequenced_read_count",  str(read_count))
                    census.set("sample_info", "indexed_sequenced_read_count",  str(read_count))
                    census.set("sample_info", "indexed_sequenced_base_count", str(base_count))
                    census.set("sample_info", "is_data_prefixed", str(0))
                    dict_data["census"] = census
                    # mapping this library on the reference
                    mapping_process = Mapper(self.meteor, dict_data,
                                             Path(output_desc.name),
                                             self.mapping_type)
                    if not mapping_process.execute():
                        raise ValueError(f"Error, TaskMainMapping failed: {library}")
            except IOError:
                logging.error("Cannot create temporary files in %s", self.meteor.tmp_dir.name)

    def launch_counting(self, workflow_ini: Path) -> None:
        """Launch meteor counter

        :param workflow_ini: A path object to workflow configuration
        """
        logging.info("Launch counting")
        # "-t", str(self.meteor.tmp_dir) + "/"
        check_call(["meteor-counter", "-i", str(self.meteor.fastq_dir) + "/",
                    "-p", str(self.meteor.ref_dir.parent.resolve()) + "/",
                    "-o", str(self.meteor.mapping_dir) + "/", "-w", workflow_ini,
                    "-c", self.counting_type, "-f"])

    def execute(self) -> bool:
        """Compute the mapping"""
        mapping_done = True
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
            logging.error("Error, no *_reference.ini file found in %s", self.meteor.ref_dir)
            sys.exit()
        try:
            census_ini_files = list(self.meteor.fastq_dir.glob("*_census_stage_0.ini"))
            assert len(census_ini_files) > 0
            #  MAPPING THIS LIBRARY ON MAIN REFERENCE
            for library in census_ini_files:
                census_ini = ConfigParser()
                with open(library, "rt", encoding="UTF-8") as cens:
                    census_ini.read_file(cens)
                    sample_info = census_ini["sample_info"]
                    # sample_file = census_ini['sample_file'] # reference
                    stage1_dir = (self.meteor.mapping_dir / sample_info["sample_name"] /
                                  f"mapping_vs_{self.meteor.ref_name}_{census_ini['sample_info']['full_sample_name']}")
                    # sample_info["full_sample_name"]
                    stage1_dir.mkdir(exist_ok=True, parents=True)
                    self.ini_data[library] = {
                        "census": census_ini,
                        "directory": stage1_dir,
                        "Stage1FileName": stage1_dir / library.name.replace("stage_0", "stage_1"),
                        "reference": ref_ini
                    }
                if not self.ini_data[library]["Stage1FileName"].exists():
                    mapping_done = False
            if not self.counting_only:
                # mapping already done and no overwriting
                if mapping_done:
                    logging.info("Mapping already done for sample: %s", sample_info["sample_name"])
                    logging.info("Skipped !")
                else:
                    self.launch_mapping()
            if not self.mapping_only:
                logging.info("Launch mapping")
                config = self.set_workflow_config(ref_ini)
                workflow_ini = self.meteor.mapping_dir / "workflow.ini"
                self.save_config(config, workflow_ini)
                self.launch_counting(workflow_ini)
            logging.info("Done ! Job finished without errors ...")
            self.meteor.tmp_dir.rmdir()
        except AssertionError:
            logging.error("Error, no *_census_stage_0.ini file found in %s", self.meteor.fastq_dir)
            sys.exit()
        return True
