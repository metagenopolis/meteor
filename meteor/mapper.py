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

from subprocess import run, Popen, PIPE
from dataclasses import dataclass
from pathlib import Path
from datetime import datetime
from packaging.version import Version, parse
from re import findall
from meteor.session import Session, Component
from time import perf_counter

# from tempfile import mkstemp
import pysam
import logging
import sys


@dataclass
class Mapper(Session):
    """Run the bowtie"""

    meteor: type[Component]
    census: dict
    fastq_list: list[str]
    mapping_type: str
    trim: int
    alignment_number: int
    counting_type: str
    identity_threshold: float

    def set_mapping_config(
        self,
        cram_file: Path,
        bowtie_version: str,
        mapping_data: list[int],
    ) -> dict:  # pragma: no cover
        """Define the census 1 configuration

        :param cmd: A string of the specific parameters
        :param cram_file: A path to the raw cram file
        :return: (Dict) A dict object with the census 1 config
        """
        config = {
            "meteor_version": self.meteor.version,
            "sample_info": self.census["census"]["sample_info"],
            "sample_file": self.census["census"]["sample_file"],
            "mapping": {
                "mapping_tool": "bowtie2",
                "mapping_tool_version": bowtie_version,
                "mapping_date": datetime.now().strftime("%Y-%m-%d"),
                "reference_name": self.census["reference"]["reference_info"][
                    "reference_name"
                ],
                "trim": str(self.trim),
                "alignment_number": self.alignment_number,
                "mapping_type": self.mapping_type,
                "identity_threshold": round(self.identity_threshold, 2),
                "total_read_count": mapping_data[0],
                "mapped_read_count": mapping_data[2] + mapping_data[3],
                "overall_alignment_rate": round(
                    (mapping_data[2] + mapping_data[3]) / mapping_data[0] * 100, 2
                ),
                "fastq_files": ",".join(self.fastq_list),
            },
            "mapping_file": {
                "bowtie_file": cram_file.name,
            },
        }
        return config

    # @profile(stream=fp)
    def execute(self) -> None:
        """Map reads"""
        cram_file = (
            self.census["directory"]
            / f"{self.census['census']['sample_info']['sample_name']}_raw.cram"
        )
        reference = (
            self.meteor.ref_dir
            / self.census["reference"]["reference_file"]["fasta_dir"]
            / self.census["reference"]["reference_file"]["fasta_filename"]
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
        bowtie_exec = run(["bowtie2", "--version"], capture_output=True)
        bowtie_version = str(bowtie_exec.stdout).split("\\n")[0].split(" ")[2]
        if bowtie_exec.returncode != 0:
            logging.error(
                "Checking bowtie2 version failed:\n%s",
                bowtie_exec.stderr.decode("utf-8"),
            )
            sys.exit(1)
        elif parse(bowtie_version) < Version("2.3.5"):
            logging.error(
                "The bowtie2 version %s is outdated for meteor. Please update bowtie2 to >=2.3.5.",
                bowtie_version,
            )
            sys.exit(1)
        # Start mapping
        start = perf_counter()
        mapping_exec = Popen(
            [
                "bowtie2",
                parameters,
                "--mm",
                "--no-unal",
                "-x",
                str(bowtie_index.resolve()),
                "-U",
                ",".join(self.fastq_list),
            ],
            stdout=PIPE,
            stderr=PIPE,
        )
        # cramfile_unsorted = Path(mkstemp(dir=self.meteor.tmp_dir)[1])
        with pysam.AlignmentFile(
            mapping_exec.stdout,
            "r",
        ) as samdesc:
            with pysam.AlignmentFile(
                str(cram_file.resolve()),
                # cramfile_unsorted,
                "wc",
                template=samdesc,
                reference_filename=str(reference.resolve()),
            ) as cram:
                for element in samdesc:
                    cram.write(element)
        # pysam.sort(
        #     "-o",
        #     str(cram_file.resolve()),
        #     "-@",
        #     str(self.meteor.threads),
        #     "-O",
        #     "cram",
        #     str(cramfile_unsorted.resolve()),
        #     catch_stdout=False,
        # )
        # pysam.index(str(cram_file.resolve()))
        # Read standard error from the process (non-blocking read)
        mapping_result = mapping_exec.stderr.read().decode("utf-8")
        mapping_exec.stderr.close()

        # Wait for the process to finish and get the exit code
        exit_code = mapping_exec.wait()

        # Check for errors and print the error output if necessary
        if exit_code != 0:
            logging.error("bowtie2 failed:\n%s" % mapping_result)
            sys.exit(1)
        try:
            mapping_log = findall(r"([0-9]+)\s+\(", mapping_result)
            assert len(mapping_log) == 4
            mapping_data = [int(i) for i in mapping_log]
        except AssertionError:
            logging.error("Could not access the mapping result from bowtie2")
            sys.exit(1)
        logging.info("Completed mapping creation in %f seconds", perf_counter() - start)
        config = self.set_mapping_config(cram_file, bowtie_version, mapping_data)
        self.save_config(config, self.census["Stage1FileName"])
