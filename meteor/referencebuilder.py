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

"""Prepare reference for meteor and index"""

import logging
import gzip
import bz2
import lzma
import sys
from subprocess import check_call, run
from pathlib import Path
from dataclasses import dataclass, field
from datetime import datetime
from packaging.version import parse
from textwrap import fill
from meteor.session import Session, Component
from typing import Iterator, Tuple, ClassVar


@dataclass
class ReferenceBuilder(Session):
    """Index the database with bowtie2 using meteor style"""

    DEFAULT_NUM_THREADS: ClassVar[int] = 1

    meteor: type[Component]
    input_fasta: Path
    fasta_dir: Path = field(default_factory=Path)
    database_dir: Path = field(default_factory=Path)
    output_annotation_file: Path = field(default_factory=Path)
    output_fasta_file: Path = field(default_factory=Path)
    # pysam_test: bool = True

    def __post_init__(self) -> None:
        # Create reference genome directory if it does not already exist
        self.meteor.ref_dir.mkdir(exist_ok=True)

        # Create subdirectories for fasta files and reference indices
        self.fasta_dir = self.meteor.ref_dir / self.meteor.ref_name / "fasta"
        self.fasta_dir.mkdir(exist_ok=True, parents=True)
        self.database_dir = self.meteor.ref_dir / self.meteor.ref_name / "database"
        self.database_dir.mkdir(exist_ok=True)

        # Read input fasta file and create new fasta file for each chromosome or contig
        self.output_annotation_file = (
            self.database_dir / f"{self.meteor.ref_name}_annotation.tsv"
        )
        self.output_fasta_file = self.fasta_dir / f"{self.meteor.ref_name}.fasta"
        self.output_index_file = self.fasta_dir / f"{self.meteor.ref_name}.dict"

        # Write configuration file
        config_ref = self.set_reference_config()
        config_path = (
            self.meteor.ref_dir
            / self.meteor.ref_name
            / f"{self.meteor.ref_name}_reference.json"
        )
        self.save_config(config_ref, config_path)

    def set_reference_config(self) -> dict:  # pragma: no cover
        """Write configuration file for reference genome"""
        config = {
            "meteor_version": self.meteor.version,
            "reference_info": {
                "reference_name": self.meteor.ref_name,
                "reference_date": datetime.now().strftime("%Y-%m-%d"),
                "database_type": "complete",
            },
            "reference_file": {
                "database_dir": "database",
                "fasta_dir": "fasta",
                "fasta_filename": self.output_fasta_file.name,
            },
        }
        return config

    def read_reference(self) -> Iterator[Tuple]:
        """Read genes by genes the catalogue including different compression format.

        :return: A generator object that iterate over each gene
        """
        header = ""
        seq = ""
        if self.input_fasta.suffix == ".gz":
            in_fasta = gzip.open(self.input_fasta, "rt")
        elif self.input_fasta.suffix == ".bz2":
            in_fasta = bz2.open(self.input_fasta, "rt")
        elif self.input_fasta.suffix == ".xz":
            in_fasta = lzma.open(self.input_fasta, "rt")
        else:
            in_fasta = open(self.input_fasta, "rt", encoding="UTF-8")
        with in_fasta:
            for line in in_fasta:
                if line.startswith(">"):
                    if len(seq) > 0:
                        yield header, len(seq), fill(seq, width=80)
                    header = line[1:].split(" ")[0].strip()
                    seq = ""
                else:
                    seq += line.strip()
            if len(seq) > 0:
                yield header, len(seq), fill(seq, width=80)

    def create_reference(self):
        """Write a new reference file for meteor with numeroted genes
        and file giving the correspondance between each gene.
        """
        with self.output_annotation_file.open(
            "wt", encoding="UTF-8"
        ) as output_annotation:
            output_annotation.write("gene_id\tgene_name\tgene_length\n")
            with self.output_fasta_file.open("wt", encoding="UTF-8") as output_fasta:
                for gene_id, (header, len_seq, seq) in enumerate(
                    self.read_reference(), start=1
                ):
                    output_annotation.write(f"{gene_id}\t{header}\t{len_seq}\n")
                    output_fasta.write(f">{gene_id}\n{seq}\n")

    def execute(self) -> None:
        """Build the database"""
        logging.info("Import %s", self.meteor.ref_name)
        # Prepare the reference for meteor
        self.create_reference()
        # Check the bowtie2 version
        bowtie_exec = run(["bowtie2", "--version"], check=False, capture_output=True)
        if bowtie_exec.returncode != 0:
            logging.error(
                "Checking bowtie2 version failed:\n%s",
                bowtie_exec.stderr.decode("utf-8"),
            )
            sys.exit(1)
        bowtie_version = bowtie_exec.stdout.decode("utf-8").split(" ")[2].split("\n")[0]
        if parse(bowtie_version) < self.meteor.MIN_BOWTIE2_VERSION:
            logging.error(
                "The bowtie2 version %s is outdated for meteor. Please update bowtie2 to >= %s.",
                bowtie_version,
                self.meteor.MIN_BOWTIE2_VERSION,
            )
            sys.exit(1)
        # Build the index with bowtie2
        check_call(
            [
                "bowtie2-build",
                "-f",
                "--threads",
                str(self.meteor.threads),
                self.output_fasta_file,
                self.fasta_dir / self.meteor.ref_name,
            ]
        )
        # Build the index with gatk
        # check_call(["gatk", "CreateSequenceDictionary", "-R", self.output_fasta_file])
