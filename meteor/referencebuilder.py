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

from subprocess import check_call
from pathlib import Path
from configparser import ConfigParser
from meteor.session import Session, Component
from dataclasses import dataclass, field
from datetime import datetime
from textwrap import fill
from typing import Type
import logging
import gzip
import bz2
import lzma

@dataclass
class ReferenceBuilder(Session):
    """Index the database with bowtie2 using meteor style"""
    meteor: Type[Component]
    input_fasta: Path
    fasta_dir: Path = field(default_factory=Path)
    database_dir: Path = field(default_factory=Path)
    output_annotation_file: Path = field(default_factory=Path)
    output_fasta_file: Path = field(default_factory=Path)

    def __post_init__(self)->None:
        # Create reference genome directory if it does not already exist
        self.meteor.ref_dir.mkdir(exist_ok=True)

        # Create subdirectories for fasta files and reference indices
        self.fasta_dir = self.meteor.ref_dir / self.meteor.ref_name / "fasta"
        self.fasta_dir.mkdir(exist_ok=True, parents=True)
        self.database_dir = self.meteor.ref_dir / self.meteor.ref_name / "database"
        self.database_dir.mkdir(exist_ok=True)

        # Read input fasta file and create new fasta file for each chromosome or contig
        self.output_annotation_file = self.database_dir / f"{self.meteor.ref_name}_lite_annotation"
        self.output_fasta_file = self.fasta_dir / f"{self.meteor.ref_name}.fasta"

        # Write configuration file
        config_ref = self.set_reference_config()
        config_path = self.meteor.ref_dir / self.meteor.ref_name / f"{self.meteor.ref_name}_reference.ini"
        self.save_config(config_ref, config_path)

    def set_reference_config(self)->ConfigParser:
        """Write configuration file for reference genome
        """
        config = ConfigParser()
        config["reference_info"] = {
            "reference_name": self.meteor.ref_name,
            "entry_type": "fragment", # Why ?
            "reference_date": datetime.now().strftime("%Y-%m-%d"),
            "database_type": "text",
            "has_lite_info": "1"
        }
        config["reference_file"] = {
            #IS_LARGE_REFERENCE_STR: 1,
            "database_dir": "database",
            "fasta_dir": "fasta",
            # is it possible to have several fasta
            "fasta_file_count": "1",
            "fasta_filename_1": self.output_fasta_file.name
        }
        config["bowtie2_index"] = {
            "is_large_reference": "0", # WTF
            "is_DNA_space_indexed": "1",
            "dna_space_bowtie_index_prefix_name_1": self.meteor.ref_name
        }
        return config

    def read_reference(self):
        """Read genes by genes the catalogue including different compression format.

        :return: A generator object that iterate over each gene
        """
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
                        yield header, fill(seq, width=80)
                    header = line.split(" ")[0].strip()
                    seq = ""
                else:
                    seq += line.strip()
            yield header, fill(seq, width=80)

    def create_reference(self):
        """Write a new reference file for meteor with numeroted genes
        and file giving the correspondance between each gene.
        """
        with self.output_annotation_file.open("wt", encoding="utf-8", newline="\n") as output_annotation:
            with self.output_fasta_file.open("wt", encoding="utf-8", newline="\n") as output_fasta:
                for gene_id, (header, seq) in enumerate(self.read_reference(), start=1):
                    output_annotation.write(f"{header}")
                    output_fasta.write(f">{gene_id}\n{seq}")

    def execute(self)->bool:
        """Build the database"""
        logging.info("Import %s", self.meteor.ref_name)
        # Prepare the reference for meteor
        self.create_reference()
        # Build the index with bowtie
        check_call(["bowtie2-build", "-f", "-t", str(self.meteor.threads),
            self.output_fasta_file, self.fasta_dir/ self.meteor.ref_name])
        return True
