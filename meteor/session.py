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

"""Define commune objects"""

from pathlib import Path
from dataclasses import dataclass, field
from typing import Protocol, Iterator, Tuple, Dict
import logging
import sys
import json
import lzma


@dataclass(kw_only=True)
class Component:
    """Set of important constant for meteor"""

    threads: int | None
    fastq_dir: Path = field(default_factory=Path)
    mapping_dir: Path = field(default_factory=Path)
    profile_dir: Path = field(default_factory=Path)
    merging_dir: Path = field(default_factory=Path)
    mapped_sample_dir: Path = field(default_factory=Path)
    strain_dir: Path = field(default_factory=Path)
    tree_dir: Path = field(default_factory=Path)
    ref_dir: Path = field(default_factory=Path)
    ref_name: str = field(default_factory=str)
    # Path given by the user
    tmp_path: Path = field(default_factory=Path)
    # Path used by meteor in the end
    tmp_dir: Path = field(default_factory=Path)
    sequence: tuple = ("_R", ".R", "_", ".", "")
    extension: tuple = (".fq", ".fastq")
    compression: tuple = (".gz", ".bz2", ".xz", "")


class Session(Protocol):
    """Class inheritating from Protocol that present shared function."""

    def check_file(self, filename: Path, expected_colnames: set[str]) -> bool:
        """Check that the expected colnames are in the file.

        :param filename: (Path) An input path object
        :param expected_colnames: (Path) An expected set of colnames
        """

        try:

            if filename.suffix == ".xz":
                header = lzma.open(filename, "rt")
            else:
                header = filename.open("rt", encoding="UTF-8")
            with header:
                real_colnames = set(header.readline().strip("\n").split("\t"))
            assert len(expected_colnames - real_colnames) == 0
        except AssertionError:
            logging.error(
                "Missing columns in %s: %s",
                filename,
                ", ".join(expected_colnames - real_colnames),
            )
            sys.exit(1)
        return True

    def save_config(self, config: Dict, config_path: Path) -> None:  # pragma: no cover
        """Save a configuration file

        :param config: A Dict object
        :param config_path: (Path) An output path object
        """
        with config_path.open("wt", encoding="utf-8") as configfile:
            json.dump(config, configfile)

    def read_json(self, input_json: Path) -> Dict:  # pragma: no cover
        """Read json file
        :param input_json: (Path) An input path to the json file
        :return: A Dict object of the json file
        """
        config = {}
        try:
            with input_json.open("rt", encoding="UTF-8") as json_data:
                config = json.load(json_data)
        except FileNotFoundError:
            logging.error("The file %s does not exist.", input_json)
            sys.exit(1)
        return config

    def get_reference_info(self, ref_dir: Path) -> Dict:
        # Get the json ref
        try:
            ref_json_file_list = list(ref_dir.glob("**/*_reference.json"))
            assert len(ref_json_file_list) == 1
            ref_json_file = ref_json_file_list[0]
            ref_json = self.read_json(ref_json_file)
        except AssertionError:
            logging.error(
                "No *_reference.json file found in %s. "
                "One *_reference.json is expected",
                ref_dir.name,
            )
            sys.exit(1)
        return ref_json

    def get_census_stage(self, mapping_dir: Path, stage: int) -> Dict:
        """Find census_stage_X.json file of a given repertory

        :param mapping_dir: A directory containing one census_stage file
        : param stage: Stage of the census file to find (census_stage_1, census_stage_2, etc)
        """
        try:
            census_json_file_list = list(
                mapping_dir.glob(f"**/*_census_stage_{stage}.json")
            )
            assert len(census_json_file_list) == 1
            census_json_file = census_json_file_list[0]
            census_json = self.read_json(census_json_file)
        except AssertionError:
            logging.error(
                "No *_census_stage_%d.json file found in %s.",
                stage,
                mapping_dir,
            )
            sys.exit(1)
        return census_json

    def update_json(
        self, config: Dict, section: str, new_fields: Dict[str, str]
    ) -> Dict:
        """Add information in the ini configuaration"""
        if section in config:
            for my_field, my_value in new_fields.items():
                if my_field in config[section]:
                    logging.error(
                        "The field %s is already present in the json file.", my_field
                    )
                    sys.exit(1)
                else:
                    config[section][my_field] = my_value
        else:
            config[section] = new_fields
        return config

    def get_sequences(self, fasta_file: Path) -> Iterator[Tuple[int, str]]:
        """Get genes sequences
        :param fasta_file: (Path) A path to fasta file
        :return: A generator providing each header and gene sequence
        """
        gene_id: int = 0
        seq: str = ""
        with fasta_file.open("rt", encoding="UTF-8") as fasta:
            for line in fasta:
                if line.startswith(">"):
                    if len(seq) > 0:
                        yield gene_id, seq
                    gene_id = int(line[1:].strip())
                    seq = ""
                else:
                    seq += line.strip().replace("\n", "")
            if len(seq) > 0:
                yield int(gene_id), seq

    def get_sequences_class(self, fasta_file: Path) -> Iterator[Tuple[str, str]]:
        """Get genes sequences
        :param fasta_file: (Path) A path to fasta file
        :return: A generator providing each header and gene sequence
        """
        gene_id: str = ""
        seq: str = ""
        with fasta_file.open("rt", encoding="UTF-8") as fasta:
            for line in fasta:
                if line.startswith(">"):
                    if len(seq) > 0:
                        yield gene_id, seq
                    gene_id = line[1:].strip()
                    seq = ""
                else:
                    seq += line.strip().replace("\n", "")
            if len(seq) > 0:
                yield gene_id, seq

    def execute(self) -> None: ...
