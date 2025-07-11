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
from typing import Protocol, Iterator, ClassVar, Union, Any
import logging
import sys
import json
import lzma
import importlib.resources
from importlib.metadata import version
from packaging.version import Version
import pandas as pd


@dataclass(kw_only=True)
class Component:
    """Set of important constant for meteor"""

    CONFIG_DATA_FILE: ClassVar[Path] = Path("data/zenodo.json")
    TEST_CATALOGUE: ClassVar[str] = "mock"
    MIN_BOWTIE2_VERSION: ClassVar[Version] = Version("2.3.5")
    MIN_FREEBAYES_VERSION: ClassVar[Version] = Version("1.3.6")
    DEFAULT_GAP_CHAR: ClassVar[str] = "N"
    DEFAULT_CORE_SIZE: ClassVar[int] = 100

    threads: int
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
    version: str = version("meteor")

    @staticmethod
    def load_catalogues_config() -> dict:
        try:
            config_data = importlib.resources.files("meteor") / str(
                Component.CONFIG_DATA_FILE
            )
            with importlib.resources.as_file(config_data) as configuration_path:
                with configuration_path.open("rt", encoding="UTF-8") as config:
                    return json.load(config)
        except FileNotFoundError:
            logging.error(
                "The file %s is missing in meteor source",
                Component.CONFIG_DATA_FILE.name,
            )
            sys.exit(1)

    @staticmethod
    def get_available_catalogues() -> list[str]:
        catalogues_config = Component.load_catalogues_config()
        available_catalogues = list(catalogues_config.keys())
        available_catalogues.remove(Component.TEST_CATALOGUE)
        return available_catalogues

    @staticmethod
    def check_catalogue(ref_json: dict) -> None:
        """Check that the catalogue is present in the config file

        :param config: A Dict object of the json file
        """
        try:
            zenodo = Component.load_catalogues_config()
            if ref_json["reference_info"]["database_type"] == "complete":
                assert (
                    ref_json["reference_info"]["reference_date"]
                    == zenodo[ref_json["reference_info"]["reference_name"]][
                        "file_info"
                    ]["reference_date"]
                )
            else:
                # we could do better
                assert (
                    ref_json["reference_info"]["reference_date"]
                    == zenodo[
                        ref_json["reference_info"]["reference_name"].replace(
                            "_" + ref_json["reference_info"]["database_type"], ""
                        )
                    ]["taxonomy_info"]["reference_date"]
                )
        except AssertionError:
            logging.warning(
                "This is not the version of the catalogue expected with "
                "Meteor %s. You might encounter crash/bug. Consider "
                "to download a new version of the catalogue with "
                "meteor download.",
                version("meteor"),
            )
            # sys.exit(1)


class Session(Protocol):
    """Class inheritating from Protocol that present shared function."""

    # def check_file(self, filename: Path, expected_colnames: set[str]) -> bool:
    #     """Check that the expected colnames are in the file.

    #     :param filename: (Path) An input path object
    #     :param expected_colnames: (Path) An expected set of colnames
    #     """
    #     try:
    #         if filename.suffix == ".xz":
    #             header = lzma.open(filename, "rt")
    #         else:
    #             header = filename.open("rt", encoding="UTF-8")
    #         with header:
    #             real_colnames = set(header.readline().strip("\n").split("\t"))
    #         assert len(expected_colnames - real_colnames) == 0
    #     except AssertionError:
    #         logging.error(
    #             "Missing columns in %s: %s",
    #             filename,
    #             ", ".join(expected_colnames - real_colnames),
    #         )
    #         sys.exit(1)
    #     return True

    def save_config(self, config: dict, config_path: Path) -> None:  # pragma: no cover
        """Save a configuration file

        :param config: A Dict object
        :param config_path: (Path) An output path object
        """
        with config_path.open("wt", encoding="utf-8") as configfile:
            json.dump(config, configfile, indent = 4)

    def read_json(self, input_json: Path) -> dict:  # pragma: no cover
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

    def get_reference_info(self, ref_dir: Path) -> dict:
        # Get the json ref

        ref_json_file_list = list(ref_dir.glob("*_reference.json"))
        if len(ref_json_file_list) == 0:
            logging.error(
                "No *_reference.json file found in %s.",
                ref_dir.name,
            )
            sys.exit(1)
        elif len(ref_json_file_list) > 1:
            logging.error(
                "Multiple *_reference.json files found in %s. "
                "One *_reference.json file is expected.",
                ref_dir.name,
            )
            sys.exit(1)

        ref_json_file = ref_json_file_list[0]
        ref_json = self.read_json(ref_json_file)

        return ref_json

    def get_census_stage(self, mapping_dir: Path, stage: int) -> dict:
        """Find census_stage_X.json file of a given directory

        :param mapping_dir: A directory containing one census_stage file
        : param stage: Stage of the census file to find (census_stage_1, census_stage_2, etc)
        """

        census_json_file_list = list(mapping_dir.glob(f"*_census_stage_{stage}.json"))

        if len(census_json_file_list) == 0:
            logging.error(
                "No *_census_stage_%d.json file found in %s.",
                stage,
                mapping_dir,
            )
            sys.exit(1)
        elif len(census_json_file_list) > 1:
            logging.error(
                "Multiple *_census_stage_%d.json files found in %s. "
                "One *_census_stage_%d.json file is expected.",
                stage,
                mapping_dir,
                stage,
            )
            sys.exit(1)

        census_json_file = census_json_file_list[0]
        census_json = self.read_json(census_json_file)

        return census_json

    def update_json(
        self, config: dict, section: str, new_fields: dict[str, Any]
    ) -> dict:
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

    def get_sequences(
        self, fasta_file: Path, use_lzma: bool = False, id_as_int: bool = False
    ) -> Iterator[tuple[Union[int, str] | None, str]]:
        """Get gene sequences
        :param fasta_file: (Path) A path to the fasta file
        :param use_lzma: (bool) Whether to use lzma for file opening
        :param id_as_int: (bool) Whether to convert gene_id to int
        :return: A generator providing each header and gene sequence
        """
        gene_id: Union[int, str] | None = None
        seq: str = ""

        if use_lzma:
            file = lzma.open(fasta_file, "rt")
        else:
            file = fasta_file.open("rt", encoding="UTF-8")

        with file as fasta:
            for line in fasta:
                if line.startswith(">"):
                    if len(seq) > 0:
                        yield gene_id, seq
                    gene_id_raw = line[1:].strip()
                    gene_id = int(gene_id_raw) if id_as_int else gene_id_raw
                    seq = ""
                else:
                    seq += line.strip().replace("\n", "")
            if len(seq) > 0:
                yield gene_id, seq

    # def get_sequences(self, fasta_file: Path) -> Iterator[tuple[int, str]]:
    #     """Get genes sequences
    #     :param fasta_file: (Path) A path to fasta file
    #     :return: A generator providing each header and gene sequence
    #     """
    #     gene_id: int = 0
    #     seq: str = ""
    #     with lzma.open(fasta_file, "rt") as fasta:
    #         for line in fasta:
    #             if line.startswith(">"):
    #                 if len(seq) > 0:
    #                     yield gene_id, seq
    #                 gene_id = int(line[1:].strip())
    #                 seq = ""
    #             else:
    #                 seq += line.strip().replace("\n", "")
    #         if len(seq) > 0:
    #             yield int(gene_id), seq

    # def get_sequences_class(self, fasta_file: Path) -> Iterator[tuple[str, str]]:
    #     """Get genes sequences
    #     :param fasta_file: (Path) A path to fasta file
    #     :return: A generator providing each header and gene sequence
    #     """
    #     gene_id: str = ""
    #     seq: str = ""
    #     with fasta_file.open("rt", encoding="UTF-8") as fasta:
    #         for line in fasta:
    #             if line.startswith(">"):
    #                 if len(seq) > 0:
    #                     yield gene_id, seq
    #                 gene_id = line[1:].strip()
    #                 seq = ""
    #             else:
    #                 seq += line.strip().replace("\n", "")
    #         if len(seq) > 0:
    #             yield gene_id, seq

    def load_data(self, file_path: Path):
        """Load data dynamically based on the file extension.
        :param file_path: The path to the file (including extension).
        :return:  pd.DataFrame: Data loaded into a pandas DataFrame.
        """
        # Choose the appropriate pandas function based on extension
        if "".join(file_path.suffixes[-2:]) in [".tsv", ".tsv.xz"]:
            return pd.read_csv(
                file_path,
                sep="\t",
                header=0,
            )  # Treat `.tsv` as a tab-separated file
        elif file_path.suffix == ".feather":
            return pd.read_feather(file_path)
        else:
            raise ValueError(f"Unsupported file extension: {file_path.suffix}")

    def execute(self) -> None: ...
