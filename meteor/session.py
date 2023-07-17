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
from configparser import ConfigParser
from dataclasses import dataclass, field
from typing import Protocol
import logging
import sys


@dataclass(kw_only=True)
class Component:
    """Set of important constant for meteor"""
    threads: int | None
    fastq_dir: Path = field(default_factory=Path)
    mapping_dir: Path = field(default_factory=Path)
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

    def save_config(self, config: ConfigParser, config_path: Path) -> None:
        """Save a configuration file

        :param config: A configparser object
        :param config_path: (Path) An output path object
        """
        with config_path.open("wt", encoding="utf-8") as configfile:
            config.write(configfile)

    def read_ini(self, input_ini: Path) -> ConfigParser:
        config = ConfigParser()
        try:
            with open(input_ini, "rt", encoding="UTF-8") as ini:
                config.read_file(ini)
        except FileNotFoundError:
            msg = f"The file {input_ini} does not exist."
            logging.error(msg)
            sys.exit()
        return config

    def get_reference_info(self, ref_dir: Path) -> ConfigParser:
        # Get the ini ref
        try:
            ref_ini_file_list = list(ref_dir.glob("**/*_reference.ini"))
            assert len(ref_ini_file_list) == 1
            ref_ini_file = ref_ini_file_list[0]
            ref_ini = self.read_ini(ref_ini_file)
        except AssertionError:
            logging.error("Error, no *_reference.ini file found in %s. "
                          "One *_reference.ini is expected", ref_dir)
            sys.exit()
        return ref_ini

    def update_ini(self, config: ConfigParser, section: str, new_fields: dict[str, str]) -> ConfigParser:
        section_suffix = 1
        while section in config.sections():
            section = section + str(section_suffix)
            section_suffix += 1
        config[section] = new_fields
        return config

    def execute(self) -> bool:
        ...
