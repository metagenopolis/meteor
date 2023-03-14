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

    def execute(self) -> bool:
        ...
