from pathlib import Path
from configparser import ConfigParser
from dataclasses import dataclass, field
from typing import Optional, Protocol

"""Start a meteor session, is too short."""

@dataclass(frozen=True, kw_only=True)
class Component:
    """Set of important constant for meteor"""
    threads: Optional [int]
    fastq_dir: Optional [Path]
    mapping_dir: Optional [Path]
    ref_dir: Optional [Path]
    ref_name: Optional [str]
    tmp_path: Optional [Path]
    tmp_dir: Path = field(default_factory=Path)
    sequence: tuple = ("R", "")
    extension: tuple = (".fq", ".fastq")
    compression: tuple = (".gz", ".bz2", ".xz")

class Session(Protocol):
    """Class inheritating from Protocol that present shared function."""

    def save_config(self, config: ConfigParser, config_path: Path)->None:
        """Save a configuration file"""
        with config_path.open("wt", encoding="utf-8") as configfile:
            config.write(configfile)

    def execute(self) -> bool:
        ...
