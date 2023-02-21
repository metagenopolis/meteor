from pathlib import Path
from configparser import ConfigParser
from dataclasses import dataclass, field
from typing import Protocol

"""Start a meteor session, is too short."""

@dataclass(kw_only=True)
class Component:
    """Set of important constant for meteor"""
    threads: int = field(default_factory=int)
    fastq_dir: Path = field(default_factory=Path)
    mapping_dir: Path = field(default_factory=Path)
    ref_dir: Path = field(default_factory=Path)
    ref_name: str = field(default_factory=str)
    tmp_path: Path = field(default_factory=Path)
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
