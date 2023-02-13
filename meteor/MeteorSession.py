from pathlib import Path
from configparser import ConfigParser
from dataclasses import dataclass
from typing import Optional, Protocol


@dataclass
class session(Protocol):
    threads: Optional[int]
    # locked_file

    def save_config(self, config: ConfigParser, config_path: Path)->None:
        """Save a configuration file"""
        with config_path.open("wt", encoding="utf-8") as configfile:
            config.write(configfile)
    
    # def acquire_lock(self, lockfile: Path):
    #     """Acquire exclusive lock file access"""
    #     locked_file_descriptor = lockfile.open("w+", encoding="UTF-8")
    #     fcntl.lockf(locked_file_descriptor, fcntl.LOCK_EX)
    #     return locked_file_descriptor

    # def release_lock(locked_file_descriptor):
    #     """Release exclusive lock file access"""
    #     locked_file_descriptor.close()

    def execute(self) -> bool:
        ...