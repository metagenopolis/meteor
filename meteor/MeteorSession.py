from pathlib import Path
from configparser import ConfigParser
from dataclasses import dataclass
from typing import Optional, Protocol

@dataclass
class session(Protocol):
    threads: Optional[int]

    def save_config(self, config: ConfigParser, config_path: Path)->None:
        with config_path.open("wt", encoding="utf-8") as configfile:
            config.write(configfile)

    def execute(self) -> bool:
        ...