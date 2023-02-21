from itertools import product
from pathlib import Path
from configparser import ConfigParser
from dataclasses import dataclass, field
from session import Session, Component
from typing import Type
import logging
import re

"""
Import and prepare fastq
"""

@dataclass
class FastqImporter(Session):
    """FastqImporter handle the fastq import
    """
    meteor: Type[Component]
    isdispatched:bool
    mask_sample_name: str
    project_name: str
    ext_r1: tuple= field(default_factory=tuple)
    ext_r2: tuple= field(default_factory=tuple)
    ext: tuple= field(default_factory=tuple)

    def __post_init__(self)->None:
        self.ext_r1 = tuple(["".join(i) for i in product(self.meteor.sequence, "1",
                            self.meteor.extension, self.meteor.compression)])
        self.ext_r2 = tuple(["".join(j) for j in product(self.meteor.sequence, "2",
                            self.meteor.extension, self.meteor.compression)])
        self.ext = tuple(["".join(i) for i in product(self.meteor.extension, 
                                                      self.meteor.compression)])
    
    def replace_ext(self, path: Path)->str:
        """Replace all fastq/compressed extension to get a fullpathname for counter
        """
        extensions = "".join(path.suffixes)
        path_str = str(path)
        for e in self.ext:
            path_str = path_str.replace(e, "")
        return path_str


    def set_fastq_config(self, sample_name:str, tag:str, fastq_file:Path,
                        full_sample_name:str)->ConfigParser:
        """Set configuration for fastq

        :param sample_name: Sample name
        :param tag: Identified tag (single or 1 or 2)
        """
        config = ConfigParser()
        config["sample_info"] = {
            "sample_name": sample_name,
            "condition_name": "NA", # What is this ?
            "project_name" : self.project_name,
            "sequencing_date": "1900-01-01", # Then it is useless
            "sequencing_device": "proton",
            "census_status": "0", # what is this ?
            "read_length": "-1", # Then it is useless
            "tag": tag,
            "full_sample_name": full_sample_name
        }
        config["sample_file"] = {
            "fastq_file": fastq_file.name,
            "is_compressed": "1"
        }
        return config

    def execute(self) -> bool:
        """Dispatch the fastq file
        """
        logging.info("Actual import")
        if self.isdispatched:
            fastq_file_list = self.meteor.fastq_dir.glob( "*/*.f*q*")
        else:
            fastq_file_list = self.meteor.fastq_dir.glob("*.f*q*")
        for fastq_file in fastq_file_list:
            logging.info(f"Import {fastq_file}")
            # Extract paired-end info
            tag = "single"
            if fastq_file.name.endswith(self.ext_r1):
                tag = "1"
            elif fastq_file.name.endswith(self.ext_r2):
                tag = "2"
            # Get rid of all possible extension
            full_sample_name = self.replace_ext(fastq_file)
            if self.isdispatched:
                sample_name = fastq_file.parent.name
            else:
                # split full sample name (in fact library/run name) in order
                # to extract sample_name according to regex mask
                full_sample_name_array = re.search(self.mask_sample_name,
                                                full_sample_name)
                if full_sample_name_array:
                    sample_name = full_sample_name_array[0]
                else:
                    # sample do not match the mask
                    continue
            if self.isdispatched:
                sample_dir = fastq_file.parent
            else:
                # create directory for the sample and move fastq file into
                sample_dir = fastq_file.parent / sample_name
                sample_dir.mkdir(exist_ok=True)
                sym_fastq = Path(sample_dir / fastq_file.name)
                sym_fastq.symlink_to(fastq_file.resolve())
            config_fastq = self.set_fastq_config(sample_name, tag, sym_fastq, 
                                                 full_sample_name)
            config_path = sample_dir / f"{full_sample_name}_census_stage_0.ini"
            self.save_config(config_fastq, config_path)
        return True
