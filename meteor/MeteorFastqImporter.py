import re
from itertools import product
from pathlib import Path
from configparser import ConfigParser
from typing import Union
import logging
from dataclasses import dataclass, field
from MeteorSession import session

@dataclass
class FastqImporter(session):
    isdispatched:bool
    fastq_dir: Path
    mask_sample_name: str
    project_name: str
    sequence: tuple = ("R", "")
    extension: tuple = (".fq", ".fastq")
    compression: tuple = (".gz", ".bz2", ".xz")
    ext_r1: tuple= field(default_factory=tuple)
    ext_r2: tuple= field(default_factory=tuple)

    def __post_init__(self)->None:
        self.ext_r1 = tuple(["".join(i) for i in product(self.sequence, "1", 
                                            self.extension, self.compression)])
        self.ext_r2 = tuple(["".join(j) for j in product(self.sequence, "2", 
                                            self.extension, self.compression)])

    def replace_ext(self, path: Union[str, Path], new_ext: str = "")->Path:
        extensions = "".join(Path(path).suffixes)
        return Path(str(path).replace(extensions, new_ext))


    def set_fastq_config(self, sample_name:str, tag:str, fastq_file:str, 
                        full_sample_name:str)->ConfigParser:
        """Set configuration for fastq

        :param ref_dir: Path object of reference directory
        :param ref_name: Name of the reference
        """
        config = ConfigParser()
        config["sample_info"] = {
            "sample_name": sample_name,
            # "condition_name": "NA", # What is this ?
            "project_name" : self.project_name,
            # "sequencing_date": "1900-01-01", # Then it is useless
            # "sequencing_device": args.techno,
            "census_status": "0", # what is this ?
            "read_length": "-1", # Then it is useless
            "tag": tag,
            "full_sample_name": full_sample_name
        }
        config["sample_file"] = {
            "fastq_file": fastq_file,
            # "is_compressed": iscompressed
        }
        return config

    def execute(self) -> bool:
        # self._decorated.execute()
        logging.info("Actual import")
        if self.isdispatched:
            fastq_file_list = self.fastq_dir.glob( "*/*.f*q*")
        else:
            fastq_file_list = self.fastq_dir.glob("*.f*q*")
        for fastq_file in fastq_file_list:
            logging.info(f"Import {fastq_file}")
            full_sample_name = fastq_file.name
            # Extract paired-end info
            tag = "single"
            if fastq_file.name.endswith(self.ext_r1):
                tag = "1"
            elif fastq_file.name.endswith(self.ext_r2):
                tag = "2"
            full_sample_name = self.replace_ext(fastq_file)
            if self.isdispatched:
                sample_name = fastq_file.parent.name
            else:
                # split full sample name (in fact library/run name) in order 
                # to extract sample_name according to regex mask
                full_sample_name_array = re.search(self.mask_sample_name, 
                                                full_sample_name.name)
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
                fastq_file.rename(sample_dir / fastq_file.name)
            config_fastq = self.set_fastq_config(sample_name, tag, 
                                            fastq_file, full_sample_name)
            config_path = sample_dir / f"{full_sample_name}_census_stage_0.ini"
            self.save_config(config_fastq, config_path)
