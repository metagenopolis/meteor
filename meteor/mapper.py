from subprocess import check_call
from dataclasses import dataclass
from pathlib import Path
from session import Session, Component
from configparser import ConfigParser
from datetime import datetime
from typing import Type
"""
Effective mapping
"""

@dataclass
class Mapper(Session):
    """
    Run the bowtie
    """
    meteor: Type[Component]
    census: dict
    fastq_reindex: Path
    mapping_type: str

    def set_mapping_config(self, cmd, sam_file):
        # census_ini = self.census["census"]
        config = ConfigParser()
        config["sample_info"] = self.census["census"]["sample_info"]
        config["sample_file"] = self.census["census"]["sample_file"]
        config["mapping"] = {
            "mapping_tool": "bowtie2",
            "mapping_tool_version": "NA",
            "mapping_date":  datetime.now().strftime("%Y-%m-%d"),
            "reference_name": self.census["reference"]["reference_info"]["reference_name"],
            "mapping_cmdline": cmd,
            "parameters": "l-1-m5",
            "mapped_read_length": "-1",
            "mapped_read_length_type": "overall",
            "mismatches": "5",
            "is_mismatches_percentage": "1",
            "matches": "10000",
            "is_local_mapping": str(int(self.mapping_type == "local")),
            "mapping_software": "Meteor",
            "mapping_software_version": "3.3",
            "processed_read_count": self.census["census"]["sample_info"]["sequenced_read_count"]
        }
        config["mapping_file"] = {
            "mapping_file_count": "1",
            "bowtie_file_1": sam_file.name,
            "mapping_file_format": "sam"
        }
        return config

    def execute(self)->bool:
        bowtie_index = self.meteor.ref_dir / self.census["reference"]["reference_info"]["reference_name"] / self.census["reference"]["reference_file"]["fasta_dir"] / self.census["reference"]["bowtie2_index"]["dna_space_bowtie_index_prefix_name_1"]
        # bowtie2 parameters
        if self.mapping_type == "local":
            parameters = f"-p {self.meteor.threads} --local --sensitive-local"
        else:
            parameters = f"-p {self.meteor.threads} --end-to-end --sensitive"
        sam_file = self.census["directory"] / f"{self.census['census']['sample_info']['full_sample_name']}_1.sam"
        # sam_file = self.census["directory"] / f"{self.census['census']['sample_info']['sample_name']}.sam"
        # execute command
        # --mm
        useless = "--trim-to 80 -k 10000"
        # check_call(["bowtie2", parameters, "--mm --no-head --no-sq --no-unal --omit-sec-seq", useless,  
        #             "-x", bowtie_index, "-U", self.fastq_reindex, "-S", sam_file])
        config = self.set_mapping_config(parameters + " "+ useless, sam_file)
        self.save_config(config, self.census["Stage1FileName"])
        return True
