from subprocess import check_call
from pathlib import Path
from configparser import ConfigParser
from MeteorSession import session
from dataclasses import dataclass, field
from datetime import datetime
from textwrap import fill
import logging

@dataclass
class ReferenceBuilder(session):
    input_fasta_file:Path
    reference_dir: Path
    ref_name: Path
    ref_dir: Path= field(default_factory=Path)
    fasta_dir: Path= field(default_factory=Path)
    database_dir: Path= field(default_factory=Path)

    def __post_init__(self)->None:
        # Create reference genome directory if it does not already exist
        self.ref_dir = self.reference_dir / self.ref_name
        self.ref_dir.mkdir(exist_ok=True)

        # Create subdirectories for fasta files and reference indices
        self.fasta_dir = self.ref_dir / "fasta"
        self.fasta_dir.mkdir(exist_ok=True)
        self.database_dir = self.ref_dir /"database"
        self.database_dir.mkdir(exist_ok=True)

        # Write configuration file
        config_ref = self.set_reference_config(self.ref_name)
        config_path = self.ref_dir / f"{self.ref_name}_reference.ini"
        self.save_config(config_ref, config_path)

        # Read input fasta file and create new fasta file for each chromosome or contig
        self.output_annotation_file = self.database_dir / f"{self.ref_name}_lite_annotation"
        self.output_fasta_file = self.fasta_dir / f"{self.ref_name}.fasta"
    
    def set_reference_config(self, ref_name:str)->ConfigParser:
        """Write configuration file for reference genome

        :param ref_name: Name of the reference
        """
        config = ConfigParser()
        config["reference_info"] = {
            "reference_name": ref_name,
            # "entry_type": "fragment", # Why ?
            "reference_date": datetime.now().strftime("%Y%m%d"),
            "database_type": "text",
            "HAS_LITE_INFO": "1"
        }
        config["reference_file"] = {
            #IS_LARGE_REFERENCE_STR: 1,
            "database_dir": "database",
            "fasta_dir": "fasta",
            #"fasta_file_count": 1,
            # is it possible to have several fasta
            "fasta_file_count": f"{ref_name}.fasta"
        }
        config["bowtie2_index"] = {
            # "is_large_reference": "1", # WTF
            "is_DNA_space_indexed": "1",
            "dna_space_bowtie_index_prefix_name_1": ref_name
        }
        return config

    def read_reference(self):
        seq = ""
        with self.input_fasta_file.open("rt", encoding="utf-8") as input_fasta:
            for line in input_fasta:
                if line.startswith(">"):
                    if len(seq) > 0:
                        yield header, fill(seq, width=80)
                    header = line.split(" ")[0].strip()
                    del(seq)
                    seq = ""
                else:
                    seq += line.strip()
            yield header, fill(seq, width=80)

    def create_reference(self):
        with self.output_annotation_file.open("wt", encoding="utf-8", newline="\n") as output_annotation:
            with self.output_fasta_file.open("wt", encoding="utf-8", newline="\n") as output_fasta:
                for gene_id, (header, seq) in enumerate(self.read_reference(), start=1):
                    output_annotation.write(f"{header}")
                    output_fasta.write(f">{gene_id}\n{seq}")

    def execute(self):
        logging.info(f"Import {self.ref_name}")
        self.create_reference()
        check_call(["bowtie2-build", '-f', '-t', str(self.threads),
            self.output_fasta_file, self.fasta_dir/ self.ref_name])