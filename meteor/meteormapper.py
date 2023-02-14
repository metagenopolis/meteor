from subprocess import check_call
from dataclasses import dataclass
from pathlib import Path
from meteorsession import Session

"""
Effective mapping
"""

@dataclass
class Mapper(Session):
    """
    Run the bowtie
    """
    census: dict
    mapping_dir: Path
    fastq_reindex: Path

    def execute(self)->bool:
        print(self.mapping_dir.name)
        print(self.census["reference_file"]["fasta_dir"])
        print(self.census["bowtie2_index"]["dna_space_bowtie_index_prefix_name_1"])
        bowtie_index = self.mapping_dir / self.census["reference_file"]["fasta_dir"] / self.census["bowtie2_index"]["dna_space_bowtie_index_prefix_name_1"]
        # bowtie2 parameters
        parameters = f"--mm -p {self.threads} {self.FMapperCmd}"
        sam_file = self.census["directory"] / f"{self.census['sample_info']['full_sample_name']}_1.sam"
        # execute command
        check_call(["bowtie2", parameters, "--no-head --no-sq --no-unal --omit-sec-seq",  "-x", bowtie_index, "-U",
                    self.fastq_reindex, "-S", sam_file.name])
        return True
