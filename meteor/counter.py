from dataclasses import dataclass, field
from tempfile import mkdtemp, NamedTemporaryFile
import tempfile
from configparser import ConfigParser
from pathlib import Path
from subprocess import check_call
from mapper import Mapper
from session import Session, Component
import gzip
import bz2
import lzma
import logging
import sys

"""
Run mapping and performs counting
"""

@dataclass
class Counter(Session):
    """Counter session map and count
    """
    meteor: Component
    counting_type: str
    mapping_type: str
    counting_only: bool
    mapping_only:bool
    ini_data: dict = field(default_factory=dict)
   
    def __post_init__(self)->None:
        self.meteor.tmp_dir = Path(mkdtemp(dir=self.meteor.tmp_path))
        self.meteor.mapping_dir.mkdir(exist_ok=True)

    def count_index_fastq(self, fastq_file: Path, output_desc: tempfile._TemporaryFileWrapper) -> tuple:
        """Count the number of bases

        :param fastq_file: Path object of the fastq_file
        :param output_desc: Output tempfile._TemporaryFileWrapper descriptor
        :return: (Tuple) A tuple giving the read count and base count
        """
        base_count = 0
        if fastq_file.suffix == ".gz":
            in_fq = gzip.open(fastq_file, "rt")
        elif fastq_file.suffix == ".bz2":
            in_fq = bz2.open(fastq_file, "rt")
        elif fastq_file.suffix == ".xz":
            in_fq = lzma.open(fastq_file, "rt")
        else:
            in_fq = open(fastq_file, "rt")
        # read input fastq line by line
        for read_count, line in enumerate(in_fq, start=1):
            output_desc.write(f"@{read_count}\n")
            # read the sequence
            read = next(in_fq)
            output_desc.write(f"{read}")
            base_count += len(read.strip())
            # pass the plus
            output_desc.write(next(in_fq))
            # pass the quality
            output_desc.write(next(in_fq))
        output_desc.flush()
        in_fq.close()
        return read_count, base_count

    def set_workflow_config(self, ref_ini)->ConfigParser:
        """Write configuration file for reference genome

        :return: (ConfigParser) A configparser object
        """
        config = ConfigParser()
        config["worksession"] = {
            "meteor.reference.dir": self.meteor.ref_dir.name,
            "meteor.db.type": "binary",
            "meteor.mapping.program": "bowtie2",
            "meteor.mapping.file.format": "sam",
            "meteor.is.cpu.percentage": "0",
            "meteor.cpu.count": str(self.meteor.threads),
            "meteor.excluded.reference.count": "0"
        }
        config["main_reference"] = {
            #IS_LARGE_REFERENCE_STR: 1,
            "meteor.reference.name": self.meteor.ref_name,
            "meteor.matches": "10000",
            "meteor.mismatches": "5",
            "meteor.is.perc.mismatches": "1",
            "meteor.bestalignment": "1",
            "meteor.mapping.prefix.name": f"mapping_vs_{ref_ini['reference_info']['reference_name']}",
            "meteor.counting.prefix.name": f"vs_{ref_ini['reference_info']['reference_name']}"
        }
        return config

    def launch_mapping(self)->None:
        """Create temporary indexed files and map aga
        """
        logging.info("Launch mapping")
        # LOOP ON EACH LIBRARY
        for library in self.ini_data:
            census = self.ini_data[library]["census"]
            # sample_info = census["sample_info"] # reference
            sample_file = census["sample_file"]
            # logging.info("Meteor Mapping task description")
            # logging.info("Sample name = " +  sample_info["sample_name"])
            # logging.info("Library name = " + sample_info["full_sample_name"])
            # logging.info("Project name = " + sample_info["project_name"])
            # logging.info("Sequencing device = " + sample_info["sequencing_device"])
            # logging.info("Workflow = " + library.name)

            # reindexing this library reads and fill FLibraryIndexerReport
            fastq_path = self.meteor.fastq_dir / sample_file["fastq_file"]
            print(fastq_path)
            try:
                with NamedTemporaryFile(mode="wt", dir = self.meteor.tmp_dir) as output_desc:
                    aReadCount, aBaseCount = self.count_index_fastq(fastq_path, output_desc)
                    census.set("sample_info", "census_status",  str(1)),
                    census.set("sample_info", "indexed_read_length", str(1))
                    census.set("sample_info", "sequenced_read_count",  str(aReadCount))
                    census.set("sample_info", "indexed_sequenced_read_count",  str(aReadCount))
                    census.set("sample_info","indexed_sequenced_base_count", str(aBaseCount))
                    census.set("sample_info", "is_data_prefixed", str(0))
                    self.ini_data[library]["census"] = census
                    # MAPPING THIS LIBRARY ON MAIN REFERENCE
                    mapping_process = Mapper(self.meteor, self.ini_data[library], 
                                             Path(output_desc.name),
                                             self.mapping_type)
                    if not mapping_process.execute():
                        raise ValueError(f"Error, TaskMainMapping failed: {library}")
            except IOError:
                logging.error(f"Cannot create temporary files in {self.meteor.tmp_dir.name}")
    
    # def launch_mapping2(self)->None:
    #     logging.info("Launch mapping")
    #     try:
    #         with NamedTemporaryFile(mode="wt", dir = self.meteor.tmp_dir) as output_desc:
    #             for library in tqdm(self.ini_data):
    #                 sample_file = self.ini_data[library]["census"]["sample_file"]
    #                 fastq_path = self.meteor.fastq_dir / sample_file["fastq_file"]
    #                 aReadCount, aBaseCount = self.count_index_fastq(fastq_path, output_desc)
    #             mapping_process = Mapper(self.meteor, self.ini_data[library], 
    #                                          Path(output_desc.name),
    #                                          self.mapping_type)
    #             if not mapping_process.execute():
    #                 raise ValueError(f"Error, TaskMainMapping failed for {library}")
    #     except IOError:
    #         logging.error(f"Cannot create temporary files in {self.meteor.tmp_dir.name}")
    
    def launch_counting(self, workflow_ini)->None:
        logging.info("Launch counting")
        # if self.force: # force overwriting former profiling results done with same parameters
        #     aparameters += " -f"
        check_call([Path(__file__).parent / "src" / "build" / "meteor-counter", 
            "-i", str(self.meteor.fastq_dir) + "/", "-p", str(self.meteor.ref_dir.parent) + "/", 
            "-o", str(self.meteor.mapping_dir) + "/", "-w", workflow_ini,
            "-c", self.counting_type])
        #"-t", str(self.meteor.tmp_dir) + "/"

    def execute(self)->bool:
        mapping_done = True
        try:
            # Get the ini ref
            ref_ini_file = list(self.meteor.ref_dir.glob("**/*_reference.ini"))
            assert len(ref_ini_file) == 1
            ref_ini_file = ref_ini_file[0]
            ref_ini = ConfigParser()
            ref_ini.read_file(open(ref_ini_file))
            self.meteor.ref_name = ref_ini["reference_info"]["reference_name"]
        except AssertionError:
            logging.error(f"Error, no *_reference.ini file found in {self.meteor.ref_dir}")
            sys.exit()
        try:
            census_ini_files = list(self.meteor.fastq_dir.glob("*_census_stage_0.ini"))
            assert len(census_ini_files) > 0
            for library in census_ini_files:
                #MAPPING THIS LIBRARY ON MAIN REFERENCE
                census_ini = ConfigParser()
                census_ini.read_file(open(library, "rt", encoding="UTF-8"))
                sample_info = census_ini["sample_info"]
                print(sample_info)
                # sample_file = census_ini['sample_file'] # reference
                stage1_dir = self.meteor.mapping_dir / sample_info["sample_name"] / f"mapping_vs_{self.meteor.ref_name}_{census_ini['sample_info']['full_sample_name']}"
                # sample_info["full_sample_name"]
                print(stage1_dir)
                stage1_dir.mkdir(exist_ok=True, parents=True)
                self.ini_data[library] = {
                    "census": census_ini,
                    "directory": stage1_dir,
                    "Stage1FileName": stage1_dir / library.name.replace("stage_0", "stage_1"),
                    "reference": ref_ini
                }
                if not self.ini_data[library]["Stage1FileName"].exists():
                    mapping_done = False
            if not self.counting_only:
                # mapping already done and no overwriting
                if mapping_done:
                    logging.info(f"Mapping already done for sample: {sample_info['sample_name']}")
                    logging.info("Skipped !")
                else:
                    self.launch_mapping()
                    # self.launch_mapping2()
            if not self.mapping_only:
                logging.info("Launch mapping")
                config = self.set_workflow_config(ref_ini)
                workflow_ini = self.meteor.mapping_dir / "workflow.ini"
                self.save_config(config, workflow_ini)
                self.launch_counting(workflow_ini)
            logging.info("Done ! Job finished without errors ...")
            self.meteor.tmp_dir.rmdir()
        except AssertionError:
            logging.error(f"Error, no *_census_stage_0.ini file found in {self.meteor.fastq_dir}")
            sys.exit()
        return True
