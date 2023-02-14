from dataclasses import dataclass, field
from tempfile import mkdtemp, NamedTemporaryFile
from configparser import ConfigParser
from pathlib import Path
from subprocess import check_call
import gzip
import bz2
import lzma
import logging
from meteormapper import Mapper
from meteorsession import Session

"""
Run mapping and performs counting
"""

@dataclass
class Counter(Session):
    """Counter session map and count
    """
    sample_dir: Path
    mapping_dir: Path
    counting_type: list
    ref_dir: Path
    tmp_path: Path
    force: bool
    counting_only: bool
    mapping_only:bool
    ini_data: dict = field(default_factory=dict)
    tmp_dir: Path = field(default_factory=Path)

    def __post_init__(self)->None:
        self.tmp_dir = Path(mkdtemp(dir=self.tmp_path))

    def count_index_fastq(self, fastq_file: Path, output_desc: NamedTemporaryFile) -> tuple:
        """Count the number of bases

        :param fastq_file: Path object of the fastq_file
        :param output_desc: Output NamedTemporaryFile descriptor
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
            output_desc.write(b"@{read_count}\n")
            # read the sequence
            read = next(in_fq)
            output_desc.write(b"{read}")
            base_count += len(read.strip())
            # pass the plus
            output_desc.write(next(in_fq).encode())
            # pass the quality
            output_desc.write(next(in_fq).encode())
        in_fq.close()
        return read_count, base_count


    def LaunchMapping(self)->None:
        """Create temporary indexed files and map
        """
        logging.info("Launch mapping")

        # LOOP ON EACH LIBRARY
        for library in self.ini_data:
            census = self.ini_data[library]["census"]
            sample_info = census["sample_info"] # reference
            aSampleFileSection = census["sample_file"]
            logging.info("Meteor Mapping task description")
            logging.info("Sample name = " +  sample_info["sample_name"])
            logging.info("Library name = " + sample_info["full_sample_name"])
            logging.info("Project name = " + sample_info["project_name"])
            # logging.info("Sequencing device = " + sample_info["sequencing_device"])
            logging.info("Workflow = " + library.name)

            # if not os.path.exists(iLibrary+".lock"):
            ## reindexing this library reads and fill FLibraryIndexerReport
            # self.CountAndReIndexReads(self.ini_data[library])
            fastq_path = self.sample_dir / aSampleFileSection["fastq_file"]
            try:
                with NamedTemporaryFile(dir = self.tmp_dir) as output_desc:
                    aReadCount, aBaseCount = self.count_index_fastq(fastq_path, output_desc)
                    # MAPPING THIS LIBRARY ON MAIN REFERENCE
                    mapping_process = Mapper(self.threads, self.ini_data[library],
                                             self.mapping_dir, Path(output_desc.name))
                    if not mapping_process.execute():
                        raise ValueError(f"Error, TaskMainMapping failed: {library}")
            except IOError:
                logging.error(f"Cannot create temporary files in {self.tmp_dir.name}")


    # def launch_counting(self)->None:
    #     logging.info("Launch counting")
    #     aparameters = " -w " + self.FMeteorJobIniFilename + " -i " + self.sample_dir + " -p " + self.project_dir + " -o " + os.path.basename(self.mapping_dir)
    #     if len(self.FCountingType) > 0:
    #         aparameters += " -c " + ",".join(self.FCountingType)
    #     if self.force: # force overwriting former profiling results done with same parameters
    #         aparameters += " -f"
    #     if self.tmp_dir: # path to directory where mapping results (SAM files) are stored. #### NEW
    #         aparameters += " -t " + self.tmp_dir.name
    #     check_call([os.path.dirname(os.path.realpath(__file__)) + os.sep + "src" + os.sep + "build" + os.sep + "meteor-counter", "-w", \
    #         self.aReferenceIniFileName, "-i", self.sample_dir, "-p", self.project_dir, "-o", os.path.basename(self.mapping_dir), \
    #         "-c", ",".join(self.FCountingType),  "-t", self.tmp_dir.name])


    def execute(self)->bool:
        mapping_done = True
        try:
            census_ini_files = list(self.sample_dir.glob("*_"+"census_stage_0.ini"))
            assert len(census_ini_files) > 0
            for library in census_ini_files:
                # Check if lock file exist
                # if os.path.exists(iLibrary+".lock") and not self.no_lock:
                #     sys.exit("Error, lock file not found: {}".format(iLibrary+".lock"))
                #MAPPING THIS LIBRARY ON MAIN REFERENCE
                census_ini = ConfigParser()
                census_ini.read_file(open(library, "rt", encoding="UTF-8"))
                sample_info = census_ini["sample_info"]
                # sample_file = census_ini['sample_file'] # reference
                stage1_dir = Path(
                    self.mapping_dir, sample_info["sample_name"],
                    f"mapping_vs_{self.ref_dir.name}",
                    sample_info["full_sample_name"])
                self.ini_data[library] = {
                    "census": census_ini,
                    "directory": stage1_dir,
                    "Stage1FileName": stage1_dir / library.name.replace("stage_0", "stage_1")
                }
                print(self.ini_data[library])
                if not self.ini_data[library]["Stage1FileName"].exists():
                    mapping_done = False
            if not self.counting_only:
                # mapping already done and no overwriting
                if mapping_done and not self.force:
                    logging.info(f"Mapping already done for sample: {sample_info['sample_name']}")
                    logging.info("Skipped !")
                else:
                    self.LaunchMapping()
            if not self.mapping_only:
                logging.info("Launch mapping")
                # self.launch_counting()
            logging.info("Done !\nJob finished without errors ...")
            self.tmp_dir.rmdir()
        except AssertionError:
            logging.error(f"Error, no *_census_stage_0.ini file found in {self.input_dir}")
