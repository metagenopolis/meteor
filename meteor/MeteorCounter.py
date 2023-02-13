from dataclasses import dataclass, field
from tempfile import TemporaryDirectory
from configparser import ConfigParser
from pathlib import Path
from subprocess import check_call
import gzip
import bz2
import lzma
import logging
from MeteorMapper import Mapper
from MeteorSession import session

@dataclass
class Counter(session):
    input_dir: Path
    mapping_dir: Path 
    counting_type: list
    ref_dir: Path
    tmp_path: Path
    force: bool
    counting_only: bool
    mapping_only:bool
    # remove_lock: bool
    # project_dir: Path = ""
    # sample_dir: Path = ""
    # FProjectName: str = ""
    # FMappingDir: str = ""
    # LibraryNames: str = ""
    # mapping_done: bool = True
    # counting_done: bool = True
    # FCountingType: list = field(default_factory=list)
    # ini_files: dict = field(default_factory=dict)
    # library_index_report: dict = field(default_factory=dict)
    # FLibraryCensusIniFileNames: list = field(default_factory=list)
    # FMainMappingCensusIniFileNames: dict = field(default_factory=dict)
    # FExcludedMappingCensusIniFileNames: list = field(default_factory=list)
    # tmp_sample_dir: TemporaryDirectory = field(default_factory=TemporaryDirectory)
    tmp_dir: TemporaryDirectory = field(default_factory=TemporaryDirectory)


    def __post_init__(self)->None:
        if self.tmp_path:
            self.tmp_dir =  TemporaryDirectory(dir=self.tmp_path)
            # self.tmp_sample_dir = TemporaryDirectory(dir=self.tmp_path)

    def count_index_fastq(self, fastq_file: str, output_file: str) -> tuple:
        aBaseCount = 0
        if fastq_file.endswith('.gz'):
            in_fq = gzip.open(fastq_file, "rt")
        elif fastq_file.endswith('.bz2'):
            in_fq = bz2.open(fastq_file, "rt")
        elif fastq_file.endswith('.xz'):
            in_fq = lzma.open(fastq_file, "rt")
        else:
            in_fq = open(fastq_file, "rt")
        # open output fastq in write mode
        with open(output_file, 'wt') as out_fq:
            # read input fastq line by line
            for aReadCount, line in enumerate(in_fq, start=1):
                out_fq.write(f"@{aReadCount}\n")
                # read the sequence
                read = next(in_fq)
                out_fq.write(f"{read}")
                aBaseCount += len(read.strip())
                # pass the plus
                out_fq.write(next(in_fq))
                # pass the quality
                out_fq.write(next(in_fq))
        in_fq.close()
        return aReadCount, aBaseCount

    # def TaskMainMapping(self, census_ini: dict) -> bool:
        # logging.info("Task main mapping ")
        # # What is Stage1FileName ?
        # aMappedCensusIniFileName = self.FMainMappingCensusIniFileNames[iLibrary]["Stage1FileName"]
        # aLibraryMappingDir = self.FMainMappingCensusIniFileNames[iLibrary]["directory"]
        # aLockedMappedCensusIniFileName = aMappedCensusIniFileName+ ".lock"
        # # remove lock file if asked
        # if self.remove_lock and os.path.exists(aLockedMappedCensusIniFileName):
        #     logging.info("Removing lock file: {}".format(aLockedMappedCensusIniFileName))
        #     os.remove(aLockedMappedCensusIniFileName)
        # remove census 1 file if asked
        # if self.force and os.path.exists(aMappedCensusIniFileName):
        #     logging.info("Removing existing census 1 file: {}".format(aMappedCensusIniFileName))
        #     os.remove(aMappedCensusIniFileName)
        # if  os.path.exists(aMappedCensusIniFileName) or os.path.exists(aLockedMappedCensusIniFileName):
        #     logging.info("Skipping library.\n{} already exists or is locked".format(aMappedCensusIniFileName))
        # create lock file in result directory
        # os.makedirs(aLibraryMappingDir, exist_ok=True)
        # open(aLockedMappedCensusIniFileName, 'a').close()

        # aOkMappingProcess = aMeteorMapper.MapRead(
        #                 aLibraryMappingDir,
        #                 main_reference["meteor.mapper.cmd"],
        #                 main_reference["meteor.matches"],
        #                 main_reference["meteor.mismatches"],
        #                 float(main_reference["meteor.is.perc.mismatches"]),
        #                 bool(main_reference["meteor.is.local.mapping"]),
        #                 main_reference["meteor.mapping.file.format"],
        #                 # int(aWorkSessionSection["meteor.is.cpu.percentage"]),
        #                 int(main_reference["meteor.cpu.count"]))
        # remove lock file
        # if os.path.exists(aLockedMappedCensusIniFileName):
        #     os.remove(aLockedMappedCensusIniFileName)
        # # Mapping failed
        # if not aOkMappingProcess:
        #     return False
        # return True

    # def LaunchMapping(self)->None:
    #     logging.info("Launch mapping")
        #~ ( build catalog reference db (bowtie2-build-l) )

        #### why not launching one bowtie process for all the libraries (see bowtie -U option in example below) ????
        #bt2_index="/projects/biodatabank/filtering_database/Homo_sapiens_2014_02_04/bowtie2/homo_sapiens.fna"
        #command_bt2="bowtie2 -q --no-head --no-sq --no-unal --omit-sec-seq --end-to-end --sensitive -k 1 -S $sample/$sample.sam -U $fastq_list -p 40 -x $bt2_index"
        # mapper = Mapper(ini_files=self.ini_files)
                            # reference_ini= reference_ini,
                            # FNGSLibraryIndexerReport=self.FLibraryIndexerReport,
                            # FMappedCensusIniFileName=census_ini["Stage1FileName"],
                            # sample_dir=self.sample_dir)
        # mapper.execute()
        # work_session = self.census_ini["worksession"] # reference
        # main_reference = self.census_ini["main_reference"] # reference
        # reference_ini = (work_session["meteor.reference.dir"] / 
        #                 main_reference["meteor.reference.name"] / 
        #                 f"{main_reference["meteor.reference.name"]}_reference.ini")
        # LOOP ON EACH LIBRARY
        # for library in self.ini_files:
        #     sample_info = self.ini_files[library]["sample_info"] # reference
        #     aSampleFileSection = self.ini_files[library]["sample_file"]
            # logging.info("Meteor Mapping task description")
            # logging.info("Sample name = " +  sample_info["sample_name"])
            # logging.info("Library name = " + sample_info["full_sample_name"])
            # logging.info("Project name = " + sample_info["project_name"])
            # logging.info("Sequencing device = " + sample_info["sequencing_device"])
            # logging.info("Workflow = " + library.name)

            # if not os.path.exists(iLibrary+".lock"):
            ### reindexing this library reads and fill FLibraryIndexerReport
            # self.CountAndReIndexReads(self.ini_files[library])
            # aInputFile = self.sample_dir / aSampleFileSection["fastq_file"]
            # aOutputFile = self.tmp_sample_dir / aSampleFileSection["fastq_file"] / '.idx'
            # aReadCount, aBaseCount = self.count_index_fastq(aInputFile, aOutputFile)
            # self.library_index_report = {
            #     "IsIndexed": 1,
            #     "IndexedfastqFilePath": aOutputFile,
            #     "IndexedFileNameExt": '.idx',
            #     "IndexedReadCount": aReadCount,
            #     "IndexedBaseCount": aBaseCount,
            #     "IndexedReadLength": -1,
            #     "MappedReadLength": -1,
            #     "MappedReadLengthType": "overall",
            # }
            #TODO library_index_report is currently useless
            # MAPPING THIS LIBRARY ON MAIN REFERENCE
            # if not self.TaskMainMapping(self.ini_files[library]):
            #     raise ValueError("Error, TaskMainMapping failed: {iLibrary}")

    def launch_counting(self)->None:
        logging.info("Launch counting")
        # does not need census_stage_0.ini
        #-w /path/to/workflow_tutorial.ini -i /path/to/sample/H1 -p /path/to/project_name -m mapping
        aparameters = " -w " + self.FMeteorJobIniFilename + " -i " + self.sample_dir + " -p " + self.project_dir + " -o " + os.path.basename(self.mapping_dir)
        if len(self.FCountingType) > 0:
            aparameters += " -c " + ",".join(self.FCountingType)
        if self.force: # force overwriting former profiling results done with same parameters
            aparameters += " -f"
        if self.tmp_dir: # path to directory where mapping results (SAM files) are stored. #### NEW
            aparameters += " -t " + self.tmp_dir.name
        check_call([os.path.dirname(os.path.realpath(__file__)) + os.sep + "src" + os.sep + "build" + os.sep + "meteor-counter", "-w", \
            self.aReferenceIniFileName, "-i", self.sample_dir, "-p", self.project_dir, "-o", os.path.basename(self.mapping_dir), \
            "-c", ",".join(self.FCountingType),  "-t", self.tmp_dir.name])


    # def launch_mapping(self)->None:
    #     # reference = self.ref_dir / 
    #     check_call(["bowtie2", f" -p {self.threads}", 
    #     "--mm --no-head --no-sq --no-unal --omit-sec-seq",  
    #     "-x", refere, "-U", 
    #     self.FNGSLibraryIndexerReport["IndexedfastqFilePath"], 
    #     "-S", FMappingOutputFile.name])

    def execute(self)->None:
        # directly from program arguments
        # TODO get it from census ini file, (FProjectName too)
        # ex: /projects/parkinson/mapping
        # self.FMeteorJobIniFile.read_file(open(workflow_ini))
        # self.FMeteorJobIniFilename = workflow_ini
        #### in PrepareWorkSpace ?
        # ex: /projects/parkinson/mapping/KCL_01
        # make directory (and parent) if needed
        #if not os.path.exists(sample_mapping_dir):
        #    os.makedirs(sample_mapping_dir)  #### maybe obsolete cause done in TaskMainMapping and/or TaskExcludedMapping
        #### in PrepareLibraryIniFileNames ?
        # # LOOP ON EACH LIBRARY
        mapping_done = True
        try:
            self.ini_files = list(self.input_dir.glob("*_"+"census_stage_0.ini"))
            assert(len(self.ini_files) > 0)
            # census_ini = ConfigParser()
            # for library in self.ini_files:
                # # Check if lock file exist
                # if os.path.exists(iLibrary+".lock") and not self.no_lock:
                #     sys.exit("Error, lock file not found: {}".format(iLibrary+".lock"))
                # MAPPING THIS LIBRARY ON MAIN REFERENCE
                # with open(library, "rt", encoding="UTF-8") as lib_desc:
                #     census_ini.read_file(lib_desc)
            #         sample_info = census_ini['sample_info'] # reference
            #         stage1_dir = Path(
            #             self.mapping_dir, sample_info["sample_name"], 
            #             f"mapping_vs_{self.ref_dir.name}", 
            #             sample_info["full_sample_name"])
            #         self.ini_files[library] = {
            #             "directory": stage1_dir,
            #             "Stage1FileName": stage1_dir / library.name.replace("stage_0", "stage_1")
            #         }
            #         if not self.ini_files[library]["Stage1FileName"].exists():
            #             mapping_done = False
            # if not self.counting_only:
                # mapping already done and no overwriting
                # if mapping_done and not self.force:
                #     logging.info(f"Mapping already done for sample: {sample_info['sample_name']}")
                #     logging.info("Skipped !")
                # else:
            
            # mapper = Mapper(self.threads, ini_files=self.ini_files)
            # mapper.execute()
                    # self.LaunchMapping()
            if not self.mapping_only:
                logging.info("Launch mapping")
                # self.launch_counting()
            logging.info("Done !\nJob finished without errors ...")
        except AssertionError:
            logging.error(f"Error, no *_census_stage_0.ini file found in {self.input_dir}")

