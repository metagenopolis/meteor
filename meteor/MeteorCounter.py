from dataclasses import dataclass, field
from tempfile import TemporaryDirectory
from configparser import ConfigParser
from pathlib import Path
from subprocess import check_call
import gzip
import bz2
import lzma
import logging


@dataclass
class MeteorSession:
    tmp_path: Path
    input_path: Path
    FMeteorJobIniFilename: str = ""
    project_dir: str = ""
    sample_dir: str = ""
    FSampleName: str = ""
    FProjectName: str = ""
    project_mapping_dir: Path = ""
    sample_mapping_dir: str = ""
    aReferenceIniFileName: str = ""
    FMappingDir: str = ""
    LibraryNames: str = ""
    FLibraryCount: int = 0
    no_lock: bool = False
    mapping_done: bool = True
    counting_done: bool = True
    FCountingTypeList: list = field(default_factory=list)
    ini_files: dict = field(default_factory=dict)
    FLibraryIndexerReport: dict = field(default_factory=dict)
    FLibraryIniFileNames: list = field(default_factory=list)
    FLibraryCensusIniFileNames: list = field(default_factory=list)
    FMainMappingCensusIniFileNames: dict = field(default_factory=dict)
    FExcludedMappingCensusIniFileNames: list = field(default_factory=list)
    tmp_sample_dir: TemporaryDirectory = field(default_factory=TemporaryDirectory)
    tmp_dir: TemporaryDirectory = field(default_factory=TemporaryDirectory)
    FMeteorJobIniFile: ConfigParser = field(default_factory=ConfigParser)

    def __post_init__(self)->None:
        if self.tmp_path:
            self.tmp_dir =  TemporaryDirectory(dir=self.tmp_path)
            self.tmp_sample_dir = TemporaryDirectory(dir=self.tmp_path)

    def CountReadAndReIndexFastqFile(self, fastq_file: str, output_file: str) -> tuple:
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
            for aReadCount, line in enumerate(in_fq):
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

    def TaskMainMapping(self, iLibrary: ConfigParser) -> bool:
        aWorkSessionSection = self.FMeteorJobIniFile["worksession"] # reference
        aReferenceSection = self.FMeteorJobIniFile["main_reference"] # reference
        self.aReferenceIniFileName = aWorkSessionSection["meteor.reference.dir"] + os.sep + \
            aReferenceSection["meteor.reference.name"] + os.sep + aReferenceSection["meteor.reference.name"] + "_reference.ini"

        logging.info("Task main mapping {}".format(iLibrary))

        # What is Stage1FileName ?
        aMappedCensusIniFileName = self.FMainMappingCensusIniFileNames[iLibrary]["Stage1FileName"]
        aLibraryMappingDir = self.FMainMappingCensusIniFileNames[iLibrary]["directory"]

        aLockedMappedCensusIniFileName = aMappedCensusIniFileName+ ".lock"

        # remove lock file if asked
        if self.no_lock and os.path.exists(aLockedMappedCensusIniFileName):
            logging.info("Removing lock file: {}".format(aLockedMappedCensusIniFileName))
            os.remove(aLockedMappedCensusIniFileName)
        # remove census 1 file if asked
        if self.force and os.path.exists(aMappedCensusIniFileName):
            logging.info("Removing existing census 1 file: {}".format(aMappedCensusIniFileName))
            os.remove(aMappedCensusIniFileName)

        if  os.path.exists(aMappedCensusIniFileName) or os.path.exists(aLockedMappedCensusIniFileName):
            logging.info("Skipping library.\n{} already exists or is locked".format(aMappedCensusIniFileName))

        # create lock file in result directory
        os.makedirs(aLibraryMappingDir, exist_ok=True)
        open(aLockedMappedCensusIniFileName, 'a').close()

        aMappingProg = aWorkSessionSection["meteor.mapping.program"]
        if aMappingProg != 'bowtie2':
            sys.exit("Error, unknown or unsupported mapping program : {}".format(aWorkSessionSection["meteor.mapping.program"]))
        aMeteorMapper = MeteorMapper(
                            FMappingProgram=aMappingProg, FLibraryCensusIniFile=self.ini_files[iLibrary],
                            aReferenceIniFileName= self.aReferenceIniFileName,
                            FNGSLibraryIndexerReport=self.FLibraryIndexerReport,
                            FMappedCensusIniFileName=aMappedCensusIniFileName,
                            sample_dir=self.sample_dir)
        aOkMappingProcess = aMeteorMapper.MapRead(
                        aLibraryMappingDir,
                        aReferenceSection["meteor.mapper.cmd"],
                        aReferenceSection["meteor.matches"],
                        aReferenceSection["meteor.mismatches"],
                        float(aReferenceSection["meteor.is.perc.mismatches"]),
                        bool(aReferenceSection["meteor.is.local.mapping"]),
                        aWorkSessionSection["meteor.mapping.file.format"],
                        # int(aWorkSessionSection["meteor.is.cpu.percentage"]),
                        int(aWorkSessionSection["meteor.cpu.count"]))
        # remove lock file
        if os.path.exists(aLockedMappedCensusIniFileName):
            os.remove(aLockedMappedCensusIniFileName)
        # Mapping failed
        if not aOkMappingProcess:
            return False
        return True

    def CountAndReIndexReads(self, aLibraryCensusIniFile: ConfigParser)->None:
        aSampleInfoSection = aLibraryCensusIniFile["sample_info"]
        aSampleFileSection = aLibraryCensusIniFile["sample_file"]

        aInputFile = self.sample_dir + os.sep + aSampleFileSection["fastq_file"]
        aOutputFile = self.tmp_sample_dir / aSampleFileSection["fastq_file"] + '.idx'
        aReadCount, aBaseCount = self.CountReadAndReIndexFastqFile(aInputFile, aOutputFile)
        if not os.path.isfile(aOutputFile):
            sys.exit("WFTFDSFDS")
        self.FLibraryIndexerReport = {
            "IsIndexed": 1,
            "IndexedfastqFilePath": aOutputFile,
            "IndexedFileNameExt": '.idx',
            "IndexedReadCount": aReadCount,
            "IndexedBaseCount": aBaseCount,
            "IndexedReadLength": -1,
            "MappedReadLength": -1,
            "MappedReadLengthType": "overall",
        }

    def LaunchMapping(self)->None:
        logging.info("Launch mapping")
        aOKToContinue = True

        #~ ( build catalog reference db (bowtie2-build-l) )

        #### why not launching one bowtie process for all the libraries (see bowtie -U option in example below) ????
        #bt2_index="/projects/biodatabank/filtering_database/Homo_sapiens_2014_02_04/bowtie2/homo_sapiens.fna"
        #command_bt2="bowtie2 -q --no-head --no-sq --no-unal --omit-sec-seq --end-to-end --sensitive -k 1 -S $sample/$sample.sam -U $fastq_list -p 40 -x $bt2_index"

        aExcludedRefCount = int(self.FMeteorJobIniFile["worksession"]["meteor.excluded.reference.count"])
        # LOOP ON EACH LIBRARY
        for iLibrary in self.ini_files:
            aLibraryCensusIniFile = self.ini_files[iLibrary]
            aSampleInfoSection = aLibraryCensusIniFile["sample_info"] # reference
            aSampleFileSection = aLibraryCensusIniFile["sample_file"] # reference

            FSampleName = aSampleInfoSection["sample_name"]
            FProjectName = aSampleInfoSection["project_name"]
            aSampleLibraryName = aSampleInfoSection["full_sample_name"]

            logging.info("Meteor Mapping task description")
            logging.info("Sample name = " + FSampleName)
            logging.info("Library name = " + aSampleLibraryName)
            logging.info("Project name = " + FProjectName)
            logging.info("Sequencing device = " + aSampleInfoSection["sequencing_device"])
            logging.info("Workflow = " + os.path.basename(iLibrary))

            if not os.path.exists(iLibrary+".lock"):
                ### reindexing this library reads and fill FLibraryIndexerReport
                self.CountAndReIndexReads(aLibraryCensusIniFile)
                # MAPPING THIS LIBRARY ON MAIN REFERENCE
                aOKToContinue = self.TaskMainMapping(iLibrary)
                if not aOKToContinue:
                    raise ValueError("Error, TaskMainMapping failed: {iLibrary}")

    def LaunchCounting(self)->None:
        logging.info("Launch counting")
        # does not need census_stage_0.ini
        #-w /path/to/workflow_tutorial.ini -i /path/to/sample/H1 -p /path/to/project_name -m mapping
        aparameters = " -w " + self.FMeteorJobIniFilename + " -i " + self.sample_dir + " -p " + self.project_dir + " -o " + os.path.basename(self.project_mapping_dir)
        if len(self.FCountingTypeList) > 0:
            aparameters += " -c " + ",".join(self.FCountingTypeList)
        if self.force: # force overwriting former profiling results done with same parameters
            aparameters += " -f"
        if self.tmp_dir: # path to directory where mapping results (SAM files) are stored. #### NEW
            aparameters += " -t " + self.tmp_dir.name
        check_call([os.path.dirname(os.path.realpath(__file__)) + os.sep + "src" + os.sep + "build" + os.sep + "meteor-counter", "-w", \
            self.aReferenceIniFileName, "-i", self.sample_dir, "-p", self.project_dir, "-o", os.path.basename(self.project_mapping_dir), \
            "-c", ",".join(self.FCountingTypeList),  "-t", self.tmp_dir.name])


import fcntl

def acquireLock():
    ''' acquire exclusive lock file access '''
    locked_file_descriptor = open('lockfile.LOCK', 'w+')
    fcntl.lockf(locked_file_descriptor, fcntl.LOCK_EX)
    return locked_file_descriptor

def releaseLock(locked_file_descriptor):
    ''' release exclusive lock file access '''
    locked_file_descriptor.close()

lock_fd = acquireLock()

# ... do stuff with exclusive access to your file(s)

releaseLock(lock_fd)


    def execute(self, mapping_dir:Path, 
                   counting_type:list, counting_only:bool,
                   mapping_only:bool)->None:
        # directly from program arguments
        # TODO get it from census ini file, (FProjectName too)
        
        self.FCountingTypeList = counting_type

        # ex: /projects/parkinson/mapping
        self.project_mapping_dir = self.project_path / mapping_dir
        # self.FMeteorJobIniFile.read_file(open(workflow_ini))
        # self.FMeteorJobIniFilename = workflow_ini
        #### in PrepareWorkSpace ?
        # ex: /projects/parkinson/mapping/KCL_01
        self.sample_mapping_dir = self.project_mapping_dir /  self.FSampleName
        # make directory (and parent) if needed
        #if not os.path.exists(sample_mapping_dir):
        #    os.makedirs(sample_mapping_dir)  #### maybe obsolete cause done in TaskMainMapping and/or TaskExcludedMapping

        #### in PrepareLibraryIniFileNames ?
        self.FLibraryIniFileNames = sorted(self.sample_dir.glob("*_"+"census_stage_0.ini"))
        print(self.FLibraryIniFileNames)
        self.FLibraryCount = len(self.FLibraryIniFileNames)
        if self.FLibraryCount == 0:
            sys.exit("Error, no census_stage_0.ini file found in {}".format(self.sample_dir))
        print("\nLaunch mapping\n")
        aMainReferenceSection = self.FMeteorJobIniFile["main_reference"] # alias
        aMappingDirPathPrefix = ""
        if aMainReferenceSection["meteor.mapping.prefix.name"] is None:
            aMappingDirPathPrefix = "mapping"+"_vs_"+aMainReferenceSection["meteor.reference.name"]+"_l"+ \
                                    "{}-{}".format(aMainReferenceSection["meteor.mapped.readlength"], "m")+aMainReferenceSection["meteor.mismatches"]+"_"
        else:
            aMappingDirPathPrefix = "{}_".format(aMainReferenceSection["meteor.mapping.prefix.name"])

        # # LOOP ON EACH LIBRARY
        for iLibrary in self.FLibraryIniFileNames:
            # Check if lock file exist
            if os.path.exists(iLibrary+".lock") and not self.no_lock:
                sys.exit("Error, lock file not found: {}".format(iLibrary+".lock"))
            # MAPPING THIS LIBRARY ON MAIN REFERENCE
            aLibraryCensusIniFile = ConfigParser()
            aLibraryCensusIniFile.read_file(open(iLibrary))
            self.ini_files[iLibrary] = aLibraryCensusIniFile

            aSampleInfo = aLibraryCensusIniFile["sample_info"] # reference
            aMappingDirPath = aMappingDirPathPrefix + aSampleInfo["full_sample_name"]
            directory = self.project_mapping_dir + os.sep + aSampleInfo["sample_name"]+os.sep+aMappingDirPath
            self.FMainMappingCensusIniFileNames[iLibrary] = {
                "directory": directory,
                "Stage1FileName": directory + os.sep + os.path.basename(iLibrary).replace("census_stage_0", "census_stage_1")
            }
            if not os.path.exists(self.FMainMappingCensusIniFileNames[iLibrary]["Stage1FileName"]):
                self.mapping_done = False

        if not counting_only:
            # mapping already done and no overwriting
            if self.mapping_done and not self.force:
                logging.info("Mapping already done for sample: #{@FSampleName}")
                logging.info("Skipped !")
            else:
                self.LaunchMapping()
        if not mapping_only:
            self.LaunchCounting()
        logging.info("Done !\nJob finished without errors ...")
