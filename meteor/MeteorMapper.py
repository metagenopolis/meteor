from subprocess import check_call
from configparser import ConfigParser
from dataclasses import dataclass, field
from pathlib import Path
from MeteorSession import session

@dataclass
class Mapper(session):
    ini_files: list
    # FMappedCensusIniFileName: str
    # aReferenceIniFileName: str
    # FSampleDir: str
    # FNGSLibraryIndexerReport: dict
    # FMappingProgram: str
    # FLibraryCensusIniFile: ConfigParser = field(default_factory=ConfigParser)
    # FReferenceIniFile: ConfigParser = field(default_factory=ConfigParser)
    # FLibraryMappingDir: Path = field(default_factory=list)
    # FTmpLibraryMappingDir:  Path = field(default_factory=list)
    # FLibraryName: str = ""
    # FMismatchesCount: int = 0
    # FMatchesCount: int = 0
    # FIsMismatchesPercentage: float = 0.0
    # FReferenceName: str = ""
    # FMappingFileFormat: str = ""
    # FParametersShortLine: str = ""
    # FMapperCmd: str = ""
    # FTmpMappingOutputFileNames: list = field(default_factory=list)
    # FIsReadyForMapping: int = 0
    # FCPUCount: int = 1
    # FIsLocalMapping: bool = False
    # FMappedReadCount: int = 0

    # def __post_init__(self)->None:
    #     self.FReferenceIniFile.read_file(open(self.aReferenceIniFileName))
    #     self.FLibraryName = self.FLibraryCensusIniFile["sample_info"]["full_sample_name"]
    #     self.FReferenceName = self.FReferenceIniFile["reference_info"]["reference_name"]

    # def Bowtie2MapRead(self)->None:
    #     aBowtieIndexList = os.path.dirname(self.aReferenceIniFileName) + os.sep + self.FReferenceIniFile["reference_file"]["fasta_dir"] + os.sep +self.FReferenceIniFile["bowtie2_index"][f"dna_space_bowtie_index_prefix_name_1"]
    #     # bowtie2 parameters
    #     aParameters = f"--mm -p {self.FCPUCount} {self.FMapperCmd}"
    #     FMappingOutputFile = self.FLibraryMappingDir / self.FLibraryName / "_1.sam"
    #     # execute command
    #     if not os.path.isfile(self.FNGSLibraryIndexerReport["IndexedfastqFilePath"]):
    #         raise TypeError("Indexer report is missing")
    #     check_call(["bowtie2", aParameters, "--no-head --no-sq --no-unal --omit-sec-seq",  "-x", aBowtieIndexList, "-U", self.FNGSLibraryIndexerReport["IndexedfastqFilePath"], "-S", FMappingOutputFile.name])

    def execute(self)->None:
        # Pairing info
        check_call(["bowtie2", f" -p {self.threads}", 
        "--mm --no-head --no-sq --no-unal --omit-sec-seq",  
        "-x", aBowtieIndexList, "-U", 
        self.FNGSLibraryIndexerReport["IndexedfastqFilePath"], 
        "-S", FMappingOutputFile.name])
        # self.FMatchesCount = aMatchesCount
        # self.FCPUCount = aCPUCount
        # # get around KCL feedback bug(?). Default FMismatchesCount is 3
        # self.FMismatchesCount = aMismatchesCount
        # # idem
        # self.FIsMismatchesPercentage = aIsMismatchesPercentage
        # self.FIsLocalMapping = aIsLocalMapping  # LOCAL
        # self.FParametersShortLine = f"l{self.FNGSLibraryIndexerReport['MappedReadLength']}-m{self.FMismatchesCount}"
        # self.FMappingFileFormat = aMappingFileFormat
        # self.FMapperCmd = aMapperCmd
        # self.FTmpLibraryMappingDir = aLibraryMappingDir
        # self.FLibraryMappingDir = aLibraryMappingDir
        # # create FLibraryMappingDir and FTmpLibraryMappingDir
        # os.makedirs(self.FLibraryMappingDir, exist_ok=True)
        # os.makedirs(self.FTmpLibraryMappingDir, exist_ok=True)
        # self.FIsReadyForMapping = 1
        # call BowtieMapRead() or Bowtie
        # self.Bowtie2MapRead()
        return True