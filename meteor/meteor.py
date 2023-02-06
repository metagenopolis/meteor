#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Meteor - A plateform for quantitative metagenomic profiling of complex ecosystems"""

__version__  = "4.3"
__copyright__ = "GPLv3"
__date__ = "2022"

#-------------------------- MODULES IMPORTATION -------------------------------#
import sys
import os
from argparse import ArgumentParser
from configparser import ConfigParser
import logging
from logging.handlers import RotatingFileHandler
import datetime
import glob
import re
import subprocess
from tempfile import TemporaryDirectory
import gzip
import bz2
import lzma
# import shutil
from dataclasses import dataclass, field
from types import NoneType
# from typing import List
# import urllib.request
# from tqdm import tqdm
# import pathlib

#---------------------------- CLASS DEFINITION --------------------------------#
class color:
    BOLD = '\033[1m'
    END = '\033[0m'

# class DownloadProgressBar(tqdm):
#     def update_to(self, b=1, bsize=1, tsize=None):
#         if tsize is not None:
#             self.total = tsize
#         self.update(b * bsize - self.n)


@dataclass
class MeteorMapper:
    FMappedCensusIniFileName: str
    aReferenceIniFileName: str
    FSampleDir: str
    FNGSLibraryIndexerReport: dict
    FMappingProgram: str
    FLibraryCensusIniFile: ConfigParser = field(default_factory=ConfigParser)
    FReferenceIniFile: ConfigParser = field(default_factory=ConfigParser)
    FLibraryMappingDir: str = ""
    FTmpLibraryMappingDir: str = ""
    FLibraryName: str = ""
    FMismatchesCount: int = 0
    FMatchesCount: int = 0
    FIsMismatchesPercentage: float = 0.0
    FReferenceName: str = ""
    FMappingFileFormat: str = ""
    FParametersShortLine: str = ""
    FMapperCmd: str = ""
    FMappingOutputFileNames: list = field(default_factory=list)
    FTmpMappingOutputFileNames: list = field(default_factory=list)
    FIsReadyForMapping: int = 0


    def __post_init__(self):
        self.FReferenceIniFile.read_file(open(self.aReferenceIniFileName))
        self.FLibraryName = self.FLibraryCensusIniFile["sample_info"]["full_sample_name"]
        self.FReferenceName = self.FReferenceIniFile['reference_info']['reference_name']

    # def FinalizeMapping(self):
    #     # Create new mapping_census_ini_file from library_census_ini_file
    #     aMappedCensusIniFile = FLibraryCensusIniFile.copy()
    #     aMappedCensusIniFile.filename = FMappedCensusIniFileName

    #     aSampleInfoSection = aMappedCensusIniFile["sample_info"]    # reference
    #     aMappingSection = aMappedCensusIniFile['mapping']           # reference
    #     aMappingFileSection = aMappedCensusIniFile["mapping_file"]  # reference

    #     aMappedCensusIniFile['mapping']["mapping_tool"] = FMappingProgram
    #     aMappedCensusIniFile['mapping']["mapping_tool_version"] = 'NA' # WTF

    #     aMappedCensusIniFile["mapping_file"]['mapping_file_count'] = len(FMappingOutputFileNames)

    #     # should iterate only once
    #     for output_names in FMappingOutputFileNames:
    #         aMappedCensusIniFile["mapping_file"]['bowtie_file'+"_{}".format(i+1)] = os.path.basename(output_names)
    #     #TODO aMappedCensusIniFile is useless for now
    #     aSampleInfoSection["indexed_read_length"] = FNGSLibraryIndexerReport["IndexedReadLength"]
    #     aSampleInfoSection["sequenced_read_count"] = FNGSLibraryIndexerReport["IndexedReadCount"]
    #     aSampleInfoSection["indexed_sequenced_read_count"] = FNGSLibraryIndexerReport["IndexedReadCount"]
    #     aSampleInfoSection["indexed_sequenced_base_count"] = FNGSLibraryIndexerReport["IndexedBaseCount"]
    #     aSampleInfoSection["database_type"] = "unknown" # WTF
    #     aSampleInfoSection["is_data_prefixed"] = 0  # to avoid too long path

    #     aSampleInfoSection["sample_indexing_software"] =  "Meteor"
    #     aSampleInfoSection["sample_indexing_software_version"] = "3.3"

    #     aSampleInfoSection["census_status"] = 1
    #     aMappingSection["mapping_date"] = datetime.date.today().strftime('%Y-%m-%d')
    #     aMappingSection['reference_name'] = FReferenceName
    #     aMappingSection["mapping_cmdline"] = FMapperCmd
    #     aMappingSection["parameters"] = FParametersShortLine;
    #     aMappingSection["mapped_read_length"] = FNGSLibraryIndexerReport["MappedReadLength"]
    #     aMappingSection["mapped_read_length_type"] = FNGSLibraryIndexerReport["MappedReadLengthType"]
    #     aMappingSection["mismatches"] = FMismatchesCount
    #     aMappingSection["is_mismatches_percentage"] = FIsMismatchesPercentage
    #     aMappingSection["matches"] = FMatchesCount
    #     aMappingSection["is_local_mapping"] = FIsLocalMapping # LOCAL

    #     aMappingSection["mapping_software"] = "Meteor"
    #     aMappingSection["mapping_software_version"] = "3.2"

    #     aMappingSection["processed_read_count"] = FNGSLibraryIndexerReport["IndexedReadCount"]

    #     aMappingFileSection["mapping_file_format"] = FMappingFileFormat
    #     with open(FMappedCensusIniFileName, 'wt') as config_file:
    #         aMappingSection.write(config_file)


    def Bowtie2MapRead(self):
        # prepare bowtie2 command
        # -q --no-head --no-sq --no-unal --omit-sec-seq --end-to-end --sensitive -k 1 -S $sample/$sample.sam -U $fastq_list -p 40 -x $bt2_index
        # references repository on IBS server: /services/projects/biodatabank/meteor-data/reference

        # aBowtieProgram = None

        # since bowtie-build 64bits
        # aIsLargeReference = sekf.FReferenceIniFile["bowtie2_index"]["is_large_reference"]
        # get around KCL feedback (inifile bug?)
        # if not str(aIsLargeReference).isdigit():
        #     aIsLargeReference = 1
        # else:
        #     aIsLargeReference = int(aIsLargeReference)
        # if aIsLargeReference == 1:
        #     aBowtieProgram = C_BOWTIE2_LARGE_MAPPER
        # else:
        #     aBowtieProgram = C_BOWTIE2_SMALL_MAPPER

        # aBowtieIndexCount should be 1
        aBowtieIndexList = os.path.dirname(self.aReferenceIniFileName) + os.sep + self.FReferenceIniFile["reference_file"]["fasta_dir"] + os.sep +self.FReferenceIniFile["bowtie2_index"][f"dna_space_bowtie_index_prefix_name_1"]
        # aBowtieIndexCount = self.FReferenceIniFile["bowtie2_index"]["index_count"]
        # aBowtieIndexCount = 1
        # aBowtieIndexList.append()
        # get around KCL feedback (inifile bug?) : loading as String instead of Fixnum => force .to_i
        # if not str(aBowtieIndexCount).isdigit():
        #     aBowtieIndexCount = 1
        # else:
        #     aBowtieIndexCount = int(aBowtieIndexCount)

        # for aIndex in range(1, aBowtieIndexCount+1):
            # aBowtieIndexList.append(FReferenceIniFile["bowtie2_index"][f"dna_space_bowtie_index_prefix_name_{aIndex}"])
        # bowtie2 parameters
        aParameters = f"--mm -p {self.FCPUCount} {self.FMapperCmd}"
        # if self.FMapperCmd is None:
            # if FIsLocalMapping == 1:
            #     aParameters += " --local --sensitive-local"
            # else:
            #     aParameters += " --end-to-end --sensitive"
            # if FMatchesCount > 0:
            #     aParameters += f" -k {self.FMatchesCount}"
            # elif FMatchesCount == 0:
            #     aParameters += " -a"
        # else:
        # aParameters += f" {FMapperCmd}"


        # aParameters += " --phred33"

        # complete command with input/output parameters, should iterate only once (because only one index)
        # for iBowtieIndex in range(aBowtieIndexCount):

        # ex: FLibraryMappingDir: /projects/parkinson/maping/KCL_01/mapping_vs_hs_9_9_metahit_Proton1-520-Parkinson_R1a_KCL_01
        #     FLibraryName + : Proton1-520-Parkinson_R1a_KCL_01
        #     iBowtieIndex+1 : 1
        FMappingOutputFileName=f"{self.FLibraryMappingDir}{os.sep}{self.FLibraryName}_1.sam"
        # aBowtie2Command = f"{aBowtieProgram} {aParameters} -x {aBowtieIndexList[iBowtieIndex]} -U {FFastqList} -S {FMappingOutputFileNames[iBowtieIndex]}"
        # execute command
        if not os.path.isfile(self.FNGSLibraryIndexerReport['IndexedfastqFilePath']):
            sys.exit("WTF HErhehreh")
        #shutil.copy(self.FNGSLibraryIndexerReport['IndexedfastqFilePath'], "/Users/aghozlan/workspace/meteor_test/mtf.fastq")
        subprocess.check_call(["bowtie2", aParameters, "--no-head --no-sq --no-unal --omit-sec-seq",  "-x", aBowtieIndexList, "-U", self.FNGSLibraryIndexerReport['IndexedfastqFilePath'], "-S", FMappingOutputFileName])


    def MapRead(self, aLibraryMappingDir, aMapperCmd, aMatchesCount, aMismatchesCount, aIsMismatchesPercentage, aIsLocalMapping, aMappingFileFormat, aIsCPUPercentage, aCPUCount):
        #aTmpDir, 
        # get around KCL feedback (inifile bug?)
        self.FMatchesCount = int(aMatchesCount)
        # FMatchesCount = C_DEFAULT_MAPPING_MATCHES_COUNT if (FMatchesCount == 0)

        # # set CPU count
        # aHostCPUCount = self.FLibraryCensusIniFile["worksession"]["meteor.cpu.count"]

        # # get around KCL feedback bug(?)
        # aCPUCount = -1 if not isinstance(aCPUCount, int) else int(aCPUCount)
        # aIsCPUPercentage = 0 if not isinstance(aIsCPUPercentage, int) else int(aIsCPUPercentage)
        # if aCPUCount == -1:
        #     FCPUCount = max(1, aHostCPUCount - 1)
        # else:
        #     FCPUCount = aCPUCount if aIsCPUPercentage == 1 else aCPUCount
        self.FCPUCount = int(aCPUCount)
        # get around KCL feedback bug(?). Default FMismatchesCount is 3
        self.FMismatchesCount = int(aMismatchesCount)
        # idem
        self.FIsMismatchesPercentage = int(aIsMismatchesPercentage)
        self.FIsLocalMapping = int(aIsLocalMapping)  # LOCAL
        self.FParametersShortLine = f"l{self.FNGSLibraryIndexerReport['MappedReadLength']}-m{self.FMismatchesCount}"
        self.FMappedReadCount = 0
        self.FMappingFileFormat = aMappingFileFormat
        self.FMapperCmd = aMapperCmd
        
        self.FTmpLibraryMappingDir = aLibraryMappingDir
        self.FLibraryMappingDir = aLibraryMappingDir
        # FTmpLibraryMappingDir: directory for temporary mapping SAM files.
        # if aTmpDir is set : "~/scratch/tmp_meteor/md5hash" + "/projects/projectname/mapping/sampleN/mapping_vs_blabla"
        # else :              "/projects/projectname/mapping/sampleN/mapping_vs_blabla"
        # if aTmpDir != "":
        #     FTmpLibraryMappingDir = aTmpDir + os.sep + aLibraryMappingDir 
        
        # create FLibraryMappingDir and FTmpLibraryMappingDir
        os.makedirs(self.FLibraryMappingDir, exist_ok=True)
        os.makedirs(self.FTmpLibraryMappingDir, exist_ok=True)
        
        self.FIsReadyForMapping = 1
        
        # call BowtieMapRead() or Bowtie
        self.Bowtie2MapRead()
        return True

@dataclass
class MeteorSession:
    logger: logging.Logger
    tmp_path: str
    FMeteorJobIniFile: ConfigParser = field(default_factory=ConfigParser)
    FMeteorJobIniFilename: str = "" 
    FProjectDir: str = ""
    FSampleDir: str = ""
    FTmpSampleDir: TemporaryDirectory = TemporaryDirectory()
    FTmpDir: TemporaryDirectory = TemporaryDirectory()
    FNoLock: bool = False
    FForce: bool = False
    FMappingDone: bool = True
    FCountingDone: bool = True
    FCountingTypeList: list = field(default_factory=list)
    FMappingDir: str = ""
    LibraryNames: str = ""
    FLibraryCount: str = ""
    FSampleName: str = ""
    FProjectName: str = ""
    FProjectMappingDir: str = ""
    FLibraryIndexerReport: dict = field(default_factory=dict)
    FLibraryIniFileNames: list = field(default_factory=list)
    FLibraryCensusIniFileNames: list = field(default_factory=list)
    FMainMappingCensusIniFileNames: dict = field(default_factory=dict)
    FExcludedMappingCensusIniFileNames: list = field(default_factory=list)

    def __post_init__(self):
        if self.tmp_path:
            self.FTmpDir.dir = self.tmp_path
            self.FTmpSampleDir.dir = self.tmp_path

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

        self.logger.info("Task main mapping {}".format(iLibrary))
        
        # What is Stage1FileName ?
        aMappedCensusIniFileName = self.FMainMappingCensusIniFileNames[iLibrary]["Stage1FileName"]
        aLibraryMappingDir = self.FMainMappingCensusIniFileNames[iLibrary]["directory"]
        
        aLockedMappedCensusIniFileName = aMappedCensusIniFileName+ ".lock"
        
        # remove lock file if asked
        if self.FNoLock and os.path.exists(aLockedMappedCensusIniFileName):
            self.logger.info("Removing lock file: {}".format(aLockedMappedCensusIniFileName))
            os.remove(aLockedMappedCensusIniFileName)
        # remove census 1 file if asked
        if self.FForce and os.path.exists(aMappedCensusIniFileName):
            self.logger.info("Removing existing census 1 file: {}".format(aMappedCensusIniFileName))
            os.remove(aMappedCensusIniFileName)
        
        if  os.path.exists(aMappedCensusIniFileName) or os.path.exists(aLockedMappedCensusIniFileName):
            self.logger.info("Skipping library.\n{} already exists or is locked".format(aMappedCensusIniFileName))

        # create lock file in result directory
        os.makedirs(aLibraryMappingDir, exist_ok=True)
        open(aLockedMappedCensusIniFileName, 'a').close()
        
        aMappingProg = aWorkSessionSection["meteor.mapping.program"]
        if (aMappingProg != 'bowtie2'): 
            sys.exit("Error, unknown or unsupported mapping program : {}".format(aWorkSessionSection["meteor.mapping.program"]))
        aMeteorMapper = MeteorMapper(
                            FMappingProgram=aMappingProg, FLibraryCensusIniFile=self.ini_files[iLibrary], 
                            aReferenceIniFileName= self.aReferenceIniFileName, 
                            FNGSLibraryIndexerReport=self.FLibraryIndexerReport, 
                            FMappedCensusIniFileName=aMappedCensusIniFileName, 
                            FSampleDir=self.FSampleDir)
        aOkMappingProcess = aMeteorMapper.MapRead(
                        #self.FTmpDir,
                        aLibraryMappingDir,
                        aReferenceSection["meteor.mapper.cmd"],
                        aReferenceSection["meteor.matches"],
                        aReferenceSection["meteor.mismatches"],
                        aReferenceSection["meteor.is.perc.mismatches"],
                        aReferenceSection["meteor.is.local.mapping"],
                        aWorkSessionSection["meteor.mapping.file.format"],
                        aWorkSessionSection["meteor.is.cpu.percentage"],
                        aWorkSessionSection["meteor.cpu.count"])
        # remove lock file
        if (os.path.exists(aLockedMappedCensusIniFileName)):
            os.remove(aLockedMappedCensusIniFileName)
                        
        if not aOkMappingProcess:
            return False
        return True

    def CountAndReIndexReads(self, aLibraryCensusIniFile: ConfigParser):
        aSampleInfoSection = aLibraryCensusIniFile["sample_info"]
        aSampleFileSection = aLibraryCensusIniFile["sample_file"]

        aInputFile = self.FSampleDir + os.sep + aSampleFileSection["fastq_file"]
        aOutputFile = self.FTmpSampleDir.name + os.sep + aSampleFileSection["fastq_file"] + '.idx'
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

    def LaunchMapping(self):
        self.logger.info("Launch mapping")
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
            
            self.logger.info("Meteor Mapping task description")
            self.logger.info("Sample name = " + FSampleName)
            self.logger.info("Library name = " + aSampleLibraryName)
            self.logger.info("Project name = " + FProjectName)
            self.logger.info("Sequencing device = " + aSampleInfoSection["sequencing_device"])
            self.logger.info("Workflow = " + os.path.basename(iLibrary))
            
            if not os.path.exists(iLibrary+".lock"):
                ### reindexing this library reads and fill FLibraryIndexerReport
                self.CountAndReIndexReads(aLibraryCensusIniFile)
                # MAPPING THIS LIBRARY ON MAIN REFERENCE
                aOKToContinue = self.TaskMainMapping(iLibrary)
                if not aOKToContinue:
                    sys.exit("Error, TaskMainMapping failed: " + iLibrary)
                
                # MAPPING THIS LIBRARY ON EACH EXCLUDED REFERENCE
                # for iExcluded in range(aExcludedRefCount):
                #     aOKToContinue = TaskExcludedMapping(iLibrary, iExcluded)
                #     if not aOKToContinue:
                #         sys.exit("Error, TaskExcludedMapping failed: " + iLibrary)


    def LaunchCounting(self):
        self.logger.info("Launch counting")
        # does not need census_stage_0.ini
        #-w /path/to/workflow_tutorial.ini -i /path/to/sample/H1 -p /path/to/project_name -m mapping
        aparameters = " -w " + self.FMeteorJobIniFilename + " -i " + self.FSampleDir + " -p " + self.FProjectDir + " -o " + os.path.basename(self.FProjectMappingDir)
        if len(self.FCountingTypeList) > 0:
            aparameters += " -c " + ",".join(self.FCountingTypeList)
        if self.FForce: # force overwriting former profiling results done with same parameters
            aparameters += " -f"
        if self.FTmpDir: # path to directory where mapping results (SAM files) are stored. #### NEW
            aparameters += " -t " + self.FTmpDir.name
        # subprocess.check_call([os.path.dirname(os.path.realpath(__file__)) + os.sep + "src" + os.sep + "build" + os.sep + "meteor-counter", aparameters])
        subprocess.check_call([os.path.dirname(os.path.realpath(__file__)) + os.sep + "src" + os.sep + "build" + os.sep + "meteor-counter", "-w", \
         self.aReferenceIniFileName, "-i", self.FSampleDir, "-p", self.FProjectDir, "-o", os.path.basename(self.FProjectMappingDir), \
         "-c", ",".join(self.FCountingTypeList),  "-t", self.FTmpDir.name])

        # ### Execute Command
        # print("\nExecuting counting command:")
        # print(aCmd)
        # aExitStatus = ExecuteCmd(aCmd)
        
        # if aExitStatus != 0:
        #     sys.stderr.write("Error in MeteorSession::LaunchCounting, counting command exited with non zero status")
        #     sys.exit(aExitStatus)
            
        # return aExitStatus


    def ProcessJob(self, workflow_ini, project_path, input_path, remove_lock, force, mapping_dir,
                   counting_type, counting_only, mapping_only): #options : Workflow, ProjectPath, InputPath, MappingBasename
        # directly from program arguments
        self.FProjectDir = project_path
        self.FSampleDir  = input_path
        # TODO get it from census ini file, (FProjectName too)
        self.FSampleName = os.path.basename(input_path)
        # self.FMappingDone = True
        # self.FCountingDone = True
        self.FCountingTypeList = counting_type
        
        # self.FForce = False
        self.FNoLock = remove_lock 
        if force:
            self.FForce = True
            self.FNoLock = True
        # ex: /projects/parkinson/mapping
        self.FProjectMappingDir = project_path + os.sep + mapping_dir

        # load information from workflow.ini ; set FMeteorJobIniFile
        # Check excluded reference
        self.LoadJobWorkflow(workflow_ini)

        #### in PrepareWorkSpace ?
        # ex: /projects/parkinson/mapping/KCL_01
        self.FSampleMappingDir = self.FProjectMappingDir + os.sep + self.FSampleName
        # make directory (and parent) if needed
        #if not os.path.exists(FSampleMappingDir):
        #    os.makedirs(FSampleMappingDir)  #### maybe obsolete cause done in TaskMainMapping and/or TaskExcludedMapping
        
        #### in PrepareLibraryIniFileNames ?
        self.FLibraryIniFileNames = sorted(glob.glob(self.FSampleDir+os.sep+"*_"+"census_stage_0.ini"))
        print(self.FLibraryIniFileNames)
        self.FLibraryCount = len(self.FLibraryIniFileNames)
        if self.FLibraryCount == 0:
            sys.exit("Error, no census_stage_0.ini file found in {}".format(self.FSampleDir))
        print("\nLaunch mapping\n")
        aMainReferenceSection = self.FMeteorJobIniFile["main_reference"] # alias
        aMappingDirPathPrefix = ""
        if aMainReferenceSection["meteor.mapping.prefix.name"] is None:
            aMappingDirPathPrefix = "mapping"+"_vs_"+aMainReferenceSection["meteor.reference.name"]+"_l"+ \
                                    "{}-{}".format(aMainReferenceSection["meteor.mapped.readlength"], "m")+aMainReferenceSection["meteor.mismatches"]+"_"
        else:
            aMappingDirPathPrefix = "{}_".format(aMainReferenceSection["meteor.mapping.prefix.name"])

        aExcludedRefCount = int(self.FMeteorJobIniFile["worksession"]["meteor.excluded.reference.count"])

        # # LOOP ON EACH LIBRARY
        self.ini_files = {}
        for iLibrary in self.FLibraryIniFileNames:
            # Check if lock file exist
            if os.path.exists(iLibrary+".lock") and not self.FNoLock:
                sys.exit("Error, lock file not found: {}".format(iLibrary+".lock"))
            # MAPPING THIS LIBRARY ON MAIN REFERENCE
            aLibraryCensusIniFile = ConfigParser()
            aLibraryCensusIniFile.read_file(open(iLibrary))
            self.ini_files[iLibrary] = aLibraryCensusIniFile

            aSampleInfo = aLibraryCensusIniFile["sample_info"] # reference
            aMappingDirPath = aMappingDirPathPrefix + aSampleInfo["full_sample_name"]
            directory = self.FProjectMappingDir + os.sep + aSampleInfo["sample_name"]+os.sep+aMappingDirPath
            self.FMainMappingCensusIniFileNames[iLibrary] = {
                "directory": directory,
                "Stage1FileName": directory + os.sep + os.path.basename(iLibrary).replace("census_stage_0", "census_stage_1")
            }
            if (not os.path.exists(self.FMainMappingCensusIniFileNames[iLibrary]["Stage1FileName"])):
                self.FMappingDone = False

        if not counting_only:
            # mapping already done and no overwriting
            if self.FMappingDone and not self.FForce:
                self.logger.info("Mapping already done for sample: #{@FSampleName}")
                self.logger.info("Skipped !")
            else: # mapping not done or we want to overwrite previous results
                # if option -t not provided, generate sample tmpdir in sample dir
                #@FTmpSampleDir = @FTmpDir.nil? ? @FSampleDir + C_PATH_SEP + aMD5 : @FTmpDir #### NEW
                #self.FTmpSampleDir = self.FSampleMappingDir + C_PATH_SEP + aMD5 : @FTmpDir #### NEW
                # if File.exists?(@FTmpSampleDir)
                #     STDERR.puts "WARNING, temporary dir #{@FTmpSampleDir} already exists !"
                # end
                # FileUtils.mkdir_p(@FTmpSampleDir)
                self.LaunchMapping()
                # remove sample tmp data if created in sample dir
                # FileUtils.rm_rf(@FTmpSampleDir) if ( @FTmpDir.nil? and File.exists?(@FTmpSampleDir) ) #### NEW
        if not mapping_only:
            self.LaunchCounting()
        self.logger.info("Done !\nJob finished without errors ...")
        #         # MAPPING THIS LIBRARY ON EACH EXCLUDED REFERENCE
        #         for iExcluded in range(aExcludedRefCount):
        #             aOKToContinue = TaskExcludedMapping(iLibrary, iExcluded)
        #             if not aOKToContinue:
        #                 sys.stderr.write("Error, TaskExcludedMapping failed: {}".format(FLibraryIniFileNames[iLibrary]))
        #                 sys.exit(1)

        #         if aSampleInfoSection[C_SEQUENCING_DEVICE_STR] == C_SOLID_DEVICE:
        #             # CGH analysis
        #             aOKToContinue = TaskCGH(iLibrary)
        #             if not aOKToContinue:
        #                 sys.stderr.write("Error, TaskCGH failed: {}".format(FLibraryIniFileNames[iLibrary]))
        #          # COPYING INI FILE TO MAPPING DIRECTORY
        #         shutil.copy(FLibraryIniFileNames[iLibrary], FProjectMappingDir)
        #         shutil.copy(FLibraryIniFileNames[iLibrary]+".lock", FProjectMappingDir)

        #     aOKToContinue = LaunchCounting()
        #     if aOKToContinue != 0:
        #         sys.stderr.write("Error in MeteorSession::ProcessJob, LaunchCounting exited with non zero status")
        #         sys.exit(aOKToContinue)

        #     # COPYING INI FILE TO PROJECT DIRECTORY
        #     shutil.copy(FMeteorJobIniFile.filename, FProjectDir)


    # def TaskExcludedMapping(self, iLibrary: ConfigParser, iExcluded) -> bool:
    #     aWorkSessionSection = FMeteorJobIniFile["worksession"]  # reference
    #     aReferenceSection = FMeteorJobIniFile["excluded_reference_" + f"{iExcluded+1}"]  # reference

    #     aReferenceIniFileName = (
    #         aWorkSessionSection["meteor.reference.dir"]
    #         + os.sep
    #         + aReferenceSection["meteor.reference.name"]
    #         + os.sep
    #         + aReferenceSection["meteor.reference.name"]
    #         + "_reference.ini"
    #     )

    #     print(f"\n######### TaskExcludedMapping({iLibrary}, {iExcluded+1}); reference name: " + aReferenceSection["meteor.reference.name"])

    #     aMappedCensusIniFileName = self.FExcludedMappingCensusIniFileNames[iLibrary][iExcluded]["Stage1FileName"]
    #     aLibraryMappingDir = self.FExcludedMappingCensusIniFileNames[iLibrary][iExcluded]["directory"]

    #     aLockedMappedCensusIniFileName = aMappedCensusIniFileName + ".lock"

    #     # remove lock file if asked
    #     if self.FNoLock and os.path.exists(aLockedMappedCensusIniFileName):
    #         print(f"\nINFO: removing lock file: {aLockedMappedCensusIniFileName}")
    #         os.remove(aLockedMappedCensusIniFileName)
    #     # remove census 1 file if asked
    #     if self.FForce and os.path.exists(aMappedCensusIniFileName):
    #         print(f"\nINFO: removing existing census 1 file: {aMappedCensusIniFileName}")
    #         os.remove(aMappedCensusIniFileName)

    #     if os.path.exists(aMappedCensusIniFileName) or os.path.exists(aLockedMappedCensusIniFileName):
    #         print(f"\nINFO: skipping library.\nINFO: {aMappedCensusIniFileName} already exists or is locked")
    #         return True

    #     # create lock file in result directory
    #     os.makedirs(aLibraryMappingDir, exist_ok=True)
    #     open(aLockedMappedCensusIniFileName, "a").close()

    #     aMappingProg = aWorkSessionSection["meteor.mapping.program"]
    #     if aMappingProg != 'bowtie2':
    #         raise ValueError(f"Error, unknown or unsupported mapping program:"+ aWorkSessionSection["meteor.mapping.program"])

    #     aMeteorMapper = MeteorMapper(
    #         aMappingProg,
    #         FLibraryIniFileNames[iLibrary],
    #         aReferenceIniFileName,
    #         FLibraryIndexerReport,
    #         aMappedCensusIniFileName
    #     )

    #     aOkMappingProcess = aMeteorMapper.MapRead(
    #         #FTmpDir,
    #         aLibraryMappingDir,
    #         aReferenceSection["meteor.mapper.cmd"],
    #         aReferenceSection["meteor.matches"],
    #         aReferenceSection["meteor.mismatches"],
    #         aReferenceSection["meteor.is.perc.mismatches"],
    #         False,
    #         aWorkSessionSection["meteor.mapping.file.format"],
    #         aWorkSessionSection["meteor.is.cpu.percentage"],
    #         aWorkSessionSection["meteor.cpu.count"]
    #     )
    #     # remove lock file
    #     if os.path.exists(aLockedMappedCensusIniFileName):
    #         os.remove(aLockedMappedCensusIniFileName)
    #     return True


    def LoadJobWorkflow(self, aWorkflowFile: str):
        if not os.path.exists(aWorkflowFile):
            sys.exit(f"Error, file {aWorkflowFile} not found")
        self.FMeteorJobIniFile.read_file(open(aWorkflowFile))
        self.FMeteorJobIniFilename = aWorkflowFile
        # check if excluded reference count is correct
        # aExcludedRefCount = FMeteorJobIniFile["worksession"]["meteor.excluded.reference.count"]
        # if not isinstance(aExcludedRefCount, int):
        #     aExcludedRefCount = False
        # else:
        #     aExcludedRefCount = int(aExcludedRefCount)
        
        # print(FMeteorJobIniFile)
        aExcludedRefCounted = len(
             [ s for s in self.FMeteorJobIniFile.sections() if s.startswith("excluded_reference_")]
        )

        # if aExcludedRefCount != aExcludedRefCounted:
        #     print(f"Error in file {aWorkflowFile}:", file=sys.stderr)
        #     print(
        #         f"meteor.excluded.reference.count is set to {aExcludedRefCount}",
        #         file=sys.stderr,
        #     )
        #     print(
        #         f"but found {aExcludedRefCounted} excluded_reference section(s) ...",
        #         file=sys.stderr,
        #     )
        #     sys.exit(1)


#-------------------------- FUNCTIONS DEFINITION ------------------------------#

def isfile(path: str): # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file
    
    :raises ArgumentTypeError: If file doesn't exist
    
    :return: (str) Path 
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return os.path.abspath(path)


def isdir(path: str): # pragma: no cover
    """Check if path is an existing file.
      
    :param path: Path to the directory

    :raises ArgumentTypeError: If directory doesn't exist

    :return: (str) Path 
    """
    if not os.path.isdir(path):
        if os.path.isfile(path):
            msg = "{0} is a file".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return os.path.abspath(path) + os.sep


# def download_url(url:str, output_path: str):
#     with DownloadProgressBar(unit='B', unit_scale=True,
#                              miniters=1, desc=url.split('/')[-1]) as t:
#         urllib.request.urlretrieve(url, filename=output_path, reporthook=t.update_to)


def get_log(path_log:str) ->logging.Logger:
    """Start logging streams

    :param path_log:  Path to the log file output

    :return: (logging.Logger) Return logging object
    """
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s')
    console_formatter = logging.Formatter('[%(levelname)s]: %(message)s')
    # Create log file
    file_handler = RotatingFileHandler(path_log, 'a', 1000000, 1)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    # Stream in the the console
    ## TO REMOVE IF daemon
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(console_formatter)
    stream_handler.setLevel(logging.INFO)
    logger.addHandler(stream_handler)
    return logger


def get_arguments()->ArgumentParser: # pragma: no cover
    """
    Meteor help and arguments
    
    No arguments
    
    Return : parser
    """
    parser = ArgumentParser(description=color.BOLD + __doc__ + color.END)
    subparsers = parser.add_subparsers(title = 'positional arguments', help="Select activity", dest= "command")
    # # Mappping commands
    download_parser = subparsers.add_parser('download',
        help='Download catalog')
    reference_parser = subparsers.add_parser('build',
        help='Index reference')
    reference_parser.add_argument("-i", dest='input_fasta_file',
        type = str, required = True, help = "Input fasta filename.")
    reference_parser.add_argument("-p", dest='reflogging.LoggerDir',
        type = str, required = True, help = "Output path of the reference repository.")
    reference_parser.add_argument("-n", dest='refName',
        type = str, required = True, help = "Name of the reference (ansi-string without space).")
    reference_parser.add_argument("-t", dest='threads', default=4,
        type = int, help = "Thread count for bowtie2 (if available).")
    fastq_parser = subparsers.add_parser('fastq',
        help='Import fastq files')
    fastq_parser.add_argument("-i", dest='fastq_dir',
        type = isdir, required = True, help = """Path to a directory containing all input FASTQ files.
FASTQ files must have the extension .FASTQ or .FQ.
For paired-ends files must be named : 
    file_R1.[fastq/fq] & file_R2.[fastq/fq] 
                       or
    file_1.[fastq/fq] & file_2.[fastq/fq]""")
    fastq_parser.add_argument("-p", dest='project_name',
        type = str, required = True, help = "Project name (ansi-string without space).")
    fastq_parser.add_argument("-m", dest='mask_sample_name',
        type = str, required = True, help = "Regular expression for extracting sample name.")
    fastq_parser.add_argument("-t", dest='techno', choices = ["illumina", "proton"],
        default="Illumina", type = str, help = "Sequencing technology (default: Illumina).")
    fastq_parser.add_argument("-c", dest='iscompressed',
        default=False, action="store_true", help = "Fastq files are compressed.")
    fastq_parser.add_argument("-d", dest='isdispatched',
        default=False, action="store_true", help = "Fastq files are already dispatched in directories.")
    mapping_parser = subparsers.add_parser('mapping',
        help='Map reads against a gene catalog')
    mapping_parser.add_argument("-w", dest="workflow_ini", type=isfile, required=True,
        help="Path to meteor configuration file, e.g. workflow.ini")
    mapping_parser.add_argument("-i", dest="input_path", type = isdir, required = True,
        help = """Path to sample directory, containing the sample sequencing metadata
        (files ending with _census_stage_0.ini)""")
    mapping_parser.add_argument("-p", dest="project_path", type=isdir, required=True,
        help="path to project directory, containing mapping and profile data (e.g. /projects/project_dir)")
    mapping_parser.add_argument("-o", dest="mapping_dir", type=str, required=True,
        help="basename of the mapping directory")
    mapping_parser.add_argument("-c", dest="counting_type", type=str, nargs='+', default=["smart_shared_reads"],
        choices = ["total_reads", "shared_reads", "smart_shared_reads", "unique_reads"],
    help="""counting type string (default smart_shared_reads.""")
    # mapping_parser.add_argument("-s", dest="strand", type=str, choices=["direct", "undirect"],
    #     help="(where direct/undirect means direct/undirect strand)")
    # mapping_parser.add_argument("-a", dest="coverage", type=str, choices=["coverage", "avg_coverage"],
    #     help="coverage or average coverage")
    mapping_parser.add_argument("-t", dest="tmp_path", type=isdir,
        help="path to the directory where temporary files (e.g. sam) are stored")
    mapping_parser.add_argument("-e", dest="clean_tmp_path", action="store_true",
        help="Remove the temporary files directory after the task counting is finished")
    mapping_parser.add_argument("-l", dest="remove_lock", action="store_true",
        help="ignore and remove any lock file")
    mapping_parser.add_argument("-f", dest="force", action="store_true",
        help="force overwriting, ignoring previous results and lock files")
    mapping_parser.add_argument("-m", dest="mapping_only", action="store_true",
        help="execute mapping only")
    mapping_parser.add_argument("-n", dest="counting_only", action="store_true",
        help="execute counting only")
    return parser.parse_args(args=None if sys.argv[2:] else ['--help'])



#==============================================================
# Main program
#==============================================================
def main()->None: # pragma: no cover
    """
    Main program function
    """
    # Let's logging
    now = datetime.datetime.now()
    path_log = "meteor_" + now.strftime("%Y%m%d_%H%M") + ".log"
    # Get arguments
    args = get_arguments()
    logger = get_log(path_log)
    # Import FASTQ
    if args.command == "fastq":
        if args.isdispatched:
            fastq_file_list = glob.glob(args.fastq_dir + "*" + os.sep + "*.f*q*")
        else:
            print(args.fastq_dir + "*.{fq,fastq}*")
            fastq_file_list = glob.glob('{}{}'.format(args.fastq_dir, "*.f*q*"))
        ext_R1 = ("1.fq", "1.fastq", "R1.fastq", "R1.fastq.gz", "R1.fq.gz")
        ext_R2 = ("2.fq", "2.fastq", "R2.fastq", "R2.fastq.gz", "R2.fq.gz")
        for fastq_file in fastq_file_list:
            # print("Import ", fastq_file)
            logger.info("Import ", fastq_file)
            full_sample_name = os.path.basename(fastq_file)
            if args.iscompressed:
                full_sample_name = ".".join(full_sample_name.split(".")[:-1])
            # Extract paired-end info
            tag = "single"
            if full_sample_name.endswith(ext_R1):
                tag = "1"
            elif full_sample_name.endswith(ext_R2):
                tag = "2"
            full_sample_name = ".".join(full_sample_name.split(".")[:-1])
            if args.isdispatched:
                sample_name = os.path.basename(os.path.dirname(fastq_file))
            else:
                try:
                    # split full sample name (in fact library/run name) in order to extract sample_name according to regex mask
                    full_sample_name_array = re.search(args.mask_sample_name, full_sample_name)
                    # risk of TypeError: 'NoneType' object is not subscriptable
                    sample_name = full_sample_name_array[0]
                except NoneType:
                    sys.exit("None matching pattern: {}".format(args.mask_sample_name))
        
            if args.isdispatched:
                sample_dir = os.path.dirname(fastq_file) + os.sep
            else:
                # create directory for the sample and move fastq file into
                sample_dir = os.path.dirname(fastq_file) + os.sep + sample_name + os.sep
                if not os.path.isdir(sample_dir):
                    os.makedirs(sample_dir)
                os.rename(fastq_file, sample_dir + os.path.basename(fastq_file))
            config = ConfigParser()
            config["sample_info"] = {
                "sample_name": sample_name,
                "condition_name": "NA", # What is this ?
                "project_name" : args.project_name,
                "sequencing_date": "1900-01-01", # Then it is useless
                "sequencing_device": args.techno,
                "census_status": 0, # what is this ?
                "read_length": -1, # Then it is useless
                "tag": tag,
                "full_sample_name": full_sample_name
            }
            config["sample_file"] = {
                "fastq_file": os.path.basename(fastq_file),
                "is_compressed": args.iscompressed
            }
            with open(sample_dir + full_sample_name + "_census_stage_0.ini", 'w') as configfile:
                config.write(configfile)
    # Import reference
    elif args.command ==  "build":
        logger.info("Import ", args.refName)
        # Create reference genome directory if it does not already exist
        ref_dir = os.path.join(args.reflogging.LoggerDir, args.refName)
        if not os.path.exists(ref_dir):
            os.makedirs(ref_dir)

        # Create subdirectories for fasta files and reference indices
        fasta_dir = os.path.join(ref_dir, "fasta")
        if not os.path.exists(fasta_dir):
            os.makedirs(fasta_dir)

        database_dir = os.path.join(ref_dir, "database")
        if not os.path.exists(database_dir):
            os.makedirs(database_dir)
        # Read input fasta file and create new fasta file for each chromosome or contig
        output_annotation_file = os.path.join(database_dir, args.refName + '_lite_annotation')
        output_fasta_file = os.path.join(fasta_dir, args.refName + '.fasta')
        with open(args.input_fasta_file, 'rt') as input_fasta:
            with open(output_annotation_file, "wt") as output_annotation:
                with open(output_fasta_file, "wt") as output_fasta:
                    gene_id = 1
                    for line in input_fasta:
                        if line.startswith(">"):
                            output_annotation.write(line[1:].split(" ")[0].strip() +"\n")
                            output_fasta.write(">" + str(gene_id) + "\n")
                            gene_id += 1
                        else:
                            output_fasta.write(line)
        # Generate configuration file for reference genome
        config = ConfigParser()
        config["reference_info"] = {
            "reference_name": args.refName,
            "entry_type": "fragment", # Why ?
            "reference_date": datetime.datetime.now().strftime("%Y%m%d"),
            "database_type": "text",
            "HAS_LITE_INFO": 1
        }

        config["reference_file"] = {
            #IS_LARGE_REFERENCE_STR: 1,
            "database_dir": "database",
            "fasta_dir": "fasta",
            #"fasta_file_count": 1,
            # is it possible to have several fasta
            "fasta_file_count": args.refName + '.fasta'
        }

        index_prefix = os.path.join(database_dir, args.refName)
        subprocess.check_call(["bowtie2-build", '-f', '-t', str(args.threads), 
            output_fasta_file, os.path.join(fasta_dir, args.refName)])
        config["bowtie2_index"] = {
            # "is_large_reference": "1", # WTF 
            "is_DNA_space_indexed": 1,
            "dna_space_bowtie_index_prefix_name_1": args.refName
        }

        # Write configuration file
        with open(os.path.join(ref_dir, args.refName + '_reference.ini'), 'wt') as config_file:
            config.write(config_file)
    elif args.command == "mapping":
        m= MeteorSession(logger, args.tmp_path)
        # counting_type = args.counting_type
        # if args.strand:
        #     counting_type = args.strand + "_" + counting_type
        # if args.coverage:
        #     counting_type = counting_type + "_" + coverage
        m.ProcessJob(args.workflow_ini,
        args.project_path,
        args.input_path,
        args.remove_lock,
        args.force,
        args.mapping_dir,
        args.counting_type,
        args.counting_only,
        args.mapping_only)


if __name__ == '__main__':
    main()