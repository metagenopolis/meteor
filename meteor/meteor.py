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
import argparse
import configparser
import logging
from logging.handlers import RotatingFileHandler
import datetime
import glob
import re
import subprocess
import hashlib
import tempfile
# from dataclasses import dataclass

#---------------------------- CLASS DEFINITION --------------------------------#
class color:
   BOLD = '\033[1m'
   END = '\033[0m'

class MeteorMapper:
    def __init__(self, aMappingProg, aLibraryCensusIniFileName, aReferenceIniFileName, 
                aLibraryIndexerReport, aMappedCensusIniFileName):
        self.FLibraryCensusIniFile = configparser.ConfigParser()
        self.FLibraryCensusIniFile.read_file(open(aLibraryCensusIniFileName))
        self.FMappedCensusIniFileName = aMappedCensusIniFileName
        self.FReferenceIniFile = configparser.ConfigParser()
        self.FReferenceIniFile.read_file(open(aReferenceIniFileName))

        self.FLibraryMappingDir = None
        self.FTmpLibraryMappingDir = None
        self.FLibraryName = self.FLibraryCensusIniFile["sample_info"]["full_sample_name"]
        self.FSampleDir = os.path.dirname(aLibraryCensusIniFileName)

        self.FNGSLibraryIndexerReport = aLibraryIndexerReport # reference to MeteorSession member
        self.FMismatchesCount = None
        self.FMatchesCount = None
        self.FIsMismatchesPercentage = None
        self.FReferenceName = self.FReferenceIniFile['reference_info']['reference_name']

        self.FQualityEncoding = self.FLibraryCensusIniFile["sample_info"]["quality_encoding"]
        self.FSequenceFileFormat = Sff[FQ]
        self.FIsColorSpaceMapping = False

        self.FMappingProgram = aMappingProg
        self.FMappingFileFormat = None
        self.FParametersShortLine = None
        self.FMapperCmd = None
        self.FMappingOutputFileNames = []
        self.FTmpMappingOutputFileNames = []
        self.FIsReadyForMapping = 0


class MeteorSession:
    def __init__(self):
        self.FMeteorJobIniFile = None
        self.FProjectDir = None
        self.FSampleDir = None
        self.FTmpSampleDir = None
        self.FTmpDir = None
        self.FNoLock = False
        self.FForce = False
        self.FMappingDone = True
        self.FCountingDone = True
        self.FCountingTypeList = None
        self.FMappingDir = None
        self.LibraryNames = None
        self.FLibraryCount = None
        self.FSampleName = None
        self.FProjectName = None
        self.FProjectMappingDir = None
        self.FLibraryIndexerReport = {}
        self.FLibraryIniFileNames = []
        self.FLibraryCensusIniFileNames = []
        self.FMainMappingCensusIniFileNames = []
        self.FExcludedMappingCensusIniFileNames = []

    def FinalizeMapping(self):
        # Create new mapping_census_ini_file from library_census_ini_file
        aMappedCensusIniFile = FLibraryCensusIniFile.copy()
        aMappedCensusIniFile.filename = FMappedCensusIniFileName

        aSampleInfoSection = aMappedCensusIniFile["sample_info"]    # reference
        aMappingSection = aMappedCensusIniFile['mapping']           # reference
        aMappingFileSection = aMappedCensusIniFile["mapping_file"]  # reference

        aMappedCensusIniFile['mapping']["mapping_tool"] = FMappingProgram
        aMappedCensusIniFile['mapping']["mapping_tool_version"] = 'NA' # WTF

        aMappedCensusIniFile["mapping_file"]['mapping_file_count'] = len(FMappingOutputFileNames)

        # should iterate only once
        for i in range(len(FMappingOutputFileNames)):
            aMappedCensusIniFile["mapping_file"]['bowtie_file'+"_{}".format(i+1)] = os.path.basename(FMappingOutputFileNames[i])

        #inherited;
        aSampleInfoSection["indexed_read_length"] = FNGSLibraryIndexerReport["IndexedReadLength"]
        aSampleInfoSection["sequenced_read_count"] = FNGSLibraryIndexerReport["IndexedReadCount"]
        aSampleInfoSection["indexed_sequenced_read_count"] = FNGSLibraryIndexerReport["IndexedReadCount"]
        aSampleInfoSection["indexed_sequenced_base_count"] = FNGSLibraryIndexerReport["IndexedBaseCount"]
        aSampleInfoSection["database_type"] = "unknown" # WTF
        aSampleInfoSection["is_data_prefixed"] = 0  # to avoid too long path

        aSampleInfoSection["sample_indexing_software"] =  "Meteor"
        aSampleInfoSection["sample_indexing_software_version"] = "3.3"

        aSampleInfoSection["census_status"] = 1
        aMappingSection["mapping_date"] = datetime.date.today().strftime('%Y-%m-%d')
        aMappingSection['reference_name'] = FReferenceName
        aMappingSection["mapping_cmdline"] = FMapperCmd
        aMappingSection["parameters"] = FParametersShortLine;
        aMappingSection["mapped_read_length"] = FNGSLibraryIndexerReport["MappedReadLength"]
        aMappingSection["mapped_read_length_type"] = FNGSLibraryIndexerReport["MappedReadLengthType"]
        aMappingSection["mismatches"] = FMismatchesCount
        aMappingSection["is_mismatches_percentage"] = FIsMismatchesPercentage
        aMappingSection["matches"] = FMatchesCount
        aMappingSection["is_local_mapping"] = FIsLocalMapping # LOCAL

        aMappingSection["mapping_software"] = "Meteor"
        aMappingSection["mapping_software_version"] = "3.2"

        aMappingSection["processed_read_count"] = FNGSLibraryIndexerReport["IndexedReadCount"]

        aMappingFileSection["mapping_file_format"] = FMappingFileFormat
        with open(FMappedCensusIniFileName, 'wt') as config_file:
            aMappingSection.write(config_file)


    def Bowtie2MapRead(self):
        # prepare bowtie2 command
        # -q --no-head --no-sq --no-unal --omit-sec-seq --end-to-end --sensitive -k 1 -S $sample/$sample.sam -U $fastq_list -p 40 -x $bt2_index
        # references repository on IBS server: /services/projects/biodatabank/meteor-data/reference

        aBowtieProgram = None

        # since bowtie-build 64bits
        aIsLargeReference = FReferenceIniFile["bowtie2_index"]["is_large_reference"]
        # get around KCL feedback (inifile bug?)
        if not str(aIsLargeReference).isdigit():
            aIsLargeReference = 1
        else:
            aIsLargeReference = int(aIsLargeReference)
        if aIsLargeReference == 1:
            aBowtieProgram = C_BOWTIE2_LARGE_MAPPER
        else:
            aBowtieProgram = C_BOWTIE2_SMALL_MAPPER

        # aBowtieIndexCount should be 1
        aBowtieIndexList = []
        aBowtieIndexCount = FReferenceIniFile["bowtie2_index"]["index_count"]
        # get around KCL feedback (inifile bug?) : loading as String instead of Fixnum => force .to_i
        if not str(aBowtieIndexCount).isdigit():
            aBowtieIndexCount = 1
        else:
            aBowtieIndexCount = int(aBowtieIndexCount)

        for aIndex in range(1, aBowtieIndexCount+1):
            aBowtieIndexList.append(FReferenceIniFile["bowtie2_index"][f"dna_space_bowtie_index_prefix_name_{aIndex}"])

        # bowtie2 parameters
        aParameters = f"--mm -p {FCPUCount}"
        if FMapperCmd is None:
            if FIsLocalMapping == 1:
                aParameters += " --local --sensitive-local"
            else:
                aParameters += " --end-to-end --sensitive"
            if FMatchesCount > 0:
                aParameters += f" -k {FMatchesCount}"
            elif FMatchesCount == 0:
                aParameters += " -a"
        else:
            aParameters += f" {FMapperCmd}"


        aParameters += " --phred33"

        # force sam format
        FMappingFileFormat = C_MAPPING_SAM_OUTPUT_STR
        aParameters += f" {C_SAM2_OUTPUT_PARAMETERS}"

        # complete command with input/output parameters, should iterate only once (because only one index)
        for iBowtieIndex in range(aBowtieIndexCount):

            # ex: FLibraryMappingDir: /projects/parkinson/maping/KCL_01/mapping_vs_hs_9_9_metahit_Proton1-520-Parkinson_R1a_KCL_01
            #     FLibraryName + : Proton1-520-Parkinson_R1a_KCL_01
            #     iBowtieIndex+1 : 1
            FMappingOutputFileNames.append(f"{FLibraryMappingDir}{os.sep}{FLibraryName}_{iBowtieIndex+1}{C_BOWTIE2_OUTPUT_FILE_SUFFIX}")
            aBowtie2Command = f"{aBowtieProgram} {aParameters} -x {aBowtieIndexList[iBowtieIndex]} -U {FFastqList} -S {FMappingOutputFileNames[iBowtieIndex]}"
            # execute command
            subprocess.check_call(["bowtie2", aParameters, "-x", aBowtieIndexList[iBowtieIndex], "-U", FFastqList, "-S", FMappingOutputFileNames[iBowtieIndex]])


    def MapRead(self, aTmpDir, aLibraryMappingDir, aMapperCmd, aMatchesCount, aMismatchesCount, aIsMismatchesPercentage, aIsLocalMapping, aMappingFileFormat, aIsCPUPercentage, aCPUCount):
        # get around KCL feedback (inifile bug?)
        if not isinstance(aMatchesCount, int) :
            FMatchesCount = -1 
        else: 
            FMatchesCount = int(aMatchesCount)
        # FMatchesCount = C_DEFAULT_MAPPING_MATCHES_COUNT if (FMatchesCount == 0)

        # set CPU count
        aHostCPUCount = cpu_count()

        # get around KCL feedback bug(?)
        aCPUCount = -1 if not isinstance(aCPUCount, int) else int(aCPUCount)
        aIsCPUPercentage = 0 if not isinstance(aIsCPUPercentage, int) else int(aIsCPUPercentage)
        if aCPUCount == -1:
            FCPUCount = max(1, aHostCPUCount - 1)
        else:
            FCPUCount = aCPUCount if aIsCPUPercentage == 1 else aCPUCount
        
        # get around KCL feedback bug(?). Default FMismatchesCount is 3
        FMismatchesCount = 3 if not isinstance(aMismatchesCount, int) else int(aMismatchesCount)
        # idem
        FIsMismatchesPercentage = 0 if not isinstance(aIsMismatchesPercentage, int) else int(aIsMismatchesPercentage)
        FIsLocalMapping = 0 if not isinstance(aIsLocalMapping, int) else int(aIsLocalMapping)  # LOCAL
        FParametersShortLine = f"l{FNGSLibraryIndexerReport['MappedReadLength']}-m{FMismatchesCount}"
        FMappedReadCount = 0
        FMappingFileFormat = aMappingFileFormat
        FMapperCmd = aMapperCmd if isinstance(aMapperCmd, str) and aMapperCmd != "" else None
        
        FTmpLibraryMappingDir = FLibraryMappingDir = aLibraryMappingDir
        # FTmpLibraryMappingDir: directory for temporary mapping SAM files.
        # if aTmpDir is set : "~/scratch/tmp_meteor/md5hash" + "/projects/projectname/mapping/sampleN/mapping_vs_blabla"
        # else :              "/projects/projectname/mapping/sampleN/mapping_vs_blabla"
        if aTmpDir is not None:
            FTmpLibraryMappingDir = aTmpDir + os.sep + aLibraryMappingDir 
        
        # create FLibraryMappingDir and FTmpLibraryMappingDir
        os.makedirs(FLibraryMappingDir, exist_ok=True)
        os.makedirs(FTmpLibraryMappingDir, exist_ok=True)
        
        FIsReadyForMapping = 1
        
        # call BowtieMapRead() or Bowtie
        Bowtie2MapRead()
        return True


    def CountReadFastqFile(self, fastq_file):
        """
        Count the number of reads and number of bases
    
        :param fastq_file: Input
    
        Return : parser
        """
        read_count = 0 
        base_Count = 0
        with open(fastq_file, 'r') as fd:
            # fastq may be gziped
            if fastq_file.endswith('.gz'):
                in_fq = gzip.GzipFile(fileobj=fd)
            else:
                in_fq = fd

            # read input fastq line by line
            for read_count,line in enumerate(in_fq):
                # read the sequence
                base_Count += len(next(in_fq).strip())
                # pass the plus
                next(in_fq)
                # pass the quality
                next(in_fq)
        return read_count, base_Count


    def ProcessJob(self, workflow_ini, project_path, input_path, tmp_path, remove_lock, force, mapping_dir,
                   counting_type, counting_only, mapping_only): #options : Workflow, ProjectPath, InputPath, MappingBasename
        # directly from program arguments
        FProjectDir = project_path
        FSampleDir  = input_path
        # TODO get it from census ini file, (FProjectName too)
        FSampleName = os.path.basename(FSampleDir)
        FMappingDone = True
        FCountingDone = True
        FCountingTypeList = [counting_type]
        
        FForce = False
        FNoLock = remove_lock 
        if force:
            FForce = True
            FNoLock = True
        # ex: /projects/parkinson/mapping
        FProjectMappingDir   = FProjectDir + os.sep + mapping_dir

        # load information from workflow.ini ; set FMeteorJobIniFile
        # Check excluded reference
        FMeteorJobIniFile = LoadJobWorkflow(workflow_ini)

        #### in PrepareWorkSpace ?
        # ex: /projects/parkinson/mapping/KCL_01
        FSampleMappingDir = FProjectMappingDir+os.sep+FSampleName
        # make directory (and parent) if needed
        #if not os.path.exists(FSampleMappingDir):
        #    os.makedirs(FSampleMappingDir)  #### maybe obsolete cause done in TaskMainMapping and/or TaskExcludedMapping
        
        #### in PrepareLibraryIniFileNames ?
        print(FSampleDir+os.sep+"*_"+"census_stage_0.ini")
        FLibraryIniFileNames = sorted(glob.glob(FSampleDir+os.sep+"*_"+"census_stage_0.ini"))
        print(FLibraryIniFileNames)
        FLibraryCount        = len(FLibraryIniFileNames)
        if FLibraryCount == 0:
            sys.stderr.write("Error, no census_stage_0.ini file found in {}".format(FSampleDir))
            sys.exit(1)
        print("\nLaunch mapping\n")
        aMainReferenceSection = FMeteorJobIniFile["main_reference"] # alias
        print(aMainReferenceSection)
        aMappingDirPathPrefix = ""
        print( aMainReferenceSection["meteor.mapping.prefix.name"])
        if aMainReferenceSection["meteor.mapping.prefix.name"] is None:
            aMappingDirPathPrefix = "mapping"+"_vs_"+aMainReferenceSection["meteor.reference.name"]+"_l"+ \
                                    "{}-{}".format(aMainReferenceSection["meteor.mapped.readlength"], "m")+aMainReferenceSection["meteor.mismatches"]+"_"
        else:
            aMappingDirPathPrefix = "{}_".format(aMainReferenceSection["meteor.mapping.prefix.name"])

        aExcludedRefCount = FMeteorJobIniFile["worksession"]["meteor.excluded.reference.count"]
        if not isinstance(aExcludedRefCount, int):
            aExcludedRefCount = 0
        else:
            aExcludedRefCount = int(aExcludedRefCount)

        aTmpHash = {}
        # # LOOP ON EACH LIBRARY
        for iLibrary in FLibraryIniFileNames:
            print(FNoLock)
            print(iLibrary)
            # Check if lock file exist
            if os.path.exists(iLibrary+".lock") and not FNoLock:
                sys.stderr.write("Error, lock file not found: {}".format(iLibrary+".lock"))
                sys.exit(1)
            # MAPPING THIS LIBRARY ON MAIN REFERENCE
            aLibraryCensusIniFile = configparser.ConfigParser()
            aLibraryCensusIniFile.read_file(open(iLibrary))
            aSampleInfo = aLibraryCensusIniFile["sample_info"] # reference
            aMappingDirPath = aMappingDirPathPrefix + aSampleInfo["full_sample_name"]

            aTmpHash["directory"] = FProjectMappingDir + os.sep + aSampleInfo["sample_name"]+os.sep+aMappingDirPath
            aTmpHash["Stage1FileName"] = aTmpHash["directory"] + os.sep + os.path.basename(iLibrary).replace("census_stage_0", "census_stage_1")
            FMainMappingCensusIniFileNames = aTmpHash
            if (not os.path.exists(aTmpHash["Stage1FileName"])):
                FMappingDone = False

        FTmpDir = tmp_path
        if not os.path.exists(FTmpDir):
            FTmpDir = tempfile.TemporaryDirectory() 
        if counting_only:
            # mapping already done and no overwriting
            if FMappingDone and not FForce:
                print("\nINFO: Mapping already done for sample: #{@FSampleName}")
                print("INFO: Skipped !")
            else: # mapping not done or we want to overwrite previous results
                # if option -t not provided, generate sample tmpdir in sample dir
                #@FTmpSampleDir = @FTmpDir.nil? ? @FSampleDir + C_PATH_SEP + aMD5 : @FTmpDir #### NEW
                # FTmpSampleDir = @FTmpDir.nil? ? @FSampleMappingDir + C_PATH_SEP + aMD5 : @FTmpDir #### NEW
                
                # if File.exists?(@FTmpSampleDir)
                #     STDERR.puts "WARNING, temporary dir #{@FTmpSampleDir} already exists !"
                # end
                # FileUtils.mkdir_p(@FTmpSampleDir)
                with tempfile.TemporaryDirectory() as tmpdirname:
                    print("le temp directory {tmpdirname} est crée")
                LaunchMapping(options)
                # remove sample tmp data if created in sample dir
                # FileUtils.rm_rf(@FTmpSampleDir) if ( @FTmpDir.nil? and File.exists?(@FTmpSampleDir) ) #### NEW
        if mapping_only:
            LaunchCounting()
        print("Done !\nJob finished without errors ...")
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


#-------------------------- FUNCTIONS DEFINITION ------------------------------#

    # def CountReadAndReIndexFastqFile(fastq_file, output_file):
    #     aInputLineCount = 0
    #     aReadCount = 0
    #     aBaseCount = 0
    #     aCpt = 0
    #     adn = None
    #     qual = None
    #     prev = None
    #     seqId = None
    #     with open(aInputFile, 'r') as fd:
    #         # fastq may be gzipped, xz or bz2
    #         if aInputFile.endswith('.gz'):
    #             in_fq = gzip.GzipFile(fileobj=fd)
    #         else:
    #             in_fq = fd

    #         # open output fastq in write mode
    #         with open(output_file, 'w') as out_fq:
    #             # read input fastq line by line
    #             for line in in_fq:
    #                 line = line.strip()
    #                 # line begins with @seqid
    #                 if re.match(r'^@(.+)$', line):
    #                     # check if the previous read is valid (quality exists and same size as adn)
    #                     aQualSize = len(qual) if qual else 0
    #                     if aQualSize > 0:
    #                         if adn and aQualSize == len(adn):
    #                             aReadCount += 1
    #                             aBaseCount += aQualSize
    #                             # then write this indexed read
    #                             out_fq.write(f"@{aReadCount}\n{adn}\n+\n{qual}\n")
    #                             adn = None
    #                     seqId = re.match(r'^@(.+)$', line).group(1)
    #                     # NB: quality line might start with @
    #                 aCpt += 1
    #                 if line == "+" or line == f"+{seqId}":
    #                     qual = ''
    #                     aCpt = 3
    #                     # previous line was adn
    #                     adn = prev
    #                 if aCpt == 4:
    #                     qual = line
    #                 prev = line

    #             # evaluate last read
    #             aQualSize = len(qual) if qual else 0
    #             if aQualSize > 0:
    #                 if adn and aQualSize == len(adn):
    #                     aReadCount += 1
    #                     aBaseCount += aQualSize
    #                     out_fq.write(f"@{aReadCount}\n{adn}\n+\n{qual}\n")

    #     return aReadCount, aBaseCount


def TaskMainMapping(iLibrary, FMeteorJobIniFile):
    aWorkSessionSection = FMeteorJobIniFile["worksession"] # reference
    aReferenceSection = FMeteorJobIniFile["main_reference"] # reference
    aReferenceIniFileName = aWorkSessionSection["meteor.reference.dir"] + os.sep + \
        aReferenceSection["meteor.reference.name"] + os.sep + aReferenceSection["meteor.reference.name"] + "_reference.ini"

    print("\n######### Task main mapping {}".format(iLibrary))
    
    # What is Stage1FileName ?
    aMappedCensusIniFileName = iLibrary["Stage1FileName"]
    aLibraryMappingDir = iLibrary["directory"]
    
    aLockedMappedCensusIniFileName = aMappedCensusIniFileName+ ".lock"
    
    # remove lock file if asked
    if (FNoLock == True and os.path.exists(aLockedMappedCensusIniFileName)):
        print("\nINFO: removing lock file: {}".format(aLockedMappedCensusIniFileName))
        os.remove(aLockedMappedCensusIniFileName)
    # remove census 1 file if asked
    if (FForce == True and os.path.exists(aMappedCensusIniFileName)):
        print("\nINFO: removing existing census 1 file: {}".format(aMappedCensusIniFileName))
        os.remove(aMappedCensusIniFileName)
    
    if ( os.path.exists(aMappedCensusIniFileName) or os.path.exists(aLockedMappedCensusIniFileName) ):
        print("\nINFO: skipping library.\nINFO: {} already exists or is locked".format(aMappedCensusIniFileName))
        return True

    # create lock file in result directory
    os.makedirs(aLibraryMappingDir, exist_ok=True)
    open(aLockedMappedCensusIniFileName, 'a').close()
    
    aMappingProg = aWorkSessionSection["meteor.mapping.program"]
    if (aMappingProg != 'bowtie2'): 
        sys.stderr.write("Error, unknown or unsupported mapping program : {}".format(aWorkSessionSection["meteor.mapping.program"]))
        sys.exit(1) # C_BAD_MAPPER_EXIT

    aMeteorMapper = MeteorMapper(aMappingProg, FLibraryIniFileNames[iLibrary], aReferenceIniFileName, FLibraryIndexerReport, aMappedCensusIniFileName)

    aOkMappingProcess = aMeteorMapper.MapRead(
                    FTmpDir,
                    aLibraryMappingDir,
                    aReferenceSection["meteor.mapper.cmd"],
                    aReferenceSection["meteor.matches"],
                    aReferenceSection["meteor.mismatches"],
                    aReferenceSection["meteor.is.perc.mismatches"],
                    aReferenceSection["meteor.is.local.mapping"],
                    aWorkSessionSection["meteor.mapping.file.format"],
                    aWorkSessionSection["meteor.is.cpu.percentage"],
                    aWorkSessionSection["meteor.cpu.count"])
                    
    if (aOkMappingProcess == True):
        # remove lock file
        if (File.exists(aLockedMappedCensusIniFileName)):
            FileUtils.rm(aLockedMappedCensusIniFileName)
        return True
    else:
        # remove lock file
        if (File.exists(aLockedMappedCensusIniFileName)):
            FileUtils.rm(aLockedMappedCensusIniFileName)
        return False

def TaskExcludedMapping(iLibrary, iExcluded):
    aWorkSessionSection = FMeteorJobIniFile["worksession"]  # reference
    aReferenceSection = FMeteorJobIniFile["excluded_reference_" + f"{iExcluded+1}"]  # reference

    aReferenceIniFileName = (
        aWorkSessionSection["meteor.reference.dir"]
        + os.sep
        + aReferenceSection["meteor.reference.name"]
        + os.sep
        + aReferenceSection["meteor.reference.name"]
        + "_reference.ini"
    )

    print(f"\n######### TaskExcludedMapping({iLibrary}, {iExcluded+1}); reference name: " + aReferenceSection["meteor.reference.name"])

    aMappedCensusIniFileName = FExcludedMappingCensusIniFileNames[iLibrary][iExcluded]["Stage1FileName"]
    aLibraryMappingDir = FExcludedMappingCensusIniFileNames[iLibrary][iExcluded]["directory"]

    aLockedMappedCensusIniFileName = aMappedCensusIniFileName + ".lock"

    # remove lock file if asked
    if FNoLock and os.path.exists(aLockedMappedCensusIniFileName):
        print(f"\nINFO: removing lock file: {aLockedMappedCensusIniFileName}")
        os.remove(aLockedMappedCensusIniFileName)
    # remove census 1 file if asked
    if FForce and os.path.exists(aMappedCensusIniFileName):
        print(f"\nINFO: removing existing census 1 file: {aMappedCensusIniFileName}")
        os.remove(aMappedCensusIniFileName)

    if os.path.exists(aMappedCensusIniFileName) or os.path.exists(aLockedMappedCensusIniFileName):
        print(f"\nINFO: skipping library.\nINFO: {aMappedCensusIniFileName} already exists or is locked")
        return True

    # create lock file in result directory
    os.makedirs(aLibraryMappingDir, exist_ok=True)
    open(aLockedMappedCensusIniFileName, "a").close()

    aMappingProg = aWorkSessionSection["meteor.mapping.program"]
    if aMappingProg != 'bowtie2':
        raise ValueError(f"Error, unknown or unsupported mapping program:"+ aWorkSessionSection["meteor.mapping.program"])

    aMeteorMapper = MeteorMapper(
        aMappingProg,
        FLibraryIniFileNames[iLibrary],
        aReferenceIniFileName,
        FLibraryIndexerReport,
        aMappedCensusIniFileName
    )

    aOkMappingProcess = aMeteorMapper.MapRead(
        FTmpDir,
        aLibraryMappingDir,
        aReferenceSection["meteor.mapper.cmd"],
        aReferenceSection["meteor.matches"],
        aReferenceSection["meteor.mismatches"],
        aReferenceSection["meteor.is.perc.mismatches"],
        False,
        aWorkSessionSection["meteor.mapping.file.format"],
        aWorkSessionSection["meteor.is.cpu.percentage"],
        aWorkSessionSection["meteor.cpu.count"]
    )
    # remove lock file
    if os.path.exists(aLockedMappedCensusIniFileName):
        os.remove(aLockedMappedCensusIniFileName)
    return True


def LoadJobWorkflow(aWorkflowFile):
    print("youpi")
    if not os.path.exists(aWorkflowFile):
        print(f"Error, file {aWorkflowFile} not found", file=sys.stderr)
        sys.exit(1)
    FMeteorJobIniFile = configparser.ConfigParser()
    FMeteorJobIniFile.read_file(open(aWorkflowFile))
    # check if excluded reference count is correct
    # aExcludedRefCount = FMeteorJobIniFile["worksession"]["meteor.excluded.reference.count"]
    # if not isinstance(aExcludedRefCount, int):
    #     aExcludedRefCount = False
    # else:
    #     aExcludedRefCount = int(aExcludedRefCount)
    
    # print(FMeteorJobIniFile)
    aExcludedRefCounted = len(
         [ s for s in FMeteorJobIniFile.sections() if s.startswith("excluded_reference_")]
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
    return FMeteorJobIniFile


def CountAndReIndexReads(aLibraryCensusIniFile):
    aSampleInfoSection = aLibraryCensusIniFile["sample_info"]
    aSampleFileSection = aLibraryCensusIniFile["sample_info"]

    FLibraryIndexerReport.clear()

    aInputFile = FSampleDir + C_PATH_SEP + aSampleFileSection["fastq_file"]
    aOutputFile = FTmpSampleDir + C_PATH_SEP + aSampleFileSection["fastq_file"] + '.idx'

    aReadCount, aBaseCount = CountReadAndReIndexFastqFile(aInputFile, aOutputFile)

    FLibraryIndexerReport = {
        "IsIndexed": 1,
        "IndexedfastqFilePath": aOutputFile,
        "IndexedFileNameExt": '.idx',
        "IndexedReadCount": aReadCount,
        "IndexedBaseCount": aBaseCount,
        "IndexedReadLength": -1,
        "MappedReadLength": -1,
        "MappedReadLengthType": "overall",
    }


def LaunchMapping():
    print("\nLaunch mapping\n")
    aOKToContinue = True

    #~ ( build catalog reference db (bowtie2-build-l) )

    #### why not launching one bowtie process for all the libraries (see bowtie -U option in example below) ????
    #bt2_index="/projects/biodatabank/filtering_database/Homo_sapiens_2014_02_04/bowtie2/homo_sapiens.fna"
    #command_bt2="bowtie2 -q --no-head --no-sq --no-unal --omit-sec-seq --end-to-end --sensitive -k 1 -S $sample/$sample.sam -U $fastq_list -p 40 -x $bt2_index"

    aExcludedRefCount = FMeteorJobIniFile["worksession"]["meteor.excluded.reference.count"]
    if not integerable(aExcludedRefCount):
        aExcludedRefCount = False
    else:
        aExcludedRefCount = int(aExcludedRefCount)
    # LOOP ON EACH LIBRARY
    for iLibrary in range(FLibraryCount):

        if (os.path.isdir(FProjectDir) and os.path.isfile(FLibraryIniFileNames[iLibrary])):

            aLibraryCensusIniFile = IniFile.load(FLibraryIniFileNames[iLibrary])
            aSampleInfoSection = aLibraryCensusIniFile["sample_info"] # reference
            aSampleFileSection = aLibraryCensusIniFile["sample_file"] # reference
            
            FSampleName = aSampleInfoSection["sample_name"]
            FProjectName = aSampleInfoSection["project_name"]
            aSampleLibraryName = aSampleInfoSection["full_sample_name"]
            
            print("\n\n######### Meteor Mapping task description")
            print("Sample name = " + FSampleName)
            print("Library name = " + aSampleLibraryName)
            print("Project name = " + FProjectName)
            print("Sequencing device = " + aSampleInfoSection["sequencing_device"])
            print("Workflow = " + os.path.basename(FMeteorJobIniFile.filename))
            
            if not os.path.exists(FLibraryIniFileNames[iLibrary]+".lock"):
                
                ### reindexing this library reads and fill FLibraryIndexerReport
                CountAndReIndexReads(aLibraryCensusIniFile)
                
                # MAPPING THIS LIBRARY ON MAIN REFERENCE
                aOKToContinue = TaskMainMapping(iLibrary)
                if not aOKToContinue:
                    sys.stderr.write("Error, TaskMainMapping failed: " + FLibraryIniFileNames[iLibrary])
                    sys.exit(1) #### TODO
                
                # MAPPING THIS LIBRARY ON EACH EXCLUDED REFER
                for iExcluded in range(aExcludedRefCount):
                    aOKToContinue = TaskExcludedMapping(iLibrary, iExcluded)
                    if not aOKToContinue:
                        sys.stderr.write("Error, TaskExcludedMapping failed: " + FLibraryIniFileNames[iLibrary])
                        sys.exit(1)
                
                if aSampleInfoSection[C_SEQUENCING_DEVICE_STR] == C_SOLID_DEVICE:
                    # MAPPING THIS LIBRARY ON MAIN REFERENCE WITH SOLID DATA PROCESSING
                    aOKToContinue = TaskMainMappingWithSolidProcessing(iLibrary)
                    if not aOKToContinue:
                        sys.stderr.write("Error, TaskMainMappingWithSolidProcessing failed: " + FLibraryIniFileNames[iLibrary])
                        sys.exit(1)
                        
                    # MAPPING THIS LIBRARY ON EACH EXCLUDED REFERENCE WITH SOLID DATA PROCESSING
                    for iExcluded in range(aExcludedRefCount):
                        aOKToContinue = TaskExcludedMappingWithSolidProcessing(iLibrary, iExcluded)
                        if not aOKToContinue:
                            sys.stderr.write("Error, TaskExcludedMappingWithSolidProcessing failed: " + FLibraryIniFileNames[iLibrary])
                            sys.exit(1)
                else:
                    # MAPPING THIS LIBRARY ON MAIN REFERENCE WITH FASTQ DATA PROCESSING
                    aOKToContinue = TaskMainMappingWithFastqProcessing(iLibrary)
                    if not aOKToContinue:
                        sys.stderr.write("Error, TaskMainMappingWithFastqProcessing failed: " + FLibraryIniFileNames[iLibrary])
                        sys.exit(1)
                        
                    # MAPPING THIS LIBRARY ON EACH EXCLUDED REFERENCE WITH FASTQ DATA PROCESSING
                    for iExcluded in range(aExcludedRefCount):
                        aOKToContinue = TaskExcludedMappingWithFastqProcessing(iLibrary, iExcluded)
                        if not aOKToContinue:
                            sys.stderr.write("Error, TaskExcludedMappingWithFastqProcessing failed: " + FLibraryIniFileNames[iLibrary])
                            sys.exit(1)
                            
                # Generate the final result file for this library
                GenerateResultFile(iLibrary)

            else:
                sys.stderr.write("Error: Cannot process library. Lock file found: " + FLibraryIniFileNames[iLibrary] + ".lock\n")
        else:
            sys.stderr.write("Error: Cannot find input files for library " + FLibraryIniFileNames[iLibrary] + "\n")
    
    # Create final summary file for the entire project
    CreateSummaryFile()
    
    return aOKToContinue


def LaunchCounting():
    print("\nLaunch counting")
    # does not need census_stage_0.ini
    #-w /path/to/workflow_tutorial.ini -i /path/to/sample/H1 -p /path/to/project_name -m mapping

    aCountingProgram = C_METEOR_COUNTER
    
    aparameters = " -w " + FMeteorJobIniFile.filename + " -i " + FSampleDir + " -p " + FProjectDir + " -o " + os.path.basename(FProjectMappingDir)
    if FCountingTypeList and not FCountingTypeList.empty():
        aparameters += " -c " + ",".join(FCountingTypeList.to_a)
    if FForce: # force overwriting former profiling results done with same parameters
        aparameters += " -f"
    if FTmpDir: # path to directory where mapping results (SAM files) are stored. #### NEW
        aparameters += " -t " + FTmpDir
    
    aCmd = aCountingProgram + aparameters

    ### Execute Command
    print("\nExecuting counting command:")
    print(aCmd)
    aExitStatus = ExecuteCmd(aCmd)
    
    if aExitStatus != 0:
        sys.stderr.write("Error in MeteorSession::LaunchCounting, counting command exited with non zero status")
        sys.exit(aExitStatus)
        
    return aExitStatus


def isfile(path): # pragma: no cover
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


def isdir(path): # pragma: no cover
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


# def get_log(path_log):
#     """
#     """
#     logger = logging.getLogger()
#     logger.setLevel(logging.DEBUG)
#     formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s')
#     # Create log file
#     file_handler = RotatingFileHandler(path_log, 'a', 1000000, 1)
#     file_handler.setFormatter(formatter)
#     logger.addHandler(file_handler)
#     # Stream in the the console
#     ## TO REMOVE IF daemon
#     stream_handler = logging.StreamHandler()
#     stream_handler.setLevel(logging.DEBUG)
#     logger.addHandler(stream_handler)
#     return logger


def get_arguments(): # pragma: no cover
    """
    Meteor help and arguments
    
    No arguments
    
    Return : parser
    """
    parser = argparse.ArgumentParser(description=color.BOLD + __doc__ + color.END)
    subparsers = parser.add_subparsers(title = 'positional arguments', help="Select activity")
    # # Mappping commands
    download_parser = subparsers.add_parser('download',
        help='Download catalog')
    reference_parser = subparsers.add_parser('build',
        help='Index reference')
    reference_parser.add_argument("-i", dest='input_fasta_file',
        type = str, required = True, help = "Input fasta filename.")
    reference_parser.add_argument("-p", dest='refRootDir',
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
    mapping_parser = subparsers.add_parser('counter',
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
    mapping_parser.add_argument("-c", dest="counting_type", type=str, nargs='+', 
        choices = ["total_reads", "shared_reads", "smart_shared_reads", "unique_reads"],
    help="""counting type string.""")
    mapping_parser.add_argument("-s", dest="strand", type=str, choices=["direct", "undirect"],
        help="(where direct/undirect means direct/undirect strand)")
    mapping_parser.add_argument("-a", dest="coverage", type=str, choices=["coverage", "avg_coverage"],
        help="coverage or average coverage")
    mapping_parser.add_argument("-t", dest="tmp_path", type=str, required=True,
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
        help="execute mapping only")
    profiler_parser = subparsers.add_parser('profiler',
        help='Map reads against a gene catalog')
    return parser.parse_args(args=None if sys.argv[2:] else ['--help'])

# class Person:
#   def __init__(mysillyobject, name, age):
#     mysillyobject.name = name
#     mysillyobject.age = age

#   def myfunc(abc):
#     print("Hello my name is " + abc.name)


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Let's logging
    now = datetime.datetime.now()
    path_log = "meteor_" + now.strftime("%Y%m%d_%H%M") + ".log"
    # Get arguments
    args = get_arguments()
    # logger = get_log(path_log)
    # Import FASTQ UGLY
    if hasattr(args, "fastq_dir"):
        if args.isdispatched:
            fastq_file_list = glob.glob(args.fastq_dir + "*" + os.sep + "*.f*q*")
        else:
            print(args.fastq_dir + "*.{fq,fastq}*")
            fastq_file_list = glob.glob('{}{}'.format(args.fastq_dir, "*.f*q*"))
        ext_R1 = ("1.fq", "1.fastq", "R1.fastq", "R1.fastq.gz", "R1.fq.gz")
        ext_R2 = ("2.fq", "2.fastq", "R2.fastq", "R2.fastq.gz", "R2.fq.gz")
        for fastq_file in fastq_file_list:
            # print("Import ", fastq_file)
            # logger.info("Import ", fastq_file)
            full_sample_name = os.path.basename(fastq_file)
            if args.iscompressed:
                full_sample_name = ".".join(full_sample_name.split(".")[:-1])
            print(full_sample_name)
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
            config = configparser.ConfigParser()
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
    # Import reference UGLY
    elif hasattr(args, "input_fasta_file"):
        # Create reference genome directory if it does not already exist
        ref_dir = os.path.join(args.refRootDir, args.refName)
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
        config = configparser.ConfigParser()
        config["reference_info"] = {
            "reference_name": args.refName,
            "entry_type": "fragment", # Why ?
            "reference_date": datetime.datetime.now().strftime("%Y%m%d"),
            "database_type": "text",
            "HAS_LITE_INFO": 1
        }

        config["reference_file"] = {
            #IS_LARGE_REFERENCE_STR: 1,
            "database_dir": "database_dir", #WTF ?
            "fasta_dir": "fasta_dir", #WTF ?
            "fasta_file_count": 1,
            # is it possible to have several fasta
            "fasta_file_count": args.refName + '.fasta'
        }

        index_prefix = os.path.join(database_dir, args.refName)
        subprocess.check_call(["bowtie2-build", '-f', '-t', str(args.threads), 
            output_fasta_file, output_fasta_file])
        config["bowtie2_index"] = {
            # "is_large_reference": "1", # WTF 
            "is_DNA_space_indexed": 1,
            "dna_space_bowtie_index_prefix_name_1": args.refName
        }

        # Write configuration file
        with open(os.path.join(ref_dir, args.refName + '_reference.ini'), 'wt') as config_file:
            config.write(config_file)
    elif hasattr(args, "input_path"):
        print("What am I supposed to do now ?")
        m= MeteorSession()
        counting_type = args.counting_type
        if args.strand:
            counting_type = args.strand + "_" + counting_type
        if args.coverage:
            counting_type = counting_type + "_" + coverage
        m.ProcessJob(args.workflow_ini,
        args.project_path,
        args.input_path,
        args.tmp_path,
        args.remove_lock,
        args.force,
        args.mapping_dir,
        counting_type,
        args.counting_only,
        args.mapping_only)
        # p1 = Person("John", 36)
        # p1.myfunc()


if __name__ == '__main__':
    main()