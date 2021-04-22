#!/usr/bin/env ruby

#
## Copyright 2017-2020, Franck Gauthier <franck.gauthier@inrae.fr>, Nicolas Pons <nicolas.pons@inrae.fr>
##
## This file is part of Meteor v3.2.
##
## Meteor v3.2 is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## Meteor v3.2 is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Meteor v3.2. If not, see <http://www.gnu.org/licenses/>.


require 'digest/md5' # for temporary dir creation
require "fileutils"
require 'optparse'
require 'rubygems'
require 'inifile'
require 'zlib'
require 'date'
require 'set'
require 'pp' # for debuging

################################################################################
#                                 Constants                                    #
################################################################################

C_IS_UNIX = (RUBY_PLATFORM =~ /mswin|mingw/).nil?
C_EXEC_EXT = C_IS_UNIX ? "" : ".exe"
C_PATH_SEP = File::SEPARATOR
C_INDEXED_EXT = '.idx'
C_LOCK_EXT = ".lock"
C_REFERENCE_INIFILE_EXT = "_reference.ini"

C_CENSUS_SEQUENCED_STR = "census_stage_0"
C_CENSUS_MAPPED_STR = "census_stage_1"
C_DEFAULT_SAMPLE_MAPPING_DIRECTORY = "mapping"

C_BOWTIE_MAPPER = 'bowtie'+C_EXEC_EXT
C_BOWTIE_LARGE_MAPPER = 'bowtie-align-l'+C_EXEC_EXT
C_BOWTIE_SMALL_MAPPER = 'bowtie-align-s'+C_EXEC_EXT
C_BOWTIE_OUTPUT_EXT = '.sam' #'.bowtie'
#~ C_LITE_OUTPUT_PARAMETERS = '--suppress 6,7'
C_SAM_OUTPUT_PARAMETERS = '--sam --sam-nohead --sam-nosq --no-unal'

C_BOWTIE2_MAPPER = 'bowtie2'+C_EXEC_EXT
C_BOWTIE2_LARGE_MAPPER  = 'bowtie2-align-l'+C_EXEC_EXT
C_BOWTIE2_SMALL_MAPPER = 'bowtie2-align-s'+C_EXEC_EXT
C_BOWTIE2_OUTPUT_EXT = '.sam'
C_METEOR_COUNTER = 'meteor-counter'+C_EXEC_EXT

### device kind
C_SOLID_DEVICE_STR = 'SOLiD'
C_SOLEXA_DEVICE_STR = 'Solexa'
C_ILLUMINA_DEVICE_STR = 'Illumina'
C_PROTON_DEVICE_STR = 'Proton'
### sequence file format : CS for csfasta + qual, FQ for fastq 
Sff = [CS=0, FQ=1] # ruby has no enum type ... use Sff[CS], Sff[FQ]

### quality encoding
C_QE_SOLID = 'solid'
C_QE_SOLEXA = 'solexa'
C_QE_ILLUMINA_13 = 'illumina 1.3+'
C_QE_ILLUMINA_15 = 'illumina 1.5+'
C_QE_ILLUMINA_18 = 'illumina 1.8+'
C_QE_SANGER = 'sanger'
C_QE_PROTON = 'proton'
C_QE_UNKNOWN = 'unknown'

C_DEFAULT_MAPPING_READ_LENGTH = 35
C_MAPPED_READ_LENGTH_CHAR = 'l'
C_MISMATCHES_COUNT_CHAR = 'm'

### census ini file
C_CENSUS_STATUS_MAPPED = 1
# sample info section
C_SAMPLE_INFO_SECTION = "sample_info"
C_SAMPLE_NAME_STR     = "sample_name"
C_SAMPLE_FULL_NAME_STR = "full_sample_name"
C_PROJECT_NAME_STR = "project_name"
C_SEQUENCING_DEVICE_STR = "sequencing_device"
C_QUALITY_ENCODING_STR = "quality_encoding"
C_READ_LENGTH_STR = "read_length"
C_METEOR_PROGRAM_NAME ="Meteor"
C_METEOR_VERSION = "3.2"
C_METEOR_DBTYPE_UNKNOWN_STR = "unknown"
C_INDEXED_READ_LENGTH_STR = "indexed_read_length"
C_SAMPLE_TOTAL_READ_COUNT_STR = "sequenced_read_count"
C_SAMPLE_INDEXED_READ_COUNT_STR = "indexed_sequenced_read_count"
C_SAMPLE_INDEXED_BASE_COUNT_STR = "indexed_sequenced_base_count"
C_DATABASE_TYPE_STR = "database_type"
C_SAMPLE_IS_DATA_PREFIXED_STR = "is_data_prefixed"
C_SAMPLE_INDEXING_SOFTWARE_STR = "sample_indexing_software"
C_SAMPLE_INDEXING_SOFTWARE_VERSION_STR = "sample_indexing_software_version"
C_CENSUS_STATUS_STR = "census_status"

# sample file section
C_SAMPLE_FILE_SECTION = "sample_file"
C_FASTQ_FILENAME_STR = "fastq_file"
C_CSFASTA_FILENAME_STR = "csfasta_file"
C_QUAL_FILENAME_STR = "qual_file"
# mapping file section
C_MAPPING_FILE_SECTION = "mapping_file"
C_MAPPING_BOWTIE_FILENAME_STR = 'bowtie_file'
C_MAPPING_FILE_FORMAT_STR = "mapping_file_format"
C_MAPPING_SAM_OUTPUT_STR = "sam"
C_SAM2_OUTPUT_PARAMETERS = '--no-head --no-sq --no-unal --omit-sec-seq'
C_MAPPING_FILE_COUNT_STR = 'mapping_file_count'
# mapping section
C_MAPPING_SECTION = 'mapping'
C_MAPPING_TOOL_VERSION_STR = "mapping_tool_version"
C_MAPPING_TOOL_STR = "mapping_tool"
C_METEOR_MAPPING_PROGRAM_NAME = "Meteor"
C_PROCESSED_READ_COUNT_STR = "processed_read_count"
C_MAPPING_DATE_STR = "mapping_date"
C_MAPPING_CMDLINE_STR = "mapping_cmdline"
C_MAPPING_PARAMETERS_STR = "parameters"
C_MAPPED_READ_LENGTH_STR = "mapped_read_length"
C_MAPPED_READ_LENGTH_TYPE_STR = "mapped_read_length_type"
C_MAPPING_MISMATCHES_STR = "mismatches"
C_MAPPING_IS_MISMATCHES_PERCENTAGE_STR = "is_mismatches_percentage"
C_MAPPING_IS_LOCAL_MAPPING_STR = "is_local_mapping"                   # LOCAL
#C_ALIGNMENT_LENGTH_CUTOFF_STR = "alignment_length_cutoff"             # LOCAL
#C_KEEP_INTERNAL_LOCAL_ALIGNMENT_STR = "keep_internal_local_alignment" # LOCAL
C_MAPPING_MATCHES_STR = "matches"
C_MAPPING_SOFTWARE_STR = "mapping_software"
C_MAPPING_SOFTWARE_VERSION_STR = "mapping_software_version"

### workflow ini file
# worksession section
C_JOB_WORKSESSION_SECTION = "worksession"
C_JOB_EXCLUDED_REFERENCE_COUNT = "meteor.excluded.reference.count"
C_JOB_REFERENCE_DIR = "meteor.reference.dir"
C_JOB_IS_CPU_PERCENTAGE = "meteor.is.cpu.percentage";
C_JOB_CPU_COUNT = "meteor.cpu.count";
C_JOB_WORK_DIR = "meteor.workdir";
C_JOB_REPOSITORY_ROOT_PATH = "meteor.repository.root.path";
C_JOB_METEOR_DB_TYPE = "meteor.db.type";
C_JOB_METEOR_ATTEMPS_COUNT = "meteor.nb.attempts";
C_JOB_INDEXED_READ_LENGTH = "meteor.indexed.readlength";
C_JOB_MAPPING_PROGRAM = "meteor.mapping.program";
C_JOB_MAPPING_DIR = "meteor.mapping.dir";
C_JOB_PROFILE_DIR = "meteor.profile.dir";
C_JOB_MAPPING_FILE_FORMAT = "meteor.mapping.file.format";
C_JOB_IS_DATA_CACHING = "meteor.datacaching";

# ref sections
C_JOB_MAIN_REFERENCE_SECTION = "main_reference"
C_JOB_EXCLUDED_REFERENCE_PREFIX_SECTION = "excluded_reference_"
C_JOB_MAPPING_PREFIX_NAME    = "meteor.mapping.prefix.name"
C_JOB_MAPPING_MISMATCHES = "meteor.mismatches"
C_JOB_MAPPING_MATCHES = "meteor.matches";
C_JOB_IS_MAPPING_MISMATCHES_PERCENTAGE = "meteor.is.perc.mismatches";
C_JOB_REFERENCE_NAME = "meteor.reference.name"
C_JOB_IS_BESTALIGNMENT = "meteor.bestalignment";
C_JOB_MAPPER_CMD = "meteor.mapper.cmd";
C_JOB_IS_LOCAL_MAPPING = "meteor.is.local.mapping";                           # LOCAL
#C_JOB_ALIGNMENT_LENGTH_CUTOFF = "meteor.alignment.length.cutoff";             # LOCAL
#C_JOB_KEEP_INTERNAL_LOCAL_ALIGNMENT = "meteor.keep.internal.local.alignment"; # LOCAL
C_JOB_MAPPED_READ_LENGTH = "meteor.mapped.readlength";
C_JOB_MAPPED_READ_LENGTH_TYPE = "meteor.mapped.readlength.type";

C_JOB_IS_ON_BOTH_STRANDS = "meteor.is.on.both.strands";
C_DEFAULT_MAPPING_MATCHES_COUNT = 10000;

C_REFERENCE_INFO_SECTION = 'reference_info'
C_REFERENCE_NAME_STR = 'reference_name'

################################################################################
#                                 methods                                      #
################################################################################

#------------------------------------------------------------------------------
# get around KCL feedback (inifile 2.0.2 issues)
def integerable(aValue)
	return true if ( aValue.is_a?(Fixnum) or aValue =~ /^ *-?\d+ *$/ )
	return false
end

#------------------------------------------------------------------------------

def cpu_count()
	require 'rbconfig'
	
	case RbConfig::CONFIG['host_os']
	when /darwin9/
		`hwprefs cpu_count`.to_i
	when /darwin/
		((`which hwprefs` != '') ? `hwprefs thread_count` : `sysctl -n hw.ncpu`).to_i
	when /linux/
		`cat /proc/cpuinfo | grep processor | wc -l`.to_i
	when /freebsd/
		`sysctl -n hw.ncpu`.to_i
	when /mswin|mingw/
		require 'win32ole'
		wmi = WIN32OLE.connect("winmgmts://")
		cpu = wmi.ExecQuery("select NumberOfCores from Win32_Processor") # TODO count hyper-threaded in this
		cpu.to_enum.first.NumberOfCores
	end
end

#------------------------------------------------------------------------------

def CheckCountingType(aCountingTypeList)

	aAcceptedTypes = Set.new [
		"total_reads", "shared_reads", "smart_shared_reads", "unique_reads",
		"direct_total_reads", "direct_shared_reads", "direct_smart_shared_reads",
		"direct_unique_reads",
		"total_reads_coverage", "unique_reads_coverage",
		"direct_total_reads_coverage", "direct_unique_reads_coverage",
		"undirect_total_reads_coverage", "undirect_unique_reads_coverage",
		"total_reads_mean_coverage", "shared_reads_mean_coverage",
		"smart_shared_reads_mean_coverage", "unique_reads_mean_coverage",
		"direct_total_reads_mean_coverage", "direct_shared_reads_mean_coverage",
		"direct_smart_shared_reads_mean_coverage",
		"direct_unique_reads_mean_coverage",
		"undirect_total_reads_mean_coverage",
		"undirect_shared_reads_mean_coverage",
		"undirect_smart_shared_reads_mean_coverage",
		"undirect_unique_reads_mean_coverage",
		"undirect_total_reads", "undirect_shared_reads",
		"undirect_smart_shared_reads", "undirect_unique_reads"]
	
	# we want all counting types
	if aCountingTypeList.nil? or aCountingTypeList.include?("all") or aCountingTypeList == aAcceptedTypes
		return Set.new ["all"]
	end
	# check the chosen types validity
	aInvalidTypes = aCountingTypeList.difference(aAcceptedTypes)
	if (not aInvalidTypes.empty?)
		raise OptionParser::InvalidArgument, "#{aInvalidTypes.to_a.join(",")} not valid for option [-C]"
	end
	# Ok, return the validated and well-ordered list (i.e. in the same order as in aAcceptedTypes)
	return aCountingTypeList.intersection aAcceptedTypes
end

#------------------------------------------------------------------------------

def CheckArgs(options)
	raise OptionParser::MissingArgument, "[-w] path to meteor configuration file." if options["Workflow"].nil?
	raise OptionParser::MissingArgument, "[-i] path to sample directory." if options["InputPath"].nil?
	raise OptionParser::MissingArgument, "[-p] path to project directory." if options["ProjectPath"].nil?
	raise OptionParser::MissingArgument, "[-o] basename of the mapping directory." if options["MappingBasename"].nil?
	options["CountingTypeList"] = CheckCountingType(options["CountingTypeList"])
	#raise OptionParser::MissingArgument, "[-t] path to the directory where temporary files (e.g. sam) are stored." if options["TmpPath"].nil?
	raise OptionParser::AmbiguousOption, "[-m] (mapping only) and [-c] (counting only) are not compatible" if (options["OnlyMapping"] and options["OnlyCounting"])

	return true
end
	
#------------------------------------------------------------------------------

def ExecuteCmd(aCmd) #### TODO

	aSuccess = system(aCmd)
	return $?.exitstatus

end


################################################################################
#                           Class definitions                                  #
################################################################################

class MeteorMapper
	
	attr_accessor :FLibraryCensusIniFile
	attr_accessor :FMappedCensusIniFileName
	attr_accessor :FReferenceIniFile
	attr_accessor :FLibraryMappingDir # /path/to/mapping/Sample/mapping_vs_catalog_lib...
	attr_accessor :FTmpLibraryMappingDir # /path/to/mapping/Sample/mapping_vs_catalog_lib...
	
	# reference to MeteorSession member (hash) FNGSLibraryIndexerReport:
	# IsIndexed, IndexedFileNameExt, IndexedReadCount, IndexedBaseCount, IndexedReadLength (=census_0 read_length), MappedReadLength (-1 or 35), MappedReadLengthType (overall or fixed)
	attr_accessor :FNGSLibraryIndexerReport
	attr_accessor :FMismatchesCount
	attr_accessor :FMatchesCount
	attr_accessor :FIsMismatchesPercentage
	attr_accessor :FReferenceName
	attr_accessor :FSequenceFileFormat
	attr_accessor :FIsColorSpaceMapping
	
	attr_accessor :FMappingProgram
	attr_accessor :FMappingFileFormat
	attr_accessor :FParametersShortLine
	attr_accessor :FMapperCmd
	attr_accessor :FMappingOutputFileNames
	attr_accessor :FTmpMappingOutputFileNames
	attr_accessor :FIsReadyForMapping
	attr_accessor :FCPUCount
	#FSequenceFileFormat
	attr_accessor :FQualityEncoding
	attr_accessor :FIsLocalMapping # LOCAL
	
	attr_accessor :FBowtieProgram
	
	# default constructor
	def initialize()
		@FLibraryCensusIniFile = nil
		@FMappedCensusIniFileName = nil
		@FReferenceIniFile = nil
		@FLibraryMappingDir = nil
		@FTmpLibraryMappingDir = nil
		@FLibraryName = nil
		@FSampleDir = nil
		
		@FNGSLibraryIndexerReport = nil # reference to MeteorSession member
		@FMismatchesCount = nil
		@FMatchesCount = nil
		@FIsMismatchesPercentage = nil
		@FReferenceName = nil
		@FSequenceFileFormat = nil
		@FIsColorSpaceMapping = false
	
		@FMappingProgram = nil
		@FMappingFileFormat = nil
		@FParametersShortLine = nil
		@FMapperCmd = nil
		@FMappingOutputFileNames = []
		@FTmpMappingOutputFileNames = []
		@FIsReadyForMapping = 0
		@FBowtieProgram = nil
		@FQualityEncoding = nil
	end
	
	#------------------------------------------------------------------------------

	# other constructor
	def initialize(aMappingProg, aLibraryCensusIniFileName, aReferenceIniFileName, aLibraryIndexerReport, aMappedCensusIniFileName)
				
		@FLibraryCensusIniFile = IniFile.load(aLibraryCensusIniFileName)
		@FMappedCensusIniFileName = aMappedCensusIniFileName
		@FReferenceIniFile = IniFile.load(aReferenceIniFileName)
		
		if @FReferenceIniFile.nil?
			STDERR.puts "Error in MeteorMapper::initialize while loading #{aReferenceIniFileName}"
			exit 1
		end

		@FLibraryMappingDir = nil
		@FTmpLibraryMappingDir = nil
		@FLibraryName = @FLibraryCensusIniFile[C_SAMPLE_INFO_SECTION][C_SAMPLE_FULL_NAME_STR]
		@FSampleDir = File.dirname(@FLibraryCensusIniFile.filename)
	
		@FNGSLibraryIndexerReport = aLibraryIndexerReport # reference to MeteorSession member
		@FMismatchesCount = nil
		@FMatchesCount = nil
		@FIsMismatchesPercentage = nil
		@FReferenceName = @FReferenceIniFile[C_REFERENCE_INFO_SECTION][C_REFERENCE_NAME_STR]
		
		@FQualityEncoding = @FLibraryCensusIniFile[C_SAMPLE_INFO_SECTION][C_QUALITY_ENCODING_STR]
		
		if (@FLibraryCensusIniFile[C_SAMPLE_INFO_SECTION][C_SEQUENCING_DEVICE_STR] == C_SOLID_DEVICE_STR)
			@FSequenceFileFormat = Sff[CS]
			@FIsColorSpaceMapping = true
			@FQualityEncoding = C_QE_SOLID if @FQualityEncoding.nil?
		else
			@FSequenceFileFormat = Sff[FQ]
			@FIsColorSpaceMapping = false
		end
		
		@FMappingProgram = aMappingProg
		@FMappingFileFormat = nil
		@FParametersShortLine = nil
		@FMapperCmd = nil
		@FMappingOutputFileNames = []
		@FTmpMappingOutputFileNames = []
		@FIsReadyForMapping = 0

	end
	
	#------------------------------------------------------------------------------
	
	def FinalizeMapping() #### TODO
		# Create new mapping_census_ini_file from library_census_ini_file
		aMappedCensusIniFile = @FLibraryCensusIniFile.clone
		aMappedCensusIniFile.filename = @FMappedCensusIniFileName

		aSampleInfoSection = aMappedCensusIniFile[C_SAMPLE_INFO_SECTION]    # reference
		aMappingSection = aMappedCensusIniFile[C_MAPPING_SECTION]           # reference
		aMappingFileSection = aMappedCensusIniFile[C_MAPPING_FILE_SECTION]  # reference

		aMappedCensusIniFile[C_MAPPING_SECTION][C_MAPPING_TOOL_STR] = @FMappingProgram
		aMappedCensusIniFile[C_MAPPING_SECTION][C_MAPPING_TOOL_VERSION_STR] = 'NA'

		aMappedCensusIniFile[C_MAPPING_FILE_SECTION][C_MAPPING_FILE_COUNT_STR] = @FMappingOutputFileNames.size

		# should iterate only once
		@FMappingOutputFileNames.size.times do |i|
			aMappedCensusIniFile[C_MAPPING_FILE_SECTION][C_MAPPING_BOWTIE_FILENAME_STR+"_#{i+1}"] = File.basename(@FMappingOutputFileNames[i])
		end

		#inherited;
		aSampleInfoSection[C_INDEXED_READ_LENGTH_STR] = @FNGSLibraryIndexerReport["IndexedReadLength"]
		aSampleInfoSection[C_SAMPLE_TOTAL_READ_COUNT_STR] = @FNGSLibraryIndexerReport["IndexedReadCount"]
		aSampleInfoSection[C_SAMPLE_INDEXED_READ_COUNT_STR] = @FNGSLibraryIndexerReport["IndexedReadCount"]
		aSampleInfoSection[C_SAMPLE_INDEXED_BASE_COUNT_STR] = @FNGSLibraryIndexerReport["IndexedBaseCount"]
		aSampleInfoSection[C_DATABASE_TYPE_STR] = C_METEOR_DBTYPE_UNKNOWN_STR # "unknown"
		aSampleInfoSection[C_SAMPLE_IS_DATA_PREFIXED_STR] = 0  # to avoid too long path

		aSampleInfoSection[C_SAMPLE_INDEXING_SOFTWARE_STR] =  C_METEOR_PROGRAM_NAME
		aSampleInfoSection[C_SAMPLE_INDEXING_SOFTWARE_VERSION_STR] = C_METEOR_VERSION

		aSampleInfoSection[C_CENSUS_STATUS_STR] = C_CENSUS_STATUS_MAPPED;
		aMappingSection[C_MAPPING_DATE_STR] = Date.today.to_s
		aMappingSection[C_REFERENCE_NAME_STR] = @FReferenceName;
		aMappingSection[C_MAPPING_CMDLINE_STR] = @FMapperCmd;
		aMappingSection[C_MAPPING_PARAMETERS_STR] = @FParametersShortLine;
		aMappingSection[C_MAPPED_READ_LENGTH_STR] = @FNGSLibraryIndexerReport["MappedReadLength"]
		aMappingSection[C_MAPPED_READ_LENGTH_TYPE_STR] = @FNGSLibraryIndexerReport["MappedReadLengthType"]
		aMappingSection[C_MAPPING_MISMATCHES_STR] = @FMismatchesCount;
		aMappingSection[C_MAPPING_IS_MISMATCHES_PERCENTAGE_STR] = @FIsMismatchesPercentage;
		aMappingSection[C_MAPPING_MATCHES_STR] = @FMatchesCount;
		aMappingSection[C_MAPPING_IS_LOCAL_MAPPING_STR] = @FIsLocalMapping # LOCAL

		aMappingSection[C_MAPPING_SOFTWARE_STR] = C_METEOR_MAPPING_PROGRAM_NAME
		aMappingSection[C_MAPPING_SOFTWARE_VERSION_STR] = C_METEOR_VERSION

		aMappingSection[C_PROCESSED_READ_COUNT_STR] = @FNGSLibraryIndexerReport["IndexedReadCount"]

		aMappingFileSection[C_MAPPING_FILE_FORMAT_STR] = @FMappingFileFormat;

		# write file to disk
		aMappedCensusIniFile.save

		#aMappedCensusIniFile = nil
	end
	
	#------------------------------------------------------------------------------
	
	def Bowtie2MapRead()
		### prepare bowtie2 command
		# -q --no-head --no-sq --no-unal --omit-sec-seq --end-to-end --sensitive -k 1 -S $sample/$sample.sam -U $fastq_list -p 40 -x $bt2_index
		# references repository on IBS server: /services/projects/biodatabank/meteor-data/reference

		aBowtieProgram = nil
		
		# since bowtie-build 64bits
		aIsLargeReference = @FReferenceIniFile["bowtie2_index"]["is_large_reference"]
		#### get around KCL feedback (inifile bug?)
		aIsLargeReference = integerable(aIsLargeReference) == false ? 1 : aIsLargeReference.to_i
		if aIsLargeReference == 1
			aBowtieProgram = C_BOWTIE2_LARGE_MAPPER
		else
			aBowtieProgram = C_BOWTIE2_SMALL_MAPPER
		end
		
		# aBowtieIndexCount should be 1
		aBowtieIndexList = []
		aBowtieIndexCount = @FReferenceIniFile["bowtie2_index"]["index_count"]
		#### get around KCL feedback (inifile bug?) : loading as String instead of Fixnum => force .to_i
		aBowtieIndexCount = integerable(aBowtieIndexCount) == false ? 1 : aBowtieIndexCount.to_i

		for aIndex in 1..aBowtieIndexCount do
			aBowtieIndexList.push(@FReferenceIniFile["bowtie2_index"]["dna_space_bowtie_index_prefix_name_#{aIndex}"])
		end
		
		# bowtie2 parameters
		aParameters = "--mm -p #{@FCPUCount}"
		if @FMapperCmd.nil?
			if @FIsLocalMapping == 1
				aParameters += " --local --sensitive-local"
			else
				aParameters += " --end-to-end --sensitive"
			end
			if @FMatchesCount > 0 
				aParameters += " -k #{@FMatchesCount}"
			elsif @FMatchesCount == 0 
				aParameters += " -a"
			end
		else
			aParameters += " #{@FMapperCmd}"
		end
		
		# add parameters for taking into account quality
		case @FQualityEncoding
			when /#{C_QE_SOLEXA}|#{C_QE_ILLUMINA_13}|#{C_QE_ILLUMINA_15}/
				aParameters += ' --phred64'
			when /#{C_QE_ILLUMINA_18}|#{C_QE_PROTON}|#{C_QE_SANGER}/
				aParameters += ' --phred33'
		end
		
		# force sam format
		@FMappingFileFormat = C_MAPPING_SAM_OUTPUT_STR
		aParameters += " "+C_SAM2_OUTPUT_PARAMETERS

		# complete command with input/output parameters, should iterate only once (because only one index)
		aBowtieIndexCount.times do |iBowtieIndex|

			# ex: @FLibraryMappingDir: /projects/parkinson/maping/KCL_01/mapping_vs_hs_9_9_metahit_Proton1-520-Parkinson_R1a_KCL_01
			#     @FLibraryName + : Proton1-520-Parkinson_R1a_KCL_01
			#     iBowtieIndex+1 : 1
			@FMappingOutputFileNames.push(@FLibraryMappingDir+C_PATH_SEP+@FLibraryName+"_#{(iBowtieIndex+1)}"+C_BOWTIE2_OUTPUT_EXT)
			@FTmpMappingOutputFileNames.push(@FTmpLibraryMappingDir+C_PATH_SEP+@FLibraryName+"_#{(iBowtieIndex+1)}"+C_BOWTIE2_OUTPUT_EXT) #### NEW

			# reindexed fastq
			#aIndexedInputFileName = @FSampleDir+C_PATH_SEP+@FLibraryCensusIniFile[C_SAMPLE_FILE_SECTION][C_FASTQ_FILENAME_STR]+C_INDEXED_EXT
			aIndexedInputFileName = @FNGSLibraryIndexerReport['IndexedfastqFilePath']
			
			aCmd = aBowtieProgram + " " + aParameters +
			 " -x #{File.dirname(@FReferenceIniFile.filename)}/#{@FReferenceIniFile["reference_file"]["fasta_dir"]}/#{aBowtieIndexList[iBowtieIndex]}" +
			 " -U #{aIndexedInputFileName}" +
			 " -S #{@FTmpMappingOutputFileNames[iBowtieIndex]}" #### NEW
			 ####" -S #{@FMappingOutputFileNames[iBowtieIndex]}"
			
			### Execute Command
			puts "Executing mapping command:"; puts aCmd	
			aExitStatus = ExecuteCmd(aCmd)
	
			if aExitStatus != 0
				STDERR.puts "Error in MeteorMapper::MapRead, mapping command exited with non zero status"
				exit aExitStatus
			end
			
		end
		
	end
	
	#------------------------------------------------------------------------------
	
	def BowtieMapRead()
		### prepare bowtie command
		# bowtie-align-l --mm -C --col-cseq -p CPUCount -v MismatchesCount -k MatchesCount --integer-quals -Q file.qual --sam --sam-nohead --sam-nosq catalogBowtieIndex -f csfastaFile output.sam
		# references repository on IBS server: /services/projects/biodatabank/meteor-data/reference
		
		# bowtie native method is faster and output.bowtie is smaller than output.sam 
		# bowtie-align-l --mm -C --col-cseq -p CPUCount -v MismatchesCount -k MatchesCount --integer-quals -Q file.qual --suppress 6,7 catalogBowtieIndex -f csfastaFile(.fa ??) output.bowtie

		aBowtieProgram = nil
		
		# since bowtie-build 64bits
		aIsLargeReference = @FReferenceIniFile["bowtie_index"]["is_large_reference"]
		#### get around KCL feedback (inifile bug?)
		aIsLargeReference = integerable(aIsLargeReference) == false ? 1 : aIsLargeReference.to_i
		if aIsLargeReference == 1
			aBowtieProgram = C_BOWTIE_LARGE_MAPPER
		else
			aBowtieProgram = C_BOWTIE_SMALL_MAPPER
		end
		
		# aBowtieIndexCount should be 1
		aBowtieIndexList = []
		aBowtieIndexCount = @FReferenceIniFile["bowtie_index"]["index_count"]
		#### get around KCL feedback (inifile bug?) : loading as String instead of Fixnum => force .to_i
		aBowtieIndexCount = integerable(aBowtieIndexCount) == false ? 1 : aBowtieIndexCount.to_i

		for aIndex in 1..aBowtieIndexCount do
			aBowtieIndexList.push(@FReferenceIniFile["bowtie_index"]["color_space_bowtie_index_prefix_name_#{aIndex}"])
		end
		
		### bowtie parameters
		# set --trim3 option if read_length > 35
		aTrim3Value = @FNGSLibraryIndexerReport['IndexedReadLength'] - @FNGSLibraryIndexerReport['MappedReadLength']
		aTrim3Option = aTrim3Value > 0 ? " --trim3 #{aTrim3Value}" : ""
		#aParameters = "--mm -f -C --col-cseq -p #{@FCPUCount}#{aTrim3Option}" # NB: --col-cseq is ignored if --sam
		aParameters = "--mm -f -C --col-cseq -p #{@FCPUCount}" # NB: --col-cseq is ignored if --sam
		if @FMapperCmd.nil? # take default parameters
			aParameters += aTrim3Option
			# end-to-end => -v mismatches
			aParameters += " -v #{@FMismatchesCount}"
			if @FMatchesCount > 0 
				aParameters += " -k #{@FMatchesCount}"
			elsif @FMatchesCount == 0
				aParameters += " -a"
			end
		else
			# might be: "-v MismatchesCount -k MatchesCount"
			#    or   : "-v MismatchesCount -a"
			aParameters += " #{@FMapperCmd}"
		end
		# add qual file (always present ?)
		aParameters += ' --integer-quals -Q ' + @FNGSLibraryIndexerReport['IndexedQualFilePath']

		if @FMappingFileFormat == C_MAPPING_SAM_OUTPUT_STR
			aParameters += " #{C_SAM_OUTPUT_PARAMETERS}";
			#~ mffLiteDefaultBowtie: #~ aParameters += " #{C_LITE_OUTPUT_PARAMETERS}";
		else
			abort("Error, mapping file format #{@FMappingFileFormat} is not supported.")
		end
		
		# complete command with input/output parameters, should iterate only once (because only one index)
		aBowtieIndexCount.times do |iBowtieIndex|

			# ex: @FLibraryMappingDir: /projects/parkinson/maping/KCL_01/mapping_vs_hs_9_9_metahit_Proton1-520-Parkinson_R1a_KCL_01
			#     @FLibraryName + : Proton1-520-Parkinson_R1a_KCL_01
			#     iBowtieIndex+1 : 1
			@FMappingOutputFileNames.push(@FLibraryMappingDir+C_PATH_SEP+@FLibraryName+"_#{(iBowtieIndex+1)}"+C_BOWTIE_OUTPUT_EXT)
			@FTmpMappingOutputFileNames.push(@FTmpLibraryMappingDir+C_PATH_SEP+@FLibraryName+"_#{(iBowtieIndex+1)}"+C_BOWTIE_OUTPUT_EXT) #### NEW

			# reindexed fastq
			#aIndexedInputFileName = @FSampleDir+C_PATH_SEP+@FLibraryCensusIniFile[C_SAMPLE_FILE_SECTION][C_FASTQ_FILENAME_STR]+C_INDEXED_EXT
			aIndexedInputFileName = @FNGSLibraryIndexerReport['IndexedCsfastaFilePath']
			
			aCmd = aBowtieProgram + " " + aParameters +
			 " #{File.dirname(@FReferenceIniFile.filename)}/#{@FReferenceIniFile["reference_file"]["fasta_dir"]}/#{aBowtieIndexList[iBowtieIndex]}" +
			 " #{aIndexedInputFileName}" +
			 " #{@FTmpMappingOutputFileNames[iBowtieIndex]}" #### NEW
			 ####" #{@FMappingOutputFileNames[iBowtieIndex]}"
			
			### Execute Command
			puts "Executing mapping command:"; puts aCmd	
			aExitStatus = ExecuteCmd(aCmd)
	
			if aExitStatus != 0
				STDERR.puts "Error in MeteorMapper::MapRead, mapping command exited with non zero status"
				exit aExitStatus
			end
			
		end
		
	end
	
	#------------------------------------------------------------------------------
	
	# return true if mapping success or already done, false if error.
	# aMappingOutputRootDir: /path/to/project_name/mapping/sample_name/
	def MapRead(
		           aTmpDir, aLibraryMappingDir, aMapperCmd,
		           aMatchesCount, aMismatchesCount,
		           aIsMismatchesPercentage,
		           aIsLocalMapping,
		           aMappingFileFormat,
		           aIsCPUPercentage,
		           aCPUCount
		         )
		#### get around KCL feedback (inifile bug?)
		@FMatchesCount = integerable(aMatchesCount) == false ? -1 : aMatchesCount.to_i
		#@FMatchesCount = C_DEFAULT_MAPPING_MATCHES_COUNT if (@FMatchesCount == 0)	

		# set CPU count
		aHostCPUCount = cpu_count()

		#### get around KCL feedback bug(?)
		aCPUCount = integerable(aCPUCount) == false ? -1 : aCPUCount.to_i
		aIsCPUPercentage = integerable(aIsCPUPercentage) == false ? 0 : aIsCPUPercentage.to_i
		if aCPUCount == -1
			@FCPUCount = [1, aHostCPUCount-1].max
		else
			@FCPUCount = aIsCPUPercentage == 1 ? [1, (aHostCPUCount*aCPUCount/100).to_int - 1].max : aCPUCount
		end
	  
		#### get around KCL feedback bug(?). Default FMismatchesCount is 3
		@FMismatchesCount = integerable(aMismatchesCount) == false ? 3 : aMismatchesCount.to_i
		#### idem
		@FIsMismatchesPercentage = integerable(aIsMismatchesPercentage) == false ? 0 : aIsMismatchesPercentage.to_i
		@FIsLocalMapping = integerable(aIsLocalMapping) == false ? 0 : aIsLocalMapping.to_i # LOCAL
		@FParametersShortLine = "#{C_MAPPED_READ_LENGTH_CHAR}#{@FNGSLibraryIndexerReport['MappedReadLength']}-#{C_MISMATCHES_COUNT_CHAR}#{@FMismatchesCount}"
		@FMappedReadCount = 0
		@FMappingFileFormat = aMappingFileFormat
		@FMapperCmd = (aMapperCmd.is_a?(String) and aMapperCmd.empty?) ? nil : aMapperCmd
		
		@FTmpLibraryMappingDir = @FLibraryMappingDir = aLibraryMappingDir
		# FTmpLibraryMappingDir: directory for temporary mapping SAM files.
		# if aTmpDir is set : "~/scratch/tmp_meteor/md5hash" + "/projects/projectname/mapping/sampleN/mapping_vs_blabla"
		# else :              "/projects/projectname/mapping/sampleN/mapping_vs_blabla"
		@FTmpLibraryMappingDir = aTmpDir + C_PATH_SEP + aLibraryMappingDir if not aTmpDir.nil? #### NEW
		
		# create FLibraryMappingDir and FTmpLibraryMappingDir
		FileUtils.mkdir_p(@FLibraryMappingDir) if not File.exists?(@FLibraryMappingDir)
		FileUtils.mkdir_p(@FTmpLibraryMappingDir) if not File.exists?(@FTmpLibraryMappingDir)
		
		@FIsReadyForMapping = 1
		
		# call BowtieMapRead() or Bowtie2MapRead() according to @FMappingProgram
		if @FMappingProgram == C_BOWTIE_MAPPER
			BowtieMapRead()
		else # C_BOWTIE2_MAPPER
			Bowtie2MapRead()
		end
		
		FinalizeMapping()

		return true

	end
	
	#------------------------------------------------------------------------------
	
end ### END class MeteorMapper


################################################################################

class MeteorSession

	### member attributs

	attr_accessor :FMeteorJobIniFile
	attr_accessor :FProjectDir
	attr_accessor :FSampleDir         # /path/to/sample/sampleN/
	attr_accessor :FTmpSampleDir
	attr_accessor :FTmpDir
	attr_accessor :FNoLock            # boolean
	attr_accessor :FForce             # boolean
	attr_accessor :FMappingDone       # boolean
	attr_accessor :FCountingDone      # boolean
	
	attr_accessor :FCountingTypeList  # array

	#~ attr_accessor :LibraryIniFileName
	attr_accessor :FSampleMappingDir  # /path/to/mapping/sampleN/

	#FRepositoryLibraryIniFileNameList
	attr_accessor :FLibraryIniFileNames
	attr_accessor :FLibraryCount

	attr_accessor :FSampleName  # from LibraryIniFile (overwritten for each libray)
	attr_accessor :FProjectName # idem

	attr_accessor :FProjectMappingDir # /path/to/mapping/

	# all census file names
	attr_accessor :FLibraryCensusIniFileNames         # array of LibraryCensusIniFile
	attr_accessor :FMainMappingCensusIniFileNames     # array of MappingCensusIniFile
	attr_accessor :FExcludedMappingCensusIniFileNames # 2D array of MappingCensusIniFiles (1st dim = library, 2nd = excluded ref)

	### member methods
  
	# Constructor
	def initialize()
		@FMeteorJobIniFile = nil
		@FProjectDir = nil
		@FSampleDir = nil
		@FTmpSampleDir = nil
		@FTmpDir = nil
		@FNoLock = false
		@FForce = false
		@FMappingDone = true
		@FCountingDone = true
		@FCountingTypeList = nil
		@FMappingDir = nil
		@LibraryNames = nil
		@FLibraryCount = nil
		@FSampleName = nil
		@FProjectName = nil
		@FProjectMappingDir = nil
		@FLibraryIndexerReport = {}
		@FLibraryIniFileNames = []
		@FLibraryCensusIniFileNames = []          # array of LibraryCensusIniFile
		@FMainMappingCensusIniFileNames = []      # array of MappingCensusIniFile
		@FExcludedMappingCensusIniFileNames = []  # 2D array of MappingCensusIniFiles
	end

	#------------------------------------------------------------------------------

	def PrepareWorkSpace(options)
		#~ CleanWorkSpace
	end
	
	#------------------------------------------------------------------------------
	#### TODO in MeteorSession :
	####       - interdire csfasta/qual + bowtie2
	####       - si csfasta/qual : 
	####         => CountAndReindexCsfastaQual() OK
	####         => @FLibraryIndexerReport["MappedReadLengthType"] = fixed OK
	####         => @FLibraryIndexerReport["MappedReadLength"] = 35 OK
	####       - faire une fonction BowtieMapRead et Bowtie2MapRead appelée par MapRead OK
	
	#------------------------------------------------------------------------------
  
	# count fastq reads
	def CountReadFastqFile(aInputFile)
		
		aLineCount = aReadCount = aBaseCount = 0
		fd = nil
		in_fq = nil
		# open input fastq in read mode
		begin
			fd = File.open(aInputFile)
			# fastq may be gziped
			in_fq = aInputFile =~ /gz$/ ? Zlib::GzipReader.new(fd) : fd
		rescue
			STDERR.puts "Error: cannot open #{aInputFile} !"
			exit 1
		end

		# read input fastq line by line
		in_fq.each do |line|
			aLineCount += 1
			aModulo = aLineCount % 4
			#if line =~ /^@/
			if aModulo == 1 # read id line
				aReadCount += 1
			else
				aBaseCount += line.size - 1 if aModulo == 0 # size - 1 because of line ending character
			end
		end
		in_fq.close if not in_fq.closed?
			
		return [aReadCount, aBaseCount]
	end
	
	#------------------------------------------------------------------------------
	
	def CountReadAndReIndexCsfastaQualFiles(aInputCsFastaFile, aInputQualFile, aOutputCsFastaFile, aOutputCsQualFile)
		
		aReadCount = aBaseCount = 0

		qual_fd = cs_fd = qual_out_fd = cs_out_fd = nil
		in_qual_fd = in_cs_fd = nil
		# open input and output files
		begin
			qual_fd = File.open(aInputQualFile)
			cs_fd = File.open(aInputCsFastaFile)
			# input files may be gzipped
			in_qual_fd = aInputQualFile =~ /gz$/ ? Zlib::GzipReader.new(qual_fd) : qual_fd
			in_cs_fd = aInputCsFastaFile =~ /gz$/ ? Zlib::GzipReader.new(cs_fd) : cs_fd
			
			qual_out_fd = File.open(aOutputCsQualFile, "w")
			cs_out_fd = File.open(aOutputCsFastaFile, "w")
		rescue Exception => e
			STDERR.puts("Error while opening input files.")
			#STDERR.puts(e.backtrace.inspect)
			abort(e.message)
		end
		
		# read each .csfasta header
		in_cs_fd.each do |header_cs|
		
			# skip potential commentaries above the csfasta header
			#~ next if header_cs[0,1] != '>'
			next if not header_cs.start_with?('>')
			# skip potential commentaries above the quality header
			loop do
				header_q = in_qual_fd.readline
				#~ break if header_q[0,1] == '>'
				break if header_q.start_with?('>')
			end		
			# abort if different (not usefull test as reads have been re-indexed)
			#abort("Error while reading solid data: csfasta and qual headers differs: #{header_cs} != #{header_q}") if header_q != header_cs
			
			cs_seq = in_cs_fd.readline.chomp			
			aReadCount += 1
			aBaseCount += cs_seq.size - 1 # -1 because leading T
			
			cs_out_fd.puts ">#{aReadCount}\n#{cs_seq}"
			qual_out_fd.puts ">#{aReadCount}\n#{in_qual_fd.readline.chomp}"
			
		end
		qual_fd.close if not qual_fd.closed?
		cs_fd.close if not cs_fd.closed?
		qual_out_fd.close if not qual_out_fd.closed?
		cs_out_fd.close if not cs_out_fd.closed?
		
		return [aReadCount, aBaseCount]
	end

	#------------------------------------------------------------------------------
	
	def CountReadAndReIndexFastqFile(aInputFile, aOutputFile)
	
		aInputLineCount = aReadCount = aBaseCount = aCpt = 0
		adn = qual = prev = seqId = nil
		fd = nil
		in_fq = nil
		# open input fastq in read mode
		begin
			fd = File.open(aInputFile)
			# fastq may be gzipped, xz or bz2
			in_fq = aInputFile =~ /gz$/ ? Zlib::GzipReader.new(fd) : fd
		rescue
			STDERR.puts "Error: cannot open #{aInputFile} !"
			exit 1
		end

		# open output fastq in write mode
		File.open(aOutputFile, "w") do |out_fq|
		
			# read input fastq line by line
			in_fq.each do |line|
				line.chomp!
				# line begins with @seqid
				if (line =~ /^@(.+)$/)
					## check if the previous read is valid (quality exists and same size as adn)
					aQualSize = (qual.nil?) ? 0 : qual.size
					if (aQualSize > 0)
						if ( (not adn.nil?) and aQualSize == adn.size)
							aReadCount += 1
							aBaseCount += aQualSize
							# then write this indexed read
							out_fq.print "@#{aReadCount}\n#{adn}\n+\n#{qual}\n"
							adn = nil
						end
					end
					seqId = $1
					## NB: quality line might start with @
				end
				aCpt += 1
				if (line == "+" or line == "+#{seqId}")
					qual = nil
					aCpt = 3
					# previous line was adn
					adn = prev
				end
				if aCpt == 4
					qual = line
				end
				prev = line
			end
			fd.close if not fd.closed?
			#in_fq.close if not in_fq.closed?
			
			# evaluate last read
			aQualSize = (qual.nil?) ? 0 : qual.size
			if (aQualSize > 0)
				if ( (not adn.nil?) and aQualSize == adn.size)
					aReadCount += 1
					aBaseCount += aQualSize
					out_fq.print "@#{aReadCount}\n#{adn}\n+\n#{qual}\n"
				end
			end
		end
		
		return [aReadCount, aBaseCount]
	end

	#------------------------------------------------------------------------------
	
	def TaskMainMapping(iLibrary)
		
		aWorkSessionSection = @FMeteorJobIniFile[C_JOB_WORKSESSION_SECTION] # reference
		aReferenceSection = @FMeteorJobIniFile[C_JOB_MAIN_REFERENCE_SECTION] # reference
		aReferenceIniFileName = aWorkSessionSection[C_JOB_REFERENCE_DIR] + C_PATH_SEP + aReferenceSection[C_JOB_REFERENCE_NAME] +
		                        C_PATH_SEP + aReferenceSection[C_JOB_REFERENCE_NAME] + C_REFERENCE_INIFILE_EXT

		puts "\n######### Task main mapping #{iLibrary}"
		
		aMappedCensusIniFileName = @FMainMappingCensusIniFileNames[iLibrary][:Stage1FileName]
		aLibraryMappingDir = @FMainMappingCensusIniFileNames[iLibrary][:directory]
		
		aLockedMappedCensusIniFileName = aMappedCensusIniFileName+C_LOCK_EXT
		
		# remove lock file if asked
		if (@FNoLock == true and File.exists?(aLockedMappedCensusIniFileName))
			puts "\nINFO: removing lock file: #{aLockedMappedCensusIniFileName}"
			FileUtils.rm(aLockedMappedCensusIniFileName)
		end
		# remove census 1 file if asked
		if (@FForce == true and File.exists?(aMappedCensusIniFileName))
			puts "\nINFO: removing existing census 1 file: #{aMappedCensusIniFileName}"
			FileUtils.rm(aMappedCensusIniFileName)
		end
		
		if ( File.exists?(aMappedCensusIniFileName) or File.exists?(aLockedMappedCensusIniFileName) )
			puts "\nINFO: skipping library.\nINFO: "+aMappedCensusIniFileName+" already exists or is locked"
			return true
		end

		# create lock file in result directory
		FileUtils.mkdir_p(aLibraryMappingDir) if not File.exists?(aLibraryMappingDir)
		FileUtils.touch(aLockedMappedCensusIniFileName)
		
		aMappingProg = aWorkSessionSection[C_JOB_MAPPING_PROGRAM]
		if (aMappingProg != C_BOWTIE2_MAPPER and aMappingProg != C_BOWTIE_MAPPER) 
			STDERR.puts "Error, unknown or unsupported mapping program : #{aWorkSessionSection[C_JOB_MAPPING_PROGRAM]}"
			exit 1 # C_BAD_MAPPER_EXIT
		end

		#~ aMeteorMapper = MeteorMapper.new(@FLibraryIniFileNames[iLibrary], aReferenceIniFileName, @FLibraryIndexerReport, aMappedCensusIniFileName)
		aMeteorMapper = MeteorMapper.new(aMappingProg, @FLibraryIniFileNames[iLibrary], aReferenceIniFileName, @FLibraryIndexerReport, aMappedCensusIniFileName)

		aOkMappingProcess = aMeteorMapper.MapRead(
					@FTmpDir,
					aLibraryMappingDir,
					aReferenceSection[C_JOB_MAPPER_CMD],
					aReferenceSection[C_JOB_MAPPING_MATCHES],
					aReferenceSection[C_JOB_MAPPING_MISMATCHES],
					aReferenceSection[C_JOB_IS_MAPPING_MISMATCHES_PERCENTAGE],
					aReferenceSection[C_JOB_IS_LOCAL_MAPPING],
					aWorkSessionSection[C_JOB_MAPPING_FILE_FORMAT],
					aWorkSessionSection[C_JOB_IS_CPU_PERCENTAGE],
					aWorkSessionSection[C_JOB_CPU_COUNT]
		)
		
		# remove census 1 lock file
		FileUtils.rm(aMappedCensusIniFileName+C_LOCK_EXT)
		
		return aOkMappingProcess
	end
	
	#------------------------------------------------------------------------------
	
	def TaskExcludedMapping(iLibrary, iExcluded)
		
		aWorkSessionSection = @FMeteorJobIniFile[C_JOB_WORKSESSION_SECTION] # reference
		aReferenceSection = @FMeteorJobIniFile[C_JOB_EXCLUDED_REFERENCE_PREFIX_SECTION+"#{iExcluded+1}"] # reference
		#puts "aWorkSessionSection:"
		#puts aWorkSessionSection
		#puts "aReferenceSection:"
		#puts aReferenceSection
		
		aReferenceIniFileName = aWorkSessionSection[C_JOB_REFERENCE_DIR] + C_PATH_SEP + aReferenceSection[C_JOB_REFERENCE_NAME] +
		                        C_PATH_SEP + aReferenceSection[C_JOB_REFERENCE_NAME] + C_REFERENCE_INIFILE_EXT
		
		puts "\n######### TaskExcludedMapping(#{iLibrary}, #{iExcluded+1}) ; reference name: "+aReferenceSection[C_JOB_REFERENCE_NAME]
		
		aMappedCensusIniFileName = @FExcludedMappingCensusIniFileNames[iLibrary][iExcluded][:Stage1FileName]
		aLibraryMappingDir = @FExcludedMappingCensusIniFileNames[iLibrary][iExcluded][:directory]
		
		aLockedMappedCensusIniFileName = aMappedCensusIniFileName+C_LOCK_EXT
		
		# remove lock file if asked
		if (@FNoLock == true and File.exists?(aLockedMappedCensusIniFileName))
			puts "\nINFO: removing lock file: #{aLockedMappedCensusIniFileName}"
			FileUtils.rm(aLockedMappedCensusIniFileName)
		end
		# remove census 1 file if asked
		if (@FForce == true and File.exists?(aMappedCensusIniFileName))
			puts "\nINFO: removing existing census 1 file: #{aMappedCensusIniFileName}"
			FileUtils.rm(aMappedCensusIniFileName)
		end
		
		if ( File.exists?(aMappedCensusIniFileName) or File.exists?(aLockedMappedCensusIniFileName) )
			puts "\nINFO: skipping library.\nINFO: "+aMappedCensusIniFileName+" already exists or is locked"
			return true
		end
		
		# create lock file in result directory
		FileUtils.mkdir_p(aLibraryMappingDir) if not File.exists?(aLibraryMappingDir)
		FileUtils.touch(aLockedMappedCensusIniFileName)

		aMappingProg = aWorkSessionSection[C_JOB_MAPPING_PROGRAM]
		if (aMappingProg != C_BOWTIE2_MAPPER and aMappingProg != C_BOWTIE_MAPPER) 
			STDERR.puts "Error, unknown or unsupported mapping program : #{aWorkSessionSection[C_JOB_MAPPING_PROGRAM]}"
			exit 1 # C_BAD_MAPPER_EXIT
		end

		aMeteorMapper = MeteorMapper.new(aMappingProg, @FLibraryIniFileNames[iLibrary], aReferenceIniFileName, @FLibraryIndexerReport, aMappedCensusIniFileName)

		aOkMappingProcess = aMeteorMapper.MapRead(
				@FTmpDir,
				aLibraryMappingDir,
				aReferenceSection[C_JOB_MAPPER_CMD],
				aReferenceSection[C_JOB_MAPPING_MATCHES],
				aReferenceSection[C_JOB_MAPPING_MISMATCHES],
				aReferenceSection[C_JOB_IS_MAPPING_MISMATCHES_PERCENTAGE],
				false, # LOCAL, aways end-to-end for excluded ref
				aWorkSessionSection[C_JOB_MAPPING_FILE_FORMAT],
				aWorkSessionSection[C_JOB_IS_CPU_PERCENTAGE],
				aWorkSessionSection[C_JOB_CPU_COUNT]
		)
		
		# remove census 1 lock file
		FileUtils.rm(aMappedCensusIniFileName+C_LOCK_EXT)
		
		return aOkMappingProcess

	end
	
	#------------------------------------------------------------------------------
	
	def IsAllMappingDone()
		return false
	end
  
	#------------------------------------------------------------------------------

	def LoadJobWorkflow(aWorkflowFile)
		if not File.exists?(aWorkflowFile)
			STDERR.puts "Error, file #{aWorkflowFile} not found"
			exit 1
		end
		@FMeteorJobIniFile = IniFile.load(aWorkflowFile)
		
		# check if excluded reference count is correct
		aExcludedRefCount = @FMeteorJobIniFile[C_JOB_WORKSESSION_SECTION][C_JOB_EXCLUDED_REFERENCE_COUNT]
		aExcludedRefCount = integerable(aExcludedRefCount) == false ? 0 : aExcludedRefCount.to_i
		aExcludedRefCounted = @FMeteorJobIniFile.sections.grep(/^#{C_JOB_EXCLUDED_REFERENCE_PREFIX_SECTION}/).size
		
		if aExcludedRefCount != aExcludedRefCounted
			STDERR.puts "Error in file #{aWorkflowFile}:"
			STDERR.puts "#{C_JOB_EXCLUDED_REFERENCE_COUNT} is set to #{aExcludedRefCount}"
			STDERR.puts "but found #{aExcludedRefCounted} excluded_reference section(s) ..."
			exit 1
		end
		
	end

	#------------------------------------------------------------------------------

	def CountAndReIndexReads(aLibraryCensusIniFile)

		aSampleInfoSection = aLibraryCensusIniFile[C_SAMPLE_INFO_SECTION] # reference
		aSampleFileSection = aLibraryCensusIniFile[C_SAMPLE_FILE_SECTION] # reference
		
		if (aSampleInfoSection[C_SEQUENCING_DEVICE_STR] != C_SOLID_DEVICE_STR) # => fastq
			@FLibraryIndexerReport.clear
			
			aInputFile = @FSampleDir+C_PATH_SEP+aSampleFileSection[C_FASTQ_FILENAME_STR]
			aOutputFile = @FTmpSampleDir+C_PATH_SEP+aSampleFileSection[C_FASTQ_FILENAME_STR]+C_INDEXED_EXT
			
			aReadCount, aBaseCount = CountReadAndReIndexFastqFile(aInputFile, aOutputFile)
			
			@FLibraryIndexerReport = {
				"IsIndexed" => 1,
				"IndexedfastqFilePath" => aOutputFile,
				"IndexedFileNameExt" => C_INDEXED_EXT,
				"IndexedReadCount" => aReadCount,
				"IndexedBaseCount" => aBaseCount,
				"IndexedReadLength" => -1,
				"MappedReadLength" => -1,
				"MappedReadLengthType" => "overall"
			}
		else # => color space data (SOLiD) => 2 files, csfasta and qual
			@FLibraryIndexerReport.clear
			
			aInputCsFastaFile = @FSampleDir+C_PATH_SEP+aSampleFileSection[C_CSFASTA_FILENAME_STR]
			aInputQualFile    = @FSampleDir+C_PATH_SEP+aSampleFileSection[C_QUAL_FILENAME_STR]
			
			aOutputCsFastaFile = @FTmpSampleDir+C_PATH_SEP+aSampleFileSection[C_CSFASTA_FILENAME_STR]+C_INDEXED_EXT
			aOutputCsQualFile  = @FTmpSampleDir+C_PATH_SEP+aSampleFileSection[C_QUAL_FILENAME_STR]+C_INDEXED_EXT
			
			aReadCount, aBaseCount = CountReadAndReIndexCsfastaQualFiles(aInputCsFastaFile, aInputQualFile, aOutputCsFastaFile, aOutputCsQualFile)
			
			@FLibraryIndexerReport = {
				"IsIndexed" => 1,
				"IndexedCsfastaFilePath" => aOutputCsFastaFile,
				"IndexedQualFilePath" => aOutputCsQualFile,
				"IndexedFileNameExt" => C_INDEXED_EXT,
				"IndexedReadCount" => aReadCount,
				"IndexedBaseCount" => aBaseCount,
				"IndexedReadLength" => aSampleInfoSection[C_READ_LENGTH_STR],
				"MappedReadLength" => C_DEFAULT_MAPPING_READ_LENGTH, # 35
				"MappedReadLengthType" => "fixed"
			}
		end
	end
	
	#------------------------------------------------------------------------------
  
	def LaunchMapping(options)
		puts "\nLaunch mapping\n"
		aOKToContinue = true

		#~ ( build catalog reference db (bowtie2-build-l) )

		#### why not launching one bowtie process for all the libraries (see bowtie -U option in example below) ????
		#bt2_index="/projects/biodatabank/filtering_database/Homo_sapiens_2014_02_04/bowtie2/homo_sapiens.fna"
		#command_bt2="bowtie2 -q --no-head --no-sq --no-unal --omit-sec-seq --end-to-end --sensitive -k 1 -S $sample/$sample.sam -U $fastq_list -p 40 -x $bt2_index"

		aExcludedRefCount = @FMeteorJobIniFile[C_JOB_WORKSESSION_SECTION][C_JOB_EXCLUDED_REFERENCE_COUNT]
		aExcludedRefCount = integerable(aExcludedRefCount) == false ? 0 : aExcludedRefCount.to_i
		
		# LOOP ON EACH LIBRARY
		@FLibraryCount.times do |iLibrary|

			if ( File.directory?(@FProjectDir) and File.file?(@FLibraryIniFileNames[iLibrary]) )

				aLibraryCensusIniFile = IniFile.load(@FLibraryIniFileNames[iLibrary])
				aSampleInfoSection = aLibraryCensusIniFile[C_SAMPLE_INFO_SECTION] # reference
				aSampleFileSection = aLibraryCensusIniFile[C_SAMPLE_FILE_SECTION] # reference
				
				@FSampleName = aSampleInfoSection[C_SAMPLE_NAME_STR]
				@FProjectName = aSampleInfoSection[C_PROJECT_NAME_STR]
				aSampleLibraryName = aSampleInfoSection[C_SAMPLE_FULL_NAME_STR]
				
				puts "\n\n######### Meteor Mapping task description"
				puts "Sample name = " + @FSampleName 
				puts "Library name = " + aSampleLibraryName 
				puts "Project name = " + @FProjectName
				puts "Sequencing device = " + aSampleInfoSection[C_SEQUENCING_DEVICE_STR] 
				puts "Workflow = " + File.basename(@FMeteorJobIniFile.filename)
				
				if not File.exists?(@FLibraryIniFileNames[iLibrary]+C_LOCK_EXT)
					
					### reindexing this library reads and fill @FLibraryIndexerReport
					CountAndReIndexReads(aLibraryCensusIniFile)
					
					# MAPPING THIS LIBRARY ON MAIN REFERENCE
					aOKToContinue = TaskMainMapping(iLibrary);
					if not aOKToContinue
						STDERR.puts "Error, TaskMainMapping failed: " + @FLibraryIniFileNames[iLibrary]
						exit 1 #### TODO
					end
					
					# MAPPING THIS LIBRARY ON EACH EXCLUDED REFERENCE
					aExcludedRefCount.times do |iExcluded|
						aOKToContinue = TaskExcludedMapping(iLibrary, iExcluded);
						if not aOKToContinue
							STDERR.puts "Error, TaskExcludedMapping failed: " + @FLibraryIniFileNames[iLibrary]
							exit 1
						end
					end
					if (aSampleInfoSection[C_SEQUENCING_DEVICE_STR] == C_SOLID_DEVICE_STR)
						puts "removing re-indexed color space files:"
						puts @FLibraryIndexerReport["IndexedCsfastaFilePath"]
						puts @FLibraryIndexerReport["IndexedQualFilePath"]
						FileUtils.rm(@FLibraryIndexerReport["IndexedCsfastaFilePath"])
						FileUtils.rm(@FLibraryIndexerReport["IndexedQualFilePath"])
					else
						puts "removing re-indexed fastq file:"
						puts @FLibraryIndexerReport["IndexedfastqFilePath"]
						FileUtils.rm(@FLibraryIndexerReport["IndexedfastqFilePath"])
					end
				else
					puts "Task already done or on going";
				end
			end
		end
	end

	#------------------------------------------------------------------------------
	
	def LaunchCounting()
		
		puts "\nLaunch counting"
		# does not need census_stage_0.ini
		#-w /path/to/workflow_tutorial.ini -i /path/to/sample/H1 -p /path/to/project_name -m mapping

		aCountingProgram = C_METEOR_COUNTER
		
		aparameters = " -w #{@FMeteorJobIniFile.filename} -i #{@FSampleDir} -p #{@FProjectDir} -o #{File.basename(@FProjectMappingDir)}"
		aparameters += " -c #{@FCountingTypeList.to_a.join(',')}" if ( (not @FCountingTypeList.nil?) and (not @FCountingTypeList.empty?) )
		aparameters += " -f" if @FForce # force overwriting former profiling results done with same parameters
		aparameters += " -t #{@FTmpDir}" if not @FTmpDir.nil? # path to directory where mapping results (SAM files) are stored. #### NEW
		
		aCmd = aCountingProgram + aparameters

		### Execute Command
		puts "\nExecuting counting command:"; puts aCmd
		aExitStatus = ExecuteCmd(aCmd)
    
		if aExitStatus != 0
			STDERR.puts "Error in MeteorSession::LaunchCounting, counting command exited with non zero status"
			exit aExitStatus
		end
		
		return aExitStatus
	end
	
	#------------------------------------------------------------------------------
	### define the MeteorSession attributes, then launch mapping, then counting

	def ProcessJob(options) #options : Workflow, ProjectPath, InputPath, MappingBasename

		# directly from program arguments
		@FProjectDir = options["ProjectPath"]
		@FSampleDir  = options["InputPath"]
		# skip if not a directory
		return if not File.directory?(@FSampleDir)
		# TODO get it from census ini file, (@FProjectName too)
		@FSampleName = File.basename(@FSampleDir)
		@FMappingDone = true
		@FCountingDone = true
		@FCountingTypeList = options["CountingTypeList"].to_a
		
		@FNoLock = false
		@FForce = false
		
		@FNoLock = true if not options["NoLock"].nil?
		if not options["Force"].nil?
			@FForce = true
			@FNoLock = true
		end

		# ex: /projects/parkinson/mapping
		@FProjectMappingDir   = @FProjectDir+C_PATH_SEP+options["MappingBasename"]

		# load information from workflow.ini ; set @FMeteorJobIniFile
		LoadJobWorkflow(options["Workflow"])

		#### in PrepareWorkSpace ?
		# ex: /projects/parkinson/mapping/KCL_01
		@FSampleMappingDir = @FProjectMappingDir+C_PATH_SEP+@FSampleName
		# make directory (and parent) if needed
		FileUtils.mkdir_p(@FSampleMappingDir) if not File.exists?(@FSampleMappingDir) #### maybe obsolete cause done in TaskMainMapping and/or TaskExcludedMapping
		
		#### in PrepareLibraryIniFileNames ?
		@FLibraryIniFileNames = Dir.glob( "#{@FSampleDir}#{C_PATH_SEP}*_"+C_CENSUS_SEQUENCED_STR+".ini" ).sort
		@FLibraryCount        = @FLibraryIniFileNames.size
		if (@FLibraryCount == 0)
			STDERR.puts "Error, no #{C_CENSUS_SEQUENCED_STR}.ini file found in #{@FSampleDir}"
			exit 1
		end

		aMainReferenceSection = @FMeteorJobIniFile[C_JOB_MAIN_REFERENCE_SECTION] # alias
		aMappingDirPathPrefix = ""

		if (aMainReferenceSection[C_JOB_MAPPING_PREFIX_NAME].nil? or aMainReferenceSection[C_JOB_MAPPING_PREFIX_NAME].empty?)
			#aMappingDirPathPrefix = C_DEFAULT_SAMPLE_MAPPING_DIRECTORY+"_vs_"+aMainReferenceSection[C_JOB_REFERENCE_NAME]+"_#{C_MAPPED_READ_LENGTH_CHAR}"+
			#						aMainReferenceSection[C_JOB_MAPPED_READ_LENGTH]+"-"+C_MISMATCHES_COUNT_CHAR+aMainReferenceSection[C_JOB_MAPPING_MISMATCHES]+"_"
			aMappingDirPathPrefix = C_DEFAULT_SAMPLE_MAPPING_DIRECTORY+"_vs_"+aMainReferenceSection[C_JOB_REFERENCE_NAME]+"_#{C_MAPPED_READ_LENGTH_CHAR}"+
									"#{aMainReferenceSection[C_JOB_MAPPED_READ_LENGTH]}-#{C_MISMATCHES_COUNT_CHAR}#{aMainReferenceSection[C_JOB_MAPPING_MISMATCHES]}_"
		else
			aMappingDirPathPrefix = "#{aMainReferenceSection[C_JOB_MAPPING_PREFIX_NAME]}_"
		end    

		aExcludedRefCount = @FMeteorJobIniFile[C_JOB_WORKSESSION_SECTION][C_JOB_EXCLUDED_REFERENCE_COUNT]
		#### get around KCL feedback bug (inifile?)
		aExcludedRefCount = integerable(aExcludedRefCount) == false ? 0 : aExcludedRefCount.to_i

		aTmpHash = {}
		#### TODO in Prepare(Main)MappingIniFileNames (and PrepareExcludedMappingIniFileNames)
		@FLibraryIniFileNames.each do |aLibraryIniFileName|

			### main reference
			aLibraryCensusIniFile = IniFile.load(aLibraryIniFileName)
			aSampleInfo = aLibraryCensusIniFile[C_SAMPLE_INFO_SECTION] # reference
			aMappingDirPath = aMappingDirPathPrefix+aSampleInfo[C_SAMPLE_FULL_NAME_STR]

			aTmpHash.clear
			aTmpHash[:directory]      = aDirectory = @FProjectMappingDir+C_PATH_SEP+aSampleInfo[C_SAMPLE_NAME_STR]+C_PATH_SEP+aMappingDirPath
			aTmpHash[:Stage1FileName] = aDirectory+C_PATH_SEP+File.basename(aLibraryIniFileName).sub(C_CENSUS_SEQUENCED_STR, C_CENSUS_MAPPED_STR)
			@FMainMappingCensusIniFileNames.push(aTmpHash.clone) # push aTmpHash value (not reference) => clone
			@FMappingDone = false if (not File.exists?(aTmpHash[:Stage1FileName]))
		
			### excuded references
			@FExcludedMappingCensusIniFileNames.push([]) # reserve space for this library
			aExcludedRefCount.times { |iExcluded|
				aExcludedReference = @FMeteorJobIniFile[C_JOB_EXCLUDED_REFERENCE_PREFIX_SECTION+"#{iExcluded+1}"];
				
				if (aExcludedReference[C_JOB_MAPPING_PREFIX_NAME].nil? or aExcludedReference[C_JOB_MAPPING_PREFIX_NAME].empty?)
					#aMappingDirPath = C_DEFAULT_SAMPLE_MAPPING_DIRECTORY+"_vs_"+aExcludedReference[C_JOB_REFERENCE_NAME]+"_#{C_MAPPED_READ_LENGTH_CHAR}"+
					#				aExcludedReference[C_JOB_MAPPED_READ_LENGTH]+"-"+C_MISMATCHES_COUNT_CHAR+aExcludedReference[C_JOB_MAPPING_MISMATCHES]+"_"+aSampleInfo[C_SAMPLE_FULL_NAME_STR]
					aMappingDirPath = C_DEFAULT_SAMPLE_MAPPING_DIRECTORY+"_vs_"+aExcludedReference[C_JOB_REFERENCE_NAME]+"_#{C_MAPPED_READ_LENGTH_CHAR}"+
									"#{aExcludedReference[C_JOB_MAPPED_READ_LENGTH]}-#{C_MISMATCHES_COUNT_CHAR}#{aExcludedReference[C_JOB_MAPPING_MISMATCHES]}_#{aSampleInfo[C_SAMPLE_FULL_NAME_STR]}"
				else
					aMappingDirPath = "#{aExcludedReference[C_JOB_MAPPING_PREFIX_NAME]}_#{aSampleInfo[C_SAMPLE_FULL_NAME_STR]}"
				end

				aTmpHash.clear
				aTmpHash[:directory]      = aDirectory = @FProjectMappingDir+C_PATH_SEP+aSampleInfo[C_SAMPLE_NAME_STR]+C_PATH_SEP+aMappingDirPath
				aTmpHash[:Stage1FileName] = aDirectory+C_PATH_SEP+File.basename(aLibraryIniFileName).sub(C_CENSUS_SEQUENCED_STR, C_CENSUS_MAPPED_STR)
				@FExcludedMappingCensusIniFileNames.last.push(aTmpHash.clone) # push aTmpHash value (not reference) => clone
				@FMappingDone = false if (not File.exists?(aTmpHash[:Stage1FileName]))
			}
			
		end
		
		@FTmpDir = options["TmpPath"]
		# generate unique path by hashing (md5) samplename path and worflow.ini path
		aMD5 = Digest::MD5.hexdigest( "#{@FSampleDir} #{options["Workflow"]}" ) #### NEW
		@FTmpDir += C_PATH_SEP + "tmp_meteor" + C_PATH_SEP + aMD5 if not @FTmpDir.nil? #### NEW
				
		### launch mapping except if we want counting only
		if options["OnlyCounting"].nil?
			# mapping already done and no overwriting
			if (@FMappingDone and @FForce == false)
				puts "\nINFO: Mapping already done for sample: #{@FSampleName}"
				puts "INFO: Skipped !"
			else # mapping not done or we want to overwrite previous results
				# if option -t not provided, generate sample tmpdir in sample dir
				#@FTmpSampleDir = @FTmpDir.nil? ? @FSampleDir + C_PATH_SEP + aMD5 : @FTmpDir #### NEW
				@FTmpSampleDir = @FTmpDir.nil? ? @FSampleMappingDir + C_PATH_SEP + aMD5 : @FTmpDir #### NEW
				
				if File.exists?(@FTmpSampleDir)
					STDERR.puts "WARNING, temporary dir #{@FTmpSampleDir} already exists !"
				end
				FileUtils.mkdir_p(@FTmpSampleDir)
				
				LaunchMapping(options)
				# remove sample tmp data if created in sample dir
				FileUtils.rm_rf(@FTmpSampleDir) if ( @FTmpDir.nil? and File.exists?(@FTmpSampleDir) ) #### NEW
			end
		end
		### launch counting except if we want mapping only
		if options["OnlyMapping"].nil?
			LaunchCounting()
			# remove temporary data if @FTmpDir exists
			FileUtils.rm_rf(@FTmpDir) if ( (not @FTmpDir.nil?) and (File.exists?(@FTmpDir)) and (not options["CleanTmpPath"].nil?) ) #### NEW
		end

		puts "Done !\nJob finished without errors ..."

	end

end ### END class MeteorSession
################################################################################


################################################################################
#                                MAIN                                          #
################################################################################

if __FILE__ == $0

	# all -w /path/to/workflow.ini -i ~/path/to/sample_dir -p /path/to/projects -o mapping_dirname
	opts = OptionParser.new
  
	options = Hash.new

	opts.on("-w", "--workflow-path WORKFLOW_PATH", String, "path to meteor configuration file, e.g. workflow.ini") { |v| options["Workflow"] = v }
	opts.on("-i", "--input-path INPUT_DIR_PATH", String, "path to sample directory, containing the sample sequencing metadata (files ending with _census_stage_0.ini)") { |v| options["InputPath"] = v }
	opts.on("-p", "--project-path PROJECT_DIR_PATH", String, "path to project directory, containing mapping and profile data (e.g. /projects/project_dir)") { |v| options["ProjectPath"] = v }
	opts.on("-o", "--mapping-dir-basename MAPPING_DIR_BASENAME", String, "basename of the mapping directory") { |v| options["MappingBasename"] = v }
	
	opts.on("-C", "--counting-type TYPE_LIST", Array,
			"counting type string",
			"\n",
			"Possible values: either \"all\" or a comma separated list of type among the following:",
			"(where direct/undirect means direct/undirect strand)",
			"\n",
			"total_reads, shared_reads, smart_shared_reads, unique_reads",
			"\n",
			"direct_total_reads, direct_shared_reads, direct_smart_shared_reads,",
			"direct_unique_reads",
			"\n",
			"total_reads_coverage, unique_reads_coverage,",
			"direct_total_reads_coverage, direct_unique_reads_coverage,",
			"undirect_total_reads_coverage, undirect_unique_reads_coverage",
			"\n",
			"total_reads_mean_coverage, shared_reads_mean_coverage,",
			"smart_shared_reads_mean_coverage, unique_reads_mean_coverage",
			"\n",
			"direct_total_reads_mean_coverage, direct_shared_reads_mean_coverage,",
			"direct_smart_shared_reads_mean_coverage,",
			"direct_unique_reads_mean_coverage",
			"\n",
			"undirect_total_reads_mean_coverage,",
			"undirect_shared_reads_mean_coverage,",
			"undirect_smart_shared_reads_mean_coverage,",
			"undirect_unique_reads_mean_coverage",
			"\n",
			"undirect_total_reads, undirect_shared_reads,",
			"undirect_smart_shared_reads, undirect_unique_reads",
			"\n"
			) { |v| options["CountingTypeList"] = v.to_set }
			
	opts.on("-t", "--tmp-path TMP_PATH", String, "path to the directory where temporary files (e.g. sam) are stored") { |v| options["TmpPath"] = v }
	opts.on("-T", "--clean-tmp-path", "remove the temporary files directory after the task counting is finished") { |v| options["CleanTmpPath"] = 1 }
	opts.on("-L", "--remove-lock", "ignore and remove any lock file") { |v| options["NoLock"] = 1 }
	opts.on("-f", "--force", "force overwriting, ignoring previous results and lock files") { |v| options["Force"] = 1 }
	opts.on("-m", "--mapping-only", "execute mapping only") { |v| options["OnlyMapping"] = 1 }
	opts.on("-c", "--counting-only", "execute counting only") { |v| options["OnlyCounting"] = 1 }

	opts.banner = "Usage: #{File.basename($0)} [options]"

	# Display global_usage if no arguments given
	if ARGV.empty?
		puts opts
		exit
	end

	begin
		opts.parse!(ARGV)
	rescue OptionParser::ParseError => e
		STDERR.puts "Error, #{e}"
		STDERR.puts
		STDERR.puts opts
		exit 1
	end

	# Check command options (raise exception + message)
	begin
		CheckArgs(options)
	rescue Exception => e
		STDERR.puts("Error while checking options.")
		#STDERR.puts(e.backtrace.inspect)
		abort(e.message)
	end
	#puts options["CountingTypeList"].to_a

	# initialize all members to nil
	aMeteorSession = MeteorSession.new
	aMeteorSession.ProcessJob(options)

end #### END MAIN ####
