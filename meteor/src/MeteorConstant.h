/*
 * Copyright 2017-2020, Franck Gauthier <franck.gauthier@inrae.fr>, Nicolas Pons <nicolas.pons@inrae.fr>
 *
 * This file is part of Meteor v3.2.
 *
 * Meteor v3.2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Meteor v3.2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Meteor v3.2. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef METEORCONSTANT_H_
#define METEORCONSTANT_H_

#include <string>
#include <limits.h>
#include "utils.h"

//using namespace std;

//// DEBUG
const int C_READ_ID = -99999;

//// unused
//const std::string C_SHA256_EXTENSION = ".sha256";
//const std::string C_MD5_EXTENSION = ".md5";

// for mapping directory parametrization
const char C_MAPPED_READ_LENGTH_CHAR = 'l';
const char C_MISMATCHES_COUNT_CHAR = 'm';

// for extended pathname (>260 characters) NB: NTFS allows longer pathnames ?
#ifdef __WINDOWS__
const std::string C_EXTENDED_PREFIX_PATHNAME = "\\?\\";
#else
const std::string C_EXTENDED_PREFIX_PATHNAME = "";
#endif

// Time format for elapsed time information
const std::string C_TIME_FORMAT = "hh:mm:ss,zzz";
//// UNUSED const std::string C_ELAPSED_TIME_STR = "Elapsed time = ";

// License file
//// UNUSED const std::string C_METEOR_IDTAG_FILENAME = "meteor.idtag";
//// UNUSED const std::string C_METEOR_LICENSE_FILENAME = "meteor.lic";

// "Official" Program name
// const std::string C_METEOR_PROGRAM_NAME = "Meteor";

// const std::string C_DEFAULT_PROJECT_SAMPLE_DIRECTORY = "sample";

const std::string C_DEFAULT_COUNTING_REPORT_FILENAME_EXTENSION = "_counting_report.csv";
const std::string C_DEFAULT_EXTENDED_COUNTING_REPORT_FILENAME_EXTENSION = "_extended_counting_report.csv";

const std::string C_METEOR_ALIENTRIMMER_PROGRAM_NAME = "alientrimmer";
const std::string C_METEOR_INTERNAL_CLEANNING_PROGRAM_NAME = "internal";

const std::string C_DEFAULT_SAMPLE_MAPPING_DIRECTORY = "mapping";
//// UNUSED const std::string C_DEFAULT_SAMPLE_PROFILES_DIRECTORY = "profiles";

const std::string C_DEFAULT_REFERENCE_INIFILE_NAME = "reference.ini";
const std::string C_DEFAULT_REFERENCE_INIFILE_NAME_EXTENSION = "_reference.ini";

//// UNUSED const std::string C_DEFAULT_LIBRARY_PROFILE_NAME = "lib";
const std::string C_DEFAULT_PROFILE_INIFILE_NAME = "gene_profile.ini";
const std::string C_DEFAULT_PROFILE_INIFILE_NAME_EXTENSION = "_gene_profile.ini";
const std::string C_DEFAULT_PROFILE_DIR = "gene_profile";

const std::string C_EXT_INIFILE = ".ini";

const int C_CENSUS_STATUS_SEQUENCED = 0;
const int C_CENSUS_STATUS_MAPPED = 1;

const std::string C_CENSUS_STATUS_SEQUENCED_STR = "census_stage_0";
const std::string C_CENSUS_STATUS_MAPPED_STR = "census_stage_1";

// Constants for counting result table (and reference table)
const std::string C_CENSUS_RESULT_TABLE_NAME = "census";

// const std::string FD_FRAGMENT_NAME = "fragment_name";
const std::string FD_FRAGMENT_ID = "id_fragment";

// const std::string C_REFERENCE_ANNOTATION_TABLE_NAME = "annotation";
const std::string C_REFERENCE_LITE_ANNOTATION_TABLE_NAME = "lite_annotation";

const std::string C_DEFAULT_CONDITION_NAME = "NA";

// misc constants
// const int C_MAX_REFERENCE_LENGTH = 20000000;
// const int C_DEFAULT_READ_LENGTH = 50;
// const int C_DEFAULT_MIN_LOCAL_QUALITY_VALUE = -1;
// const int C_DEFAULT_MAX_LOCAL_QUALITY_VALUE = 100;

// const std::string C_PREFIX_IDX = "IDX";

/// /////////////////////////////////////////////////////////////////////////////

// --------------- constant for counted read type
const std::string C_TYPE_ALL_READ_COUNT = "all";
const std::string C_TYPE_READ_COUNT = "total_reads";
const std::string C_TYPE_SHARED_READ_COUNT = "shared_reads";
const std::string C_TYPE_SMART_SHARED_READ_COUNT = "smart_shared_reads";
const std::string C_TYPE_UNIQUE_READ_COUNT = "unique_reads";
const std::string C_TYPE_DIRECT_READ_COUNT = "direct_total_reads";
const std::string C_TYPE_DIRECT_SHARED_READ_COUNT = "direct_shared_reads";
const std::string C_TYPE_DIRECT_SMART_SHARED_READ_COUNT = "direct_smart_shared_reads";
const std::string C_TYPE_DIRECT_UNIQUE_READ_COUNT = "direct_unique_reads";

const std::string C_TYPE_READ_COVERAGE = "total_reads_coverage";
const std::string C_TYPE_UNIQUE_READ_COVERAGE = "unique_reads_coverage";
const std::string C_TYPE_DIRECT_READ_COVERAGE = "direct_total_reads_coverage";
const std::string C_TYPE_DIRECT_UNIQUE_READ_COVERAGE = "direct_unique_reads_coverage";
const std::string C_TYPE_UNDIRECT_READ_COVERAGE = "undirect_total_reads_coverage";
const std::string C_TYPE_UNDIRECT_UNIQUE_READ_COVERAGE = "undirect_unique_reads_coverage";

const std::string C_TYPE_READ_MEAN_COVERAGE = "total_reads_mean_coverage";
const std::string C_TYPE_SHARED_READ_MEAN_COVERAGE = "shared_reads_mean_coverage";
const std::string C_TYPE_SMART_SHARED_READ_MEAN_COVERAGE = "smart_shared_reads_mean_coverage";
const std::string C_TYPE_UNIQUE_READ_MEAN_COVERAGE = "unique_reads_mean_coverage";
const std::string C_TYPE_DIRECT_READ_MEAN_COVERAGE = "direct_total_reads_mean_coverage";
const std::string C_TYPE_DIRECT_SHARED_READ_MEAN_COVERAGE = "direct_shared_reads_mean_coverage";
const std::string C_TYPE_DIRECT_SMART_SHARED_READ_MEAN_COVERAGE = "direct_smart_shared_reads_mean_coverage";
const std::string C_TYPE_DIRECT_UNIQUE_READ_MEAN_COVERAGE = "direct_unique_reads_mean_coverage";
const std::string C_TYPE_UNDIRECT_READ_MEAN_COVERAGE = "undirect_total_reads_mean_coverage";
const std::string C_TYPE_UNDIRECT_SHARED_READ_MEAN_COVERAGE = "undirect_shared_reads_mean_coverage";
const std::string C_TYPE_UNDIRECT_SMART_SHARED_READ_MEAN_COVERAGE = "undirect_smart_shared_reads_mean_coverage";
const std::string C_TYPE_UNDIRECT_UNIQUE_READ_MEAN_COVERAGE = "undirect_unique_reads_mean_coverage";

const std::string C_TYPE_UNDIRECT_READ_COUNT = "undirect_total_reads";
const std::string C_TYPE_UNDIRECT_SHARED_READ_COUNT = "undirect_shared_reads";
const std::string C_TYPE_UNDIRECT_SMART_SHARED_READ_COUNT = "undirect_smart_shared_reads";
const std::string C_TYPE_UNDIRECT_UNIQUE_READ_COUNT = "undirect_unique_reads";

const std::string C_TYPE_UNKNOWN_READ_COUNT = "unknown";
// ------------------------------------------------------

const std::string C_NORMALIZATION_NO = "NO";

const std::string C_NORMALIZATION_NONE = "none";
const std::string C_NORMALIZATION_LENGTH = "frequency";
const std::string C_NORMALIZATION_FREQUENCY = "length";
const std::string C_NORMALIZATION_LENGTH_AND_FREQUENCY = "lengthfreq";
const std::string C_NORMALIZATION_RPKM = "rpkm";
const std::string C_NORMALIZATION_UNKNOWN = "unknown";

/// /////////////////////////////////////////////////////////////////////////////

// constant for device Kind
const std::string C_SOLID_DEVICE_STR = "SOLiD";
const std::string C_SOLEXA_DEVICE_STR = "Solexa";
const std::string C_ILLUMINA_DEVICE_STR = "Illumina";
const std::string C_PROTON_DEVICE_STR = "Proton";

/// /////////////////////////////////////////////////////////////////////////////
// Parameters for mapper output format

const std::string C_MAPPING_CSFASTA_OUTPUT_STR = "csfasta";
const std::string C_MAPPING_DEFAULT_BOWTIE_OUTPUT_STR = "default_bowtie";
const std::string C_MAPPING_LITE_DEFAULT_BOWTIE_OUTPUT_STR = "lite_default_bowtie";
const std::string C_MAPPING_SAM_OUTPUT_STR = "sam";
const std::string C_MAPPING_BAM_OUTPUT_STR = "bam";
const std::string C_MAPPING_UNKNOWN_OUTPUT_STR = "unknown";

const int C_DEFAULT_MAPPING_MATCHES_COUNT = 10000;
const int C_DEFAULT_MAPPING_MISMATCHES_COUNT = 3;
const int C_DEFAULT_MAPPING_READ_LENGTH = 35;

const int C_DEFAULT_LOCAL_QUALITY_CUTOFF = 9;
const int C_DEFAULT_LOCAL_QUALITY_COUNT_CUTOFF = 1;
const int C_DEFAULT_MEAN_QUALITY_CUTOFF = 20;

/// /////////////////////////////////////////////////////////////////////////////
// Constants for meteor database type

const std::string C_METEOR_DBTYPE_ADT_STR = "advantagedb";
const std::string C_METEOR_DBTYPE_PARSTREAM_STR = "parstream";
const std::string C_METEOR_DBTYPE_GIGABASE_STR = "gigabase";
const std::string C_METEOR_DBTYPE_BINARY_STR = "binary";
const std::string C_METEOR_DBTYPE_UNKNOWN_STR = "unknown";

/// /////////////////////////////////////////////////////////////////////////////
/// Constats for mapping position offset and orientation
/// see bowtie documentation: http://bowtie-bio.sourceforge.net/manual.shtml
// const std::string C_0_BASED_MAPPING_OFFSET_TYPE = "0-based";
// const std::string C_1_BASED_MAPPING_OFFSET_TYPE = "1-based";

/// /////////////////////////////////////////////////////////////////////////////
// Constants for mapped read length type
const std::string C_FIXED_MAPPED_READ_LENGTH_TYPE = "fixed";
const std::string C_OVERALL_MAPPED_READ_LENGTH_TYPE = "overall";
const std::string C_MAX_MAPPED_READ_LENGTH_TYPE = "max";
const std::string C_MIN_MAPPED_READ_LENGTH_TYPE = "min";

/// /////////////////////////////////////////////////////////////////////////////
// Constants for census and profile ini files

// census file
const std::string C_DATE_FORMAT = "yyyy-mm-dd";
const std::string C_SAMPLE_INFO_SECTION = "sample_info";
const std::string C_SAMPLE_NAME_STR = "sample_name";
// const std::string C_OLD_SAMPLE_NAME_STR = "old_sample_name";
// const std::string C_USUAL_SAMPLE_NAME_STR = "usual_sample_name";
const std::string C_CONDITION_NAME_STR = "condition_name";
const std::string C_PROJECT_NAME_STR = "project_name";
const std::string C_SEQUENCING_DATE_STR = "sequencing_date";
const std::string C_SEQUENCING_DEVICE_STR = "sequencing_device";
const std::string C_LIBRARY_TAG_STR = "tag";
// const std::string C_READ_LENGTH_STR = "read_length";
// const std::string C_READ_LENGTH_TYPE_STR = "read_length_type";
// const std::string C_INDEXED_READ_LENGTH_STR = "indexed_read_length";
// const std::string C_MAX_QUALITY_VALUE_STR = "max_quality_value";
// const std::string C_LOW_QUALITY_THRESHOLD_STR = "low_quality_threshold";
// const std::string C_QUALITY_ENCODING_STR = "quality_encoding";
// const std::string C_CENSUS_STATUS_STR = "census_status";
const std::string C_SAMPLE_FULL_NAME_STR = "full_sample_name";
const std::string C_DATABASE_TYPE_STR = "database_type";
const std::string C_SAMPLE_INDEXING_SOFTWARE_STR = "sample_indexing_software";
const std::string C_SAMPLE_INDEXING_SOFTWARE_VERSION_STR = "sample_indexing_software_version";
const std::string C_SAMPLE_INDEXING_SOFTWARE_DIGEST_STR = "sample_indexing_software";
const std::string C_SAMPLE_READ_CLEANING_METHOD_STR = "read_cleaning_method";
const std::string C_SAMPLE_READ_CLEANING_PARAMETERS_STR = "read_cleaning_parameters";
const std::string C_SAMPLE_RAW_TOTAL_READ_COUNT_STR = "sequenced_raw_read_count";
const std::string C_SAMPLE_TOTAL_READ_COUNT_STR = "sequenced_read_count";
const std::string C_SAMPLE_HOST_READ_COUNT_STR = "host_read_count";
const std::string C_SAMPLE_INDEXED_READ_COUNT_STR = "indexed_sequenced_read_count";
const std::string C_SAMPLE_INDEXED_BASE_COUNT_STR = "indexed_sequenced_base_count";
const std::string C_SAMPLE_FILE_SECTION = "sample_file";
const std::string C_SAMPLE_FILE_IS_COMPRESSED_STR = "is_compressed";
const std::string C_CSFASTA_FILENAME_STR = "csfasta_file";
// const std::string C_OLD_CSFASTA_FILENAME_STR = "old_csfasta_file";
const std::string C_QUAL_FILENAME_STR = "qual_file";
// const std::string C_OLD_QUAL_FILENAME_STR = "old_qual_file";
const std::string C_FASTQ_FILENAME_STR = "fastq_file";
// const std::string C_MD5_STR = "md5";

const std::string C_MAPPING_SECTION = "mapping";
// const std::string C_MAPPING_DATE_STR = "mapping_date";
const std::string C_REFERENCE_NAME_STR = "reference_name";
const std::string C_MAPPING_TOOL_STR = "mapping_tool";
const std::string C_MAPPING_TOOL_VERSION_STR = "mapping_tool_version";
const std::string C_MAPPING_TOOL_DIGEST_STR = "mapping_tool";
const std::string C_MAPPING_SOFTWARE_STR = "mapping_software";
const std::string C_MAPPING_SOFTWARE_VERSION_STR = "mapping_software_version";
const std::string C_MAPPING_SOFTWARE_DIGEST_STR = "mapping_software";
const std::string C_MAPPING_CMDLINE_STR = "mapping_cmdline";
const std::string C_MAPPED_READ_LENGTH_STR = "mapped_read_length";
const std::string C_MAPPED_READ_LENGTH_TYPE_STR = "mapped_read_length_type";
const std::string C_MAPPING_MISMATCHES_STR = "mismatches";
const std::string C_MAPPING_IS_MISMATCHES_PERCENTAGE_STR = "is_mismatches_percentage";
const std::string C_MAPPING_IS_LOCAL_MAPPING_STR = "is_local_mapping";
const std::string C_MAPPING_MATCHES_STR = "matches";
const std::string C_KEEP_ONLY_BEST_MATCHES_STR = "keep_only_best_matches";
const std::string C_MAPPING_FILE_FORMAT_STR = "mapping_file_format";
//const std::string C_READ_COUNT_STR = "read_count";
const std::string C_PROCESSED_READ_COUNT_STR = "processed_read_count";
//const std::string C_MAPPED_READ_COUNT_STR = "mapped_read_count";
//const std::string C_UNIQUE_MAPPED_READ_COUNT_STR = "unique_mapped_read_count";

const std::string C_MAPPING_FILE_SECTION = "mapping_file";
const std::string C_MAPPING_BOWTIE_FILENAME_STR = "bowtie_file";
const std::string C_MAPPING_FILE_COUNT_STR = "mapping_file_count";

const std::string C_COUNTING_SECTION = "counting";
//const std::string C_COUNTING_DATE_STR = "counting_date";

const std::string C_ALIGNMENT_LENGTH_CUTOFF_STR = "alignment_length_cutoff";
const std::string C_SOFT_CLIPPING_LENGTH_CUTOFF_STR = "soft_clipping_length_cutoff";
const std::string C_KEEP_INTERNAL_LOCAL_ALIGNMENT_STR = "keep_internal_local_alignment";
const std::string C_EXCLUDED_REFERENCE_COUNT_STR = "excluded_reference_count";
const std::string C_PREFIX_EXCLUDED_REFERENCE_NAME_STR = "excluded_reference_name_";
const std::string C_COUNTED_READ_COUNT_STR = "counted_read_count";
const std::string C_UNIQUE_COUNTED_READ_COUNT_STR = "unique_counted_read_count";
const std::string C_REJECTED_READ_COUNT_STR = "rejected_read_count";
const std::string C_EXCLUDED_REFERENCE_REJECTED_READ_COUNT = "excluded_reference_rejected_read_count";
const std::string C_QUALITY_REJECTED_READ_COUNT_STR = "quality_rejected_read_count";
const std::string C_NOT_COUNTED_READ_COUNT_STR = "not_counted_read_count";
const std::string C_NOT_COUNTED_CLEAN_READ_COUNT_STR = "not_counted_clean_read_count";

const std::string C_INTERSECT_COUNT_PREFIX_STR = "intersect_count_";

// profile ini file (profile constants common with census ones are commented)
////const std::string C_DATE_FORMAT = "yyyy-mm-dd";
const std::string C_PROFILE_INFO_SECTION = "profile_info";
const std::string C_METEOR_CONFIG_PATH_STR = "meteor_config_path";
const std::string C_PROFILE_NAME_STR = "profile_name";
////const std::string C_PROJECT_NAME_STR = "project_name";
////const std::string C_REFERENCE_NAME_STR = "reference_name";
const std::string C_COUNTED_READ_TYPE_STR = "counted_read_type";
const std::string C_NORMALIZATION_STR = "normalization";
const std::string C_PROFILE_DATE_STR = "profile_date";
const std::string C_PROFILE_FILENAME_STR = "profile_filename";
////const std::string C_DATABASE_TYPE_STR = "database_type";
////const std::string C_EXCLUDED_REFERENCE_COUNT_STR = "excluded_reference_count";
////const std::string C_PREFIX_EXCLUDED_REFERENCE_NAME_STR = "excluded_reference_name_";

////const std::string C_SAMPLE_READ_CLEANING_METHOD_STR = "read_cleaning_method";
////const std::string C_SAMPLE_READ_CLEANING_PARAMETERS_STR = "read_cleaning_parameters";
////const std::string C_SAMPLE_TOTAL_READ_COUNT_STR = "sequenced_read_count";
////const std::string C_SAMPLE_RAW_TOTAL_READ_COUNT_STR = "sequenced_raw_read_count";
////const std::string C_SAMPLE_INDEXED_READ_COUNT_STR = "indexed_sequenced_read_count";
////const std::string C_SAMPLE_INDEXED_BASE_COUNT_STR = "indexed_sequenced_base_count";

////const std::string C_COUNTED_READ_COUNT_STR = "counted_read_count";
////const std::string C_UNIQUE_COUNTED_READ_COUNT_STR = "unique_counted_read_count";
////const std::string C_REJECTED_READ_COUNT_STR = "rejected_read_count";
////const std::string C_EXCLUDED_REFERENCE_REJECTED_READ_COUNT = "excluded_reference_rejected_read_count";
////const std::string C_QUALITY_REJECTED_READ_COUNT_STR = "quality_rejected_read_count";
////const std::string C_NOT_COUNTED_READ_COUNT_STR = "not_counted_read_count";
////const std::string C_NOT_COUNTED_CLEAN_READ_COUNT_STR = "not_counted_clean_read_count";

////const std::string C_INTERSECT_COUNT_PREFIX_STR = "intersect_count_";

////const std::string C_MAPPING_SOFTWARE_STR = "mapping_software";
////const std::string C_MAPPED_READ_LENGTH_STR = "mapped_read_length";
////const std::string C_MAPPED_READ_LENGTH_TYPE_STR = "mapped_read_length_type";
////const std::string C_MAPPING_MISMATCHES_STR = "mismatches";
////const std::string C_MAPPING_IS_MISMATCHES_PERCENTAGE_STR = "is_mismatches_percentage";
////const std::string C_MAPPING_MATCHES_STR = "matches";
////const std::string C_KEEP_ONLY_BEST_MATCHES_STR = "keep_only_best_matches";
////const std::string C_MAPPING_IS_LOCAL_MAPPING_STR = "is_local_mapping"; //// LOCAL
////const std::string C_ALIGNMENT_LENGTH_CUTOFF_STR = "alignment_length_cutoff";
////const std::string C_KEEP_INTERNAL_LOCAL_ALIGNMENT_STR = "keep_internal_local_alignment";
////const std::string C_MAPPING_CMDLINE_STR = "mapping_cmdline";
const std::string C_PREFIX_PROFILED_SAMPLE_SECTION = "profiled_sample";
////const std::string C_SAMPLE_NAME_STR = "sample_name";
////const std::string C_SAMPLE_FULL_NAME_STR = "full_sample_name";
////const std::string C_USUAL_SAMPLE_NAME_STR = "usual_sample_name";
////const std::string C_CONDITION_NAME_STR = "condition_name";
////const std::string C_SEQUENCING_DATE_STR = "sequencing_date";

const std::string C_PROFILE_DIRECTORY = "directory";

/// /////////////////////////////////////////////////////////////////////////////
const int C_DEFAULT_ALIGNMENT_LENGTH_CUTOFF = 45;
const int C_DEFAULT_SOFT_CLIPPING_LENGTH_CUTOFF = INT_MAX;
/// /////////////////////////////////////////////////////////////////////////////

/// /////////////////////////////////////////////////////////////////////////////
/// Constants for digest algorithm //// unused
//const std::string C_SHA256_DIGEST_ALGORITHM = "sha256";
//const std::string C_MD5_DIGEST_ALGORITHM = "md5";
//const std::string C_NONE_DIGEST_ALGORITHM = "none";

enum TMappedReadLengthType {
	mrltFixed, mrltMin, mrltMax, mrltOverall
};

enum TMappingPositionOffSetType {
	mpot0based, mpot1based
};

enum TMappingPositionOrientationType {
	mporitLeftMost, mporitReadStart
};

enum TReadLengthType {
	rltFixed, rltVariable, rltUnknown
};

enum TReferenceEntryType {
	retContig, retFragment, retUnknown
};

enum TQualityEncoding {
	qeSOLiD,
	qeSolexa,
	qeIllumina13,
	qeIllumina15,
	qeIllumina18,
	qeProton,
	qeSanger,
	qeUnknown
};
enum TQualityEncodingParam {
	qepMinAscii, qepMaxAscii, qepMinQuality, qepMaxQuality, qepOffSet
};

enum TSequenceFileFormat {
	sffCSFastaAndQual, sffFastQ, sffXSQ, sffUnknown
};

enum TMappingFileFormat {
	mffCSFasta,
	mffSAM,
	mffBAM,
	mffDefaultBowtie,
	mffLiteDefaultBowtie,
	mffUnknown
};

/// counted read type enumeration.
// first part defines the column number in the input files (census.dat) first column is zero.
// second part are calculated from the others.
// ORDER IS VERY IMPORTANT as it corresponds to the order of the elements in C_CENSUS_TABLE_FILE_HEADER_ARRAY (see below).

enum TCountedReadType {
	crtTotalReads = 0,
	crtSharedReads,
	crtSmartSharedReads,
	crtUniqueReads,
	crtDirectTotalReads,
	crtDirectSharedReads,
	crtDirectSmartSharedReads,
	crtDirectUniqueReads,
	crtTotalReadsCoverage,
	crtUniqueReadsCoverage,
	crtDirectTotalReadsCoverage,
	crtDirectUniqueReadsCoverage,
	crtUnDirectTotalReadsCoverage,
	crtUnDirectUniqueReadsCoverage,
	crtTotalReadsMeanCoverage,
	crtSharedReadsMeanCoverage,
	crtSmartSharedReadsMeanCoverage,
	crtUniqueReadsMeanCoverage,
	crtDirectTotalReadsMeanCoverage,
	crtDirectSharedReadsMeanCoverage,
	crtDirectSmartSharedReadsMeanCoverage,
	crtDirectUniqueReadsMeanCoverage,
	crtUndirectTotalReadsMeanCoverage,
	crtUndirectSharedReadsMeanCoverage,
	crtUndirectSmartSharedReadsMeanCoverage,
	crtUndirectUniqueReadsMeanCoverage,

	crtUndirectTotalReads,
	crtUndirectSharedReads,
	crtUndirectSmartSharedReads,
	crtUndirectUniqueReads,
	
	crtUnknownReads,
	crtUnknownReadsMeanCoverage
};

static const char* const C_CENSUS_TABLE_FILE_HEADER_ARRAY[30] = {
		"TotalReadCount",
		"SharedReadCount",
		"SmartSharedReadCount",
		"UniqueReadCount",
		"DirectTotalReadCount",
		"DirectSharedReadCount",
		"DirectSmartSharedReadCount",
		"DirectUniqueReadCount",
		"TotalReadCoverage",
		"UniqueReadCoverage",
		"DirectTotalReadCoverage",
		"DirectUniqueReadCoverage",
		"UnDirectTotalReadCoverage",
		"UnDirectUniqueReadCoverage",
		"TotalReadMeanCoverage",
		"SharedReadMeanCoverage",
		"SmartSharedReadMeanCoverage",
		"UniqueReadMeanCoverage",
		"DirectTotalReadMeanCoverage",
		"DirectSharedReadMeanCoverage",
		"DirectSmartSharedReadMeanCoverage",
		"DirectUniqueReadMeanCoverage",
		"UnDirectTotalReadMeanCoveragev",
		"UnDirectSharedReadMeanCoverage",
		"UnDirectSmartSharedReadMeanCoverage",
		"UnDirectUniqueReadMeanCoverage",
		"UnDirectTotalReadCount",
		"UnDirectSharedReadCount",
		"UnDirectSmartSharedReadCount",
		"UnDirectUniqueReadCount"};

enum TMeteorDBType {
	mdbADT, mdbParstream, mdbGigabase, mdbBinary, mdbUnknown
};

enum TMeteorProfileNormalizationType {
	mpnNone, mpnFrequency, mpnLength, mpnLengthAndFrequency, mpnRPKM, mpnUnknown
};

typedef struct {
	int RawSequencedReads;
	int HostSequencedReads;
	int SequencedReads;
	int IndexedReads;
	// total reads count after cleaning (on rejected reference and quality)
	int TotalCleanCountedReads;
	// unique reads count after cleaning (on rejected reference and quality)
	int UniqueCleanCountedReads;
	// all rejected reads (on rejected reference and quality)
	int RejectedReads;
	// rejected reads only on rejected reference
	int RejectedReferenceReads;
	// rejected reads only on quality
	int RejectedQualityReads;
	// all other reads not mapped on all reference
	int NotMappedReads;
	// all other reads not mapped on all reference after cleaning on quality
	int NotMappedCleanReads;
} TCountingStatistics;

enum TReadCleaningMethodType {
	rcmNone, rcmAlienTrimmer, rcmInternal
};

//enum TDigestAlgorithm {daMD5, daSHA256, daNone}; //// unused

#endif /* METEORCONSTANT_H_ */
