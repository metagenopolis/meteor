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

#ifndef TMETEORJOBINIFILE_H_
#define TMETEORJOBINIFILE_H_

#include <string>
#include <vector>
#include "MeteorConstant.h"
#include "IniFileManager.h"


//using namespace std;

//const std::string C_DEFAULT_REPOSITORY_DRIVE = "Y:";

const std::string C_JOB_WORKSESSION_SECTION = "worksession";
const std::string C_JOB_WORK_DIR = "meteor.workdir";
const std::string C_JOB_REFERENCE_DIR = "meteor.reference.dir";
const std::string C_JOB_REPOSITORY_ROOT_PATH = "meteor.repository.root.path";
const std::string C_JOB_REPOSITORY_DRIVE = "meteor.repository.drive";
const std::string C_JOB_IS_LOCAL_COPY = "meteor.localcopy";
const std::string C_JOB_METEOR_DB_TYPE = "meteor.db.type";
const std::string C_JOB_METEOR_ATTEMPS_COUNT = "meteor.nb.attempts";
const std::string C_JOB_INDEXED_READ_LENGTH = "meteor.indexed.readlength";
const std::string C_JOB_MAPPING_PROGRAM = "meteor.mapping.program";
const std::string C_JOB_MAPPING_DIR = "meteor.mapping.dir";
const std::string C_JOB_PROFILE_DIR = "meteor.profile.dir";
const std::string C_JOB_MAPPING_FILE_FORMAT = "meteor.mapping.file.format";
const std::string C_JOB_IS_DATA_CACHING = "meteor.datacaching";
const std::string C_JOB_IS_CPU_PERCENTAGE = "meteor.is.cpu.percentage";
const std::string C_JOB_CPU_COUNT = "meteor.cpu.count";
const std::string C_JOB_EXCLUDED_REFERENCE_COUNT = "meteor.excluded.reference.count";
const std::string C_JOB_CHECKSUM = "meteor.checksum";

const std::string C_JOB_READ_CLEANING_PROGRAM = "meteor.read.cleaning.program";
const std::string C_JOB_READ_CLEANING_CMD = "meteor.read.cleaning.cmd";

const std::string C_JOB_MAIN_REFERENCE_SECTION = "main_reference";
const std::string C_JOB_EXCLUDED_REFERENCE_PREFIX_SECTION = "excluded_reference_";
const std::string C_JOB_REFERENCE_NAME = "meteor.reference.name";
const std::string C_JOB_MAPPING_MATCHES = "meteor.matches";
const std::string C_JOB_MAPPING_MISMATCHES = "meteor.mismatches";
const std::string C_JOB_IS_MAPPING_MISMATCHES_PERCENTAGE = "meteor.is.perc.mismatches";
const std::string C_JOB_MAPPED_READ_LENGTH = "meteor.mapped.readlength";
const std::string C_JOB_MAPPED_READ_LENGTH_TYPE = "meteor.mapped.readlength.type";
const std::string C_JOB_LOCAL_QUALITY_CUTOFF = "meteor.local.quality.cutoff";
const std::string C_JOB_LOCAL_QUALITY_COUNT_CUTOFF = "meteor.local.quality.count.cutoff";
const std::string C_JOB_MEAN_QUALITY_CUTOFF = "meteor.mean.quality.cutoff";
const std::string C_JOB_IS_ON_BOTH_STRANDS = "meteor.is.on.both.strands";
const std::string C_JOB_IS_BESTALIGNMENT = "meteor.bestalignment";
const std::string C_JOB_MAPPING_PREFIX_NAME = "meteor.mapping.prefix.name";
const std::string C_JOB_COUNTING_PREFIX_NAME = "meteor.counting.prefix.name";
const std::string C_JOB_MAPPER_CMD = "meteor.mapper.cmd";

const std::string C_JOB_IS_LOCAL_MAPPING = "meteor.is.local.mapping";
const std::string C_JOB_ALIGNMENT_LENGTH_CUTOFF = "meteor.alignment.length.cutoff";
const std::string C_JOB_SOFT_CLIPPING_LENGTH_CUTOFF = "meteor.soft.clipping.length.cutoff";
const std::string C_JOB_KEEP_INTERNAL_LOCAL_ALIGNMENT = "meteor.keep.internal.local.alignment";

const std::string C_JOB_LIBRARY_EXPORT_OPTIONS_SECTION = "library_export_options";
const std::string C_JOB_MAPPING_EXPORT_OPTIONS_SECTION = "mapping_export_options";
const std::string C_JOB_EXPORT_FORMAT = "format";
const std::string C_JOB_EXPORT_OUTPUT_DIR = "output.dir";

const std::string C_JOB_EMAIL_INFO_SECTION = "email_info";
const std::string C_JOB_EMAIL_LOGIN = "login";
const std::string C_JOB_EMAIL_PASSWORD = "password";
const std::string C_JOB_EMAIL_HOST = "host";
const std::string C_JOB_EMAIL_PORT = "port";
const std::string C_JOB_EMAIL_SENDER = "sender";

typedef struct {
    bool IsExcludedReference;
    std::string ReferenceName;
    int MappingMatches;
    int MappingMismatches;
    bool IsMappingMismatchesPercentage;
    int MappedReadLength;
    TMappedReadLengthType MappedReadLengthType;
    int LocalQualityCutoff;
    int LocalQualityCountCutoff;
    int MeanQualityCutoff;
    bool IsBestAlignment;
    bool IsOnBothStrand;
    bool IsLocalAlignment; //// LOCAL
    bool KeepInternalLocalAlignment; //// LOCAL
    int AlignmentLengthCutoff; //// LOCAL
    int SoftClippingLengthCutoff; //// LOCAL
    std::string MappingPrefixName;
    std::string CountingPrefixName;
    std::string MapperCmd;
} TMappingReferenceProperties;

typedef std::vector<TMappingReferenceProperties> TMappingReferencePropertiesArray;

class TMeteorJobIniFile : public IniFileManager {
public:
	TMeteorJobIniFile();
	TMeteorJobIniFile(const std::string & filename);
	virtual ~TMeteorJobIniFile();

	int GetExcludedReferenceCount() const;
	void SetExcludedReferenceCount(int Value);

	void GetReferenceProperties(const std::string & aReferenceSection, TMappingReferenceProperties & aMappingReferenceProperties) const;
	void GetExcludedReferenceProperties(const int aExcludedReferenceIndex, TMappingReferenceProperties & aMappingReferenceProperties) const;
	void GetAllExcludedReferenceProperties(TMappingReferencePropertiesArray & aMappingReferencePropertiesArray) const;

	void GetMainReferenceProperties(TMappingReferenceProperties & aMappingReferenceProperties) const;

	std::string GetReferenceDir() const;
	void GetExcludedReferenceNameList(TStrings & aExcludedReferenceNameList) const;
	bool GetIsDataCaching() const;
	TMeteorDBType GetMeteorDBType() const;
	std::string GetMappingProgram() const;


};

#endif /* TMETEORJOBINIFILE_H_ */
