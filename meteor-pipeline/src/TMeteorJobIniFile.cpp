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

#include "TMeteorJobIniFile.h"

using namespace std;

TMeteorJobIniFile::TMeteorJobIniFile():IniFileManager() {;}

TMeteorJobIniFile::TMeteorJobIniFile(const string & filename):IniFileManager(filename) {
	if (getSectionsNumber() == 0){
		cerr << "Error, no section found in " << filename << endl;
		exit(1);
	}
}

TMeteorJobIniFile::~TMeteorJobIniFile() {;}

// ------------------------------------------------------------------------------

int TMeteorJobIniFile::GetExcludedReferenceCount() const{
	return ReadInteger(C_JOB_WORKSESSION_SECTION, C_JOB_EXCLUDED_REFERENCE_COUNT, 0);
}

// ------------------------------------------------------------------------------

void TMeteorJobIniFile::SetExcludedReferenceCount(int Value){
	WriteInteger(C_JOB_WORKSESSION_SECTION, C_JOB_EXCLUDED_REFERENCE_COUNT, Value);
}

// ------------------------------------------------------------------------------

void TMeteorJobIniFile::GetReferenceProperties(const string & aReferenceSection, TMappingReferenceProperties & aMappingReferenceProperties) const
{
  aMappingReferenceProperties.ReferenceName = ReadString(aReferenceSection, C_JOB_REFERENCE_NAME, "");
  aMappingReferenceProperties.MappingMatches = ReadInteger(aReferenceSection, C_JOB_MAPPING_MATCHES, C_DEFAULT_MAPPING_MATCHES_COUNT);
  aMappingReferenceProperties.MappingMismatches = ReadInteger(aReferenceSection, C_JOB_MAPPING_MISMATCHES, C_DEFAULT_MAPPING_MISMATCHES_COUNT);
  aMappingReferenceProperties.IsMappingMismatchesPercentage = (ReadInteger(aReferenceSection, C_JOB_IS_MAPPING_MISMATCHES_PERCENTAGE, 0) == 1);
  aMappingReferenceProperties.MappedReadLength = ReadInteger(aReferenceSection, C_JOB_MAPPED_READ_LENGTH, C_DEFAULT_MAPPING_READ_LENGTH);
  aMappingReferenceProperties.LocalQualityCutoff = ReadInteger(aReferenceSection, C_JOB_LOCAL_QUALITY_CUTOFF, C_DEFAULT_LOCAL_QUALITY_CUTOFF);
  aMappingReferenceProperties.LocalQualityCountCutoff = ReadInteger(aReferenceSection, C_JOB_LOCAL_QUALITY_COUNT_CUTOFF, C_DEFAULT_LOCAL_QUALITY_COUNT_CUTOFF);
  aMappingReferenceProperties.MeanQualityCutoff = ReadInteger(aReferenceSection, C_JOB_MEAN_QUALITY_CUTOFF, C_DEFAULT_MEAN_QUALITY_CUTOFF);
  aMappingReferenceProperties.IsBestAlignment = (ReadInteger(aReferenceSection, C_JOB_IS_BESTALIGNMENT, 1) == 1);
  aMappingReferenceProperties.IsOnBothStrand = (ReadInteger(aReferenceSection, C_JOB_IS_ON_BOTH_STRANDS, 1) == 1);

  aMappingReferenceProperties.IsLocalAlignment = (ReadInteger(aReferenceSection, C_JOB_IS_LOCAL_MAPPING, 0) == 1); //// LOCAL
  aMappingReferenceProperties.KeepInternalLocalAlignment = (ReadInteger(aReferenceSection, C_JOB_KEEP_INTERNAL_LOCAL_ALIGNMENT, 0) == 1);
  aMappingReferenceProperties.AlignmentLengthCutoff = ReadInteger(aReferenceSection, C_JOB_ALIGNMENT_LENGTH_CUTOFF, C_DEFAULT_ALIGNMENT_LENGTH_CUTOFF);
  aMappingReferenceProperties.SoftClippingLengthCutoff = ReadInteger(aReferenceSection, C_JOB_SOFT_CLIPPING_LENGTH_CUTOFF, C_DEFAULT_SOFT_CLIPPING_LENGTH_CUTOFF);

  string aMappedReadLengthTypeStr = ReadString(aReferenceSection, C_JOB_MAPPED_READ_LENGTH_TYPE, C_OVERALL_MAPPED_READ_LENGTH_TYPE);

  if (aMappedReadLengthTypeStr == C_FIXED_MAPPED_READ_LENGTH_TYPE)
    aMappingReferenceProperties.MappedReadLengthType = mrltFixed;
  else if (aMappedReadLengthTypeStr == C_MAX_MAPPED_READ_LENGTH_TYPE)
    aMappingReferenceProperties.MappedReadLengthType = mrltMax;
  else if (aMappedReadLengthTypeStr == C_MIN_MAPPED_READ_LENGTH_TYPE)
    aMappingReferenceProperties.MappedReadLengthType = mrltMin;
  else if (aMappedReadLengthTypeStr == C_OVERALL_MAPPED_READ_LENGTH_TYPE)
    aMappingReferenceProperties.MappedReadLengthType = mrltOverall;
  else
    aMappingReferenceProperties.MappedReadLengthType = mrltOverall;

  //// trim the strings ????
  aMappingReferenceProperties.MappingPrefixName = ReadString(aReferenceSection, C_JOB_MAPPING_PREFIX_NAME, "");
  aMappingReferenceProperties.CountingPrefixName = ReadString(aReferenceSection, C_JOB_COUNTING_PREFIX_NAME, "");
  aMappingReferenceProperties.MapperCmd = ReadString(aReferenceSection, C_JOB_MAPPER_CMD, "");
}

// ------------------------------------------------------------------------------

void TMeteorJobIniFile::GetAllExcludedReferenceProperties(TMappingReferencePropertiesArray & aMappingReferencePropertiesArray) const
{
  aMappingReferencePropertiesArray.resize(GetExcludedReferenceCount());

  for ( int aExcludedReferenceIndex = 0; aExcludedReferenceIndex < GetExcludedReferenceCount(); aExcludedReferenceIndex++)
  {
    GetReferenceProperties(C_JOB_EXCLUDED_REFERENCE_PREFIX_SECTION + numberToString(aExcludedReferenceIndex + 1), aMappingReferencePropertiesArray.at(aExcludedReferenceIndex)); //// AT
    aMappingReferencePropertiesArray.at(aExcludedReferenceIndex).IsExcludedReference = true; //// AT
  }
}

// ------------------------------------------------------------------------------

void TMeteorJobIniFile::GetMainReferenceProperties(TMappingReferenceProperties & aMappingReferenceProperties) const
{
  GetReferenceProperties(C_JOB_MAIN_REFERENCE_SECTION, aMappingReferenceProperties);
  aMappingReferenceProperties.IsExcludedReference = false;
}

// ------------------------------------------------------------------------------

void TMeteorJobIniFile::GetExcludedReferenceProperties(const int aExcludedReferenceIndex, TMappingReferenceProperties & aMappingReferenceProperties) const
{
  GetReferenceProperties(C_JOB_EXCLUDED_REFERENCE_PREFIX_SECTION + numberToString(aExcludedReferenceIndex), aMappingReferenceProperties);
  aMappingReferenceProperties.IsExcludedReference = true;
}

// ------------------------------------------------------------------------------

string TMeteorJobIniFile::GetReferenceDir() const
{
	return ReadString(C_JOB_WORKSESSION_SECTION, C_JOB_REFERENCE_DIR, "");
}

// ------------------------------------------------------------------------------

void TMeteorJobIniFile::GetExcludedReferenceNameList(TStrings & aExcludedReferenceNameList) const
{
  aExcludedReferenceNameList.clear();
  for (int i = 1; i <= GetExcludedReferenceCount(); i++)
  {
    aExcludedReferenceNameList.push_back(ReadString(C_JOB_EXCLUDED_REFERENCE_PREFIX_SECTION + numberToString(i), C_JOB_REFERENCE_NAME, ""));
  }
}

// ------------------------------------------------------------------------------

bool TMeteorJobIniFile::GetIsDataCaching() const
{
  return ReadInteger(C_JOB_WORKSESSION_SECTION, C_JOB_IS_DATA_CACHING, 0) == 1;
}

// ------------------------------------------------------------------------------

TMeteorDBType TMeteorJobIniFile::GetMeteorDBType() const

{
  string aMeteorDBTypeStr;
  aMeteorDBTypeStr = ReadString(C_JOB_WORKSESSION_SECTION, C_JOB_METEOR_DB_TYPE, C_METEOR_DBTYPE_ADT_STR);

  if (aMeteorDBTypeStr == C_METEOR_DBTYPE_ADT_STR) return mdbADT;
  if (aMeteorDBTypeStr == C_METEOR_DBTYPE_PARSTREAM_STR) return mdbParstream;
  if (aMeteorDBTypeStr == C_METEOR_DBTYPE_GIGABASE_STR) return mdbGigabase;
  if (aMeteorDBTypeStr == C_METEOR_DBTYPE_BINARY_STR) return mdbBinary;
  return mdbUnknown;
}

// ------------------------------------------------------------------------------

string TMeteorJobIniFile::GetMappingProgram() const {
	return ReadString(C_JOB_WORKSESSION_SECTION, C_JOB_MAPPING_PROGRAM, "");
}
