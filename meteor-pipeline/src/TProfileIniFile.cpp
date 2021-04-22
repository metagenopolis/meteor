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

#include "TProfileIniFile.h"

using namespace std;

TProfileIniFile::TProfileIniFile(const string & filename, bool aIsReadOnly):IniFileManager(){
	FReadOnly = aIsReadOnly;

	// do not load existing data if open in write mode (overwrite)
	if (FReadOnly) loadIniData(filename);
	else iniFilename = filename;
}

// --------------------------------------------------------------------

TProfileIniFile::~TProfileIniFile() {
	// with no argument, printIniData() overwrite IniFileManager::iniFilename (inherited).
	if (!FReadOnly) printIniData();
}

// --------------------------------------------------------------------

//void TProfileIniFile::SetProfileDate(const TDateTime & Value){ ////
//	string aStrDateTime;
//
//	aStrDateTime = FormatDateTime(C_DATE_FORMAT, Value);
//	WriteString(C_PROFILE_INFO_SECTION, C_PROFILE_DATE_STR, aStrDateTime);
//}

void TProfileIniFile::SetProfileDate(const string & Value){
	WriteString(C_PROFILE_INFO_SECTION, C_PROFILE_DATE_STR, Value);
}

// --------------------------------------------------------------------

void TProfileIniFile::SetReferenceName(const string & Value){
	WriteString(C_PROFILE_INFO_SECTION, C_REFERENCE_NAME_STR, Value);
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetProfileFileName(const string & Value){
	WriteString(C_PROFILE_INFO_SECTION, C_PROFILE_FILENAME_STR, Value);
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetProfileName(const string & Value){
	WriteString(C_PROFILE_INFO_SECTION, C_PROFILE_NAME_STR, Value);
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetMeteorConfigPath(const string & Value){
	WriteString(C_PROFILE_INFO_SECTION, C_METEOR_CONFIG_PATH_STR, Value);
}

// --------------------------------------------------------------------

void TProfileIniFile::SetProjectName(const string & Value){
	WriteString(C_PROFILE_INFO_SECTION, C_PROJECT_NAME_STR, Value);
}

// --------------------------------------------------------------------

void TProfileIniFile::SetDatabaseType(const TMeteorDBType & Value){

	string aMeteorDBTypeStr;

	switch(Value) {
		case mdbADT:
			aMeteorDBTypeStr = C_METEOR_DBTYPE_ADT_STR; break;
		case mdbParstream:
			aMeteorDBTypeStr = C_METEOR_DBTYPE_PARSTREAM_STR; break;
		case mdbGigabase:
			aMeteorDBTypeStr = C_METEOR_DBTYPE_GIGABASE_STR; break;
		case mdbBinary:
			aMeteorDBTypeStr = C_METEOR_DBTYPE_BINARY_STR; break;
		case mdbUnknown:
			aMeteorDBTypeStr = C_METEOR_DBTYPE_UNKNOWN_STR; break;
	}
	WriteString(C_PROFILE_INFO_SECTION, C_DATABASE_TYPE_STR, aMeteorDBTypeStr);
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetCountedReadType(const TCountedReadType & Value)
{
	switch(Value) {
		case crtTotalReads:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_READ_COUNT); break;
		case crtSharedReads:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_SHARED_READ_COUNT); break;
		case crtUniqueReads:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_UNIQUE_READ_COUNT); break;
		case crtDirectTotalReads:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_DIRECT_READ_COUNT); break;
		case crtDirectSharedReads:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_DIRECT_SHARED_READ_COUNT); break;
		case crtDirectUniqueReads:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_DIRECT_UNIQUE_READ_COUNT); break;
		case crtUndirectTotalReads:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_UNDIRECT_READ_COUNT); break;
		case crtUndirectSharedReads:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_UNDIRECT_SHARED_READ_COUNT); break;
		case crtUndirectUniqueReads:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_UNDIRECT_UNIQUE_READ_COUNT); break;
		case crtSmartSharedReads:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_SMART_SHARED_READ_COUNT); break;
		case crtDirectSmartSharedReads:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_DIRECT_SMART_SHARED_READ_COUNT); break;
		case crtUndirectSmartSharedReads:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_UNDIRECT_SMART_SHARED_READ_COUNT); break;
		case crtTotalReadsMeanCoverage:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_READ_MEAN_COVERAGE); break;
		case crtDirectTotalReadsMeanCoverage:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_DIRECT_READ_MEAN_COVERAGE); break;
		case crtUndirectTotalReadsMeanCoverage:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_UNDIRECT_READ_MEAN_COVERAGE); break;
		case crtSharedReadsMeanCoverage:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_SHARED_READ_MEAN_COVERAGE); break;
		case crtDirectSharedReadsMeanCoverage:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_DIRECT_SHARED_READ_MEAN_COVERAGE); break;
		case crtUndirectSharedReadsMeanCoverage:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_UNDIRECT_SHARED_READ_MEAN_COVERAGE); break;
		case crtSmartSharedReadsMeanCoverage:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_SMART_SHARED_READ_MEAN_COVERAGE); break;
		case crtDirectSmartSharedReadsMeanCoverage:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_DIRECT_SMART_SHARED_READ_MEAN_COVERAGE); break;
		case crtUndirectSmartSharedReadsMeanCoverage:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_UNDIRECT_SMART_SHARED_READ_MEAN_COVERAGE); break;
		case crtUniqueReadsMeanCoverage:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_UNIQUE_READ_MEAN_COVERAGE); break;
		case crtDirectUniqueReadsMeanCoverage:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_DIRECT_UNIQUE_READ_MEAN_COVERAGE); break;
		case crtUndirectUniqueReadsMeanCoverage:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_UNDIRECT_UNIQUE_READ_MEAN_COVERAGE); break;
		case crtDirectTotalReadsCoverage:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_DIRECT_READ_COVERAGE); break;
		case crtDirectUniqueReadsCoverage:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_DIRECT_UNIQUE_READ_COVERAGE); break;
		case crtTotalReadsCoverage:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_READ_COVERAGE); break;
		case crtUnDirectTotalReadsCoverage:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_UNDIRECT_READ_COVERAGE); break;
		case crtUnDirectUniqueReadsCoverage:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_UNDIRECT_UNIQUE_READ_COVERAGE); break;
		case crtUniqueReadsCoverage:
		  WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_UNIQUE_READ_COVERAGE); break;
		default:
			WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, C_TYPE_UNKNOWN_READ_COUNT);
	}
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetCountedReadType(const string & Value){
	WriteString(C_PROFILE_INFO_SECTION, C_COUNTED_READ_TYPE_STR, Value);
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetNormalizationType(const TMeteorProfileNormalizationType & Value){
	switch(Value) {
		case mpnNone:
		  WriteString(C_PROFILE_INFO_SECTION, C_NORMALIZATION_STR, C_NORMALIZATION_NONE); break;
		case mpnFrequency:
		  WriteString(C_PROFILE_INFO_SECTION, C_NORMALIZATION_STR, C_NORMALIZATION_FREQUENCY); break;
		case mpnLength:
		  WriteString(C_PROFILE_INFO_SECTION, C_NORMALIZATION_STR, C_NORMALIZATION_LENGTH); break;
		case mpnLengthAndFrequency:
		  WriteString(C_PROFILE_INFO_SECTION, C_NORMALIZATION_STR, C_NORMALIZATION_LENGTH_AND_FREQUENCY); break;
		case mpnRPKM:
		  WriteString(C_PROFILE_INFO_SECTION, C_NORMALIZATION_STR, C_NORMALIZATION_RPKM); break;
		case mpnUnknown:
		  WriteString(C_PROFILE_INFO_SECTION, C_NORMALIZATION_STR, C_NORMALIZATION_UNKNOWN); break;
	}
}

// --------------------------------------------------------------------

void TProfileIniFile::SetConditionName(const string & Value){
	WriteString(C_PROFILE_INFO_SECTION, C_CONDITION_NAME_STR, Value);
}

// --------------------------------------------------------------------

void TProfileIniFile::SetReadCleaningMethod(const TReadCleaningMethodType & Value){
	switch(Value) {
		case rcmNone:
			WriteString(C_PROFILE_INFO_SECTION, C_SAMPLE_READ_CLEANING_METHOD_STR, "");
			break;
		case rcmAlienTrimmer:
			WriteString(C_PROFILE_INFO_SECTION, C_SAMPLE_READ_CLEANING_METHOD_STR, C_METEOR_ALIENTRIMMER_PROGRAM_NAME);
			break;
		case rcmInternal:
			WriteString(C_PROFILE_INFO_SECTION, C_SAMPLE_READ_CLEANING_METHOD_STR, C_METEOR_INTERNAL_CLEANNING_PROGRAM_NAME);
			break;
	}
}

// --------------------------------------------------------------------

void TProfileIniFile::SetReadCleaningParameters(const string & Value){
	WriteString(C_PROFILE_INFO_SECTION, C_SAMPLE_READ_CLEANING_PARAMETERS_STR, Value);
}

// --------------------------------------------------------------------

void TProfileIniFile::SetSampleName(const string & Value){
	WriteString(C_PROFILE_INFO_SECTION, C_SAMPLE_NAME_STR, Value);
}

// --------------------------------------------------------------------

void TProfileIniFile::SetProfileDirectory(const string & Value){
	WriteString(C_PROFILE_INFO_SECTION, C_PROFILE_DIRECTORY, Value);
}

// --------------------------------------------------------------------

int TProfileIniFile::AddProfiledSample(const string & aSampleName, const string & aSampleFullName, const string & aConditionName, const string & aSequencingDate)
{
  ////TStrings aSectionList;
  ////ReadSections(aSectionList);
  ////aProfiledSampleIndex = aSectionList.Count - 1;
  size_t aProfiledSampleIndex = getSectionsNumber() - 1;
  // -1 for eliminate profile_info section

  string aSection = C_PREFIX_PROFILED_SAMPLE_SECTION + "_" + numberToString(aProfiledSampleIndex);

  WriteString(aSection, C_SAMPLE_NAME_STR, aSampleName);
  WriteString(aSection, C_SAMPLE_FULL_NAME_STR, aSampleFullName);
  WriteString(aSection, C_CONDITION_NAME_STR, aConditionName);
  WriteString(aSection, C_SEQUENCING_DATE_STR, aSequencingDate);

  ////aSectionList.Free;
  ////aSectionList = nil;

  return aProfiledSampleIndex;
}

// ------------------------------------------------------------------------------

string TProfileIniFile::GetSampleName() const
{
  return ReadString(C_PROFILE_INFO_SECTION, C_SAMPLE_NAME_STR, "");
}

// ------------------------------------------------------------------------------

string TProfileIniFile::GetConditionName() const
{
	//// pourquoi faire un test empty ????
	  string Result = ReadString(C_PROFILE_INFO_SECTION, C_CONDITION_NAME_STR, C_DEFAULT_CONDITION_NAME);
	  if (Result.empty()) return C_DEFAULT_CONDITION_NAME;
	  return Result;
}

// ------------------------------------------------------------------------------

bool TProfileIniFile::GetIsBestAlignment() const
{
	return ReadBool(C_PROFILE_INFO_SECTION, C_KEEP_ONLY_BEST_MATCHES_STR, true);
}

// ------------------------------------------------------------------------------

string TProfileIniFile::GetProfiledSampleFullName(int Index) const
{
  return ReadString(C_PREFIX_PROFILED_SAMPLE_SECTION + '_' + numberToString(Index), C_SAMPLE_FULL_NAME_STR, "");
}

// ------------------------------------------------------------------------------

string TProfileIniFile::GetReferenceName() const {
	return ReadString(C_PROFILE_INFO_SECTION, C_REFERENCE_NAME_STR, "");
}
int TProfileIniFile::GetMappingMatchesCount() const {
	return ReadInteger(C_PROFILE_INFO_SECTION, C_MAPPING_MATCHES_STR, -1);
}
int TProfileIniFile::GetMappingMismatchesCount() const {
	return ReadInteger(C_PROFILE_INFO_SECTION, C_MAPPING_MISMATCHES_STR, -1);
}
string TProfileIniFile::GetMappingToolName() const
{
	return ReadString(C_PROFILE_INFO_SECTION, C_MAPPING_SOFTWARE_STR, "");
}
int TProfileIniFile::GetMappedReadLength() const {
	return ReadInteger(C_PROFILE_INFO_SECTION, C_MAPPED_READ_LENGTH_STR, -1);
}
bool TProfileIniFile::GetIsMappingMismatchesPercentage() const {
	return ReadBool(C_PROFILE_INFO_SECTION, C_MAPPING_IS_MISMATCHES_PERCENTAGE_STR, false);
}

// ------------------------------------------------------------------------------

string TProfileIniFile::GetProfiledSequencingDate(int Index) const
{
  return ReadString(C_PREFIX_PROFILED_SAMPLE_SECTION + "_" + numberToString(Index), C_SEQUENCING_DATE_STR, "0");
}

// ------------------------------------------------------------------------------

TMappedReadLengthType TProfileIniFile::GetMappedReadLengthType() const {
	string aMappedReadLengthTypeStr = ReadString(C_PROFILE_INFO_SECTION, C_MAPPED_READ_LENGTH_TYPE_STR, C_OVERALL_MAPPED_READ_LENGTH_TYPE);

	if (aMappedReadLengthTypeStr == C_FIXED_MAPPED_READ_LENGTH_TYPE) return mrltFixed;
	if (aMappedReadLengthTypeStr == C_MAX_MAPPED_READ_LENGTH_TYPE) return mrltMax;
	if (aMappedReadLengthTypeStr == C_MIN_MAPPED_READ_LENGTH_TYPE) return mrltMin;
	if (aMappedReadLengthTypeStr == C_OVERALL_MAPPED_READ_LENGTH_TYPE) return mrltOverall;
	return mrltOverall;
}

string TProfileIniFile::GetMappingCmdLine() const {
	return ReadString(C_PROFILE_INFO_SECTION, C_MAPPING_CMDLINE_STR, "");
}

// ------------------------------------------------------------------------------

void TProfileIniFile::GetExcludedReferenceNameList(TStrings & aExcludedReferenceNameList) const
{
	int excludedRefCount = ReadInteger(C_PROFILE_INFO_SECTION, C_EXCLUDED_REFERENCE_COUNT_STR, 0);
	aExcludedReferenceNameList.clear();
	aExcludedReferenceNameList.resize(excludedRefCount);

	for (int i = 0; i < excludedRefCount; i++) {
//		aExcludedReferenceNameList[i] = ReadString(C_PROFILE_INFO_SECTION, C_PREFIX_EXCLUDED_REFERENCE_NAME_STR + to_string(i+1), ""); //// c++11
		aExcludedReferenceNameList.at(i) = ReadString(C_PROFILE_INFO_SECTION, C_PREFIX_EXCLUDED_REFERENCE_NAME_STR + numberToString(i+1), "");
	}
}

// ------------------------------------------------------------------------------

int TProfileIniFile::GetExcludedReferenceReadCount(const int aExcludedReferenceIndex) const
{
  return ReadInteger(C_PROFILE_INFO_SECTION, C_EXCLUDED_REFERENCE_REJECTED_READ_COUNT + "_" + numberToString(aExcludedReferenceIndex), 0);
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetCountingStatistics(const int aProfiledSampleIndex, const TCountingStatistics & aCountingStatistics)
{
  string aSection = C_PREFIX_PROFILED_SAMPLE_SECTION + "_" + numberToString(aProfiledSampleIndex);

  WriteInteger(aSection, C_SAMPLE_RAW_TOTAL_READ_COUNT_STR, aCountingStatistics.RawSequencedReads);
  WriteInteger(aSection, C_SAMPLE_HOST_READ_COUNT_STR, aCountingStatistics.HostSequencedReads);
  WriteInteger(aSection, C_SAMPLE_TOTAL_READ_COUNT_STR, aCountingStatistics.SequencedReads);
  WriteInteger(aSection, C_SAMPLE_INDEXED_READ_COUNT_STR, aCountingStatistics.IndexedReads);
  WriteInteger(aSection, C_COUNTED_READ_COUNT_STR, aCountingStatistics.TotalCleanCountedReads);
  WriteInteger(aSection, C_UNIQUE_COUNTED_READ_COUNT_STR, aCountingStatistics.UniqueCleanCountedReads);
  WriteInteger(aSection, C_REJECTED_READ_COUNT_STR, aCountingStatistics.RejectedReads);
  WriteInteger(aSection, C_EXCLUDED_REFERENCE_REJECTED_READ_COUNT, aCountingStatistics.RejectedReferenceReads);
  WriteInteger(aSection, C_QUALITY_REJECTED_READ_COUNT_STR, aCountingStatistics.RejectedQualityReads);
  WriteInteger(aSection, C_NOT_COUNTED_READ_COUNT_STR, aCountingStatistics.NotMappedReads);
  WriteInteger(aSection, C_NOT_COUNTED_CLEAN_READ_COUNT_STR, aCountingStatistics.NotMappedCleanReads);
  }

// ------------------------------------------------------------------------------

int TProfileIniFile::GetExcludedReferenceReadCount(const int aProfiledSampleIndex, const int aExcludedReferenceIndex) const
{
  return ReadInteger(C_PREFIX_PROFILED_SAMPLE_SECTION + "_" + numberToString(aProfiledSampleIndex), C_EXCLUDED_REFERENCE_REJECTED_READ_COUNT + "_" + numberToString(aExcludedReferenceIndex), 0);
}

// ------------------------------------------------------------------------------

TCountingStatistics TProfileIniFile::GetCountingStatistics() const
{
	TCountingStatistics aCountingStatistics;

	aCountingStatistics.RawSequencedReads = GetSequencedRawReadCount();
	aCountingStatistics.HostSequencedReads = GetSequencedHostReadCount();
	aCountingStatistics.SequencedReads = ReadInteger(C_PROFILE_INFO_SECTION, C_SAMPLE_TOTAL_READ_COUNT_STR, 0);
	aCountingStatistics.IndexedReads = ReadInteger(C_PROFILE_INFO_SECTION, C_SAMPLE_INDEXED_READ_COUNT_STR, 0);
	aCountingStatistics.TotalCleanCountedReads = ReadInteger(C_PROFILE_INFO_SECTION, C_COUNTED_READ_COUNT_STR, 0);
	aCountingStatistics.UniqueCleanCountedReads = ReadInteger(C_PROFILE_INFO_SECTION, C_UNIQUE_COUNTED_READ_COUNT_STR, 0);
	aCountingStatistics.RejectedReads = ReadInteger(C_PROFILE_INFO_SECTION, C_REJECTED_READ_COUNT_STR, 0);
	aCountingStatistics.RejectedReferenceReads = ReadInteger(C_PROFILE_INFO_SECTION, C_EXCLUDED_REFERENCE_REJECTED_READ_COUNT, 0);
	aCountingStatistics.RejectedQualityReads = ReadInteger(C_PROFILE_INFO_SECTION, C_QUALITY_REJECTED_READ_COUNT_STR, 0);
	aCountingStatistics.NotMappedReads = ReadInteger(C_PROFILE_INFO_SECTION, C_NOT_COUNTED_READ_COUNT_STR, 0);
	aCountingStatistics.NotMappedCleanReads = ReadInteger(C_PROFILE_INFO_SECTION, C_NOT_COUNTED_CLEAN_READ_COUNT_STR, 0);

	return aCountingStatistics;
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetCountingStatistics(const TCountingStatistics & aCountingStatistics)
{
  WriteInteger(C_PROFILE_INFO_SECTION, C_SAMPLE_RAW_TOTAL_READ_COUNT_STR, aCountingStatistics.RawSequencedReads);
  WriteInteger(C_PROFILE_INFO_SECTION, C_SAMPLE_HOST_READ_COUNT_STR, aCountingStatistics.HostSequencedReads);
  WriteInteger(C_PROFILE_INFO_SECTION, C_SAMPLE_TOTAL_READ_COUNT_STR, aCountingStatistics.SequencedReads);
  WriteInteger(C_PROFILE_INFO_SECTION, C_SAMPLE_INDEXED_READ_COUNT_STR, aCountingStatistics.IndexedReads);
  WriteInteger(C_PROFILE_INFO_SECTION, C_COUNTED_READ_COUNT_STR, aCountingStatistics.TotalCleanCountedReads);
  WriteInteger(C_PROFILE_INFO_SECTION, C_UNIQUE_COUNTED_READ_COUNT_STR, aCountingStatistics.UniqueCleanCountedReads);
  WriteInteger(C_PROFILE_INFO_SECTION, C_REJECTED_READ_COUNT_STR, aCountingStatistics.RejectedReads);
  WriteInteger(C_PROFILE_INFO_SECTION, C_EXCLUDED_REFERENCE_REJECTED_READ_COUNT, aCountingStatistics.RejectedReferenceReads);
  WriteInteger(C_PROFILE_INFO_SECTION, C_QUALITY_REJECTED_READ_COUNT_STR, aCountingStatistics.RejectedQualityReads);
  WriteInteger(C_PROFILE_INFO_SECTION, C_NOT_COUNTED_READ_COUNT_STR, aCountingStatistics.NotMappedReads);
  WriteInteger(C_PROFILE_INFO_SECTION, C_NOT_COUNTED_CLEAN_READ_COUNT_STR, aCountingStatistics.NotMappedCleanReads);
}

// ------------------------------------------------------------------------------

TCountingStatistics TProfileIniFile::GetCountingStatistics(const int aProfiledSampleIndex) const
{
	TCountingStatistics aCountingStatistics;

	aCountingStatistics.RawSequencedReads = ReadInteger(C_PREFIX_PROFILED_SAMPLE_SECTION + "_" + numberToString(aProfiledSampleIndex), C_SAMPLE_RAW_TOTAL_READ_COUNT_STR, 0);
	aCountingStatistics.HostSequencedReads = ReadInteger(C_PREFIX_PROFILED_SAMPLE_SECTION + "_" + numberToString(aProfiledSampleIndex), C_SAMPLE_HOST_READ_COUNT_STR, 0);
	aCountingStatistics.SequencedReads = ReadInteger(C_PREFIX_PROFILED_SAMPLE_SECTION + "_" + numberToString(aProfiledSampleIndex), C_SAMPLE_TOTAL_READ_COUNT_STR, 0);
	aCountingStatistics.IndexedReads = ReadInteger(C_PREFIX_PROFILED_SAMPLE_SECTION + "_" + numberToString(aProfiledSampleIndex), C_SAMPLE_INDEXED_READ_COUNT_STR, 0);
	aCountingStatistics.TotalCleanCountedReads = ReadInteger(C_PREFIX_PROFILED_SAMPLE_SECTION + "_" + numberToString(aProfiledSampleIndex), C_COUNTED_READ_COUNT_STR, 0);
	aCountingStatistics.UniqueCleanCountedReads = ReadInteger(C_PREFIX_PROFILED_SAMPLE_SECTION + "_" + numberToString(aProfiledSampleIndex), C_UNIQUE_COUNTED_READ_COUNT_STR, 0);
	aCountingStatistics.RejectedReads = ReadInteger(C_PREFIX_PROFILED_SAMPLE_SECTION + "_" + numberToString(aProfiledSampleIndex), C_REJECTED_READ_COUNT_STR, 0);
	aCountingStatistics.RejectedReferenceReads = ReadInteger(C_PREFIX_PROFILED_SAMPLE_SECTION + "_" + numberToString(aProfiledSampleIndex), C_EXCLUDED_REFERENCE_REJECTED_READ_COUNT, 0);
	aCountingStatistics.RejectedQualityReads = ReadInteger(C_PREFIX_PROFILED_SAMPLE_SECTION + "_" + numberToString(aProfiledSampleIndex), C_QUALITY_REJECTED_READ_COUNT_STR, 0);
	aCountingStatistics. NotMappedReads = ReadInteger(C_PREFIX_PROFILED_SAMPLE_SECTION + "_" + numberToString(aProfiledSampleIndex), C_NOT_COUNTED_READ_COUNT_STR, 0);
	aCountingStatistics.NotMappedCleanReads = ReadInteger(C_PREFIX_PROFILED_SAMPLE_SECTION + "_" + numberToString(aProfiledSampleIndex), C_NOT_COUNTED_CLEAN_READ_COUNT_STR, 0);

	return aCountingStatistics;
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetExcludedReferenceCountingStatistics(const TStrings & aExcludedReferenceNameList, const vector<int> & aExcludedReferenceReadCountArray)
{
  WriteInteger(C_PROFILE_INFO_SECTION, C_EXCLUDED_REFERENCE_COUNT_STR, aExcludedReferenceNameList.size());
  for(size_t i = 0; i < aExcludedReferenceNameList.size(); i++)
  {
    WriteString(C_PROFILE_INFO_SECTION, C_PREFIX_EXCLUDED_REFERENCE_NAME_STR + numberToString(i+1), aExcludedReferenceNameList[i]);
    WriteInteger(C_PROFILE_INFO_SECTION, C_EXCLUDED_REFERENCE_REJECTED_READ_COUNT + "_" + numberToString(i+1), aExcludedReferenceReadCountArray[i]);
  }
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetExcludedReferenceCountingStatistics(const int aProfiledSampleIndex, const TStrings & aExcludedReferenceNameList, const vector<int> & aExcludedReferenceReadCountArray)
{
  string aSection = C_PREFIX_PROFILED_SAMPLE_SECTION + "_" + numberToString(aProfiledSampleIndex);

  WriteInteger(aSection, C_EXCLUDED_REFERENCE_COUNT_STR, aExcludedReferenceNameList.size());
  for (size_t i = 0; i < aExcludedReferenceNameList.size(); i++)
  {
    WriteString(aSection, C_PREFIX_EXCLUDED_REFERENCE_NAME_STR + numberToString(i+1), aExcludedReferenceNameList[i]);
    WriteInteger(aSection, C_EXCLUDED_REFERENCE_REJECTED_READ_COUNT + "_" + numberToString(i+1), aExcludedReferenceReadCountArray[i]);
  }
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetExcludedReferenceNameList(const TStrings & aExcludedReferenceNameList){
  WriteInteger(C_PROFILE_INFO_SECTION, C_EXCLUDED_REFERENCE_COUNT_STR, aExcludedReferenceNameList.size());
  for (size_t i = 0; i < aExcludedReferenceNameList.size(); i++)
  {
    WriteString(C_PROFILE_INFO_SECTION, C_PREFIX_EXCLUDED_REFERENCE_NAME_STR + numberToString(i + 1), aExcludedReferenceNameList[i]);
  }
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetIndexedSequencedBaseCount(const int64 Value)
{
  WriteString(C_PROFILE_INFO_SECTION, C_SAMPLE_INDEXED_BASE_COUNT_STR, numberToString(Value));
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetIndexedSequencedReadCount(const int Value)
{
  WriteInteger(C_PROFILE_INFO_SECTION, C_SAMPLE_INDEXED_READ_COUNT_STR, Value);
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetIntersectCount(const int aProfiledSampleIndex, const int aIntersectDescrIndex, const int aCountValue)
{
  WriteInteger(C_PREFIX_PROFILED_SAMPLE_SECTION + "_" + numberToString(aProfiledSampleIndex), C_INTERSECT_COUNT_PREFIX_STR + numberToString(aIntersectDescrIndex), aCountValue);
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetIntersectCount(const int aIntersectDescrIndex, const int aCountValue)
{
  WriteInteger(C_PROFILE_INFO_SECTION, C_INTERSECT_COUNT_PREFIX_STR + numberToString(aIntersectDescrIndex), aCountValue);
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetIsBestAlignment(const bool Value)
{
  WriteBool(C_PROFILE_INFO_SECTION, C_KEEP_ONLY_BEST_MATCHES_STR, Value);
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetIsMappingMismatchesPercentage(const bool Value)
{
  WriteBool(C_PROFILE_INFO_SECTION, C_MAPPING_IS_MISMATCHES_PERCENTAGE_STR, Value);
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetMappedReadLength(const int Value)
{
  WriteInteger(C_PROFILE_INFO_SECTION, C_MAPPED_READ_LENGTH_STR, Value);
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetMappedReadLengthType(const TMappedReadLengthType Value)
{
  switch(Value) {
    case mrltFixed:
      WriteString(C_PROFILE_INFO_SECTION, C_MAPPED_READ_LENGTH_TYPE_STR, C_FIXED_MAPPED_READ_LENGTH_TYPE);
      break;
    case mrltMin:
      WriteString(C_PROFILE_INFO_SECTION, C_MAPPED_READ_LENGTH_TYPE_STR, C_MIN_MAPPED_READ_LENGTH_TYPE);
      break;
    case mrltMax:
      WriteString(C_PROFILE_INFO_SECTION, C_MAPPED_READ_LENGTH_TYPE_STR, C_MAX_MAPPED_READ_LENGTH_TYPE);
      break;
    case mrltOverall:
      WriteString(C_PROFILE_INFO_SECTION, C_MAPPED_READ_LENGTH_TYPE_STR, C_OVERALL_MAPPED_READ_LENGTH_TYPE);
      break;
  }
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetMappingCmdLine(const string & Value)
{
  WriteString(C_PROFILE_INFO_SECTION, C_MAPPING_CMDLINE_STR, Value);
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetMappingMatchesCount(const int Value)
{
  WriteInteger(C_PROFILE_INFO_SECTION, C_MAPPING_MATCHES_STR, Value);
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetMappingMismatchesCount(const int Value)
{
  WriteInteger(C_PROFILE_INFO_SECTION, C_MAPPING_MISMATCHES_STR, Value);
}

// ------------------------------------------------------------------------------
//// LOCAL
int TProfileIniFile::GetAlignmentLengthCutoff() const
{
  if (GetIsLocalMapping())
    return ReadInteger(C_PROFILE_INFO_SECTION, C_ALIGNMENT_LENGTH_CUTOFF_STR, C_DEFAULT_ALIGNMENT_LENGTH_CUTOFF);
  return -1;
}

// ------------------------------------------------------------------------------
//// LOCAL
int TProfileIniFile::GetSoftClippingLengthCutoff() const
{
  if (GetIsLocalMapping())
    return ReadInteger(C_PROFILE_INFO_SECTION, C_SOFT_CLIPPING_LENGTH_CUTOFF_STR, C_DEFAULT_SOFT_CLIPPING_LENGTH_CUTOFF);
  return -1;
}

// ------------------------------------------------------------------------------

bool TProfileIniFile::GetIsLocalMapping() const
{
  return ReadBool(C_PROFILE_INFO_SECTION, C_MAPPING_IS_LOCAL_MAPPING_STR, false);
}

// ------------------------------------------------------------------------------

bool TProfileIniFile::GetKeepInternalLocalAlignment() const
{
  if (GetIsLocalMapping())
    return ReadBool(C_PROFILE_INFO_SECTION, C_KEEP_INTERNAL_LOCAL_ALIGNMENT_STR, false);

  return false;
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetAlignmentLengthCutoff(const int Value)
{
  WriteInteger(C_PROFILE_INFO_SECTION, C_ALIGNMENT_LENGTH_CUTOFF_STR, Value);
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetSoftClippingLengthCutoff(const int Value)
{
  WriteInteger(C_PROFILE_INFO_SECTION, C_SOFT_CLIPPING_LENGTH_CUTOFF_STR, Value);
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetIsLocalMapping(const bool Value)
{
  WriteBool(C_PROFILE_INFO_SECTION, C_MAPPING_IS_LOCAL_MAPPING_STR, Value);
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetKeepInternalLocalAlignment(const bool Value)
{
  WriteBool(C_PROFILE_INFO_SECTION, C_KEEP_INTERNAL_LOCAL_ALIGNMENT_STR, Value);
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetMappingToolName(const string & Value)
{
  WriteString(C_PROFILE_INFO_SECTION, C_MAPPING_SOFTWARE_STR, Value);
}

// ------------------------------------------------------------------------------

string TProfileIniFile::GetProjectName() const
{
  return ReadString(C_PROFILE_INFO_SECTION, C_PROJECT_NAME_STR, "");
}

// ------------------------------------------------------------------------------

TReadCleaningMethodType TProfileIniFile::GetReadCleaningMethod() const
{
  string aReadCleaningMethodStr = ReadString(C_PROFILE_INFO_SECTION, C_SAMPLE_READ_CLEANING_METHOD_STR, "");
  if (aReadCleaningMethodStr == C_METEOR_ALIENTRIMMER_PROGRAM_NAME) return rcmAlienTrimmer;
  return rcmNone;
}

// ------------------------------------------------------------------------------

string TProfileIniFile::GetReadCleaningParameters() const
{
  return ReadString(C_PROFILE_INFO_SECTION, C_SAMPLE_READ_CLEANING_PARAMETERS_STR, "");
}

// ------------------------------------------------------------------------------

int TProfileIniFile::GetSequencedRawReadCount() const
{
  int Result = ReadInteger(C_PROFILE_INFO_SECTION, C_SAMPLE_RAW_TOTAL_READ_COUNT_STR, -1);
  return Result == -1 ? GetSequencedReadCount() : Result;
}

// ------------------------------------------------------------------------------

int TProfileIniFile::GetSequencedHostReadCount() const
{
  return ReadInteger(C_PROFILE_INFO_SECTION, C_SAMPLE_HOST_READ_COUNT_STR, -1);
}

// ------------------------------------------------------------------------------

int TProfileIniFile::GetSequencedReadCount() const
{
  return ReadInteger(C_PROFILE_INFO_SECTION, C_SAMPLE_TOTAL_READ_COUNT_STR, -1);
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetSequencedRawReadCount(const int Value)
{
  WriteInteger(C_PROFILE_INFO_SECTION, C_SAMPLE_RAW_TOTAL_READ_COUNT_STR, Value);
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetSequencedHostReadCount(const int Value)
{
  WriteInteger(C_PROFILE_INFO_SECTION, C_SAMPLE_HOST_READ_COUNT_STR, Value);
}

// ------------------------------------------------------------------------------

void TProfileIniFile::SetSequencedReadCount(const int Value)
{
  WriteInteger(C_PROFILE_INFO_SECTION, C_SAMPLE_TOTAL_READ_COUNT_STR, Value);
}

// --------------------------------------------------------------------

int TProfileIniFile::GetProfiledSampleCount() const
{
	// -1 because profile_info section is not counted
	return getSectionsNumber()-1;
}

// --------------------------------------------------------------------

bool TProfileIniFile::GetProfiledSample(const int aProfileSampleIndex, map<string, string> & aMapSection) const
{
	aMapSection.clear();
	return getSection(C_PREFIX_PROFILED_SAMPLE_SECTION + "_" + numberToString(aProfileSampleIndex), aMapSection);
}

// ------------------------------------------------------------------------------

int TProfileIniFile::AddProfiledSample(const map<string, string> & aMapSection){

  int aProfiledSampleIndex = GetProfiledSampleCount();

  setSection(C_PREFIX_PROFILED_SAMPLE_SECTION + "_" + numberToString(aProfiledSampleIndex), aMapSection);

  return aProfiledSampleIndex;
}

