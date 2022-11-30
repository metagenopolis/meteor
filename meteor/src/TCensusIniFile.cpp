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

#include <vector>
#include <string>
#include <cstdlib>
#include "utils.h"
#include "TCensusIniFile.h"

using namespace std;

TCensusIniFile::TCensusIniFile(const string & filename, bool aIsReadOnly):IniFileManager() {
	FReadOnly = aIsReadOnly;
	// do not load existing data if open in write mode (overwrite)
	if (FReadOnly) loadIniData(filename);
	else iniFilename = filename;
}

// ------------------------------------------------------------------------------

TCensusIniFile::~TCensusIniFile() {
	// with no argument, printIniData() overwrite IniFileManager::iniFilename (inherited).
	if (!FReadOnly)
		printIniData();
}

// ------------------------------------------------------------------------------

string TCensusIniFile::GetConditionName() const
{
  string Result = ReadString(C_SAMPLE_INFO_SECTION, C_CONDITION_NAME_STR, C_DEFAULT_CONDITION_NAME);
  return Result.empty() ? C_DEFAULT_CONDITION_NAME : Result;
}

// ------------------------------------------------------------------------------

string TCensusIniFile::GetFullSampleName() const
{
  return ReadString(C_SAMPLE_INFO_SECTION, C_SAMPLE_FULL_NAME_STR, "");
}

// ------------------------------------------------------------------------------

int64 TCensusIniFile::GetIndexedSequenceBaseCount() const
{
  return atoll(ReadString(C_SAMPLE_INFO_SECTION, C_SAMPLE_INDEXED_BASE_COUNT_STR, "-1").c_str());
}

// ------------------------------------------------------------------------------

int TCensusIniFile::GetIndexedSequencedReadCount() const
{
  return ReadInteger(C_SAMPLE_INFO_SECTION, C_SAMPLE_INDEXED_READ_COUNT_STR, -1);
}

// ------------------------------------------------------------------------------

string TCensusIniFile::GetProjectName() const
{
  return ReadString(C_SAMPLE_INFO_SECTION, C_PROJECT_NAME_STR, "");
}

// ------------------------------------------------------------------------------

TReadCleaningMethodType TCensusIniFile::GetReadCleaningMethod() const
{
  string aReadCleaningMethodStr = ReadString(C_SAMPLE_INFO_SECTION, C_SAMPLE_READ_CLEANING_METHOD_STR, "");
  return aReadCleaningMethodStr == C_METEOR_ALIENTRIMMER_PROGRAM_NAME ? rcmAlienTrimmer : rcmNone;
}

// ------------------------------------------------------------------------------

string TCensusIniFile::GetReadCleaningParameters() const
{
  return ReadString(C_SAMPLE_INFO_SECTION, C_SAMPLE_READ_CLEANING_PARAMETERS_STR, "");
}

// ------------------------------------------------------------------------------

string TCensusIniFile::GetSampleName() const
{
  return ReadString(C_SAMPLE_INFO_SECTION, C_SAMPLE_NAME_STR, "");
}

// ------------------------------------------------------------------------------

int TCensusIniFile::GetSequencedRawReadCount() const
{
  int Result = ReadInteger(C_SAMPLE_INFO_SECTION, C_SAMPLE_RAW_TOTAL_READ_COUNT_STR, -1);
  return Result == -1 ? GetSequencedReadCount() : Result;
}

// -------------------------------------------------------------------------------

int TCensusIniFile::GetSequencedReadCount() const
{
  return ReadInteger(C_SAMPLE_INFO_SECTION, C_SAMPLE_TOTAL_READ_COUNT_STR, -1);
}

// -------------------------------------------------------------------------------

int TCensusIniFile::GetSequencedHostReadCount() const
{
  return ReadInteger(C_SAMPLE_INFO_SECTION, C_SAMPLE_HOST_READ_COUNT_STR, -1);
}

// ------------------------------------------------------------------------------

TMappingFileFormat TCensusIniFile::GetMappingFileFormat() const {


	  TMappingFileFormat aMappingFileFormat;
	  string aMappingFileFormatStr = ReadString(C_MAPPING_FILE_SECTION, C_MAPPING_FILE_FORMAT_STR, C_MAPPING_UNKNOWN_OUTPUT_STR);

	  if(aMappingFileFormatStr == C_MAPPING_DEFAULT_BOWTIE_OUTPUT_STR)
	    aMappingFileFormat = mffDefaultBowtie;
	  else if(aMappingFileFormatStr == C_MAPPING_LITE_DEFAULT_BOWTIE_OUTPUT_STR)
	    aMappingFileFormat = mffLiteDefaultBowtie;
	  else if(aMappingFileFormatStr == C_MAPPING_SAM_OUTPUT_STR)
	    aMappingFileFormat = mffSAM;
	  else if(aMappingFileFormatStr == C_MAPPING_BAM_OUTPUT_STR)
	    aMappingFileFormat = mffBAM;
	  else if(aMappingFileFormatStr == C_MAPPING_CSFASTA_OUTPUT_STR)
	    aMappingFileFormat = mffCSFasta;
	  else
	    aMappingFileFormat = mffUnknown;

	  // Test on mapping tool name (old version of census_ini file)
	  if (aMappingFileFormat == mffUnknown)
	  {
		string aMappingTool = ReadString(C_MAPPING_SECTION, C_MAPPING_TOOL_STR, "");
	    if(aMappingTool == "corona_lite_v2")
	      aMappingFileFormat = mffCSFasta;
	    else if(aMappingTool ==  "bowtie")
	      aMappingFileFormat = mffDefaultBowtie;
	    else
	      aMappingFileFormat = mffUnknown;
	  }
	  return aMappingFileFormat;
}

// ------------------------------------------------------------------------------

void TCensusIniFile::SetMappingFileFormat(const TMappingFileFormat & Value)
{
	//// passer directement la string C_MAPPING_XXXX_OUTPUT_STR a la fonction ????
	switch(Value) {
		case mffCSFasta:
			WriteString(C_MAPPING_FILE_SECTION, C_MAPPING_FILE_FORMAT_STR, C_MAPPING_CSFASTA_OUTPUT_STR);
			break;
		case mffSAM:
			WriteString(C_MAPPING_FILE_SECTION, C_MAPPING_FILE_FORMAT_STR, C_MAPPING_SAM_OUTPUT_STR);
			break;
		case mffBAM:
			WriteString(C_MAPPING_FILE_SECTION, C_MAPPING_FILE_FORMAT_STR, C_MAPPING_BAM_OUTPUT_STR);
			break;
		case mffDefaultBowtie:
			WriteString(C_MAPPING_FILE_SECTION, C_MAPPING_FILE_FORMAT_STR, C_MAPPING_DEFAULT_BOWTIE_OUTPUT_STR);
			break;
		case mffLiteDefaultBowtie:
			WriteString(C_MAPPING_FILE_SECTION, C_MAPPING_FILE_FORMAT_STR, C_MAPPING_LITE_DEFAULT_BOWTIE_OUTPUT_STR);
			break;
		case mffUnknown:
			WriteString(C_MAPPING_FILE_SECTION, C_MAPPING_FILE_FORMAT_STR, C_MAPPING_UNKNOWN_OUTPUT_STR);
			break;
	}
}

// ------------------------------------------------------------------------------

string TCensusIniFile::GetSequencingDate() const {
	  return ReadString(C_SAMPLE_INFO_SECTION, C_SEQUENCING_DATE_STR, "0");
}

// ------------------------------------------------------------------------------

bool TCensusIniFile::GetIsLocalMapping() const {
	return ReadBool(C_MAPPING_SECTION, C_MAPPING_IS_LOCAL_MAPPING_STR, false);
}

// ------------------------------------------------------------------------------

void TCensusIniFile::SetIsLocalMapping(const bool Value){
	WriteBool(C_MAPPING_SECTION, C_MAPPING_IS_LOCAL_MAPPING_STR, Value);
}

// ------------------------------------------------------------------------------

int TCensusIniFile::GetAlignmentLengthCutoff() const{
	if (GetIsLocalMapping()) return ReadInteger(C_COUNTING_SECTION, C_ALIGNMENT_LENGTH_CUTOFF_STR, C_DEFAULT_ALIGNMENT_LENGTH_CUTOFF);
	return -1;
}

// ------------------------------------------------------------------------------

int TCensusIniFile::GetSoftClippingLengthCutoff() const{
	if (GetIsLocalMapping()) return ReadInteger(C_COUNTING_SECTION, C_SOFT_CLIPPING_LENGTH_CUTOFF_STR, C_DEFAULT_SOFT_CLIPPING_LENGTH_CUTOFF);
	return -1;
}

// ------------------------------------------------------------------------------

bool TCensusIniFile::GetKeepInternalLocalAlignment() const{
	  if (GetIsLocalMapping()) return ReadBool(C_COUNTING_SECTION, C_KEEP_INTERNAL_LOCAL_ALIGNMENT_STR, false);
	  return false;
}

// ------------------------------------------------------------------------------

void TCensusIniFile::SetAlignmentLengthCutoff(const int Value){
	WriteInteger(C_COUNTING_SECTION, C_ALIGNMENT_LENGTH_CUTOFF_STR, Value);
}

// ------------------------------------------------------------------------------

void TCensusIniFile::SetSoftClippingLengthCutoff(const int Value){
	WriteInteger(C_COUNTING_SECTION, C_SOFT_CLIPPING_LENGTH_CUTOFF_STR, Value);
}

// ------------------------------------------------------------------------------

void TCensusIniFile::SetKeepInternalLocalAlignment(const bool Value){
	WriteBool(C_COUNTING_SECTION, C_KEEP_INTERNAL_LOCAL_ALIGNMENT_STR, Value);
}

// ------------------------------------------------------------------------------

TCountingStatistics TCensusIniFile::GetCountingStatistics() const
{
	TCountingStatistics Result;

	Result.TotalCleanCountedReads  = ReadInteger(C_COUNTING_SECTION, C_COUNTED_READ_COUNT_STR, 0);
	Result.UniqueCleanCountedReads = ReadInteger(C_COUNTING_SECTION, C_UNIQUE_COUNTED_READ_COUNT_STR, 0);
	Result.RejectedReads           = ReadInteger(C_COUNTING_SECTION, C_REJECTED_READ_COUNT_STR, 0);
	Result.RejectedReferenceReads  = ReadInteger(C_COUNTING_SECTION, C_EXCLUDED_REFERENCE_REJECTED_READ_COUNT, 0);
	Result.RejectedQualityReads    = ReadInteger(C_COUNTING_SECTION, C_QUALITY_REJECTED_READ_COUNT_STR, 0);
	Result.NotMappedReads          = ReadInteger(C_COUNTING_SECTION, C_NOT_COUNTED_READ_COUNT_STR, 0);
	Result.NotMappedCleanReads     = ReadInteger(C_COUNTING_SECTION, C_NOT_COUNTED_CLEAN_READ_COUNT_STR, 0);

	return Result;
}

// ------------------------------------------------------------------------------

void TCensusIniFile::GetExcludedReferenceNameList(TStrings & aExcludedReferenceNameList) const {

	int excludedRefCount = ReadInteger(C_COUNTING_SECTION, C_EXCLUDED_REFERENCE_COUNT_STR, 0);
	aExcludedReferenceNameList.clear();
	aExcludedReferenceNameList.resize(excludedRefCount);

	for (int i = 0; i < excludedRefCount; i++) {
//		aExcludedReferenceNameList[i] = ReadString(C_COUNTING_SECTION, C_PREFIX_EXCLUDED_REFERENCE_NAME_STR + to_string(i+1), ""); //// c++11
		aExcludedReferenceNameList.at(i) = ReadString(C_COUNTING_SECTION, C_PREFIX_EXCLUDED_REFERENCE_NAME_STR + numberToString(i+1), ""); //// AT
	}
}

// ------------------------------------------------------------------------------

void TCensusIniFile::GetMappingBowtieFileNameList(TStrings & aMappingBowtieFileNameList){

	int aMappingFileCount;
	string aMappingBowtieFileName;

	aMappingBowtieFileNameList.clear();

	aMappingFileCount = ReadInteger(C_MAPPING_FILE_SECTION, C_MAPPING_FILE_COUNT_STR, 0);

	if (aMappingFileCount == 0) {

		aMappingBowtieFileName = ReadString(C_MAPPING_FILE_SECTION, C_MAPPING_BOWTIE_FILENAME_STR, "");

		if (!aMappingBowtieFileName.empty())
			aMappingBowtieFileNameList.push_back(aMappingBowtieFileName);
	}
	else {
		aMappingBowtieFileNameList.reserve(aMappingFileCount);
		for (int i=1; i<= aMappingFileCount; i++) {
//			aMappingBowtieFileName = ReadString(C_MAPPING_FILE_SECTION, C_MAPPING_BOWTIE_FILENAME_STR + "_" + to_string(i), ""); //// C++11
			aMappingBowtieFileName = ReadString(C_MAPPING_FILE_SECTION, C_MAPPING_BOWTIE_FILENAME_STR + "_" + numberToString(i), "");
			if (!aMappingBowtieFileName.empty())
				aMappingBowtieFileNameList.push_back(aMappingBowtieFileName);
		}
	}
}

// ------------------------------------------------------------------------------

string TCensusIniFile::GetMappingTool() const {
  return ReadString(C_MAPPING_SECTION, C_MAPPING_TOOL_STR, "");
}

// ------------------------------------------------------------------------------
//// ajouter un membre privé TSampleFileInfo samplefileInfo initialisé dans le constructeur et simplement faire return samplefileInfo ????
TSampleFileInfo TCensusIniFile::GetSampleFileInfo() const {

TSampleFileInfo Result;
  switch(GetSequenceFileFormat()) {
    case sffCSFastaAndQual:
        Result.Exists = true;
        Result.Tag = ReadString(C_SAMPLE_INFO_SECTION, C_LIBRARY_TAG_STR, "");
        Result.Name = ReadString(C_SAMPLE_INFO_SECTION, C_SAMPLE_FULL_NAME_STR, "");
        Result.IsCompressed = GetIsCompressedSampleFile();
        Result.SequenceFileFormat = sffCSFastaAndQual;
        Result.CSFastaFileName = ReadString(C_SAMPLE_FILE_SECTION, C_CSFASTA_FILENAME_STR, "");
        Result.QualFastaFileName = ReadString(C_SAMPLE_FILE_SECTION, C_QUAL_FILENAME_STR, "");
        Result.FastQFileName = "";
        break;
    case sffFastQ:
        Result.Exists = true;
        Result.Tag = ReadString(C_SAMPLE_INFO_SECTION, C_LIBRARY_TAG_STR, "");
        Result.Name = ReadString(C_SAMPLE_INFO_SECTION, C_SAMPLE_FULL_NAME_STR, "");
        Result.IsCompressed = GetIsCompressedSampleFile();
        Result.SequenceFileFormat = sffFastQ;
        Result.CSFastaFileName = "";
        Result.QualFastaFileName = "";
        Result.FastQFileName = ReadString(C_SAMPLE_FILE_SECTION, C_FASTQ_FILENAME_STR, "");
        break;
    default: break;
  }
  return Result;
}

// ------------------------------------------------------------------------------
//// voir commentaire pour TCensusIniFile::GetSampleFileInfo() ????
TSequenceFileFormat TCensusIniFile::GetSequenceFileFormat() const
{
	string aSequencingDevice = ReadString(C_SAMPLE_INFO_SECTION, C_SEQUENCING_DEVICE_STR, "");
	if (aSequencingDevice == C_SOLID_DEVICE_STR) return sffCSFastaAndQual;
	if (aSequencingDevice == C_PROTON_DEVICE_STR || aSequencingDevice == C_ILLUMINA_DEVICE_STR || aSequencingDevice == C_SOLEXA_DEVICE_STR)
		return sffFastQ;
	else return sffUnknown;
}

// ------------------------------------------------------------------------------

bool TCensusIniFile::GetIsCompressedSampleFile() const {
	return ReadBool(C_SAMPLE_FILE_SECTION, C_SAMPLE_FILE_IS_COMPRESSED_STR, false);
}

// ------------------------------------------------------------------------------

void TCensusIniFile::SetIsCompressedSampleFile(const bool & value) {
	WriteBool(C_SAMPLE_FILE_SECTION, C_SAMPLE_FILE_IS_COMPRESSED_STR, value);
}

// ------------------------------------------------------------------------------

int TCensusIniFile::GetProcessedReadCount() const {
	return ReadInteger(C_MAPPING_SECTION, C_PROCESSED_READ_COUNT_STR, -1);
}

// ------------------------------------------------------------------------------

void TCensusIniFile::SetProcessedReadCount(const int & value) {
	WriteInteger(C_MAPPING_SECTION, C_PROCESSED_READ_COUNT_STR, value);
}

// ------------------------------------------------------------------------------

bool TCensusIniFile::GetReadOnly() const { return FReadOnly;}
void TCensusIniFile::SetReadOnly(bool Value) { FReadOnly=Value; }

// ------------------------------------------------------------------------------

string TCensusIniFile::GetSequencingDevice() const {
  return ReadString(C_SAMPLE_INFO_SECTION, C_SEQUENCING_DEVICE_STR, "");
}

// ------------------------------------------------------------------------------

//// plus tard ????
//void TCensusIniFile::SetLibraryIndexingSoftware(
//		const string & aSoftwareName, const string & aSoftwareVersion,
//		const string & aSoftwareDigest, const TDigestAlgorithm & aDigestAlgorithm)
//{
//	string aDigestAttribute;
//	switch(aDigestAlgorithm) {
//		case daMD5:
//			aDigestAttribute = C_SAMPLE_INDEXING_SOFTWARE_DIGEST_STR + "_" + C_MD5_DIGEST_ALGORITHM;
//			break;
//		case daSHA256:
//			aDigestAttribute = C_SAMPLE_INDEXING_SOFTWARE_DIGEST_STR + "_" + C_SHA256_DIGEST_ALGORITHM;
//			break;
//		case daNone:
//			break;
//	}
//	WriteString(C_SAMPLE_INFO_SECTION, C_SAMPLE_INDEXING_SOFTWARE_STR, aSoftwareName);
//	WriteString(C_SAMPLE_INFO_SECTION, C_SAMPLE_INDEXING_SOFTWARE_VERSION_STR, aSoftwareVersion);
//
//	if (aDigestAlgorithm != daNone)
//		WriteString(C_SAMPLE_INFO_SECTION, aDigestAttribute, aSoftwareDigest);
//}

// ------------------------------------------------------------------------------
//// plus tard ????
//void TCensusIniFile::SetMappingBowtieFileNameList(const TStrings & aMappingBowtieFileNameList);
//var
//  int i;
//{
//  WriteInteger(C_MAPPING_FILE_SECTION, C_MAPPING_FILE_COUNT_STR,
//    aMappingBowtieFileNameList.Count);
//  for i = 0 to aMappingBowtieFileNameList.Count - 1 do
//    WriteString(C_MAPPING_FILE_SECTION, C_MAPPING_BOWTIE_FILENAME_STR + '_' +
//      IntToStr(i + 1), ExtractFileName(aMappingBowtieFileNameList.Strings[i]));
//}

// ------------------------------------------------------------------------------
//// plus tard ????
//void TCensusIniFile::SetMappingSoftware(const string & aSoftwareName, const string & aSoftwareVersion, const string & aSoftwareDigest; const TDigestAlgorithm & aDigestAlgorithm);
//var
//  string aDigestAttribute;
//{
//  switch(aDigestAlgorithm) {
//    case daMD5:
//      aDigestAttribute = C_MAPPING_SOFTWARE_DIGEST_STR + '_' +
//        C_MD5_DIGEST_ALGORITHM;
//    case daSHA256:
//      aDigestAttribute = C_MAPPING_SOFTWARE_DIGEST_STR + '_' +
//        C_SHA256_DIGEST_ALGORITHM;
//    case daNone:
//      ;
//  }
//  WriteString(C_MAPPING_SECTION, C_MAPPING_SOFTWARE_STR, aSoftwareName);
//  WriteString(C_MAPPING_SECTION, C_MAPPING_SOFTWARE_VERSION_STR,
//    aSoftwareVersion);
//  if aDigestAlgorithm != daNone then
//    WriteString(C_MAPPING_SECTION, aDigestAttribute, aSoftwareDigest);
//}

// ------------------------------------------------------------------------------
//// plus tard ????
//void TCensusIniFile::SetMappingTool(const string & aToolName, const string & aToolVersion, const string & aToolDigest; const TDigestAlgorithm & aDigestAlgorithm);
//var
//  string aDigestAttribute;
//{
//  switch(aDigestAlgorithm) {
//    case daMD5:
//      aDigestAttribute = C_MAPPING_TOOL_DIGEST_STR + '_' +
//        C_MD5_DIGEST_ALGORITHM;
//    case daSHA256:
//      aDigestAttribute = C_MAPPING_TOOL_DIGEST_STR + '_' +
//        C_SHA256_DIGEST_ALGORITHM;
//    case daNone:
//      ;
//  }
//  WriteString(C_MAPPING_SECTION, C_MAPPING_TOOL_STR, aToolName);
//  WriteString(C_MAPPING_SECTION, C_MAPPING_TOOL_VERSION_STR, aToolVersion);
//  if aDigestAlgorithm != daNone then
//    WriteString(C_MAPPING_SECTION, aDigestAttribute, aToolDigest);
//}
