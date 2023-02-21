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

#include "TMappingDataFile.h"
#include "TProfileIniFile.h"
#include "TMeteorCounter.h"
#include "TNGSCensusData.h"
#include "TNGSCensusFile.h"
#include "TProfileIniFile.h"
#include "TReferenceIniFile.h"
#include "utils.h"
#include <cmath>

using namespace std;

TMeteorCounter::TMeteorCounter(const vector<TCountedReadType> & aCountedReadTypeList):FMappingDataList(),FCountedReadTypeList(aCountedReadTypeList) {

}

TMeteorCounter::~TMeteorCounter() {

	for (size_t i = 0; i<FMappingDataList.size(); i++){ delete FMappingDataList[i]; FMappingDataList[i]=NULL; }
	////FReferenceAnnotationDatabase.Disconnect; //// unused
	////FreeAndNil(FReferenceAnnotationDatabase);

}

void TMeteorCounter::AddMappingLibrary(
		const string & aTmpDir,
		const string & aMainMappingCensusIniFileName,
		const double aMaxDistance,
		const bool aIsRelativeDistance,
		const bool aIsLocalAlignment, const bool aKeepInternalLocalAlignment,
		const int aAlignmentLengthCutoff, const int aSoftClippingLengthCutoff,
		const TStrings & aExcludedMappingCensusIniFileNameList,
		const TMappingReferencePropertiesArray & aExcludedMappingReferencePropertiesArray)
{
	double aExcludedMaxDistance;
	TStrings aMappingFileNameList;

	TCensusIniFile aMappingCensusIniFile(aMainMappingCensusIniFileName, true);

	aMappingCensusIniFile.GetMappingBowtieFileNameList(aMappingFileNameList);
	if (aMappingFileNameList.empty()) {
		cerr << "Error, file: "<< __FILE__ << ", line: " << __LINE__ <<", no mapping filename found in "<< aMainMappingCensusIniFileName << endl;
		exit(1);
	}

////  splited sam files SHOULD NOT HAPPEN (test bowtie2-build-l ok with big files)
//	  if (aMappingFileNameList.size() > 1 ) // use in memory dataset
//	  {}
	  //else // use file
	  //{

	    TMappingDataFile* aMappingData = new TMappingDataFile(
	    		aMappingCensusIniFile.GetProcessedReadCount(),
				aMaxDistance, aIsRelativeDistance,
				aIsLocalAlignment, aKeepInternalLocalAlignment, aAlignmentLengthCutoff, aSoftClippingLengthCutoff,
				aMappingCensusIniFile.GetMappingFileFormat());

	    // MS Windows bug when path is over MAX_PATH constant = 255 (see delphi code)

		string aDirName = ExtractFileDir(aMainMappingCensusIniFileName);
	    // the user has set a tmp directory path
	    // if (!aTmpDir.empty()) {
	    //     aDirName = aTmpDir + C_PATH_SEP + aDirName;
	    // }
      // MODIFIED
	    string aOriginalFileName = aDirName + C_PATH_SEP + aMappingFileNameList[0];
      // string aOriginalFileName = aDirName + C_PATH_SEP;
      if (aMappingFileNameList.size() != 1 ) cerr << "je fais chier" << endl; 
      cout << "size " << aMappingFileNameList.size() << endl;
      cout << "shit " << aOriginalFileName << endl;
      // string aOriginalFileName = aMappingFileNameList[0];
	    // string aOriginalFileName = ExtractFileDir(aMainMappingCensusIniFileName) + C_PATH_SEP + aMappingFileNameList[0];

		// aOriginalFileName (SAM file) is not read yet. First we read excluded mapping files (see the for loop below).
	    aMappingData->AddMainMapping(aOriginalFileName);
		//cerr << "aOriginalFileName (SAM): " <<aOriginalFileName << endl;

	    /*
	      aMappingData.AddMainMapping(ExtractFileDir(aMainMappingCensusIniFileName) + C_PATH_SEP + aMappingFileNameList[iMappingFile]);
	    */

	  //FreeAndNil(aMappingCensusIniFile);

	  // add excluded mapping. Each excluded MappingFileName (SAM file) is read once via AddSAMExcludedMapping
	  for (size_t iExcludedMapping = 0; iExcludedMapping < aExcludedMappingCensusIniFileNameList.size(); iExcludedMapping++)
	  {
	    aExcludedMaxDistance = (double)aExcludedMappingReferencePropertiesArray.at(iExcludedMapping).MappingMismatches; //// AT
	    if (aExcludedMappingReferencePropertiesArray.at(iExcludedMapping).IsMappingMismatchesPercentage)
	      aExcludedMaxDistance /= 100.0;
	    // relative distance

	      TCensusIniFile aMappingCensusIniFile(aExcludedMappingCensusIniFileNameList[iExcludedMapping], true);

	    aMappingFileNameList.clear();
	    aMappingCensusIniFile.GetMappingBowtieFileNameList(aMappingFileNameList);

	    for (size_t iMappingFile = 0; iMappingFile < aMappingFileNameList.size(); iMappingFile++) // should iterate only once
	    {
	      //// skip Windows bug when path is over MAX_PATH constant = 255 (see delphi code)
	      aDirName = ExtractFileDir(aExcludedMappingCensusIniFileNameList[iExcludedMapping]);
	      if (!aTmpDir.empty()) {
			aDirName = aTmpDir + C_PATH_SEP + aDirName;
		  }
	      aOriginalFileName = aDirName + C_PATH_SEP + aMappingFileNameList.at(iMappingFile); //// AT
	      //aOriginalFileName = ExtractFileDir(aExcludedMappingCensusIniFileNameList[iExcludedMapping]) + C_PATH_SEP + aMappingFileNameList.at(iMappingFile); //// AT

	      aMappingData->AddExcludedMapping(
	    		  aOriginalFileName, iExcludedMapping + 1, aExcludedMaxDistance,
				  aExcludedMappingReferencePropertiesArray.at(iExcludedMapping).IsMappingMismatchesPercentage,
				  aExcludedMappingReferencePropertiesArray.at(iExcludedMapping).IsOnBothStrand ); //// AT
		  // aOriginalFileName (SAM file) is read once via AddSAMExcludedMapping

	      /*
	        aMappingData->AddExcludedMapping(
	        	aOriginalFileName, iExcludedMapping+1, aExcludedMaxDistance,
	        	aExcludedMappingReferencePropertiesArray[iExcludedMapping].IsMappingMismatchesPercentage );
	      */
	  	 }
	  }
	  // Add MappingData to objectList for counting
	  FMappingDataList.push_back(aMappingData); // address memory assignation. Will be deleted in destructor ~TMeteorCounter()

	  ////aMappingData = NULL; //// not useful as end of method is reached.
}

// ------------------------------------------------------------------------------

void TMeteorCounter::GetCountingStatistics(TCountingStatistics & aCountingStatistics) const
{
	// total reads count after cleaning (on rejected reference and quality)
	aCountingStatistics.TotalCleanCountedReads = 0;
	// unique reads count after cleaning (on rejected reference and quality)
	aCountingStatistics.UniqueCleanCountedReads = 0;
	// all rejected reads (on rejected reference and quality)
	aCountingStatistics.RejectedReads = 0;
	// rejected reads only on rejected reference
	aCountingStatistics.RejectedReferenceReads = 0;
	// rejected reads only on quality
	aCountingStatistics.RejectedQualityReads = 0;
	// all other reads not mapped on all reference
	aCountingStatistics.NotMappedReads = 0;
	// all other reads not mapped on all reference after cleaning on quality
	aCountingStatistics.NotMappedCleanReads = 0;

	for (size_t iMappingData = 0; iMappingData < FMappingDataList.size(); iMappingData++){
		FMappingDataList[iMappingData]->IncCountingStatistics(aCountingStatistics);
	}
}

// ------------------------------------------------------------------------------

void TMeteorCounter::GetIntersectionCountArray(vector<int> & aIntersectionCountArray, const int aReferenceCount) {

//	    initialize array with size = 2^N (N is the number of used references (main and excluded)
//	    2^N because intersections are represented in bit combinations

	aIntersectionCountArray.assign((size_t)pow(2, aReferenceCount), 0); // regular cast is safe with numeric type

	for (size_t i = 0; i < FMappingDataList.size(); i++)
		FMappingDataList[i]->IncIntersectionCount(aIntersectionCountArray); //// FMappingDataList[i] must be previously allocated with new
}

// ------------------------------------------------------------------------------

int TMeteorCounter::LoadReferenceAnnotationDatabase(const string & aReferenceIniFileName)
{
  TReferenceIniFile aReferenceIniFile(aReferenceIniFileName);

  int Result = 0;

//  switch(aReferenceIniFile.GetDatabaseType()) { //// TODO
//    case mdbADT:
//    	  //// FReferenceAnnotationDatabase is actually TReferenceAnnotationDatabase (see later if a pointer is needed => polymorphisme)
//        TAdsReferenceAnnotationDatabase FReferenceAnnotationDatabase;
        if (aReferenceIniFile.GetHasLiteInfo())
        {
          FReferenceAnnotationDatabase.OpenLiteReferenceAnnotationTable(
        		  ExtractFileDir(aReferenceIniFileName) + C_PATH_SEP + aReferenceIniFile.GetReferenceDatabaseDirectory()
				  + C_PATH_SEP + aReferenceIniFile.GetReferenceName() + '_' + C_REFERENCE_LITE_ANNOTATION_TABLE_NAME);
          Result = FReferenceAnnotationDatabase.GetLiteReferenceAnnotationTable().GetFragmentCount();
        }
        else {
        	cerr << "Error in TMeteorCounter::LoadReferenceAnnotationDatabase: lite_annotation is missing" << endl;
        	exit(1);
        }
//        break;
//    case mdbParstream:
//    	 break;
//    case mdbGigabase:
//    	 break;
//    case mdbBinary:
//    	 break;
//    case mdbUnknown:
//      ;
//  }
  return Result;
}

// ------------------------------------------------------------------------------

void TMeteorCounter::ProcessCounting(
		const string & aCountingPrefixName,
		const string & aNGSMappingDirectory,
		const string & aReferenceIniFileName,
		const string & aProfileIniFileName,
		const TStrings & aLibraryCensusIniFileNameList,
		const TStrings & aExcludedReferenceNameList,
		const TMappingReferenceProperties & aMainMappingReferenceProperties,
		const bool aOKToCacheData, const TMeteorDBType & aMeteorDBType)
{
  TNGSCensusData* aNGSCensusTable = NULL; // base class
  string aNGSCensusTableName;
  string aCensusDatabaseDirectory;

  // declared in u_MeteorConstant
  TCountingStatistics aCountingStatistics;
  TCountingStatistics aLibraryCountingStatistics;

  vector<int> aExcludedReferenceReadCountArray;
  vector<int> aLibraryExcludedReferenceReadCountArray;

  vector<int> aIntersectionCountArray;
  vector<int> aLibraryIntersectionCountArray;

  int aSequencedReadCount, aSequencedRawReadCount, aSequencedHostReadCount, aIndexedReadCount;
  int64 aIndexedBaseCount;
  string aMappingTool;

//  TDateTime t1, t2;

  int aTotalFragmentCount;

  FCountingPrefixName = aCountingPrefixName;

  // Census data preparation

  // aCensusDatabaseDirectory = aNGSMappingDirectory + C_PATH_SEP + aCountingPrefixName;

  aCensusDatabaseDirectory = ExtractFileDir(aProfileIniFileName) + C_PATH_SEP + aCountingPrefixName;
  if ( access(aCensusDatabaseDirectory.c_str(), F_OK) != 0 ) {// test if file or directory exists
	  if ( my_mkdir(aCensusDatabaseDirectory.c_str()) != 0) {
		  cerr << "Error, could not create "<< aCensusDatabaseDirectory << endl;
		  exit(1);
	  }
  }

  // Prepare Census table
  aNGSCensusTableName = C_CENSUS_RESULT_TABLE_NAME;

  switch(aMeteorDBType) {
    case mdbADT:
      break;
    case mdbParstream:
      break;
    case mdbGigabase:
      break;
    case mdbBinary:
    	aNGSCensusTable = new TNGSCensusFile(aCensusDatabaseDirectory, aNGSCensusTableName); // TNGSCensusFile is a TNGSCensusData child class
        break;
    case mdbUnknown: ;
  }
  aTotalFragmentCount = LoadReferenceAnnotationDatabase(aReferenceIniFileName);

  TCountedFragmentList FCountedFragmentList(aTotalFragmentCount);

  // PASS 1
  for (size_t i = 0; i < FMappingDataList.size(); i++) {// FMappingDataList.size() corresponds to the number of libraries
    FMappingDataList[i]->CountFragment(FCountedFragmentList, FReferenceAnnotationDatabase); //// TMappingData. Main mapping file (SAM) is read here
  }

  // PASS 2 for updating smart shared counting
  for (size_t i = 0; i < FMappingDataList.size(); i++)
  {
    FMappingDataList[i]->UpdateSmartSharedCount(FCountedFragmentList, FReferenceAnnotationDatabase); //// TMappingData. Main mapping file (SAM) is also read here
  }
  FCountedFragmentList.Finalize();

  aNGSCensusTable->WriteCountedFragmentList(FCountedFragmentList, FCountedReadTypeList, 't');

  // Calculate statistics on read counting according to quality and excluded reference
  GetCountingStatistics(aCountingStatistics);
  // Retrieve read counts for excluded reference
  aExcludedReferenceReadCountArray.resize(aExcludedReferenceNameList.size());
  for (size_t iMappingData = 0; iMappingData < FMappingDataList.size(); iMappingData++)
    for (size_t i = 0 ; i < aExcludedReferenceNameList.size(); i++)
      aExcludedReferenceReadCountArray[i] += FMappingDataList[iMappingData]->GetExcludedReadCount(i);

  // Create new profile_ini_file and update
  TProfileIniFile aProfileIniFile(aProfileIniFileName);

  // initialization with info contained in the first library info file
  TCensusIniFile aLibraryCensusIniFile(C_EXTENDED_PREFIX_PATHNAME + aLibraryCensusIniFileNameList[0], true);

	aProfileIniFile.SetSampleName(aLibraryCensusIniFile.GetSampleName());
	aProfileIniFile.SetConditionName(aLibraryCensusIniFile.GetConditionName());
	aProfileIniFile.SetProjectName(aLibraryCensusIniFile.GetProjectName());
	aProfileIniFile.SetDatabaseType(aMeteorDBType);
	aProfileIniFile.SetProfileDirectory(FCountingPrefixName);
	aProfileIniFile.SetProfileDate(str_now()); ////
	aProfileIniFile.SetReferenceName(aMainMappingReferenceProperties.ReferenceName);
	aProfileIniFile.SetReadCleaningMethod(aLibraryCensusIniFile.GetReadCleaningMethod());

    aMappingTool = aLibraryCensusIniFile.GetMappingTool();

    ////FreeAndNil(aLibraryCensusIniFile);

    aSequencedRawReadCount = 0;
    aSequencedHostReadCount = 0;
    aSequencedReadCount = 0;
    aIndexedReadCount = 0;
    aIndexedBaseCount = 0;

    for (size_t i = 0; i < aLibraryCensusIniFileNameList.size(); i++)
    {
      TCensusIniFile aLibraryCensusIniFile(C_EXTENDED_PREFIX_PATHNAME + aLibraryCensusIniFileNameList[i], true);
      aSequencedRawReadCount += aLibraryCensusIniFile.GetSequencedRawReadCount();
      aSequencedHostReadCount += aLibraryCensusIniFile.GetSequencedHostReadCount();
      aSequencedReadCount += aLibraryCensusIniFile.GetSequencedReadCount();
      aIndexedReadCount += aLibraryCensusIniFile.GetIndexedSequencedReadCount();
      aIndexedBaseCount += aLibraryCensusIniFile.GetIndexedSequenceBaseCount();

      aProfileIniFile.AddProfiledSample(aLibraryCensusIniFile.GetSampleName(),
                        aLibraryCensusIniFile.GetFullSampleName(),
                        aLibraryCensusIniFile.GetConditionName(),
                        aLibraryCensusIniFile.GetSequencingDate());

    	aLibraryCountingStatistics.RawSequencedReads = aLibraryCensusIniFile.GetSequencedRawReadCount();
    	aLibraryCountingStatistics.HostSequencedReads = aLibraryCensusIniFile.GetSequencedHostReadCount();
    	aLibraryCountingStatistics.SequencedReads = aLibraryCensusIniFile.GetSequencedReadCount();
    	aLibraryCountingStatistics.IndexedReads = aLibraryCensusIniFile.GetIndexedSequencedReadCount();
    	aLibraryCountingStatistics.TotalCleanCountedReads = 0;
    	aLibraryCountingStatistics.UniqueCleanCountedReads = 0;
    	aLibraryCountingStatistics.RejectedReads = 0;
    	aLibraryCountingStatistics.RejectedReferenceReads = 0;
    	aLibraryCountingStatistics.RejectedQualityReads = 0;
    	aLibraryCountingStatistics.NotMappedReads = 0;
    	aLibraryCountingStatistics.NotMappedCleanReads = 0;

      FMappingDataList[i]->IncCountingStatistics(aLibraryCountingStatistics); //// as TMappingData

      aProfileIniFile.SetCountingStatistics(i, aLibraryCountingStatistics);

      // Retrieve read counts for excluded reference
      aLibraryExcludedReferenceReadCountArray.resize(aExcludedReferenceNameList.size());
      for (size_t iReference = 0 ; iReference < aExcludedReferenceNameList.size();  iReference++)
        aLibraryExcludedReferenceReadCountArray[iReference] = FMappingDataList[i]->GetExcludedReadCount(iReference); //// as TMappingData

      aProfileIniFile.SetExcludedReferenceCountingStatistics(i, aExcludedReferenceNameList, aLibraryExcludedReferenceReadCountArray);

//      // Count intersection
//        initialize array with size = 2^N (N is the number of used references (main and excluded)
//        2^N because intersections are represented in bit combinations

      ////aIntersectionCountArray = nil;
      aIntersectionCountArray.resize( (int)pow(2, aExcludedReferenceNameList.size()+1) );
      // +1 with main
      for ( size_t iIntersection = 0; iIntersection < aIntersectionCountArray.size(); iIntersection++)
        aIntersectionCountArray[iIntersection] = 0;

      FMappingDataList[i]->IncIntersectionCount(aIntersectionCountArray); //// as TMappingData

      for ( size_t iIntersection = 0; iIntersection < aIntersectionCountArray.size(); iIntersection++)
    	  aProfileIniFile.SetIntersectCount(i, iIntersection, aIntersectionCountArray[iIntersection]);
    }

    aProfileIniFile.SetSequencedRawReadCount(aSequencedRawReadCount);
    aProfileIniFile.SetSequencedHostReadCount(aSequencedHostReadCount);
    aProfileIniFile.SetSequencedReadCount(aSequencedReadCount);
    aProfileIniFile.SetIndexedSequencedReadCount(aIndexedReadCount);
    aProfileIniFile.SetIndexedSequencedBaseCount(aIndexedBaseCount);

    aCountingStatistics.RawSequencedReads = aSequencedRawReadCount;
    aCountingStatistics.HostSequencedReads = aSequencedHostReadCount;
    aCountingStatistics.SequencedReads = aSequencedReadCount;
    aCountingStatistics.IndexedReads = aIndexedReadCount;

    aProfileIniFile.SetMappedReadLength(aMainMappingReferenceProperties.MappedReadLength); ////
    aProfileIniFile.SetMappedReadLengthType(aMainMappingReferenceProperties.MappedReadLengthType);
    aProfileIniFile.SetMappingToolName(aMappingTool);
    aProfileIniFile.SetMappingCmdLine(aMainMappingReferenceProperties.MapperCmd);
    aProfileIniFile.SetIsMappingMismatchesPercentage(aMainMappingReferenceProperties.IsMappingMismatchesPercentage);
    aProfileIniFile.SetMappingMismatchesCount(aMainMappingReferenceProperties.MappingMismatches);
    aProfileIniFile.SetMappingMatchesCount(aMainMappingReferenceProperties.MappingMatches);
    aProfileIniFile.SetIsBestAlignment(aMainMappingReferenceProperties.IsBestAlignment);

    aProfileIniFile.SetIsLocalMapping(aMainMappingReferenceProperties.IsLocalAlignment);
    if (aMainMappingReferenceProperties.IsLocalAlignment)
    {
      aProfileIniFile.SetKeepInternalLocalAlignment(aMainMappingReferenceProperties.KeepInternalLocalAlignment);
      aProfileIniFile.SetAlignmentLengthCutoff(aMainMappingReferenceProperties.AlignmentLengthCutoff);
      aProfileIniFile.SetSoftClippingLengthCutoff(aMainMappingReferenceProperties.SoftClippingLengthCutoff);
    }

    aProfileIniFile.SetCountingStatistics(aCountingStatistics);
    aProfileIniFile.SetExcludedReferenceCountingStatistics(aExcludedReferenceNameList, aExcludedReferenceReadCountArray);

    // Add intersection count info in census_stage ini file
    GetIntersectionCountArray(aIntersectionCountArray, 1 + aExcludedReferenceNameList.size());
    // 1+ because main reference too

    for ( size_t i = 0; i < aIntersectionCountArray.size(); i++)
    	aProfileIniFile.SetIntersectCount(i, aIntersectionCountArray[i]);
    //// for (i = Low(aIntersectionCountArray) to High(aIntersectionCountArray))

  if (aNGSCensusTable != NULL) delete aNGSCensusTable; //// free memory here ????

   // Important:
   // Create SHA256 hash file after closing all files (by FreeAndNil)

//
//  // Create SHA256/MD5 hash for profile ini file
//  switch(DigestAlgorithm) {
//    case daMD5:
//      CreateMD5File(aProfileIniFileName);
//    case daSHA256:
//      CreateSHA256File(aProfileIniFileName);
//    case daNone:
//      ;
//  }
//  // Create SHA256/MD5 hash for files in census_result directory
//  switch(DigestAlgorithm) {
//    case daMD5:
//      CreateMD5FileFromDirectory(aCensusDatabaseDirectory);
//    case daSHA256:
//      CreateSHA256FileFromDirectory(aCensusDatabaseDirectory);
//    case daNone:
//      ;
//  }
}

