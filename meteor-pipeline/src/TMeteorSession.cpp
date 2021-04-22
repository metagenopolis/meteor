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

#include "TMeteorCounter.h"
#include "TMeteorSession.h"
#include "TCensusIniFile.h"
#include "TProfileIniFile.h"
#include "utils.h"
//#include <algorithm>

using namespace std;

TMeteorSession::TMeteorSession() {FLibraryCount=0;}

TMeteorSession::~TMeteorSession() {}

// ------------------------------------------------------------------------------

bool TMeteorSession::IsProfileDone()
{
  string aRepositoryProfileIniFileName;
  TMappingReferenceProperties aMappingReferenceProperties;

  bool Result = false;

  if (FLibraryCount > 0)
  {
    TCensusIniFile aLibraryIniFile(FRepositoryLibraryIniFileNameList.at(0), true); //// AT
//    aLibraryIniFile.SetReadOnly(true); // done in constructor

    FMeteorJobIniFile.GetMainReferenceProperties(aMappingReferenceProperties);

    aRepositoryProfileIniFileName = FRepositoryMappingDir + C_PATH_SEP + aLibraryIniFile.GetSampleName();
    // /projects/parkinson/mapping/KCL_01/KCL_01

    if (! aMappingReferenceProperties.CountingPrefixName.empty()) {
      aRepositoryProfileIniFileName += "_" + aMappingReferenceProperties.CountingPrefixName + C_DEFAULT_PROFILE_INIFILE_NAME_EXTENSION;
      // /projects/parkinson/mapping/KCL_01/KCL_01_vs_hs9_9_metahit_id95_rmHost_gene_profile.ini
    }
    else aRepositoryProfileIniFileName += C_DEFAULT_PROFILE_INIFILE_NAME_EXTENSION;
         // /project/parkinson/mapping/KCL_01/KCL_01_gene_profile.ini

    if (access(aRepositoryProfileIniFileName.c_str(), F_OK) == 0)
    {
      TProfileIniFile aProfileIniFile(aRepositoryProfileIniFileName, true); // readonly = true
      Result = aProfileIniFile.GetProfiledSampleCount() == FLibraryCount;
    }
  }
  return Result;
}

// ------------------------------------------------------------------------------

void TMeteorSession::PrepareWorkSpace()
{
	FWorkProjectDir = FRepositoryProjectDir;

	FWorkLibraryDir = FRepositorySampleDir; // /projects/vegan/sample/VG_000
	FWorkMappingDir = FRepositoryProjectDir + C_PATH_SEP + FMappingOutputDirBasename + replaceStr(FRepositorySampleDir, ExtractFileDir(FRepositorySampleDir), "");
	//                /projects/vegan       +     /      +     mapping       +    /VG_000
	// NB: ExtractFileDir("/projects/vegan/sample/VG_000") returns /projects/vegan/sample (ending path separator is not included in the resulting substring)

	FRepositoryMappingDir = FWorkMappingDir;

    // Create output directory in repository. e.g. /projects/parkinson/mapping/KCL_01
    int aDirCreateStatus = 0; //succes

	// ForceDirectories do nothing if path already exists.
	//aDirCreateStatus = ForceDirectories(C_EXTENDED_PREFIX_PATHNAME + FRepositoryProjectDir + C_PATH_SEP + FMappingOutputDirBasename + replaceStr(FRepositorySampleDir, ExtractFileDir(FRepositorySampleDir), ""));
	aDirCreateStatus = ForceDirectories(C_EXTENDED_PREFIX_PATHNAME + FWorkMappingDir);

	if (aDirCreateStatus != 0) exit(1);
}

// ------------------------------------------------------------------------------

void TMeteorSession::PrepareLibraryCensusIniFileNameList()
{
	const string C_CENSUS_STAGE_0_SUFFIX = "census_stage_0.ini";
	//const string C_CENSUS_STAGE_0_MASK = "*_" + C_CENSUS_STAGE_0_SUFFIX;

	TStrings aRepositoryLibraryIniFileNameArray;
	string aRepositoryLibraryIniFileName;

	FRepositoryLibraryIniFileNameList.clear();

	getFiles(FRepositorySampleDir, C_CENSUS_STAGE_0_SUFFIX, FRepositoryLibraryIniFileNameList);

	FLibraryCount = FRepositoryLibraryIniFileNameList.size();
}

// ------------------------------------------------------------------------------

void TMeteorSession::PrepareMainMappingCensusIniFileNameList()
{
  string aMappingDir, aDirectory;
  TMappingReferenceProperties aMappingReferenceProperties;

  string aRepositoryLibraryIniFileName, aRepositoryLibraryIniFileBaseName;
  FMainMappingCensusIniFileNameList.resize(FLibraryCount);
  ////FMainMappingDataCopyStatus.resize(FLibraryCount); //// TODO ???? or not

  FMeteorJobIniFile.GetMainReferenceProperties(aMappingReferenceProperties);

  for (int iLibrary = 0; iLibrary<FLibraryCount; iLibrary++)
  {
    aRepositoryLibraryIniFileName = FRepositoryLibraryIniFileNameList.at(iLibrary); //// AT
    TCensusIniFile aLibraryCensusIniFile(aRepositoryLibraryIniFileName, true);
//    aLibraryCensusIniFile.SetReadOnly(true); // done in constructor

    aRepositoryLibraryIniFileBaseName = BaseName(aRepositoryLibraryIniFileName);

    if (aMappingReferenceProperties.MappingPrefixName.empty())
    {
      aMappingDir = C_DEFAULT_SAMPLE_MAPPING_DIRECTORY + "_vs_" + aMappingReferenceProperties.ReferenceName + "_" + C_MAPPED_READ_LENGTH_CHAR +
    		  numberToString(aMappingReferenceProperties.MappedReadLength) + "-" + C_MISMATCHES_COUNT_CHAR +
			  numberToString(aMappingReferenceProperties.MappingMismatches) + "_" + aLibraryCensusIniFile.GetFullSampleName();
    }
    else aMappingDir = aMappingReferenceProperties.MappingPrefixName + "_" + aLibraryCensusIniFile.GetFullSampleName();

    // aMappingCensusIniFileNameWork is a reference (alias) to FMainMappingCensusIniFileNameList[iLibrary].FileLocalization[flWprk];
    TMappingCensusIniFileName & aMappingCensusIniFileNameWork = FMainMappingCensusIniFileNameList.at(iLibrary); //// AT

    aMappingCensusIniFileNameWork.Directory = aDirectory = FWorkProjectDir + C_PATH_SEP + FMappingOutputDirBasename + C_PATH_SEP + aLibraryCensusIniFile.GetSampleName() + C_PATH_SEP + aMappingDir;
    aMappingCensusIniFileNameWork.CensusStage1IniFileName = aDirectory + C_PATH_SEP + replaceStr(aRepositoryLibraryIniFileBaseName, C_CENSUS_STATUS_SEQUENCED_STR, C_CENSUS_STATUS_MAPPED_STR);

  }
}

// ------------------------------------------------------------------------------

void TMeteorSession::PrepareExcludedMappingCensusIniFileNameList()
{
  string aMappingDir, aDirectory;
  TMappingReferenceProperties aMappingReferenceProperties;
//  TCensusIniFile aLibraryCensusIniFile;
  string aRepositoryLibraryIniFileName, aRepositoryLibraryIniFileBaseName;

  //SetLength(FExcludedMappingCensusIniFileNameList, FLibraryCount, FMeteorJobIniFile.ExcludedReferenceCount);
  //// WARNING the 2Dvector resizing method below is safe only because FExcludedMappingCensusIniFileNameList is initially empty.
  //// => TODO: initialize this 2D vector in a new TMeteorSession constructor (with adequat parameters)
  FExcludedMappingCensusIniFileNameList.resize( FLibraryCount, vector<TMappingCensusIniFileName>(FMeteorJobIniFile.GetExcludedReferenceCount()) );

//  SetLength(FExcludedMappingDataCopyStatus, FLibraryCount, FMeteorJobIniFile.ExcludedReferenceCount);

  for (int iLibrary = 0; iLibrary < FLibraryCount; iLibrary++)
  {
	aRepositoryLibraryIniFileName = FRepositoryLibraryIniFileNameList.at(iLibrary);
    TCensusIniFile aLibraryCensusIniFile(aRepositoryLibraryIniFileName, true);
//    aLibraryCensusIniFile.SetReadOnly(true); // done in constructor

	aRepositoryLibraryIniFileBaseName = BaseName(aRepositoryLibraryIniFileName);

    for (int iExcludedReference = 0; iExcludedReference < FMeteorJobIniFile.GetExcludedReferenceCount(); iExcludedReference++)
    {
      FMeteorJobIniFile.GetExcludedReferenceProperties(iExcludedReference + 1, aMappingReferenceProperties);

      if (aMappingReferenceProperties.MappingPrefixName.empty())
      {
        aMappingDir = C_DEFAULT_SAMPLE_MAPPING_DIRECTORY + "_vs_" + aMappingReferenceProperties.ReferenceName + "_" + C_MAPPED_READ_LENGTH_CHAR + numberToString(aMappingReferenceProperties.MappedReadLength) + "-"
		+ C_MISMATCHES_COUNT_CHAR + numberToString(aMappingReferenceProperties.MappingMismatches) + "_" + aLibraryCensusIniFile.GetFullSampleName();
      }
      else aMappingDir = aMappingReferenceProperties.MappingPrefixName + "_" + aLibraryCensusIniFile.GetFullSampleName();

      TMappingCensusIniFileName & aMappingCensusIniFileNameWork = FExcludedMappingCensusIniFileNameList.at(iLibrary).at(iExcludedReference); //// AT

   	  aMappingCensusIniFileNameWork.Directory = aDirectory = FWorkProjectDir + C_PATH_SEP + FMappingOutputDirBasename + C_PATH_SEP + aLibraryCensusIniFile.GetSampleName() + C_PATH_SEP + aMappingDir;
   	  aMappingCensusIniFileNameWork.CensusStage1IniFileName = aDirectory + C_PATH_SEP + replaceStr(aRepositoryLibraryIniFileBaseName, C_CENSUS_STATUS_SEQUENCED_STR, C_CENSUS_STATUS_MAPPED_STR);
    }
  }
}

// ------------------------------------------------------------------------------

bool TMeteorSession::TaskMainCounting(const vector<TCountedReadType> & aCountedReadTypeList)
{
//  TMeteorCounter aMeteorCounter;
  TMappingReferenceProperties aMappingReferenceProperties;
  string aReferenceIniFileName;
  TStrings aMainMappingCensusIniFileNameList;

  TMappingReferencePropertiesArray aExcludedMappingReferencePropertiesArray;
  double aMaxDistance;

  string aProfileIniFileName;
  string aProfilePrefixDir;

  TStrings aExcludedReferenceNameList;
  bool Result = false;
  bool aOKToProcessTask = true;

  cout << "######### Task main counting" << endl;

  if (aOKToProcessTask)
  {
    TStrings aExcludedMappingCensusIniFileNameList;

    // get mapping and counting properties for main and excluded reference
    FMeteorJobIniFile.GetAllExcludedReferenceProperties(aExcludedMappingReferencePropertiesArray);
    FMeteorJobIniFile.GetMainReferenceProperties(aMappingReferenceProperties);

    aMaxDistance = (double)aMappingReferenceProperties.MappingMismatches;
    if (aMappingReferenceProperties.IsMappingMismatchesPercentage)
      aMaxDistance = aMaxDistance / 100.0; // relative distance

    TMeteorCounter aMeteorCounter(aCountedReadTypeList);

    aReferenceIniFileName = FMeteorJobIniFile.GetReferenceDir() + C_PATH_SEP + aMappingReferenceProperties.ReferenceName +
      C_PATH_SEP + aMappingReferenceProperties.ReferenceName + C_DEFAULT_REFERENCE_INIFILE_NAME_EXTENSION;

    TStrings aMainMappingCensusIniFileNameList;

    for (int iLibrary = 0; iLibrary<FLibraryCount; iLibrary++)
    {
      aMainMappingCensusIniFileNameList.push_back(FMainMappingCensusIniFileNameList.at(iLibrary).CensusStage1IniFileName); //// AT
      aExcludedMappingCensusIniFileNameList.clear();
      for (int iExcludedReference = 0; iExcludedReference<FMeteorJobIniFile.GetExcludedReferenceCount(); iExcludedReference++)
      {
        aExcludedMappingCensusIniFileNameList.push_back(C_EXTENDED_PREFIX_PATHNAME + FExcludedMappingCensusIniFileNameList.at(iLibrary).at(iExcludedReference).CensusStage1IniFileName); //// AT
      }
//      cerr <<"path: "<< C_EXTENDED_PREFIX_PATHNAME + FMainMappingCensusIniFileNameList.at(iLibrary).FileLocalization[flWork].CensusStage1IniFileName<<endl; //// AT
      
      // AddMappingLibrary(): set main mapping file name and excluded mapping file names (one SAM per excluded mapping file) but we read only the excluded mapping files.
      aMeteorCounter.AddMappingLibrary(
    		  FTmpDir,
    		  C_EXTENDED_PREFIX_PATHNAME + FMainMappingCensusIniFileNameList.at(iLibrary).CensusStage1IniFileName,
    		  aMaxDistance, aMappingReferenceProperties.IsMappingMismatchesPercentage,
		      aMappingReferenceProperties.IsLocalAlignment, /**** LOCAL*/
		      aMappingReferenceProperties.KeepInternalLocalAlignment,
		      aMappingReferenceProperties.AlignmentLengthCutoff,
			  aMappingReferenceProperties.SoftClippingLengthCutoff,
			  aExcludedMappingCensusIniFileNameList,
			  aExcludedMappingReferencePropertiesArray); //// AT
    }
    TStrings aExcludedReferenceNameList;
    FMeteorJobIniFile.GetExcludedReferenceNameList(aExcludedReferenceNameList);

    // Process counting
    if ( ! aMappingReferenceProperties.CountingPrefixName.empty() )
    {
      aProfileIniFileName = FWorkMappingDir + C_PATH_SEP + FSampleName + '_' + aMappingReferenceProperties.CountingPrefixName + C_DEFAULT_PROFILE_INIFILE_NAME_EXTENSION;
      aProfilePrefixDir = FSampleName + '_' + aMappingReferenceProperties.CountingPrefixName + '_' + C_DEFAULT_PROFILE_DIR;
    }
    else
    {
      aProfileIniFileName = FWorkMappingDir + C_PATH_SEP + FSampleName + C_DEFAULT_PROFILE_INIFILE_NAME_EXTENSION;
      aProfilePrefixDir = FSampleName + '_' + C_DEFAULT_PROFILE_DIR;
    }

    // Here we read the main mapping file (SAM) twice //// TODO try reading only once
    aMeteorCounter.ProcessCounting(
    		aProfilePrefixDir, C_EXTENDED_PREFIX_PATHNAME + FWorkMappingDir,
			aReferenceIniFileName, aProfileIniFileName,
			aMainMappingCensusIniFileNameList, aExcludedReferenceNameList,
			aMappingReferenceProperties, FMeteorJobIniFile.GetIsDataCaching(),
			FMeteorJobIniFile.GetMeteorDBType());

//    FreeAndNil(aMeteorCounter);
//    FreeAndNil(aExcludedMappingCensusIniFileNameList);
//    FreeAndNil(aMainMappingCensusIniFileNameList);
//    FreeAndNil(aExcludedReferenceNameList);

    Result = true;
  }
  else
  {
    cerr << "Unable to process this task because error in file copy" << endl;
    Result = false;
  }
  return Result;
}

// ------------------------------------------------------------------------------

// -c aMeteorJobIniFileName -i aSampleDirectory -p aProjectDirectory -o aMappingDirectory -m => true -l => true
void TMeteorSession::ProcessJobOnDir(
		const string & aTmpDir,
		const string & aMeteorJobIniFileName,
		const string & aSampleDirectory,
		const string & aProjectDirectory,
		const string & aMappingDirBasename,
		const vector<TCountedReadType> & aCountedReadTypeList,
		const bool aOverwrite)
{
//  int iExcludedReference;
  bool aOKToContinue=true;

  FMeteorJobIniFile.loadIniData(aMeteorJobIniFileName); // FMeteorJobIniFile was implicitely initialized via TMeteorSession default constructor.
  FRepositoryProjectDir = aProjectDirectory;
  FRepositorySampleDir = aSampleDirectory;
  FRepositoryLibraryIniFileName = "";
  FMappingOutputDirBasename = aMappingDirBasename;
  FTmpDir = aTmpDir;

  //if (ConnectRepository)
  if (true)
  {
    // Check if global counting is done on these files
//    aJobStartTime := Now;
//    writeln('Job started on ', FormatDateTime('c', aJobStartTime));
    cout << "######### workspace preparation" << endl;

    PrepareWorkSpace();
    PrepareLibraryCensusIniFileNameList();
    PrepareMainMappingCensusIniFileNameList();
    PrepareExcludedMappingCensusIniFileNameList();

    // Check if profile is already built

    if (FLibraryCount > 0)
    {
      // double test for generate intermediate library profile (after check)
      if (!IsProfileDone() || aOverwrite)//|| aOKToGenerateLibraryProfile)
      {
        for (int iLibrary = 0; iLibrary < FLibraryCount; iLibrary++)
        {
          // Test file path for project dir and library file (access with F_OK returns 0 if file or directory exists)
          if ( access(FRepositoryLibraryIniFileNameList.at(iLibrary).c_str(), F_OK) == 0 && access(FRepositoryProjectDir.c_str(), F_OK) == 0) //// 1st test should imply 2nd one ????
          {
            TCensusIniFile aLibraryIniFile(FRepositoryLibraryIniFileNameList.at(iLibrary), true);
//            aLibraryIniFile.SetReadOnly(true); // done in constructor

            cout << "######### Meteor task description" << endl;
            cout << "Sample name = " + aLibraryIniFile.GetSampleName() << endl;
            cout << "Library name = " + aLibraryIniFile.GetFullSampleName() << endl;
            cout << "Project name = " + aLibraryIniFile.GetProjectName() << endl;
            cout << "Sequencing device = " + aLibraryIniFile.GetSequencingDevice() << endl;
            cout << "Workflow = " + BaseName(aMeteorJobIniFileName) << endl;

            FSampleName = aLibraryIniFile.GetSampleName();
            FProjectName = aLibraryIniFile.GetProjectName();
          }
          else {
        	  cerr << "Warning, skipping not found "+FRepositoryLibraryIniFileNameList.at(iLibrary) << endl;
        	  aOKToContinue = false; //// are census_stage_0.ini files really needed for the MainCountingTask ????
        	  continue;
          }
        }
        // Now counting data
        if (aOKToContinue){
          //if (!IsProfileDone()) {
        	  aOKToContinue = TaskMainCounting(aCountedReadTypeList);
        	  if (!aOKToContinue) { cerr << "Error, TaskMainCounting failed ..." << endl; exit(1); }
          //}
          //else cout << "Profile already done" << endl;
        }
      }
      else cout << "Skipping already done profiling (counting) for sample: " << FWorkMappingDir << endl;
    }
    else { cerr << "Meteor job aborted because no library in the current directory" << endl; exit(1); }

    cout << "\nJob finished !" << endl;

  }
}
