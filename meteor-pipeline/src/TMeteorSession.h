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

#ifndef TMETEORSESSION_H_
#define TMETEORSESSION_H_

#include "TMeteorJobIniFile.h"
#include "MeteorConstant.h"
#include <vector>
//#include "TNGSLibraryCleaner.h"

const std::string C_LOCK_FILE_EXTENSION = ".lock";

typedef struct {
  std::string Directory;
  std::string CensusStage1IniFileName;

} TMappingCensusIniFileName;

class TMeteorSession {
private:

  TMeteorJobIniFile FMeteorJobIniFile;

  // repository data (on the data server)
  std::string FRepositoryProjectDir;
  std::string FRepositorySampleDir;
  std::string FRepositoryMappingDir;
  std::string FRepositoryLibraryIniFileName;

  // list of all census_stage_0.ini in one sample directory
  // e.g. /projects/project_name/sample/sample1/*_census_stage_0.ini
  TStrings FRepositoryLibraryIniFileNameList;

//  bool FIsLocalCopy;
  // meteor working directory
  std::string FWorkProjectDir;
  std::string FWorkSampleDir;
  std::string FWorkMappingDir;
  std::string FWorkLibraryDir;

  int FLibraryCount;

  std::vector<TMappingCensusIniFileName> FMainMappingCensusIniFileNameList;
  std::vector< std::vector<TMappingCensusIniFileName> > FExcludedMappingCensusIniFileNameList;
  //  FExcludedMappingCensusIniFileNameList: array of array of TMappingCensusIniFileNameList;

  std::string FSampleName;
  std::string FProjectName;
  std::string FMappingOutputDirBasename;
  std::string FTmpDir; // path to the directory where temporary files (e.g. sam) are stored

  void PrepareWorkSpace();
  void PrepareLibraryCensusIniFileNameList();
  void PrepareMainMappingCensusIniFileNameList();
  void PrepareExcludedMappingCensusIniFileNameList();

  bool TaskMainCounting(const std::vector<TCountedReadType> & aCountedReadTypeList);

  //bool IsAllMappingDone(const int IndexLib) const;
  bool IsProfileDone();

public:
	TMeteorSession();
	virtual ~TMeteorSession();

//    property LicenseInfo: TLicenseInfo read FLicenseInfo;
    void ProcessJobOnDir(
			const std::string & aTmpDir,
    		const std::string & aMeteorJobIniFileName,
			const std::string & aSampleDirectory,
			const std::string & aProjectDirectory,
			const std::string & aMappingDirectory,
			const std::vector<TCountedReadType> & aCountedReadTypeList,
            const bool aOverwrite);

};




#endif /* TMETEORSESSION_H_ */
