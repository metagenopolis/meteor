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

#ifndef TMETEORPROFILER_H_
#define TMETEORPROFILER_H_

#include "utils.h"
#include "MeteorConstant.h"
#include "TMeteorJobIniFile.h"
#include "TCensusIniFile.h"
#include "TProfileIniFile.h"
#include "TProfileTable.h"
#include "TReferenceAnnotationDatabase.h"

//using namespace std;

const std::string C_GENE_PROFILE_SUFFIX = "_gene_profile";
const std::string C_CENSUS_STAGE_1_SUFFIX = "_census_stage_1.ini";

const std::string C_INIFILENAME = "inifilename";
const std::string C_COUNTING_PARAMETERS_ID = "id_counting_parameters";
const std::string C_COUNTED_SAMPLE_ID = "id_counted_sample";
const std::string IDX_COUNTING = "idx_counting";

typedef struct {
	std::string SampleName;
	std::string ConditionName;
	std::string SampleRunName;
	std::string SequencingDate;
	std::string ReferenceName;
	std::string MappingTool;
	int Mismatches;
	int Matches;
	int MappedReadLength;
	int SequencedReadCount;
	int HostReadCount;
	int IndexedReadCount;
	int MappedReadCount;
	int UniqueMappedReadCount;
	int CountedReadCount;
	int UniqueCountedReadCount;
	int RejectedReadCount;
	int ExcludedReferenceReadCount;
	int NotAssignedReadCount;
	int NotAssignedCleadReadCount;
	std::vector<int> ExcludedReferenceReadCountList;
} TCountingStatisticsItem;

typedef struct {
	std::string ResultFileName;
	std::string IniFileName;
} TCountingResult;

typedef std::vector< std::vector<double> > TCountedSamplesTable;

class TMeteorProfiler {

  private:

    TMeteorJobIniFile FMeteorJobIniFile;
    std::string FProjectDirectory;
    std::string FProfileDirectory;

    TReferenceAnnotationDatabase FReferenceAnnotationDatabase;

    TStrings FExcludedReferenceNameList;
    TMappingReferenceProperties FMappingMainReferenceProperties;

    // for comparing excluded reference for current mapping with job config
    TStrings FCurrentExcludedReferenceNameList;

    // stores each sample counting result (census.dat and corresponding ini file)
    std::vector<std::string> FSampleProfileFileNames;
    std::vector<std::string> FSampleProfileIniFileNames;
    std::vector<std::string> FSampleNames;

    // stores each sample statistics
    std::vector<TCountingStatisticsItem> FCountingStatisticsArray;

    // stores each sample extended statistics. Contains counting information for every run
    std::vector<TCountingStatisticsItem> FExtendedCountingStatisticsArray;

    TCountedReadType FCountedReadType;

    int FRunSampleCount;
    int FSampleCount;
    bool FIsHostRemoved; // true if host reads have been removed from sequenced data, and .
//    bool FIsUsualSampleName;
    bool FOkToCreateMatrix;

    void LoadCountingInfo(const std::string & aProfileIniFileName);
    bool CheckProfileIniFile(const TProfileIniFile & aProfileIniFile);
    void FillExtendedCountingStatisticsArray(const TProfileIniFile & aProfileIniFile);

  public:

	TMeteorProfiler(const std::string & aMeteorJobIniFileName, const std::string & aProjectDirectory, const std::string & aProfileDirectory, const TCountedReadType aCountedReadType, const std::vector<std::string> & aInputData, const bool aOkToCreateMatrix);
	virtual ~TMeteorProfiler();

    void CreateProfile(const std::string & aProfileName);
    void ExportCountingStatistics(const std::string & aFileName);
    void ExportExtendedCountingStatistics(const std::string & aFileName);

};

#endif /* TMETEORPROFILER_H_ */
