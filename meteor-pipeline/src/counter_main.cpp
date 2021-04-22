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

#include "tclap/CmdLine.h" // argument parser
//#include "TSAMReader.h"
//#include "TSAMDefinition.h"
#include "TMeteorCounter.h"
#include "TMeteorSession.h"
#include "MeteorConstant.h"
#include "utils.h"
#include <algorithm>
#include <iostream>
#include <sstream>

using namespace std;

// ------------------------------------------------------------------------------

TCountedReadType CheckProfileOutputType(const string & aCountedReadTypeStr)
{
  if (aCountedReadTypeStr == C_TYPE_READ_COUNT)
    return  crtTotalReads;
  if (aCountedReadTypeStr == C_TYPE_SHARED_READ_COUNT)
    return crtSharedReads;
  if (aCountedReadTypeStr == C_TYPE_SMART_SHARED_READ_COUNT)
    return crtSmartSharedReads;
  if (aCountedReadTypeStr == C_TYPE_UNIQUE_READ_COUNT)
    return crtUniqueReads;
  if (aCountedReadTypeStr == C_TYPE_DIRECT_READ_COUNT)
    return crtDirectTotalReads;
  if (aCountedReadTypeStr == C_TYPE_DIRECT_SHARED_READ_COUNT)
    return crtDirectSharedReads;
  if (aCountedReadTypeStr == C_TYPE_DIRECT_SMART_SHARED_READ_COUNT)
    return crtDirectSmartSharedReads;
  if (aCountedReadTypeStr == C_TYPE_DIRECT_UNIQUE_READ_COUNT)
    return crtDirectUniqueReads;

  if (aCountedReadTypeStr == C_TYPE_READ_COVERAGE) // total_reads_coverage
    return crtTotalReadsCoverage;
  if (aCountedReadTypeStr == C_TYPE_UNIQUE_READ_COVERAGE) // unique_reads_coverage
    return crtUniqueReadsCoverage;
  if (aCountedReadTypeStr == C_TYPE_DIRECT_READ_COVERAGE) // direct_total_reads_coverage
	return crtDirectTotalReadsCoverage;
  if (aCountedReadTypeStr == C_TYPE_DIRECT_UNIQUE_READ_COVERAGE) // direct_unique_reads_coverage
	return crtDirectUniqueReadsCoverage;
  if (aCountedReadTypeStr == C_TYPE_UNDIRECT_READ_COVERAGE) // undirect_total_reads_coverage
	return crtUnDirectTotalReadsCoverage;
  if (aCountedReadTypeStr == C_TYPE_UNDIRECT_UNIQUE_READ_COVERAGE) // undirect_unique_reads_coverage
	return crtUnDirectUniqueReadsCoverage;

  if (aCountedReadTypeStr == C_TYPE_READ_MEAN_COVERAGE)
    return crtTotalReadsMeanCoverage;
  if (aCountedReadTypeStr == C_TYPE_SHARED_READ_MEAN_COVERAGE)
    return crtSharedReadsMeanCoverage;
  if (aCountedReadTypeStr == C_TYPE_SMART_SHARED_READ_MEAN_COVERAGE)
    return crtSmartSharedReadsMeanCoverage;
  if (aCountedReadTypeStr == C_TYPE_UNIQUE_READ_MEAN_COVERAGE)
    return crtUniqueReadsMeanCoverage;
  if (aCountedReadTypeStr == C_TYPE_DIRECT_READ_MEAN_COVERAGE)
    return crtDirectTotalReadsMeanCoverage;
  if (aCountedReadTypeStr == C_TYPE_DIRECT_SHARED_READ_MEAN_COVERAGE)
    return crtDirectSharedReadsMeanCoverage;
  if (aCountedReadTypeStr == C_TYPE_DIRECT_SMART_SHARED_READ_MEAN_COVERAGE)
    return crtDirectSmartSharedReadsMeanCoverage;
  if (aCountedReadTypeStr == C_TYPE_DIRECT_UNIQUE_READ_MEAN_COVERAGE)
    return crtDirectUniqueReadsMeanCoverage;
  if (aCountedReadTypeStr == C_TYPE_UNDIRECT_READ_MEAN_COVERAGE)
    return crtUndirectTotalReadsMeanCoverage;
  if (aCountedReadTypeStr == C_TYPE_UNDIRECT_SHARED_READ_MEAN_COVERAGE)
    return crtUndirectSharedReadsMeanCoverage;
  if (aCountedReadTypeStr == C_TYPE_UNDIRECT_SMART_SHARED_READ_MEAN_COVERAGE)
    return crtUndirectSmartSharedReadsMeanCoverage;
  if (aCountedReadTypeStr == C_TYPE_UNDIRECT_UNIQUE_READ_MEAN_COVERAGE)
    return crtUndirectUniqueReadsMeanCoverage;

  if (aCountedReadTypeStr == C_TYPE_UNDIRECT_READ_COUNT)
    return crtUndirectTotalReads;
  if (aCountedReadTypeStr == C_TYPE_UNDIRECT_SHARED_READ_COUNT)
    return crtUndirectSharedReads;
  if (aCountedReadTypeStr == C_TYPE_UNDIRECT_SMART_SHARED_READ_COUNT)
    return crtUndirectSmartSharedReads;
  if (aCountedReadTypeStr == C_TYPE_UNDIRECT_UNIQUE_READ_COUNT)
    return crtUndirectUniqueReads;

  // Error in command line
  return crtUnknownReads;
}

// ------------------------------------------------------------------------------

vector<TCountedReadType> CheckCountingTypes(const vector<string> & aCountingTypeList){

	vector<TCountedReadType> aCountedReadTypeList(aCountingTypeList.size());
	string badTypes("");

	for (size_t i = 0; i<aCountingTypeList.size(); i++){
		// Check the wanted counting type
		TCountedReadType aCountedReadType = CheckProfileOutputType(aCountingTypeList[i]);
		if ( aCountedReadType < crtUnknownReads )
			aCountedReadTypeList[i] = aCountedReadType;
		else { badTypes += aCountingTypeList[i] + '\n'; }
	}
	if (!badTypes.empty()) {
		// raise exception
		throw badTypes;
	}
	sort( aCountedReadTypeList.begin(), aCountedReadTypeList.end() );
	return aCountedReadTypeList;
}

// ------------------------------------------------------------------------------

int MeteorCounterSession(int argc, char**argv){

	//          -p projectDir       /projects/parkinson
	//          -w workflowPath     /path/to/workflow.ini
	//// TODO ?? : -m sampleMappingDir /projects/parkinson/mapping/KCL_01)
	try {
		TCLAP::CmdLine cmd("Command description message", ' ', "v3.2");

		// ValueArg object: option expecting a value (--option value)
		// -w, --workflow, description, required?, default_value, value_type
		TCLAP::ValueArg<string> workflowArg("w","workflow","Path to the meteor configuration ini file",true,"","string");
		cmd.add( workflowArg );

		string countTypeDesc("");

//		countTypeDesc += "type of counting. A comma separated list of type among:\n";
		countTypeDesc += "type of counting.\nPossible values: \""+C_TYPE_ALL_READ_COUNT+"\" (default) or a comma separated list of type among the following (where direct/undirect means direct/undirect strand):\n";
		countTypeDesc += "    "+C_TYPE_READ_COUNT+", "+C_TYPE_SHARED_READ_COUNT+", "+C_TYPE_SMART_SHARED_READ_COUNT+", "+C_TYPE_UNIQUE_READ_COUNT+"\n";
		countTypeDesc += "    "+C_TYPE_DIRECT_READ_COUNT+", "+C_TYPE_DIRECT_SHARED_READ_COUNT+", "+C_TYPE_DIRECT_SMART_SHARED_READ_COUNT+", "+C_TYPE_DIRECT_UNIQUE_READ_COUNT+"\n";
//		countTypeDesc += "\n";
		countTypeDesc += "    "+C_TYPE_READ_COVERAGE+", "+C_TYPE_UNIQUE_READ_COVERAGE+", "+C_TYPE_DIRECT_READ_COVERAGE+", "+C_TYPE_DIRECT_UNIQUE_READ_COVERAGE+", "+C_TYPE_UNDIRECT_READ_COVERAGE+", "+C_TYPE_UNDIRECT_UNIQUE_READ_COVERAGE+"\n";
//		countTypeDesc += "\n";
		countTypeDesc += "    "+C_TYPE_READ_MEAN_COVERAGE+", "+C_TYPE_SHARED_READ_MEAN_COVERAGE+", "+C_TYPE_SMART_SHARED_READ_MEAN_COVERAGE+", "+C_TYPE_UNIQUE_READ_MEAN_COVERAGE+"\n";
		countTypeDesc += "    "+C_TYPE_DIRECT_READ_MEAN_COVERAGE+", "+C_TYPE_DIRECT_SHARED_READ_MEAN_COVERAGE+", "+C_TYPE_DIRECT_SMART_SHARED_READ_MEAN_COVERAGE+", "+C_TYPE_DIRECT_UNIQUE_READ_MEAN_COVERAGE+"\n";
		countTypeDesc += "    "+C_TYPE_UNDIRECT_READ_MEAN_COVERAGE+", "+C_TYPE_UNDIRECT_SHARED_READ_MEAN_COVERAGE+", "+C_TYPE_UNDIRECT_SMART_SHARED_READ_MEAN_COVERAGE+", "+C_TYPE_UNDIRECT_UNIQUE_READ_MEAN_COVERAGE+"\n";
		countTypeDesc += "    "+C_TYPE_UNDIRECT_READ_COUNT+", "+C_TYPE_UNDIRECT_SHARED_READ_COUNT+", "+C_TYPE_UNDIRECT_SMART_SHARED_READ_COUNT+", "+C_TYPE_UNDIRECT_UNIQUE_READ_COUNT+"\n";

		TCLAP::ValueArg<string> typeArg("c","counting-type", countTypeDesc,false,"all","string");
		cmd.add( typeArg );

		TCLAP::ValueArg<string> inputArg("i","input","path to sample directory, containing the sample sequencing metadata (files ending with _census_stage_0.ini)",true,"","string");
		cmd.add( inputArg );
		TCLAP::ValueArg<string> projectArg("p","project","path to project directory, containing mapping and profile data (e.g. /projects/project_dir)",true,"","string");
		cmd.add( projectArg );
		TCLAP::ValueArg<string> mappingArg("o","mapping-basename","basename of the input mapping directory",true,"","string");
		cmd.add( mappingArg );
		TCLAP::SwitchArg forceOverwriteSwitch("f","force-overwrite","force overwriting former profiling results done with same parameters");
		cmd.add( forceOverwriteSwitch );
		TCLAP::ValueArg<string> tmpArg("t","tmp-path","path to the directory where temporary files (e.g. sam) are stored",false,"","string"); //// NEW
		cmd.add( tmpArg );

		// SwitchArg object: flag option without value (--flag)
		// -r, --reverse, option description, [referent CmdLine object,] default_value=false
		//TCLAP::SwitchArg flagSwitch("f","flag","flag description", [cmd,]);

		// parse argv
		cmd.parse( argc, argv );
		// Get the value parsed by each arg.
		string aMeteorJobIniFileName = workflowArg.getValue();
		string aCountingTypeListStr  = typeArg.getValue();
		string aInputData            = inputArg.getValue();
		string aProjectDirectory     = projectArg.getValue();
		string aMappingDirBasename   = mappingArg.getValue();
		bool aOverwrite              = forceOverwriteSwitch.getValue();
		string aTmpDir               = tmpArg.getValue(); // may be empty //// NEW

		vector<string> aCountingTypeList;
		vector<TCountedReadType> aCountedReadTypeList;

		if(!aCountingTypeListStr.empty() && aCountingTypeListStr != "all"){
			splitString(aCountingTypeList, aCountingTypeListStr, ',');

			try { aCountedReadTypeList = CheckCountingTypes(aCountingTypeList);	}
			catch(const string & e) {
				cerr << "Error, [-c] unrecognized counting type(s) : \n" << e << "See complete USAGE with option -h" << endl;
				exit(1);
			}
		}
//		exit(0);

		//bool aOKToGenerateLibraryProfile = false; //obsolete

		TMeteorSession aMeteorSession;
		
		//// TODO: get absolute path of aProjectDirectory returning a std::string
		//         on unix/linux: use char* realpath(const char *path, char *resolved_path); (stdlib.h + limits.h)
		//         on Windows :   use char* _fullpath(char *absPath, const char *relPath, size_t maxLength)
		//                        see https://docs.microsoft.com/en-us/cpp/c-runtime-library/reference/fullpath-wfullpath?view=msvc-160

		aMeteorSession.ProcessJobOnDir(aTmpDir, aMeteorJobIniFileName, aInputData, aProjectDirectory, aMappingDirBasename, aCountedReadTypeList, aOverwrite);
	}
	catch (TCLAP::ArgException &e) { cerr << "error: " << e.error() << " for arg " << e.argId() << endl; exit(1); }

	return 0;

}

// ------------------------------------------------------------------------------

int main(int argc, char**argv){

	MeteorCounterSession(argc, argv);
}
