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
#include "TMeteorSession.h"
#include "TMeteorProfiler.h"
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

int MeteorProfilerSession(int argc, char**argv){
	try {
		TCLAP::CmdLine cmd("Command example:\nmeteor-profiler -p /path/to/project -f /path/to/profiles -t smart_shared_reads -w /path/to/workflow.ini -o output_prefix /path/to/input_files.txt", ' ', "v3.2");
		//cmd.ignoreUnmatched(true); // Tell the command line to ignore any unmatched args.

		TCLAP::SwitchArg noMatrixSwitch("M","no-matrix","Do not create the matrix gathering all sample counting results");
		cmd.add( noMatrixSwitch );
		
		// ValueArg object: option expecting a value (--option value)
		// -w, --workflow, description, required?, default_value, value_type
		TCLAP::ValueArg<string> workflowArg("w","workflow","Path to the meteor configuration ini file",true,"","string");
		cmd.add( workflowArg );

		string countTypeDesc("");
//		countTypeDesc += "type of counting. A comma separated list of type among:\n";
		countTypeDesc += "type of counting.\npossible value (where direct/undirect means direct/undirect strand):\n";
		countTypeDesc += "    "+C_TYPE_READ_COUNT+", "+C_TYPE_SHARED_READ_COUNT+", "+C_TYPE_SMART_SHARED_READ_COUNT+", "+C_TYPE_UNIQUE_READ_COUNT+"\n";
		countTypeDesc += "    "+C_TYPE_DIRECT_READ_COUNT+", "+C_TYPE_DIRECT_SHARED_READ_COUNT+", "+C_TYPE_DIRECT_SMART_SHARED_READ_COUNT+", "+C_TYPE_DIRECT_UNIQUE_READ_COUNT+"\n";
//		countTypeDesc += "\n";
		countTypeDesc += "    "+C_TYPE_READ_COVERAGE+", "+C_TYPE_UNIQUE_READ_COVERAGE+", "+C_TYPE_DIRECT_READ_COVERAGE+", "+C_TYPE_DIRECT_UNIQUE_READ_COVERAGE+", "+C_TYPE_UNDIRECT_READ_COVERAGE+", "+C_TYPE_UNDIRECT_UNIQUE_READ_COVERAGE+"\n";
//		countTypeDesc += "\n";
		countTypeDesc += "    "+C_TYPE_READ_MEAN_COVERAGE+", "+C_TYPE_SHARED_READ_MEAN_COVERAGE+", "+C_TYPE_SMART_SHARED_READ_MEAN_COVERAGE+", "+C_TYPE_UNIQUE_READ_MEAN_COVERAGE+"\n";
		countTypeDesc += "    "+C_TYPE_DIRECT_READ_MEAN_COVERAGE+", "+C_TYPE_DIRECT_SHARED_READ_MEAN_COVERAGE+", "+C_TYPE_DIRECT_SMART_SHARED_READ_MEAN_COVERAGE+", "+C_TYPE_DIRECT_UNIQUE_READ_MEAN_COVERAGE+"\n";
		countTypeDesc += "    "+C_TYPE_UNDIRECT_READ_MEAN_COVERAGE+", "+C_TYPE_UNDIRECT_SHARED_READ_MEAN_COVERAGE+", "+C_TYPE_UNDIRECT_SMART_SHARED_READ_MEAN_COVERAGE+", "+C_TYPE_UNDIRECT_UNIQUE_READ_MEAN_COVERAGE+"\n";
		countTypeDesc += "    "+C_TYPE_UNDIRECT_READ_COUNT+", "+C_TYPE_UNDIRECT_SHARED_READ_COUNT+", "+C_TYPE_UNDIRECT_SMART_SHARED_READ_COUNT+", "+C_TYPE_UNDIRECT_UNIQUE_READ_COUNT+"\n";

		TCLAP::ValueArg<string> typeArg("t","type", countTypeDesc,true,"","string");
		cmd.add( typeArg );
		TCLAP::ValueArg<string> outputPrefixArg("o","output-prefix","output files prefix (stored in the profile directory, see option -f",true,"","string");
		cmd.add( outputPrefixArg );
		TCLAP::ValueArg<string>profileArg("f","profile-path","Path to profile directory, containing the output data",true,"","string");
		cmd.add( profileArg );
		TCLAP::ValueArg<string>projectArg("p","project-path","Path to project directory",true,"","string");
		cmd.add( projectArg );

		// input files. UnlabeledMultiArg HAS TO BE THE LAST ARGUMENT ADDED to the object cmd.
		string inputDesc("");
		inputDesc += "The list of input files (paths of files census.dat created by meteor-counter), or a unique text file ending with .txt countaining this list.\n";
		inputDesc += "    On unix-like OS, you can easily create this list in a text file with the command find:\n";
		inputDesc += "    find <starting_point_path> -name census.dat > input_files.txt";
		TCLAP::UnlabeledMultiArg<string> inputArg("input-file-list",inputDesc,true,"","string");
		cmd.add( inputArg );

		// parse argv
		cmd.parse( argc, argv );

		////TODO: check counting type
//		CheckProfileOutputType();

		// Get the value parsed by each arg.
		string aMeteorJobIniFileName = workflowArg.getValue();
		string aProjectDirectory     = projectArg.getValue();
		string aCountedReadTypeStr   = typeArg.getValue();
		string aOutputTableName      = outputPrefixArg.getValue();
		string aProfileOutputDir     = profileArg.getValue();
		vector<string> aInputData    = inputArg.getValue();
		bool aOkToCreateMatrix       = !noMatrixSwitch.getValue();

		// Check the wanted counting type
		TCountedReadType aCountedReadType = CheckProfileOutputType(aCountedReadTypeStr);
		if (aCountedReadType == crtUnknownReads){
			cerr << "Error, unrecognized counting type: " << aCountedReadTypeStr << endl;
			exit(1);
		}
		// prepare input files list
		TMeteorProfiler aMeteorProfiler(aMeteorJobIniFileName, aProjectDirectory, aProfileOutputDir, aCountedReadType, aInputData, aOkToCreateMatrix);
		aMeteorProfiler.CreateProfile(aOutputTableName);
		//aMeteorProfiler.CreateProfile(aOutputTableName, aCountedReadType);
	}
	catch (TCLAP::ArgException &e) { cerr << "error: " << e.error() << " for arg " << e.argId() << endl; exit(1); }

	return 0;
}

int main(int argc, char**argv){

	MeteorProfilerSession(argc, argv);
}
