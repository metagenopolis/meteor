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

#include <cstdlib>  // exit()
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "MeteorConstant.h"
#include "TProfileTable.h"

using namespace std;

//// to do in constructor TProfileTable. allocate and initialize private members a2DArray and aProgressionArray
TProfileTable::TProfileTable(const vector<string> & aFileNames, const TCountedReadType aCountedReadType){

	F2DArray = NULL;
	FMaxRowsNumber = 100000;
	FLineNumber = 0;
	FCountedReadType = aCountedReadType;
	string aCountedReadTypeStr(C_CENSUS_TABLE_FILE_HEADER_ARRAY[aCountedReadType]);
	FSampleCount = aFileNames.size();

	FInputFileInfoArray.resize(FSampleCount);

	// collect informations from each input file header (column index to consider, number of columns)
	string line;
	for(size_t iSample = 0; iSample < FSampleCount; iSample++){

		TInputFileInfo & aFileInfo = FInputFileInfoArray[iSample]; // ref
		string & aFileName = aFileInfo.FileName = aFileNames[iSample];
		aFileInfo.Progression = 0L;
		int & iWantedCol = aFileInfo.WantedColId = -1;
		int & iOptionalCol = aFileInfo.OptionalColId = 0;
		aFileInfo.NbOfColumns = 0;

		// open in read mode
		ifstream aFile(aFileName.c_str());
		/// open file failed
		if (!aFile) {
			cerr<<"Error in "<<__FILE__<<", line "<<__LINE__<<", cannot open file: "<<aFileName<<endl;
			exit(1);
		}

		getline(aFile, line);
		vector<string> aHeaderVect;
		splitString(aHeaderVect, line, '\t');

		int & NbOfColumns = aFileInfo.NbOfColumns = aHeaderVect.size() - 2;

		iWantedCol = find(aHeaderVect.begin(), aHeaderVect.end(), aCountedReadTypeStr) - aHeaderVect.begin() - 2;

		// counted type not found
		if(iWantedCol == NbOfColumns) {

			// case where we want one of the undirect counting type (crtUndirectTotalReads, crtUndirectSharedReads, crtUndirectSmartSharedReads, crtUndirectUniqueReads)
			// => let see if it can be infered ( undirect = total - direct )
			if (aCountedReadType == crtUndirectTotalReads) {
				iWantedCol = find(aHeaderVect.begin(), aHeaderVect.end(), C_CENSUS_TABLE_FILE_HEADER_ARRAY[crtTotalReads]) - aHeaderVect.begin() - 2;
				iOptionalCol = find(aHeaderVect.begin(), aHeaderVect.end(), C_CENSUS_TABLE_FILE_HEADER_ARRAY[crtDirectTotalReads]) - aHeaderVect.begin() - 2;
			}
			else if (aCountedReadType == crtUndirectSharedReads) {
				iWantedCol = find(aHeaderVect.begin(), aHeaderVect.end(), C_CENSUS_TABLE_FILE_HEADER_ARRAY[crtSharedReads]) - aHeaderVect.begin() - 2;
				iOptionalCol = find(aHeaderVect.begin(), aHeaderVect.end(), C_CENSUS_TABLE_FILE_HEADER_ARRAY[crtDirectSharedReads]) - aHeaderVect.begin() - 2;
			}
			else if (aCountedReadType == crtUndirectSmartSharedReads) {
				iWantedCol = find(aHeaderVect.begin(), aHeaderVect.end(), C_CENSUS_TABLE_FILE_HEADER_ARRAY[crtSmartSharedReads]) - aHeaderVect.begin() - 2;
				iOptionalCol = find(aHeaderVect.begin(), aHeaderVect.end(), C_CENSUS_TABLE_FILE_HEADER_ARRAY[crtDirectSmartSharedReads]) - aHeaderVect.begin() - 2;
			}
			else if (aCountedReadType == crtUndirectUniqueReads) {
				iWantedCol = find(aHeaderVect.begin(), aHeaderVect.end(), C_CENSUS_TABLE_FILE_HEADER_ARRAY[crtUniqueReads]) - aHeaderVect.begin() - 2;
				iOptionalCol = find(aHeaderVect.begin(), aHeaderVect.end(), C_CENSUS_TABLE_FILE_HEADER_ARRAY[crtDirectUniqueReads]) - aHeaderVect.begin() - 2;
			}
			// counted type missing and cannot be computed
			if (iWantedCol == NbOfColumns || iOptionalCol == NbOfColumns) {
				cerr << "Error, counting type " << aCountedReadTypeStr << " not found in file: " << aFileName << endl;
				exit(1);
			}
		}
		aFile.close();
	}

	// initialize the 2D buffer
	F2DArray = new double*[FMaxRowsNumber];       // allocate an array of pointers (each will point on a row, see the loop below)
	F2DArray[0] = new double[FMaxRowsNumber * FSampleCount]; // allocate the array containing the data (double)

	for (int i=1; i<FMaxRowsNumber; i++){                    // assign each pointer to the right address in the data array (pointed by a2DArray[0])
		F2DArray[i] = F2DArray[i-1] + FSampleCount;
	}
}

void TProfileTable::Reset2DArray(){ /* really usefull ?*/}
TProfileTable::~TProfileTable(){
	// release memory
	delete[] F2DArray[0];
	delete[] F2DArray;
}

int TProfileTable::GetMaxRowsNumber() const { return FMaxRowsNumber; }

size_t TProfileTable::GetSampleCount() const { return FSampleCount; }

int TProfileTable::GetLineNumber() const { return FLineNumber; }

/// TProfileTable::MergeColumns read the file(s) and returns the number of lines read:
//  exactly FMaxRowsNumber, except for the last block (0 < n <= FMaxRowsNumber).
//  finaly 0 when the end of the file(s) is reached.
int TProfileTable::MergeColumns(){

	int aLineNumber = 0;
	string line;
	string aValueStr;
	double aValue;

	// loop on each file
	for (size_t iSample = 0; iSample<FInputFileInfoArray.size(); iSample++){


		TInputFileInfo & aFileInfo = FInputFileInfoArray[iSample]; // reference
		string & aFileName = aFileInfo.FileName;
		long & aProgression = aFileInfo.Progression;
		int & iWantedCol = aFileInfo.WantedColId;
		int & iOptionalCol = aFileInfo.OptionalColId;
		int aLastColumnId = aFileInfo.NbOfColumns - 1;
		// nb of bytes to skip if the fragment is not counted (FragmentSize == -1): nb of columns x 2 characters ('0'+'\t')
		int aZeroBytesOffset = aFileInfo.NbOfColumns * 2;

		aLineNumber = 0;

		ifstream aFile(aFileName.c_str());       // read only stream on the file

		/// open file failed
		if (!aFile) {
			cerr<<"Error in "<<__FILE__<<", line "<<__LINE__<<", cannot open file: "<<aFileName<<endl;
			exit(1);
		}

		/// open file success

		// skip the header
		if (aProgression == 0) { getline(aFile, line); }

		// jump to the byte where we left off, during the previous reading of this file.
		else aFile.seekg(aProgression);

		// and continue reading the next FMaxRowsNumber lines.
		while (aLineNumber < FMaxRowsNumber){

			string aFragmentId, aFragmentSize;

			// read the 1st field (fragmentId),
			//if no success, we have reached the end of file (eof) => stop read
			if (!getline(aFile, aFragmentId, '\t')) break; // eof
			// else read the second field (fragmentSize)
			getline(aFile, aFragmentSize, '\t');

			// if fragment size == -1 => following fields are composed of zeros: "0\t ... 0\n"
			// => don't read the rest of the line and jump to next line (skip aZeroBytesOffset characters)
			if(aFragmentSize == "-1") {
				//if(aFragmentSize[0] == '-') { //// TO TEST
				// skip 52 characters from current position to jump to next line (without reading data ???)
				aFile.seekg( aZeroBytesOffset, aFile.cur );
				// place a "0" in F2DArray[aFragmentId][iSample]
				F2DArray[aLineNumber][iSample] = 0.0;
			}
			// else fragment is counted => evaluate the rest of the line
			else {
				// special case: the last column is wanted
				if (iWantedCol == aLastColumnId){
					if (!getline(aFile, line)) break; // eof
					// get the value after the last '\t'
					size_t tabPos = line.rfind('\t');
					//// TO TEST: line.find_last_of('\t');
					//// TO TEST: F2DArray[aLineNumber][iSample] = atof(&line[tabPos+1]);
					F2DArray[aLineNumber][iSample] = atof(line.substr(tabPos+1).c_str());
				}
				else {
					// read the rest of the line column by column till iWantedCol is reached
					int iValue = 0;
					while(iValue <= iWantedCol ) {
						if (!getline(aFile, aValueStr, '\t')) break; // eof
						iValue++;
					}
					if (aFile.eof()) break; // eof

					// column is reached and corresponding value is placed in F2DArray[aFragmentId][iSample]
					if(iOptionalCol == 0) {
						//~ if(aValueStr == "0") F2DArray[aLineNumber][iSample] = 0.0;
						//~ else F2DArray[aLineNumber][iSample] = atof(aValueStr.c_str());
						F2DArray[aLineNumber][iSample] = atof(aValueStr.c_str());
					}
					else{
						/// we need the data of another column to compute the final value (undirect type).
						// keep current one (total type) in tmp value
						aValue = atof(aValueStr.c_str());
						//  => move to the optional column (direct type)
						while (iValue <= iOptionalCol) { getline(aFile, aValueStr, '\t'); iValue++; }
						// place the difference in F2DArray[aFragmentId][iSample]
						F2DArray[aLineNumber][iSample] =  aValue - atof(aValueStr.c_str());
					}
					// read the rest of the line until '\n'
					if (!getline(aFile, line)) break; // eof
				}
			}
			//cout << aFragmentId << '\t' << aFragmentSize << '\t' << F2DArray[aLineNumber][iSample] << endl;
			aLineNumber++;

		} // END of while (aLineNumber < FMaxRowsNumber)

		aProgression = aFile.tellg();
		aFile.close();

	} // END of loop on input files

	FLineNumber += aLineNumber;

	return aLineNumber;
}

void TProfileTable::Write2DArray(int aReadLinesNb, ofstream & aOutputFile) const {

	// get the first Fragment id of the current block.
	int iFragment = FLineNumber - aReadLinesNb;

	// for each row (fragment) of this block
	for (int i = 0; i < aReadLinesNb; i++) {

		aOutputFile << ++iFragment;

		// write the value of each sample
		for (size_t j = 0; j < FSampleCount; j++){
			if(F2DArray[i][j] == 0.0) aOutputFile << "\t0";
			else aOutputFile << '\t' << F2DArray[i][j];
//			fmt::print(std::ofstream & file, "{}\n", doublePrecisionNum);
//			fmt::print(aOutputFile, "{:.8f}\n", F2DArray[i][j]);
//			fmt::print(aOutputFile, "{:.15g}\n", F2DArray[i][j]);
//			fmt::print(aOutputFile, "{:.15G}\n", F2DArray[i][j]);
		}
		aOutputFile << '\n';
	}
}
