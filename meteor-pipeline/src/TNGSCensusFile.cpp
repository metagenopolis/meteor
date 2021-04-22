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

#include "TNGSCensusFile.h"
#include "TCountedFragment.h"
#include <fstream>
#include <sstream>

using namespace std;

TNGSCensusFile::TNGSCensusFile(const string & aNGSCensusTablePathName, const string & aNGSCensusTableName):TNGSCensusData()
{
  FNGSCensusTableName     = aNGSCensusTableName;
  FNGSCensusTablePathName = aNGSCensusTablePathName;
  FCensusTableFileName    = FNGSCensusTablePathName + C_PATH_SEP + FNGSCensusTableName + C_CENSUS_TABLE_FILE_EXT;
  //FCensusProfileFileName  = FNGSCensusTablePathName + C_PATH_SEP + FNGSCensusTableName + C_CENSUS_PROFILE_TABLE_FILE_EXT; //// not used
}

TNGSCensusFile::~TNGSCensusFile() {

}

void TNGSCensusFile::WriteCountedFragmentList(TCountedFragmentList & aCountedFragmentList, const vector<TCountedReadType> & aCountedReadTypeList, char aWriteMode)
{
	if (aWriteMode == 'b')
		WriteCountedFragmentListBinary(aCountedFragmentList);
	else {
		if (aCountedReadTypeList.empty()) WriteCountedFragmentListAsciiText(aCountedFragmentList); // write all type (2 + 30 columns)
		else WriteCountedFragmentListAsciiText(aCountedFragmentList, aCountedReadTypeList); // write chosen type(s) only (2 + N columns)
	}
}

void TNGSCensusFile::WriteCountedFragmentListAsciiText(TCountedFragmentList & aCountedFragmentList, const vector<TCountedReadType> & aCountedReadTypeList){

	ofstream aCensusFile(FCensusTableFileName.c_str(), ios::out);

	if(aCensusFile){

		// write the header
		aCensusFile<<"FragmentID\tFragmentSize";
		for( size_t iType = 0; iType < aCountedReadTypeList.size(); iType++){
			aCensusFile << '\t' << C_CENSUS_TABLE_FILE_HEADER_ARRAY[aCountedReadTypeList[iType]];
		}
		aCensusFile << '\n';

		//// begin at 1 because id_fragment=0 does not exist (see TCountedFragmentList constructor: size is set to aFragmentCount+1).
		for (int aFragmentID = 1; aFragmentID <= aCountedFragmentList.GetFragmentCount(); aFragmentID++)
		{
			const TCountedFragment * aCountedFragment = aCountedFragmentList.GetItemByFragmentID(aFragmentID);

			if (aCountedFragment == NULL) {
				// write data for uncounted fragment: aFragmentID, -1, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
				aCensusFile<<aFragmentID<<"\t-1";
				for(size_t i=0; i<aCountedReadTypeList.size(); i++) aCensusFile << "\t0";
				aCensusFile << '\n';
			}
			else {
				aCensusFile << aFragmentID << '\t' << aCountedFragment->GetFragmentSize();

				// write the selected columns (stored in aCountedReadTypeList)
				for( size_t iType = 0; iType < aCountedReadTypeList.size(); iType++){

					aCensusFile.precision(15);

					switch(aCountedReadTypeList[iType]) {
						case crtTotalReads:
							aCensusFile  << '\t' << aCountedFragment->GetTotalReadCount(); break;
						case crtSharedReads:
							aCensusFile << '\t' << aCountedFragment->GetSharedReadCount(); break;
						case crtSmartSharedReads:
							aCensusFile << '\t' << aCountedFragment->GetSmartSharedReadCount(); break;
						case crtUniqueReads:
							aCensusFile << '\t' << aCountedFragment->GetUniqueReadCount(); break;
						case crtDirectTotalReads:
							aCensusFile << '\t' << aCountedFragment->GetDirectTotalReadCount(); break;
						case crtDirectSharedReads:
							aCensusFile << '\t' << aCountedFragment->GetDirectSharedReadCount(); break;
						case crtDirectSmartSharedReads:
							aCensusFile << '\t' << aCountedFragment->GetDirectSmartSharedReadCount(); break;
						case crtDirectUniqueReads:
							aCensusFile << '\t' << aCountedFragment->GetDirectUniqueReadCount(); break;
						case crtTotalReadsCoverage:
							aCensusFile << '\t' << aCountedFragment->GetTotalReadCoverage(); break;
						case crtUniqueReadsCoverage:
							aCensusFile << '\t' << aCountedFragment->GetUniqueReadCoverage(); break;
						case crtDirectTotalReadsCoverage:
							aCensusFile << '\t' << aCountedFragment->GetDirectTotalReadCoverage(); break;
						case crtDirectUniqueReadsCoverage:
							aCensusFile << '\t' << aCountedFragment->GetDirectUniqueReadCoverage(); break;
						case crtUnDirectTotalReadsCoverage:
							aCensusFile << '\t' << aCountedFragment->GetUnDirectTotalReadCoverage(); break;
						case crtUnDirectUniqueReadsCoverage:
							aCensusFile << '\t' << aCountedFragment->GetUnDirectUniqueReadCoverage(); break;
						case crtTotalReadsMeanCoverage:
							aCensusFile << '\t' << aCountedFragment->GetTotalReadMeanCoverage(); break;
						case crtSharedReadsMeanCoverage:
							aCensusFile << '\t' << aCountedFragment->GetSharedReadMeanCoverage(); break;
						case crtSmartSharedReadsMeanCoverage:
							aCensusFile << '\t' << aCountedFragment->GetSmartSharedReadMeanCoverage(); break;
						case crtUniqueReadsMeanCoverage:
							aCensusFile << '\t' << aCountedFragment->GetUniqueReadMeanCoverage(); break;
						case crtDirectTotalReadsMeanCoverage:
							aCensusFile << '\t' << aCountedFragment->GetDirectTotalReadMeanCoverage(); break;
						case crtDirectSharedReadsMeanCoverage:
							aCensusFile << '\t' << aCountedFragment->GetDirectSharedReadMeanCoverage(); break;
						case crtDirectSmartSharedReadsMeanCoverage:
							aCensusFile << '\t' << aCountedFragment->GetDirectSmartSharedReadMeanCoverage(); break;
						case crtDirectUniqueReadsMeanCoverage:
							aCensusFile << '\t' << aCountedFragment->GetDirectUniqueReadMeanCoverage(); break;
						case crtUndirectTotalReadsMeanCoverage:
							aCensusFile << '\t' << aCountedFragment->GetUnDirectTotalReadMeanCoverage(); break;
						case crtUndirectSharedReadsMeanCoverage:
							aCensusFile << '\t' << aCountedFragment->GetUnDirectSharedReadMeanCoverage(); break;
						case crtUndirectSmartSharedReadsMeanCoverage:
							aCensusFile << '\t' << aCountedFragment->GetUnDirectSmartSharedReadMeanCoverage(); break;
						case crtUndirectUniqueReadsMeanCoverage:
							aCensusFile << '\t' << aCountedFragment->GetUnDirectUniqueReadMeanCoverage(); break;

						case crtUndirectTotalReads:
							aCensusFile << '\t' << aCountedFragment->GetUnDirectTotalReadCount(); break;
						case crtUndirectSharedReads:
							aCensusFile << '\t' << aCountedFragment->GetUnDirectSharedReadCount(); break;
						case crtUndirectSmartSharedReads:
							aCensusFile << '\t' << aCountedFragment->GetUnDirectSmartSharedReadCount(); break;
						case crtUndirectUniqueReads:
							aCensusFile << '\t' << aCountedFragment->GetUnDirectUniqueReadCount(); break;
					}
				}
				aCensusFile << '\n';
			}
		}
		aCensusFile.close();
	}
	else {
		cerr << "Error, opening "<<FCensusTableFileName<<" failed."<< endl;
		exit(1);
	}
}

void TNGSCensusFile::WriteCountedFragmentListAsciiText(TCountedFragmentList & aCountedFragmentList){

	ofstream aCensusFile(FCensusTableFileName.c_str(), ios::out);

	if(aCensusFile){

		// write the header
		aCensusFile<<C_CENSUS_TABLE_FILE_HEADER<<endl;

		//// begin at 1 because id_fragment=0 does not exist (see TCountedFragmentList constructor: size is set to aFragmentCount+1).
		for (int aFragmentID = 1; aFragmentID <= aCountedFragmentList.GetFragmentCount(); aFragmentID++)
		{
			const TCountedFragment * aCountedFragment = aCountedFragmentList.GetItemByFragmentID(aFragmentID);

			if (aCountedFragment == NULL) {
				// write data for uncounted fragment: aFragmentID, -1, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
				aCensusFile<<aFragmentID<<"\t-1";
				for(int i=0; i<C_CENSUS_TABLE_FILE_COL_NUMBER; i++) aCensusFile << "\t0";
				aCensusFile << '\n';
			}
			else {
				aCensusFile.precision(15);
				aCensusFile << aFragmentID << '\t';
				aCensusFile << aCountedFragment->GetFragmentSize() << '\t';
				aCensusFile << aCountedFragment->GetTotalReadCount() << '\t';
				aCensusFile << aCountedFragment->GetSharedReadCount() << '\t';
				aCensusFile << aCountedFragment->GetSmartSharedReadCount() << '\t';
				aCensusFile << aCountedFragment->GetUniqueReadCount() << '\t';
				aCensusFile << aCountedFragment->GetDirectTotalReadCount() << '\t';
				aCensusFile << aCountedFragment->GetDirectSharedReadCount() << '\t';
				aCensusFile << aCountedFragment->GetDirectSmartSharedReadCount() << '\t';
				aCensusFile << aCountedFragment->GetDirectUniqueReadCount() << '\t';
				aCensusFile << aCountedFragment->GetTotalReadCoverage() << '\t';
				aCensusFile << aCountedFragment->GetUniqueReadCoverage() << '\t';
				aCensusFile << aCountedFragment->GetDirectTotalReadCoverage() << '\t';
				aCensusFile << aCountedFragment->GetDirectUniqueReadCoverage() << '\t';
				aCensusFile << aCountedFragment->GetUnDirectTotalReadCoverage() << '\t';
				aCensusFile << aCountedFragment->GetUnDirectUniqueReadCoverage() << '\t';
				aCensusFile << aCountedFragment->GetTotalReadMeanCoverage() << '\t';
				aCensusFile << aCountedFragment->GetSharedReadMeanCoverage() << '\t';
				aCensusFile << aCountedFragment->GetSmartSharedReadMeanCoverage() << '\t';
				aCensusFile << aCountedFragment->GetUniqueReadMeanCoverage() << '\t';
				aCensusFile << aCountedFragment->GetDirectTotalReadMeanCoverage() << '\t';
				aCensusFile << aCountedFragment->GetDirectSharedReadMeanCoverage() << '\t';
				aCensusFile << aCountedFragment->GetDirectSmartSharedReadMeanCoverage() << '\t';
				aCensusFile << aCountedFragment->GetDirectUniqueReadMeanCoverage() << '\t';
				aCensusFile << aCountedFragment->GetUnDirectTotalReadMeanCoverage() << '\t';
				aCensusFile << aCountedFragment->GetUnDirectSharedReadMeanCoverage() << '\t';
				aCensusFile << aCountedFragment->GetUnDirectSmartSharedReadMeanCoverage() << '\t';
				aCensusFile << aCountedFragment->GetUnDirectUniqueReadMeanCoverage() << '\t';

				aCensusFile << aCountedFragment->GetUnDirectTotalReadCount() << '\t';
				aCensusFile << aCountedFragment->GetUnDirectSharedReadCount() << '\t';
				aCensusFile << aCountedFragment->GetUnDirectSmartSharedReadCount() << '\t';
				aCensusFile << aCountedFragment->GetUnDirectUniqueReadCount() << '\n';
			}
		}
		aCensusFile.close();
	}
	else {
		cerr << "Error, opening "<<FCensusTableFileName<<" failed."<< endl;
		exit(1);
	}
}

void TNGSCensusFile::WriteCountedFragmentListBinary(TCountedFragmentList & aCountedFragmentList)
{
  TNGSCensusItem aNGSCensusItem;

  // open stream in write/binary mode
  ofstream aCensusFile(FCensusTableFileName.c_str(), ios::out | ios::binary);

  if(aCensusFile){
	  //// begin at 1 because id_fragment=0 does not exist
	  //// (see TCountedFragmentList constructor: size is set to aFragmentCount+1)
	  for (int aFragmentID = 1; aFragmentID <= aCountedFragmentList.GetFragmentCount(); aFragmentID++)
	  {
		const TCountedFragment * aCountedFragment = aCountedFragmentList.GetItemByFragmentID(aFragmentID);
		////const TCountedFragment * aCountedFragment = aCountedFragmentList.GetItemByFragmentID[aFragmentID];

		if (aCountedFragment == NULL)
		{
		  aNGSCensusItem.FragmentID = aFragmentID;
		  // fragment size unknown in the data structure but not needed for the next task because count=0
		  // --> set to -1
		  aNGSCensusItem.FragmentSize = -1;
		  aNGSCensusItem.TotalReadCount = 0;
		  aNGSCensusItem.SharedReadCount = 0.0;
		  aNGSCensusItem.SmartSharedReadCount = 0.0;
		  aNGSCensusItem.UniqueReadCount = 0;
		  aNGSCensusItem.DirectTotalReadCount = 0;
		  aNGSCensusItem.DirectSharedReadCount = 0.0;
		  aNGSCensusItem.DirectSmartSharedReadCount = 0.0;
		  aNGSCensusItem.DirectUniqueReadCount = 0;
		  aNGSCensusItem.TotalReadCoverage = 0.0;
		  aNGSCensusItem.UniqueReadCoverage = 0.0;
		  aNGSCensusItem.DirectTotalReadCoverage = 0.0;
		  aNGSCensusItem.DirectUniqueReadCoverage = 0.0;
		  aNGSCensusItem.UnDirectTotalReadCoverage = 0.0;
		  aNGSCensusItem.UnDirectUniqueReadCoverage = 0.0;
		  aNGSCensusItem.TotalReadMeanCoverage = 0.0;
		  aNGSCensusItem.SharedReadMeanCoverage = 0.0;
		  aNGSCensusItem.SmartSharedReadMeanCoverage = 0.0;
		  aNGSCensusItem.UniqueReadMeanCoverage = 0.0;
		  aNGSCensusItem.DirectTotalReadMeanCoverage = 0.0;
		  aNGSCensusItem.DirectSharedReadMeanCoverage = 0.0;
		  aNGSCensusItem.DirectSmartSharedReadMeanCoverage = 0.0;
		  aNGSCensusItem.DirectUniqueReadMeanCoverage = 0.0;
		  aNGSCensusItem.UnDirectTotalReadMeanCoverage = 0.0;
		  aNGSCensusItem.UnDirectSharedReadMeanCoverage = 0.0;
		  aNGSCensusItem.UnDirectSmartSharedReadMeanCoverage = 0.0;
		  aNGSCensusItem.UnDirectUniqueReadMeanCoverage = 0.0;

		  aNGSCensusItem.UnDirectTotalReadCount = 0;
		  aNGSCensusItem.UnDirectSharedReadCount = 0.0;
		  aNGSCensusItem.UnDirectSmartSharedReadCount = 0.0;
		  aNGSCensusItem.UnDirectUniqueReadCount = 0;
		  //aNGSCensusItem = {aFragmentID, -1, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; //// c++11
		}
		else
		{
		  aNGSCensusItem.FragmentID = aFragmentID;
		  aNGSCensusItem.FragmentSize = aCountedFragment->GetFragmentSize();
		  aNGSCensusItem.TotalReadCount = aCountedFragment->GetTotalReadCount();
		  aNGSCensusItem.SharedReadCount = aCountedFragment->GetSharedReadCount();
		  aNGSCensusItem.SmartSharedReadCount = aCountedFragment->GetSmartSharedReadCount();
		  aNGSCensusItem.UniqueReadCount = aCountedFragment->GetUniqueReadCount();
		  aNGSCensusItem.DirectTotalReadCount = aCountedFragment->GetDirectTotalReadCount();
		  aNGSCensusItem.DirectSharedReadCount = aCountedFragment->GetDirectSharedReadCount();
		  aNGSCensusItem.DirectSmartSharedReadCount = aCountedFragment->GetDirectSmartSharedReadCount();
		  aNGSCensusItem.DirectUniqueReadCount = aCountedFragment->GetDirectUniqueReadCount();
		  aNGSCensusItem.TotalReadCoverage = aCountedFragment->GetTotalReadCoverage();
		  aNGSCensusItem.UniqueReadCoverage = aCountedFragment->GetUniqueReadCoverage();
		  aNGSCensusItem.DirectTotalReadCoverage = aCountedFragment->GetDirectTotalReadCoverage();
		  aNGSCensusItem.DirectUniqueReadCoverage = aCountedFragment->GetDirectUniqueReadCoverage();
		  aNGSCensusItem.UnDirectTotalReadCoverage = aCountedFragment->GetUnDirectTotalReadCoverage();
		  aNGSCensusItem.UnDirectUniqueReadCoverage = aCountedFragment->GetUnDirectUniqueReadCoverage();
		  aNGSCensusItem.TotalReadMeanCoverage = aCountedFragment->GetTotalReadMeanCoverage();
		  aNGSCensusItem.SharedReadMeanCoverage = aCountedFragment->GetSharedReadMeanCoverage();
		  aNGSCensusItem.SmartSharedReadMeanCoverage = aCountedFragment->GetSmartSharedReadMeanCoverage();
		  aNGSCensusItem.UniqueReadMeanCoverage = aCountedFragment->GetUniqueReadMeanCoverage();
		  aNGSCensusItem.DirectTotalReadMeanCoverage = aCountedFragment->GetDirectTotalReadMeanCoverage();
		  aNGSCensusItem.DirectSharedReadMeanCoverage = aCountedFragment->GetDirectSharedReadMeanCoverage();
		  aNGSCensusItem.DirectSmartSharedReadMeanCoverage = aCountedFragment->GetDirectSmartSharedReadMeanCoverage();
		  aNGSCensusItem.DirectUniqueReadMeanCoverage = aCountedFragment->GetDirectUniqueReadMeanCoverage();
		  aNGSCensusItem.UnDirectTotalReadMeanCoverage = aCountedFragment->GetUnDirectTotalReadMeanCoverage();
		  aNGSCensusItem.UnDirectSharedReadMeanCoverage = aCountedFragment->GetUnDirectSharedReadMeanCoverage();
		  aNGSCensusItem.UnDirectSmartSharedReadMeanCoverage = aCountedFragment->GetUnDirectSmartSharedReadMeanCoverage();
		  aNGSCensusItem.UnDirectUniqueReadMeanCoverage = aCountedFragment->GetUnDirectUniqueReadMeanCoverage();

		  aNGSCensusItem.UnDirectTotalReadCount = aCountedFragment->GetUnDirectTotalReadCount();
		  aNGSCensusItem.UnDirectSharedReadCount = aCountedFragment->GetUnDirectSharedReadCount();
		  aNGSCensusItem.UnDirectSmartSharedReadCount = aCountedFragment->GetUnDirectSmartSharedReadCount();
		  aNGSCensusItem.UnDirectUniqueReadCount = aCountedFragment->GetUnDirectUniqueReadCount();
		}
		// push on file
		// write counting stats (binary mode)
		aCensusFile.write((char*)&aNGSCensusItem, sizeof(aNGSCensusItem));
	  }
	  aCensusFile.close();
  }
  else {
	cerr << "Error, opening "<<FCensusTableFileName<<" failed."<< endl;
	exit(1);
  }
}
