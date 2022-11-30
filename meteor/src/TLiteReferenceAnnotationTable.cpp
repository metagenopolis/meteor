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

#include "TLiteReferenceAnnotationTable.h"
#include <fstream>
#include "utils.h"

using namespace std;

// ------------------------------------------------------------------------------

TLiteReferenceAnnotationTable::TLiteReferenceAnnotationTable() { ; }

// ------------------------------------------------------------------------------

TLiteReferenceAnnotationTable::~TLiteReferenceAnnotationTable(){ ; }

// ------------------------------------------------------------------------------

int TLiteReferenceAnnotationTable::GetFragmentCount() const
{
  // FTable size has been set to number of fragments + 1 because id_fragment starts from 1
  return FTable.size() - 1; // -1
}

// ------------------------------------------------------------------------------

int TLiteReferenceAnnotationTable::GetFragmentSize(const int  aFragmentID) const
{
  return FTable.at(aFragmentID);
}

// ------------------------------------------------------------------------------

void TLiteReferenceAnnotationTable::OpenReferenceAnnotationTable(const string & aReferenceAnnotationTableFileName)
{
  FTable.clear();
  string line;

  ifstream aLitefile(aReferenceAnnotationTableFileName.c_str()); // read only stream on the file

  // 1- initialize table size
  int aFragmentCount = 0;

  if (aLitefile){
	  while (getline(aLitefile, line)) aFragmentCount++;
	  FTable.resize(aFragmentCount + 1); // +1 because id_fragment starts from 1

	  // 2- fill table
	  aLitefile.clear(); aLitefile.seekg(0, aLitefile.beg); // rewind to the begining of the file
	  for ( string aFragmentId, aFragmentSize; getline(aLitefile, aFragmentId, '\t') && std::getline(aLitefile, aFragmentSize); )
	  {
		  FTable.at(atoi(aFragmentId.c_str())) = atoi(aFragmentSize.c_str());
	  }
  }
  else {
	  cerr<<"Error in "<<__FILE__<<", line "<<__LINE__<<", cannot open file: "<<aReferenceAnnotationTableFileName<<endl;
	  exit(1);
  }
  aLitefile.close();
}

// ------------------------------------------------------------------------------

//void TLiteReferenceAnnotationTable::SaveTable(const string & aReferenceAnnotationTableFileName)
//{
////  FTempTable.SaveToFile(aReferenceAnnotationTableFileName); //// TODO ????
//}


