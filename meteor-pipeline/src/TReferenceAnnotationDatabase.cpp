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

#include "TReferenceAnnotationDatabase.h"

using namespace std;

TReferenceAnnotationDatabase::TReferenceAnnotationDatabase() {
	  FMeteorDBType = mdbUnknown;
	  FIsLiteInfoActive = false;

}

TReferenceAnnotationDatabase::~TReferenceAnnotationDatabase() {
	// TODO Auto-generated destructor stub
}

// ------------------------------------------------------------------------------

bool TReferenceAnnotationDatabase::GetIsLiteInfoActive(){
	return FIsLiteInfoActive;
}

// ------------------------------------------------------------------------------

TLiteReferenceAnnotationTable & TReferenceAnnotationDatabase::GetLiteReferenceAnnotationTable()
{
	return FLiteReferenceAnnotationTable;
}

// ------------------------------------------------------------------------------

void TReferenceAnnotationDatabase::OpenLiteReferenceAnnotationTable(const string & aLiteReferenceAnnotationTableFileName)
{
  FLiteReferenceAnnotationTable.OpenReferenceAnnotationTable(aLiteReferenceAnnotationTableFileName);
  FIsLiteInfoActive = true;
}

// ------------------------------------------------------------------------------

bool TReferenceAnnotationDatabase::GetFragmentLiteInfo(const int aFragmentID, int & aFragmentSize, int & aStartLocation, int & aEndLocation, unsigned char & aIsOnPlusStrand) const
{
  if (FIsLiteInfoActive)
  {
    aFragmentSize = FLiteReferenceAnnotationTable.GetFragmentSize(aFragmentID);
    aStartLocation = 1;
    aEndLocation = aFragmentSize;
    aIsOnPlusStrand = 1;
    return true;
  }
//  else
//  {
//    FReferenceAnnotationTable.FilterFragment(aFragmentID, false);
//    aFragmentSize = ReferenceAnnotationTable.DataSet.Fields
//      [ReferenceAnnotationTable.FragmentSizeFieldNo].AsInteger;
//    aStartLocation = ReferenceAnnotationTable.DataSet.Fields
//      [ReferenceAnnotationTable.FragmentStartLocationFieldNo].AsInteger;
//    aEndLocation = ReferenceAnnotationTable.DataSet.Fields
//      [ReferenceAnnotationTable.FragmentEndLocationFieldNo].AsInteger;
//    aIsOnPlusStrand = ReferenceAnnotationTable.DataSet.Fields
//      [ReferenceAnnotationTable.FragmentStrandFieldNo].AsInteger;
//    return true;
//  }
  return false;
}

// ------------------------------------------------------------------------------

int TReferenceAnnotationDatabase::GetFragmentCount() const
{
  if (FIsLiteInfoActive) return FLiteReferenceAnnotationTable.GetFragmentCount();
//  else
//    FReferenceAnnotationTable.DataSet.RecordCount;
  return false;
}

// ------------------------------------------------------------------------------

void TReferenceAnnotationDatabase::FinalizeBuilding(){

	  //FLiteReferenceAnnotationTable.SaveTable(FReferenceDirectory + C_PATH_SEP + FReferenceName + "_" + C_REFERENCE_LITE_ANNOTATION_TABLE_NAME);

}
