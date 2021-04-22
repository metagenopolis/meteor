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

#include "TMappingData.h"
#include "TSAMReader.h"
#include <cstdlib>

using namespace std;

//TMappingData::TMappingData() {
//	// TODO Auto-generated constructor stub

//}

TMappingData::~TMappingData() {
}

// ------------------------------------------------------------------------------

TMappingData::TMappingData(
		const int aMappedReadCount, const double aMaxDistance, const bool aIsRelativeDistance,
		const bool aIsLocalAlignment, const bool aKeepInternalLocalAlignment,
		const int aAlignmentLengthCutoff, const int aSoftClippingLengthCutoff, const TMappingFileFormat & aMappingFileFormat)
{
  FMappedReadCount = aMappedReadCount;
  FMaxDistance = aMaxDistance;
  FIsRelativeDistance = aIsRelativeDistance;

  FIsLocalAlignment = aIsLocalAlignment; //// LOCAL
  FKeepInternalLocalAlignment = aKeepInternalLocalAlignment; //// LOCAL
  FAlignmentLengthCutoff = aAlignmentLengthCutoff; //// LOCAL
  FSoftClippingLengthCutoff = aSoftClippingLengthCutoff; //// LOCAL

  FMappingFileFormat = aMappingFileFormat;

  CreateMappingIntersectionArray(aMappedReadCount); // aMappedReadCount is known from field "processed_read_count" in census_stage_1.ini file
}

// ------------------------------------------------------------------------------

void TMappingData::AddBowtieExcludedMapping(
		const string & aFileName, const int & aIndexExcludedMapping,
		const double & aExcludedMaxDistance, const bool aIsRelativeMaxDistance, const bool aIsOnBothStrand)
{
//  int aReadID;
//  double aDistance;
//  TIntersectionStatus aCurrentIntersectionStatus;
//  /*
//    by default IntersectionStatus = 0 (no mapping)
//    intersection = 1 --> mapping on master reference
//    In bit representation Intersection = 00001
//    If a read is mapped also on an other reference for example the first
//    in IntersectNGSMappingDatabaseList (with iMapping=0) then
//    we move 1 bit in left (shl) the IntersectionStatus=1 and so on:
//    Here (in this example):
//    aCurrentIntersectionStatus = 1 shl 1 (iMapping=0 + 1) = 2
//    bit representation is 00010
//  */
//  aCurrentIntersectionStatus = C_MAIN_MAPPING_INTERSECTION_STATUS << aIndexExcludedMapping;
//
//  ifstream aBowtieFile(aFileName.c_str());
//
//  if (aBowtieFile)
//  {
//    TBowtieReader aBowtieReader(FMappingFileFormat == mffLiteDefaultBowtie); //// TODO
//    TBowtieMatch aBowtieMatch; //// TODO
//
//    string line;
//    while (getline(aBowtieFile, line)) // ending '\n' is removed
//    {
//      aBowtieReader.GetMatch(line, aBowtieMatch);
//      // test if read accepted on both strand or if not if is mapped on + strand
//      if (aIsOnBothStrand || (!aIsOnBothStrand && aBowtieMatch.Strand == 1))
//      {
//        aReadID = atoi(aBowtieMatch.ReadName.c_str());
//        aDistance = aBowtieMatch.EditDistance;
//        if (aIsRelativeMaxDistance)
//          aDistance = aDistance / aBowtieMatch.AlignmentLength;
//
//
//        {
//          if ((aDistance <= aExcludedMaxDistance))
//          // in the template
//          {
//            if (aDistance < FMappingIntersectionArray[aReadID].MinDistance)
//            // better alignement than on the main reference
//            {
//              /*
//                Here IntersectionStatus is associated with IntersectionStatus for
//                the current reference with the 'OR' bit-to-bit operateur
//                For example, if current read is mapped on master reference:
//                IntersectionStatus = 1
//                If the same read is mapped on the other reference with iMapping=0 then
//                IntersectionStatus = IntersectionStatus or aCurrentIntersectionStatus
//                IntersectionStatus = 1 or 2
//                IntersectionStatus = 3
//                In bit representation : 00011
//              */
//            	FMappingIntersectionArray[aReadID].IntersectionStatus = FMappingIntersectionArray[aReadID].IntersectionStatus | aCurrentIntersectionStatus;
//              FExcludedReadCountArray[aIndexExcludedMapping - 1]++;
//              FMappingIntersectionArray[aReadID].IsRejected = true;
//            }
//          }
//        }
//      }
//    }
//    aBowtieFile.close();
//  }
}

// ------------------------------------------------------------------------------

void TMappingData::AddExcludedMapping(
		const string & aFileName,  const int aIndexExcludedMapping,	const double aExcludedMaxDistance,
		const bool aIsRelativeMaxDistance, const bool aIsOnBothStrand)
{
	FExcludedReadCountArraySize++;
	FExcludedReadCountArray.push_back(0);

	switch(FMappingFileFormat){
		case mffCSFasta: break;
		case mffSAM:
			AddSAMExcludedMapping(aFileName, aIndexExcludedMapping, aExcludedMaxDistance, aIsRelativeMaxDistance, aIsOnBothStrand);
			break;
		case mffBAM: break;
		case mffDefaultBowtie:
			AddBowtieExcludedMapping(aFileName, aIndexExcludedMapping, aExcludedMaxDistance, aIsRelativeMaxDistance, aIsOnBothStrand);
			break;
		case mffLiteDefaultBowtie:
			AddBowtieExcludedMapping(aFileName, aIndexExcludedMapping, aExcludedMaxDistance, aIsRelativeMaxDistance, aIsOnBothStrand);
			break;
		case mffUnknown: break;
	}
}

// ------------------------------------------------------------------------------

void TMappingData::AddSAMExcludedMapping(
		const string & aSAMFileName, const int & aIndexExcludedMapping,
		const double & aExcludedMaxDistance, const bool aIsRelativeMaxDistance, const bool aIsOnBothStrand){

//  int64 aFileSize;
  int aReadID;
  double aDistance;
  TIntersectionStatus aCurrentIntersectionStatus;

  /*
    by default IntersectionStatus = 0 (no mapping)
    intersection = 1 --> mapping on master reference
    In bit representation Intersection = 00001
    If a read is mapped also on an other reference for example the first
    in IntersectNGSMappingDatabaseList (with iMapping=0) then
    we move 1 bit in left (shl) the IntersectionStatus=1 and so on:
    Here (in this example):
    aCurrentIntersectionStatus = 1 shl 1 (iMapping=0 + 1) = 2
    bit representation is 00010
  */
  aCurrentIntersectionStatus = C_MAIN_MAPPING_INTERSECTION_STATUS << aIndexExcludedMapping; // << is left shift (shl delphi)

  ifstream aSAMFile(aSAMFileName.c_str());

  if (aSAMFile)
  {
	TSAMReader aSAMReader;
	TSAMMatch aSAMMatch;
	string line;
	while (getline(aSAMFile, line)) // ending '\n' is removed
    {
      aSAMReader.GetMatch(line, aSAMMatch);
      if (! aSAMMatch.Unmapped)
      {
        // test if read accepted on both strand or if not if is mapped on + strand
        if ( aIsOnBothStrand || ( (!aIsOnBothStrand) && (aSAMMatch.Strand == 1) ) )
        //if ( aIsOnBothStrand || aSAMMatch.Strand == 1 ) //// this test should be equivalent
        {
          aReadID = atoi(aSAMMatch.QName.c_str());
          aDistance = aSAMMatch.EditDistance;

          if (aIsRelativeMaxDistance)
            aDistance = aDistance / (double)aSAMMatch.AlignmentLength;

		  if (aDistance <= aExcludedMaxDistance)
		  // in the template
		  {
			if (aDistance < FMappingIntersectionArray.at(aReadID).MinDistance) //// AT
		    {
		      // record the best match min distance for this read among all excuded references
		      FMappingIntersectionArray.at(aReadID).MinDistance = aDistance;
			  /*
			    Here IntersectionStatus is associated with IntersectionStatus for
			    the current reference with the 'OR' bit-to-bit operateur
			    For example, if current read is mapped on master reference:
			    IntersectionStatus = 1
			    If the same read is mapped on the other reference with iMapping=0 then
			    IntersectionStatus = IntersectionStatus or aCurrentIntersectionStatus
			    IntersectionStatus = 1 or 2
			    IntersectionStatus = 3
			    In bit representation : 00011
			  */
		      FMappingIntersectionArray.at(aReadID).IntersectionStatus = FMappingIntersectionArray.at(aReadID).IntersectionStatus | aCurrentIntersectionStatus;
			  ////FMappingIntersectionArray.at(aReadID).IntersectionStatus |= aCurrentIntersectionStatus;
			  FExcludedReadCountArray.at(aIndexExcludedMapping - 1)++; //// AT
			  FMappingIntersectionArray.at(aReadID).IsRejected = true; //// AT
		    }
		  }
        }
      }
    }
    aSAMFile.close();
  }
}

// ------------------------------------------------------------------------------

void TMappingData::IncCountingStatistics(TCountingStatistics & aCountingStatistics) const {

	for ( size_t i = 1; i < FMappingIntersectionArray.size(); i++ ){

		if (FMappingIntersectionArray[i].IntersectionStatus == 0) // no mapping
		{
			aCountingStatistics.NotMappedReads++;
			aCountingStatistics.NotMappedCleanReads++;
		}
		else if ( FMappingIntersectionArray[i].IntersectionStatus == C_MAIN_MAPPING_INTERSECTION_STATUS ) {

			if (!FMappingIntersectionArray[i].IsRejected) {

				aCountingStatistics.TotalCleanCountedReads++;
				if ( FMappingIntersectionArray[i].MatchesCount == 1 )
					aCountingStatistics.UniqueCleanCountedReads++;
			}
			else {
				aCountingStatistics.RejectedReads++;
				aCountingStatistics.RejectedQualityReads++;
			}
		}
		else {
			aCountingStatistics.RejectedReads++;
			aCountingStatistics.RejectedReferenceReads++;
		}
	}
}

// ------------------------------------------------------------------------------

void TMappingData::IncIntersectionCount(vector<int> & aIntersectionCountArray) {

	//for (size_t i = 1 ; High(FMappingIntersectionArray; i++)
	for (size_t i = 1 ; i < FMappingIntersectionArray.size(); i++) //// i = 1 cause first readID is 1.
		aIntersectionCountArray.at(FMappingIntersectionArray[i].IntersectionStatus)++; //// AT
}

// ------------------------------------------------------------------------------

int TMappingData::GetExcludedReadCount(const int aIndexExcludedMapping) const
{
  return FExcludedReadCountArray.at(aIndexExcludedMapping); //// AT
}

// ------------------------------------------------------------------------------

void TMappingData::CreateMappingIntersectionArray(const int aMappedReadCount){
	//// Y a quoi dans FMappingIntersectionArray[0] ????

	FMappingIntersectionArray.resize(aMappedReadCount + 1);

	for (int i = 1; i <= aMappedReadCount; i++)
	{
		FMappingIntersectionArray[i].MatchesCount = 0;
		FMappingIntersectionArray[i].MinDistance = 1000.0;
		FMappingIntersectionArray[i].IsRejected = true;
		FMappingIntersectionArray[i].IntersectionStatus = C_NO_MAPPING_INTERSECTION_STATUS;

		// SEARCH_BUG
//		MinDistance := 1000.0;
//		MatchesCount := 1;
//		IsRejected := FALSE;
//		IntersectionStatus := C_MAIN_MAPPING_INTERSECTION_STATUS;
	}
}

// ------------------------------------------------------------------------------
