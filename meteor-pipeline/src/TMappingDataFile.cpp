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

#include "TMappingDataFile.h"

using namespace std;

TMappingDataFile::TMappingDataFile(
		const int aMappedReadCount, const double aMaxDistance, const bool aIsRelativeDistance,
		const bool aIsLocalAlignment, const bool aKeepInternalLocalAlignment,
		const int aAlignmentLengthCutoff, const int aSoftClippingLengthCutoff, const TMappingFileFormat & aMappingFileFormat)
		:TMappingData(aMappedReadCount, aMaxDistance, aIsRelativeDistance, aIsLocalAlignment, aKeepInternalLocalAlignment, aAlignmentLengthCutoff, aSoftClippingLengthCutoff, aMappingFileFormat)
{
	FReferenceAnnotationDatabase = NULL;
}

// ------------------------------------------------------------------------------

TMappingDataFile::~TMappingDataFile() {
	// TODO Auto-generated destructor stub
}

// ------------------------------------------------------------------------------

void TMappingDataFile::AddMainMapping(const string & aFileName)
{
  FMainMappingFileName = aFileName;
}

// ------------------------------------------------------------------------------

void TMappingDataFile::CountFragmentFromBowtieFile(const TCountedFragmentList & aCountedFragmentList){;} //// TODO

// ------------------------------------------------------------------------------

void TMappingDataFile::CountFragmentFromSAMFile_old(TCountedFragmentList & aCountedFragmentList)
{
// 11 ALIGNMENT SECTION MANDATORY FIELDS
//	QNAME
//	FLAG
//	RNAME
//	POS
//	MAPQ
//	CIGAR
//	RNEXT
//	PNEXT
//	TLEN
//	SEQ
//	QUAL

	TSAMReader aSAMReader;
	TSAMMatch aSAMMatch;
	int aReadID, aCurrentReadID = -1;
	//double aDistance; //// unused

	////  char aSAMFileBuf[65535]; //64K buffer, used with alternative std::getline function

	TReadMatchInfoList aReadMatchInfoList;

	ifstream aSAMFile(FMainMappingFileName.c_str()); // read only stream on the file
	string line;

	if(aSAMFile){

		while (getline(aSAMFile, line)) { // ending '\n' is removed

			// store current alignment (match) properties in aSAMMatch.
			aSAMReader.GetMatch(line, aSAMMatch);

			if (!aSAMMatch.Unmapped)
			{
				aReadID = atoi(aSAMMatch.QName.c_str());
				//// NB: In SAM specification (june 2017), QNAME field is a string including any non-blank characters (regex: [!-?A-~]{1,254})
				//     But in meteor context, it is composed of digits only

				if (aCurrentReadID != aReadID)
				{
					if (aCurrentReadID != -1)
					{
//						if (aCurrentReadID == C_READ_ID){ //// DEBUG
//							cerr <<"Read "<<aCurrentReadID <<" ; aReadMatchInfoList.Count : "<< aReadMatchInfoList.Count() << endl;
//						}
						// Fill MappingIntersectionArray for statistics on aCurrentReadID (aReadID's changed)
						SetMappingIntersectionArray(aReadMatchInfoList);
						// and now count...
						CountSAMReadList(aReadMatchInfoList, aCountedFragmentList);
					}
					// and finally reinitialize for next read block
					aReadMatchInfoList.Initialize(aReadID);
					aCurrentReadID = aReadID;
				}

				/// continue to fill ReadMatchInfoList

				// add an empty TReadMatchInfoItem and return its reference (IsRejected is set to true).
				TReadMatchInfoItem & aReadMatchInfoItem = aReadMatchInfoList.Add();
				// 1-based offset into the forward reference strand where leftmost character of the alignment occurs
				aReadMatchInfoItem.ReadLocation = aSAMMatch.Pos;
				aReadMatchInfoItem.ReadStrand = (aSAMMatch.Strand == 1);
				aReadMatchInfoItem.ReadMismatches = aSAMMatch.Mismatches;
				aReadMatchInfoItem.ReadEditDistance = aSAMMatch.EditDistance;
				aReadMatchInfoItem.ReadEntryID = atoi(aSAMMatch.RName.c_str());
				aReadMatchInfoItem.ReadMappingLength = aSAMMatch.AlignmentLength;
				aReadMatchInfoItem.ReadSAMFlag = aSAMMatch.Flag;
				aReadMatchInfoItem.ReadSAMCigar = aSAMMatch.Cigar;
				aReadMatchInfoItem.ReadSAMOpt = aSAMMatch.Opt;

				aReadMatchInfoItem.ReadDistance = aSAMMatch.EditDistance;
				if (FIsRelativeDistance) aReadMatchInfoItem.ReadDistance /= (double)aSAMMatch.AlignmentLength;
			}
		}
		// Process last read block
		if (aCurrentReadID != -1)
		{
			// Fill MappingIntersectionArray for statistics
			SetMappingIntersectionArray(aReadMatchInfoList);
			// and now count...
			CountSAMReadList(aReadMatchInfoList, aCountedFragmentList);
		}
		aSAMFile.close();
	}
	else {
		cerr<<"Error in "<<__FILE__<<", line "<<__LINE__<<", cannot open file: "<<FMainMappingFileName<<endl;
		exit(1);
	}
}

// ------------------------------------------------------------------------------

void TMappingDataFile::CountFragmentFromSAMFile(TCountedFragmentList & aCountedFragmentList)
{
// 11 ALIGNMENT SECTION MANDATORY FIELDS
//	QNAME
//	FLAG
//	RNAME
//	POS
//	MAPQ
//	CIGAR
//	RNEXT
//	PNEXT
//	TLEN
//	SEQ
//	QUAL

	TSAMReader aSAMReader;
	TSAMMatch aSAMMatch;
	int aReadID, aCurrentReadID = -1;

	bool aOKToProcessRead;

	TLocalMatchInfo aLocalMatchInfo;

	////  char aSAMFileBuf[65535]; //64K buffer, used with alternative std::getline function

	TReadMatchInfoList aReadMatchInfoList;

	ifstream aSAMFile(FMainMappingFileName.c_str()); // read only stream on the file
	string line;

	if(aSAMFile){

		while (getline(aSAMFile, line)) { // ending '\n' is removed

			// store current alignment (match) properties in aSAMMatch.
			aSAMReader.GetMatch(line, aSAMMatch);

			//if (aSAMMatch.RName == "1") { cerr << endl << "line: "<< line << endl << "fragment 1 unmatched ? " << aSAMMatch.Unmapped << endl; }

			if (!aSAMMatch.Unmapped)
			{
				aReadID = atoi(aSAMMatch.QName.c_str());
				//// NB: In SAM specification (june 2017), QNAME field is a string including any non-blank characters (regex: [!-?A-~]{1,254})
				//     But in meteor context, it is composed of digits only

				if (aCurrentReadID != aReadID)
				{
					if (aCurrentReadID != -1)
					{
//						if (aCurrentReadID == C_READ_ID){ //// DEBUG (other read to test 1954315)
//							cerr <<"Read 5059252 ; aReadMatchInfoList.Count : "<< aReadMatchInfoList.Count() << endl;
//						}
						// Fill MappingIntersectionArray for statistics on aCurrentReadID (aReadID's changed)
						SetMappingIntersectionArray(aReadMatchInfoList);
						// and now count...
						CountSAMReadList(aReadMatchInfoList, aCountedFragmentList);
					}
					// and finally reinitialize for next read block
					aReadMatchInfoList.Initialize(aReadID);
					aCurrentReadID = aReadID;
				}

				/// continue to fill ReadMatchInfoList
		        // 1- check if local alignment is acceptable
		        // retrieve fragment (reference) length
		        int aFragmentSize = FReferenceAnnotationDatabase->GetLiteReferenceAnnotationTable().GetFragmentSize(atoi(aSAMMatch.RName.c_str()));

		        if (FIsLocalAlignment){
		          aLocalMatchInfo = IsAcceptedLocalMatch(aFragmentSize, aSAMMatch.Cigar, aSAMMatch.Pos, FAlignmentLengthCutoff, FSoftClippingLengthCutoff, FKeepInternalLocalAlignment);
		          aOKToProcessRead = aLocalMatchInfo.IsAccepted;
		        }
		        else {
				  aLocalMatchInfo.IsAccepted = true;
		          aLocalMatchInfo.NewReadStartLocation = aSAMMatch.Pos;
		          aLocalMatchInfo.NewReadEndLocation = aSAMMatch.Pos + aSAMMatch.AlignmentLength - 1;
		          aLocalMatchInfo.SupplementaryMismatchesCount = 0;
		          aOKToProcessRead = true;
		        }
		        if (aOKToProcessRead){
					// add an empty TReadMatchInfoItem and return its reference (IsRejected is set to true).
					TReadMatchInfoItem & aReadMatchInfoItem = aReadMatchInfoList.Add();
					// 1-based offset into the forward reference strand where leftmost character of the alignment occurs
					// Warning : take into account correction brought by clipped (local) alignment
					aReadMatchInfoItem.ReadLocation = aLocalMatchInfo.NewReadStartLocation;
					aReadMatchInfoItem.ReadStrand = (aSAMMatch.Strand == 1);
					aReadMatchInfoItem.ReadMismatches = aSAMMatch.Mismatches;
					aReadMatchInfoItem.ReadEditDistance = aSAMMatch.EditDistance + aLocalMatchInfo.SupplementaryMismatchesCount;
					aReadMatchInfoItem.ReadEntryID = atoi(aSAMMatch.RName.c_str());
					aReadMatchInfoItem.ReadMappingLength = aSAMMatch.AlignmentLength + aLocalMatchInfo.SupplementaryMismatchesCount;
					aReadMatchInfoItem.ReadSAMFlag = aSAMMatch.Flag;
					aReadMatchInfoItem.ReadSAMCigar = aSAMMatch.Cigar;
					aReadMatchInfoItem.ReadSAMOpt = aSAMMatch.Opt;

					aReadMatchInfoItem.ReadDistance = (double)aSAMMatch.EditDistance + aLocalMatchInfo.SupplementaryMismatchesCount;
					if (FIsRelativeDistance) aReadMatchInfoItem.ReadDistance /= (double)(aSAMMatch.AlignmentLength + aLocalMatchInfo.SupplementaryMismatchesCount);
		        }
			}
		}
		// Process last read block
		if (aCurrentReadID != -1)
		{
			// Fill MappingIntersectionArray for statistics
			SetMappingIntersectionArray(aReadMatchInfoList);
			// and now count...
			CountSAMReadList(aReadMatchInfoList, aCountedFragmentList);
		}
		aSAMFile.close();
	}
	else {
		cerr<<"Error in "<<__FILE__<<", line "<<__LINE__<<", cannot open file: "<<FMainMappingFileName<<endl;
		exit(1);
	}
}

// ------------------------------------------------------------------------------

void TMappingDataFile::CountFragment(TCountedFragmentList & aCountedFragmentList, TReferenceAnnotationDatabase & aReferenceAnnotationDatabase)
{
  // Reference annotation dataset preparation
	FReferenceAnnotationDatabase = &aReferenceAnnotationDatabase; // address assignation

	if(!FReferenceAnnotationDatabase->GetIsLiteInfoActive()) { //// TODO later
//		FReferenceAnnotationDatabase.ReferenceAnnotationTable.CacheGeneRecords(1000);
//		FReferenceAnnotationDatabase.ReferenceAnnotationTable.DataSet.First;
	}
	switch(FMappingFileFormat){
		case mffSAM:
			CountFragmentFromSAMFile(aCountedFragmentList); break;
		case mffDefaultBowtie:
//			CountFragmentFromBowtieFile(aCountedFragmentList); break;
		case mffLiteDefaultBowtie:
//			CountFragmentFromBowtieFile(aCountedFragmentList); break;
		default: break;
	}
}

// ------------------------------------------------------------------------------
//// TODO !!!! pourquoi aReferenceAnnotationDatabase et pas membre FReferenceAnnotationDatabase ????
void TMappingDataFile::UpdateSmartSharedCountFromSAMFile(TCountedFragmentList & aCountedFragmentList, TReferenceAnnotationDatabase & aReferenceAnnotationDatabase)
{
  int aReadID, aCurrentReadID;

  //double aDistance; // not used

  ifstream aSAMFile(FMainMappingFileName.c_str());
  string s;

  TReadMatchInfoList aReadMatchInfoList;

  int aFragmentSize;
  bool aOKToProcessRead;

  TLocalMatchInfo aLocalMatchInfo;

  if (aSAMFile)
  {
    TSAMReader aSAMReader;
    TSAMMatch aSAMMatch;
    string line;

    aCurrentReadID = -1;

    while (getline(aSAMFile, line)) // ending '\n' is removed
    {
      aSAMReader.GetMatch(line, aSAMMatch);

      if (not aSAMMatch.Unmapped)
      {
        aReadID = atoi(aSAMMatch.QName.c_str());

        if (aCurrentReadID != aReadID)
        {
          if (aCurrentReadID != -1)
          {
            TagReadMatchInfoList(aReadMatchInfoList);

            // and now count...
            UpdateSmartSharedCountSAMReadList(aReadMatchInfoList, aCountedFragmentList);
          }
          // and finally reinitialize for next read block
          aReadMatchInfoList.Initialize(aReadID);
          aCurrentReadID = aReadID;
        }

        // continue to fill ReadMatchInfoList

        // 1- check if local alignment is acceptable
        // retrieve fragment (reference) length
        aFragmentSize = aReferenceAnnotationDatabase.GetLiteReferenceAnnotationTable().GetFragmentSize(atoi(aSAMMatch.RName.c_str()));

        if (FIsLocalAlignment) {
          aLocalMatchInfo = IsAcceptedLocalMatch(aFragmentSize, aSAMMatch.Cigar, aSAMMatch.Pos, FAlignmentLengthCutoff, FSoftClippingLengthCutoff, FKeepInternalLocalAlignment);
          aOKToProcessRead = aLocalMatchInfo.IsAccepted;
        }
        else {
		  aLocalMatchInfo.IsAccepted = true;
          aLocalMatchInfo.NewReadStartLocation = aSAMMatch.Pos;
          aLocalMatchInfo.NewReadEndLocation = aSAMMatch.Pos + aSAMMatch.AlignmentLength - 1;
          aLocalMatchInfo.SupplementaryMismatchesCount = 0;
          aOKToProcessRead = true;
        }

        if (aOKToProcessRead) {
            TReadMatchInfoItem & aReadMatchInfoItem = aReadMatchInfoList.Add();

			// 1-based offset into the forward reference strand where leftmost character of the alignment occurs
			aReadMatchInfoItem.ReadLocation = aLocalMatchInfo.NewReadStartLocation;
			aReadMatchInfoItem.ReadStrand = (aSAMMatch.Strand == 1);
			aReadMatchInfoItem.ReadMismatches = aSAMMatch.Mismatches;
			aReadMatchInfoItem.ReadEditDistance = aSAMMatch.EditDistance + aLocalMatchInfo.SupplementaryMismatchesCount;
			aReadMatchInfoItem.ReadEntryID = atoi(aSAMMatch.RName.c_str());
			aReadMatchInfoItem.ReadMappingLength = aSAMMatch.AlignmentLength + aLocalMatchInfo.SupplementaryMismatchesCount;
			aReadMatchInfoItem.ReadSAMFlag = aSAMMatch.Flag;
			aReadMatchInfoItem.ReadSAMCigar = aSAMMatch.Cigar;
			aReadMatchInfoItem.ReadSAMOpt = aSAMMatch.Opt;

			aReadMatchInfoItem.ReadDistance = (double)aSAMMatch.EditDistance + aLocalMatchInfo.SupplementaryMismatchesCount; //// CAST
			if (FIsRelativeDistance)
			  aReadMatchInfoItem.ReadDistance /= (double)(aSAMMatch.AlignmentLength + aLocalMatchInfo.SupplementaryMismatchesCount);

        }
      }
    }
    // Process last read block
    if (aCurrentReadID != -1)
    {
      TagReadMatchInfoList(aReadMatchInfoList);

      // and now count...
      UpdateSmartSharedCountSAMReadList(aReadMatchInfoList, aCountedFragmentList);
    }
    aSAMFile.close();
  }
//  FreeAndNil(aReadMatchInfoList);
}

// ------------------------------------------------------------------------------

void TMappingDataFile::UpdateSmartSharedCountSAMReadList(TReadMatchInfoList & aReadMatchInfoList, TCountedFragmentList & aCountedFragmentList)
{
  TCountedFragment * aCountedFragment = NULL;
  TCigarTypeArray aCigarArray;

  // update smart shared count with current read info if list != empty
  // 1- calculate new fraction of each smart shared count
  // = unique / sum(unique)
  long long aTotalUniqueCount = 0;
  int aAcceptedMatchCount = 0;
  size_t iMatch;

  for (iMatch = 0 ; iMatch < aReadMatchInfoList.Count(); iMatch++)
  {
	 TReadMatchInfoItem & aReadMatchInfoItem = aReadMatchInfoList.RefItem(iMatch);

    if (! aReadMatchInfoItem.IsRejected)
    {
      //aTotalUniqueCount += aCountedFragmentList.GetItemByFragmentID[aReadMatchInfoList.Items[iMatch].ReadEntryID].UniqueReadCount;
      aTotalUniqueCount += (aCountedFragmentList.GetItemByFragmentID(aReadMatchInfoItem.ReadEntryID))->GetUniqueReadCount();
      aAcceptedMatchCount++;
    }
  }

  // 2- distribute reads with smart shared count
  for (iMatch = 0 ; iMatch < aReadMatchInfoList.Count(); iMatch++)
  {
	TReadMatchInfoItem & aReadMatchInfoItem = aReadMatchInfoList.RefItem(iMatch);

    if (! aReadMatchInfoItem.IsRejected)
    {
      aCountedFragment = aCountedFragmentList.GetItemByFragmentID(aReadMatchInfoItem.ReadEntryID);

      // Extract CIGAR string of current SAM entry
      CIGARToArrayOfByte(aReadMatchInfoItem.ReadSAMCigar, aCigarArray);

      if (aTotalUniqueCount > 0)
        aCountedFragment->UpdateSmartSharedCount(aCigarArray, aReadMatchInfoItem.ReadStrand, aCountedFragment->GetUniqueReadCount() / (double)aTotalUniqueCount);
      else
        aCountedFragment->UpdateSmartSharedCount(aCigarArray, aReadMatchInfoItem.ReadStrand, 1 / (double)aAcceptedMatchCount);
    }
  }
}

// ------------------------------------------------------------------------------

void TMappingDataFile::UpdateSmartSharedCountFromBowtieFile(const TCountedFragmentList & aCountedFragmentList){;} //// TODO

// ------------------------------------------------------------------------------

void TMappingDataFile::UpdateSmartSharedCount(TCountedFragmentList & aCountedFragmentList, TReferenceAnnotationDatabase & aReferenceAnnotationDatabase)
{
	switch(FMappingFileFormat){
		case mffSAM:
			UpdateSmartSharedCountFromSAMFile(aCountedFragmentList, aReferenceAnnotationDatabase); break;
		case mffDefaultBowtie:
			UpdateSmartSharedCountFromBowtieFile(aCountedFragmentList); break;
		case mffLiteDefaultBowtie:
			UpdateSmartSharedCountFromBowtieFile(aCountedFragmentList); break;
		default: break;
	}
}

// ------------------------------------------------------------------------------

void TMappingDataFile::SetMappingIntersectionArray(TReadMatchInfoList & aReadMatchInfoList)
{
	TMappingIntersectionItem & aMappingIntersectionItem = FMappingIntersectionArray.at(aReadMatchInfoList.ReadID()); //// AT

	  // not yet counted or already counted
	  if (aMappingIntersectionItem.IntersectionStatus <= C_MAIN_MAPPING_INTERSECTION_STATUS) ////
	  { ////
	    // 1- tag kept alignment
	    for (size_t iMatch = 0; iMatch < aReadMatchInfoList.Count(); iMatch++)
	    {
	      TReadMatchInfoItem & aReadMatchInfoItem = aReadMatchInfoList.RefItem(iMatch);

	      if (aReadMatchInfoItem.ReadDistance <= FMaxDistance) // FMaxDistance = meteor.mismatches(workflow)/100
	      // in the template
	      {
	        if (aReadMatchInfoItem.ReadDistance < aMappingIntersectionItem.MinDistance)
	        {
	          // Reinitialization
	        	aMappingIntersectionItem.MatchesCount = 1;
	        	aMappingIntersectionItem.MinDistance = aReadMatchInfoItem.ReadDistance;
	        	aMappingIntersectionItem.IsRejected = false;
	        	aMappingIntersectionItem.IntersectionStatus = C_MAIN_MAPPING_INTERSECTION_STATUS;
	        }
	        else if (aReadMatchInfoItem.ReadDistance == aMappingIntersectionItem.MinDistance)
	        {
	        	aMappingIntersectionItem.MatchesCount++;
	        	aMappingIntersectionItem.IsRejected = false;
	        	aMappingIntersectionItem.IntersectionStatus = C_MAIN_MAPPING_INTERSECTION_STATUS;
	        }
	      }
	    }
	  } ////
	  TagReadMatchInfoList(aReadMatchInfoList);
}

void TMappingDataFile::CountSAMReadList(TReadMatchInfoList & aReadMatchInfoList, TCountedFragmentList & aCountedFragmentList)
{
	  TCigarTypeArray aCigarArray;
	  int aFragmentSize, aFragmentStartLocation, aFragmentEndLocation;
	  unsigned char aFragmentIsOnPlusStrand;

	  for (size_t iMatch = 0; iMatch < aReadMatchInfoList.Count(); iMatch++)
	  {
		TReadMatchInfoItem & aReadMatchInfoItem = aReadMatchInfoList.RefItem(iMatch);
	    if (! aReadMatchInfoItem.IsRejected)
	    {
	      TCountedFragment * aCountedFragment = aCountedFragmentList.GetItemByFragmentID(aReadMatchInfoItem.ReadEntryID);

	      if (aCountedFragment == NULL)
	      {
	        // Searching for next corresponding fragment in annotation table
	        FReferenceAnnotationDatabase->GetFragmentLiteInfo(aReadMatchInfoItem.ReadEntryID, aFragmentSize, aFragmentStartLocation, aFragmentEndLocation, aFragmentIsOnPlusStrand);

	        aCountedFragment = &(aCountedFragmentList.AddCountedFragment(aReadMatchInfoItem.ReadEntryID, aFragmentStartLocation, aFragmentEndLocation, aFragmentSize, aFragmentIsOnPlusStrand));
	      }

	      /*
	        Definition of the read start and end location compared to the 0-based location of the contig -1 for 1-based to 0-based conversion
	        * ReadStartLocation = Location - 1
	        * ReadEndLocation = ReadStartLocation + AlignmentLength - 1
	        * ReadEndLocation = Location + AlignmentLength - 2
	      */
	      // Extract CIGAR string of current SAM entry. Ex: 3M2I2D => aCigarArray = [ctMatch,ctMatch,ctMatch, ctInsertion,ctInsertion, ctDeletion,ctDeletion]
	      CIGARToArrayOfByte(aReadMatchInfoItem.ReadSAMCigar, aCigarArray);

	      aCountedFragment->AddRead(
	    		  aReadMatchInfoItem.ReadLocation - 1,
				  aReadMatchInfoItem.ReadLocation + aReadMatchInfoItem.ReadMappingLength - 2,
				  aReadMatchInfoItem.ReadStrand, (int)(FMappingIntersectionArray.at(aReadMatchInfoList.ReadID()).MatchesCount), aCigarArray); //// AT
	    }
	  }
}

// ------------------------------------------------------------------------------

void TMappingDataFile::TagReadMatchInfoList(TReadMatchInfoList & aReadMatchInfoList)
{
	TMappingIntersectionItem & aMappingIntersectionItem = FMappingIntersectionArray.at(aReadMatchInfoList.ReadID()); //// AT

    for (size_t iMatch = 0 ; iMatch < aReadMatchInfoList.Count(); iMatch++)
    {
      TReadMatchInfoItem & aReadMatchInfoItem = aReadMatchInfoList.RefItem(iMatch);

      if ( (!aMappingIntersectionItem.IsRejected) && (aReadMatchInfoItem.ReadDistance <= aMappingIntersectionItem.MinDistance) )
        aReadMatchInfoItem.IsRejected = false;
      else
        aReadMatchInfoItem.IsRejected = true;
    }
}

