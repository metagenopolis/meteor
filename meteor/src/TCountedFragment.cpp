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

#include "TCountedFragment.h"

TCountedFragment::TCountedFragment(
	const long long aFragmentID,
	const long long aFragmentStartLocation,
	const long long aFragmentEndLocation,
	const long long aFragmentSize,
	const unsigned char aFragmentStrand):FFragmentProfile(aFragmentSize, C_NO_READ)
{
	FFragmentID = aFragmentID;
	FFragmentStartLocation = aFragmentStartLocation;
	FFragmentEndLocation = aFragmentEndLocation;
	FFragmentSize = aFragmentSize;
	FFragmentStrand = aFragmentStrand == 1;

	FTotalReadCount = 0;
	FSharedReadCount = 0.0;
	FUniqueReadCount = 0;
	FDirectTotalReadCount = 0;
	FDirectSharedReadCount = 0.0;
	FDirectUniqueReadCount = 0;
	FTotalReadSumCoverage = 0;
	FSharedReadSumCoverage = 0.0;
	FUniqueReadSumCoverage = 0;
	FDirectTotalReadSumCoverage = 0;
	FDirectSharedReadSumCoverage = 0.0;
	FDirectUniqueReadSumCoverage = 0;

	FSmartSharedReadCount = 0.0;
	FDirectSmartSharedReadCount = 0.0;
	FSmartSharedReadSumCoverage = 0.0;
	FDirectSmartSharedReadSumCoverage = 0.0;

	FTotalReadMeanCoverage = 0.0;
	FSharedReadMeanCoverage = 0.0;
	FSmartSharedReadMeanCoverage = 0.0;
	FUniqueReadMeanCoverage = 0.0;
	FDirectTotalReadMeanCoverage = 0.0;
	FDirectSharedReadMeanCoverage = 0.0;
	FDirectSmartSharedReadMeanCoverage = 0.0;
	FDirectUniqueReadMeanCoverage = 0.0;

	FTotalReadCoverage = 0.0;
	FUniqueReadCoverage = 0.0;
	FDirectTotalReadCoverage = 0.0;
	FDirectUniqueReadCoverage = 0.0;
	FUnDirectTotalReadCoverage = 0.0;
	FUnDirectUniqueReadCoverage = 0.0;

	FUnDirectTotalReadSumCoverage = 0;
	FUnDirectSharedReadSumCoverage = 0.0;
	FUnDirectSmartSharedReadSumCoverage = 0.0;
	FUnDirectUniqueReadSumCoverage = 0;
	FUnDirectTotalReadMeanCoverage = 0.0;
	FUnDirectSharedReadMeanCoverage = 0.0;
	FUnDirectSmartSharedReadMeanCoverage = 0.0;
	FUnDirectUniqueReadMeanCoverage = 0.0;
}
TCountedFragment::TCountedFragment()/*:FFragmentProfile()*/ ////
{
	  FFragmentID = -1;
	  FFragmentStartLocation = -1;
	  FFragmentEndLocation = -1;
	  FFragmentSize = -1;
	  FFragmentStrand = true;
	  FTotalReadCount = 0;
	  FSharedReadCount = 0.0;
	  FUniqueReadCount = 0;
	  FDirectTotalReadCount = 0;
	  FDirectSharedReadCount = 0.0;
	  FDirectUniqueReadCount = 0;
	  FTotalReadSumCoverage = 0;
	  FSharedReadSumCoverage = 0.0;
	  FUniqueReadSumCoverage = 0;
	  FDirectTotalReadSumCoverage = 0;
	  FDirectSharedReadSumCoverage = 0.0;
	  FDirectUniqueReadSumCoverage = 0;

	  FSmartSharedReadCount = 0.0;
	  FDirectSmartSharedReadCount = 0.0;
	  FSmartSharedReadSumCoverage = 0.0;
	  FDirectSmartSharedReadSumCoverage = 0.0;

	  FTotalReadMeanCoverage = 0.0;
	  FSharedReadMeanCoverage = 0.0;
	  FSmartSharedReadMeanCoverage = 0.0;
	  FUniqueReadMeanCoverage = 0.0;
	  FDirectTotalReadMeanCoverage = 0.0;
	  FDirectSharedReadMeanCoverage = 0.0;
	  FDirectSmartSharedReadMeanCoverage = 0.0;
	  FDirectUniqueReadMeanCoverage = 0.0;

	  FTotalReadCoverage = 0.0;
	  FUniqueReadCoverage = 0.0;
	  FDirectTotalReadCoverage = 0.0;
	  FDirectUniqueReadCoverage = 0.0;
	  FUnDirectTotalReadCoverage = 0.0;
	  FUnDirectUniqueReadCoverage = 0.0;

	  FUnDirectTotalReadSumCoverage = 0.0;
	  FUnDirectSharedReadSumCoverage = 0.0;
	  FUnDirectSmartSharedReadSumCoverage = 0.0;
	  FUnDirectUniqueReadSumCoverage = 0;
	  FUnDirectTotalReadMeanCoverage = 0.0;
	  FUnDirectSharedReadMeanCoverage = 0.0;
	  FUnDirectSmartSharedReadMeanCoverage = 0.0;
	  FUnDirectUniqueReadMeanCoverage = 0.0;
}

TCountedFragment::~TCountedFragment() {
	// TODO Auto-generated destructor stub
}

long long TCountedFragment::GetFragmentID() const {return FFragmentID;}
long long TCountedFragment::GetFragmentSize() const {return FFragmentSize;}
long long TCountedFragment::GetFragmentStartLocation() const {return FFragmentStartLocation;}
long long TCountedFragment::GetFragmentEndLocation() const {return FFragmentEndLocation;}
long long TCountedFragment::GetTotalReadCount() const {return FTotalReadCount;}
double TCountedFragment::GetSharedReadCount() const {return FSharedReadCount;}
double TCountedFragment::GetSmartSharedReadCount() const {return FSmartSharedReadCount;}
long long TCountedFragment::GetUniqueReadCount() const {return FUniqueReadCount;}
long long TCountedFragment::GetDirectTotalReadCount() const {return FDirectTotalReadCount;}
double TCountedFragment::GetDirectSharedReadCount() const {return FDirectSharedReadCount;}
double TCountedFragment::GetDirectSmartSharedReadCount() const {return FDirectSmartSharedReadCount;}
long long TCountedFragment::GetDirectUniqueReadCount() const {return FDirectUniqueReadCount;}

long long TCountedFragment::GetTotalReadSumCoverage() const {return FTotalReadSumCoverage;}
double TCountedFragment::GetSharedReadSumCoverage() const {return FSharedReadSumCoverage;}
double TCountedFragment::GetSmartSharedReadSumCoverage() const {return FSmartSharedReadSumCoverage;}
long long TCountedFragment::GetUniqueReadSumCoverage() const {return FUniqueReadSumCoverage;}
long long TCountedFragment::GetDirectTotalReadSumCoverage() const {return FDirectTotalReadSumCoverage;}
double TCountedFragment::GetDirectSharedReadSumCoverage() const {return FDirectSharedReadSumCoverage;}
double TCountedFragment::GetDirectSmartSharedReadSumCoverage() const {return FDirectSmartSharedReadSumCoverage;}
long long TCountedFragment::GetDirectUniqueReadSumCoverage() const {return FDirectUniqueReadSumCoverage;}

double TCountedFragment::GetTotalReadCoverage() const {return FTotalReadCoverage;}
double TCountedFragment::GetUniqueReadCoverage() const {return FUniqueReadCoverage;}
double TCountedFragment::GetDirectTotalReadCoverage() const {return FDirectTotalReadCoverage;}
double TCountedFragment::GetDirectUniqueReadCoverage() const {return FDirectUniqueReadCoverage;}
double TCountedFragment::GetUnDirectTotalReadCoverage() const {return FUnDirectTotalReadCoverage;}
double TCountedFragment::GetUnDirectUniqueReadCoverage() const {return FUnDirectUniqueReadCoverage;}

double TCountedFragment::GetTotalReadMeanCoverage() const {return FTotalReadMeanCoverage;}
double TCountedFragment::GetSharedReadMeanCoverage() const {return FSharedReadMeanCoverage;}
double TCountedFragment::GetSmartSharedReadMeanCoverage() const {return FSmartSharedReadMeanCoverage;}
double TCountedFragment::GetUniqueReadMeanCoverage() const {return FUniqueReadMeanCoverage;}
double TCountedFragment::GetDirectTotalReadMeanCoverage() const {return FDirectTotalReadMeanCoverage;}
double TCountedFragment::GetDirectSharedReadMeanCoverage() const {return FDirectSharedReadMeanCoverage;}
double TCountedFragment::GetDirectSmartSharedReadMeanCoverage() const {return FDirectSmartSharedReadMeanCoverage;}
double TCountedFragment::GetDirectUniqueReadMeanCoverage() const {return FDirectUniqueReadMeanCoverage;}

long long TCountedFragment::GetUnDirectTotalReadSumCoverage() const {return FUnDirectTotalReadSumCoverage;}
double TCountedFragment::GetUnDirectSharedReadSumCoverage() const {return FUnDirectSharedReadSumCoverage;}
double TCountedFragment::GetUnDirectSmartSharedReadSumCoverage() const {return FUnDirectSmartSharedReadSumCoverage;}
long long TCountedFragment::GetUnDirectUniqueReadSumCoverage() const {return FUnDirectUniqueReadSumCoverage;}

double TCountedFragment::GetUnDirectTotalReadMeanCoverage() const {return FUnDirectTotalReadMeanCoverage;}
double TCountedFragment::GetUnDirectSharedReadMeanCoverage() const {return FUnDirectSharedReadMeanCoverage;}
double TCountedFragment::GetUnDirectSmartSharedReadMeanCoverage() const {return FUnDirectSmartSharedReadMeanCoverage;}
double TCountedFragment::GetUnDirectUniqueReadMeanCoverage() const {return FUnDirectUniqueReadMeanCoverage;}

long long TCountedFragment::GetUnDirectTotalReadCount() const {return FTotalReadCount - FDirectTotalReadCount;}
double TCountedFragment::GetUnDirectSharedReadCount() const {return FSharedReadCount - FDirectSharedReadCount;}
double TCountedFragment::GetUnDirectSmartSharedReadCount() const {return FSmartSharedReadCount - FDirectSmartSharedReadCount;}
long long TCountedFragment::GetUnDirectUniqueReadCount() const {return FUniqueReadCount - FDirectUniqueReadCount;}

const TByteArray & TCountedFragment::GetFragmentProfile() const {return FFragmentProfile;}

// ------------------------------------------------------------------------------

void TCountedFragment::UpdateSmartSharedCount(const TCigarTypeArray & aCigarArray, const bool aReadStrand, const double aSmartSharedCount)
{

  int aMatchCount;
  // number of aligned position (including mismatches)

  aMatchCount = 0;
  for (size_t i = 0; i < aCigarArray.size(); i++)
	if (aCigarArray[i] == ctMatch) aMatchCount++;

  FSmartSharedReadCount = FSmartSharedReadCount + aSmartSharedCount;
  FSmartSharedReadSumCoverage = FSmartSharedReadSumCoverage + aSmartSharedCount	* aMatchCount;

  if (aReadStrand == FFragmentStrand)
  {
	FDirectSmartSharedReadCount = FDirectSmartSharedReadCount + aSmartSharedCount;
	FDirectSmartSharedReadSumCoverage = FDirectSmartSharedReadSumCoverage + aSmartSharedCount * aMatchCount;
  }
  else
	FUnDirectSmartSharedReadSumCoverage = FUnDirectSmartSharedReadSumCoverage + aSmartSharedCount * aMatchCount;
}

// ------------------------------------------------------------------------------

void TCountedFragment::AddRead(
	const int aReadStartLocation,
	const int aReadEndLocation,
	const bool aReadStrand,
	const int aReadMatchesCount,
	const TCigarTypeArray & aCigarArray)
{
	  int aFragmentPos;
	  double aTmpSharedReadCount;
	  //char aStateLocation; //// not used

	  // number of aligned position (including mismatches)

	  int aMatchCount = 0;
	  for (size_t i = 0; i < aCigarArray.size(); i++)
	    if (aCigarArray[i] == ctMatch) aMatchCount++;

	  //int aMappedReadLength = aReadEndLocation - aReadStartLocation + 1; //// not used
	  aTmpSharedReadCount = 1 / (double)aReadMatchesCount;
	  FTotalReadCount++;
	  FTotalReadSumCoverage += aMatchCount;
	  FSharedReadCount = FSharedReadCount + aTmpSharedReadCount;
	  FSharedReadSumCoverage = FSharedReadSumCoverage + aTmpSharedReadCount * aMatchCount;

	  if (aReadStrand == FFragmentStrand)
	  {
	    FDirectTotalReadCount++;
	    FDirectTotalReadSumCoverage += aMatchCount;
	    FDirectSharedReadCount = FDirectSharedReadCount + aTmpSharedReadCount;
	    FDirectSharedReadSumCoverage = FDirectSharedReadSumCoverage + aTmpSharedReadCount * aMatchCount;

	    if (aReadMatchesCount > 1)
	    {
	      aFragmentPos = aReadStartLocation;

	      for (size_t aReadPos = 0; aReadPos < aCigarArray.size(); aReadPos++)
	      {
	        switch(aCigarArray[aReadPos]) {
	          case ctUnknown:
	            aFragmentPos++; break;
	          case ctMatch:
	              FFragmentProfile.at(aFragmentPos) = FFragmentProfile.at(aFragmentPos) | C_MULTIPLE_DIRECT_READ; //// AT
	              aFragmentPos++; break;
	          case ctInsertion:
	        	  break; // nothing
	          case ctDeletion:
	              aFragmentPos++; break;
	          case ctSoftClipping:
	        	  break; // nothing
	        	default:
	        		break;
	        }
	      }
	    }
	    else
	    {
	      FUniqueReadCount++;
	      FUniqueReadSumCoverage += aMatchCount;
	      FDirectUniqueReadCount++;
	      FDirectUniqueReadSumCoverage += aMatchCount;

	      aFragmentPos = aReadStartLocation;

	      for (size_t aReadPos = 0; aReadPos < aCigarArray.size(); aReadPos++)
	      {
	        switch(aCigarArray[aReadPos]) {
	          case ctUnknown:
	            aFragmentPos++; break;
	          case ctMatch:
	              FFragmentProfile.at(aFragmentPos) = FFragmentProfile.at(aFragmentPos) | C_UNIQUE_DIRECT_READ; //// AT
	              aFragmentPos++; break;
	          case ctInsertion:
	        	  break; // nothing
	          case ctDeletion:
	              aFragmentPos++; break;
	          case ctSoftClipping:
	          	  break; // nothing
	          default:
	        		break;
	        }
	      }
	    }
	  }
	  else
	  {
	    FUnDirectTotalReadSumCoverage += aMatchCount;
	    FUnDirectSharedReadSumCoverage = FUnDirectSharedReadSumCoverage + aTmpSharedReadCount * aMatchCount;
	    if (aReadMatchesCount > 1)
	    {
	      aFragmentPos = aReadStartLocation;

	      for (size_t aReadPos = 0; aReadPos < aCigarArray.size(); aReadPos++)
	      {
	        switch(aCigarArray[aReadPos]) {
	          case ctUnknown:
	            aFragmentPos++; break;
	          case ctMatch:
	              FFragmentProfile.at(aFragmentPos) = FFragmentProfile.at(aFragmentPos) | C_MULTIPLE_UNDIRECT_READ; //// AT
	              aFragmentPos++; break;
	          case ctInsertion:
	        	  break; // nothing
	          case ctDeletion:
	              aFragmentPos++; break;
	          case ctSoftClipping:
	              break; // nothing
	          default:
	        		break;
	        }
	      }
	    }
	    else
	    {
	      FUniqueReadCount++;
	      FUniqueReadSumCoverage += aMatchCount;
	      FUnDirectUniqueReadSumCoverage += aMatchCount;

	      aFragmentPos = aReadStartLocation;

	      for (size_t aReadPos = 0; aReadPos < aCigarArray.size(); aReadPos++)
	      {
	        switch(aCigarArray[aReadPos]) {
	          case ctUnknown:
	            aFragmentPos++; break;
	          case ctMatch:
	            {
	              FFragmentProfile.at(aFragmentPos) = FFragmentProfile.at(aFragmentPos) | C_UNIQUE_UNDIRECT_READ; //// AT
	              aFragmentPos++; break;
	            }
	          case ctInsertion:
	        	  break; // nothing
	          case ctDeletion:
	              aFragmentPos++; break;
	          case ctSoftClipping:
	          	  break; // nothing
	        	default:
	        		break;
	        }
	      }
	    }
	  }
}

// ------------------------------------------------------------------------------

void TCountedFragment::Finalize() {

	  int aTotalReadCoverage=0, aUniqueReadCoverage=0;
	  int aDirectTotalReadCoverage=0, aDirectUniqueReadCoverage=0;
	  int aUnDirectTotalReadCoverage=0, aUnDirectUniqueReadCoverage=0;

	  bool aIsMultipleDirectPosition, aIsMultipleUndirectPosition, aIsUniqueDirectPosition, aIsUniqueUndirectPosition;

	  unsigned char aProfileValue;

	  for (int i = 0; i < FFragmentSize; i++)
	  {
	    aProfileValue = FFragmentProfile.at(i); //// AT

	    if (aProfileValue != C_NO_READ)
	    {
	      aIsMultipleDirectPosition = (aProfileValue == (aProfileValue | C_MULTIPLE_DIRECT_READ)); //// bitwise or a verifier
	      aIsMultipleUndirectPosition = (aProfileValue == (aProfileValue | C_MULTIPLE_UNDIRECT_READ));
	      aIsUniqueDirectPosition = (aProfileValue == (aProfileValue | C_UNIQUE_DIRECT_READ));
	      aIsUniqueUndirectPosition = (aProfileValue == (aProfileValue | C_UNIQUE_UNDIRECT_READ));

	      aTotalReadCoverage++;

	      if (aIsUniqueDirectPosition || aIsUniqueUndirectPosition) aUniqueReadCoverage++;
	      if (aIsMultipleDirectPosition) aDirectTotalReadCoverage++;
	      if (aIsMultipleUndirectPosition) aUnDirectTotalReadCoverage++;
	      if (aIsUniqueDirectPosition) aDirectUniqueReadCoverage++;
	      if (aIsUniqueUndirectPosition) aUnDirectUniqueReadCoverage++;
	    }
	  }

	  // reajust: total=multiple+unique
	  aDirectTotalReadCoverage += aDirectUniqueReadCoverage;
	  aUnDirectTotalReadCoverage += aUnDirectUniqueReadCoverage;

	  FTotalReadCoverage = aTotalReadCoverage / (double)FFragmentSize;
	  FUniqueReadCoverage = aUniqueReadCoverage / (double)FFragmentSize;
	  FDirectTotalReadCoverage = aDirectTotalReadCoverage / (double)FFragmentSize;
	  FDirectUniqueReadCoverage = aDirectUniqueReadCoverage / (double)FFragmentSize;
	  FUnDirectTotalReadCoverage = aUnDirectTotalReadCoverage / (double)FFragmentSize;
	  FUnDirectUniqueReadCoverage = aUnDirectUniqueReadCoverage / (double)FFragmentSize;

	  FTotalReadMeanCoverage = FTotalReadSumCoverage / (double)FFragmentSize;
	  FSharedReadMeanCoverage = FSharedReadSumCoverage / (double)FFragmentSize;
	  FSmartSharedReadMeanCoverage = FSmartSharedReadSumCoverage / (double)FFragmentSize;
	  FUniqueReadMeanCoverage = FUniqueReadSumCoverage / (double)FFragmentSize;

	  FDirectTotalReadMeanCoverage = FDirectTotalReadSumCoverage / (double)FFragmentSize;
	  FDirectSharedReadMeanCoverage = FDirectSharedReadSumCoverage / (double)FFragmentSize;
	  FDirectSmartSharedReadMeanCoverage = FDirectSmartSharedReadSumCoverage / (double)FFragmentSize;
	  FDirectUniqueReadMeanCoverage = FDirectUniqueReadSumCoverage / (double)FFragmentSize;

	  FUnDirectTotalReadMeanCoverage = FUnDirectTotalReadSumCoverage / (double)FFragmentSize;
	  FUnDirectSharedReadMeanCoverage = FUnDirectSharedReadSumCoverage /(double)FFragmentSize;
	  FUnDirectSmartSharedReadMeanCoverage = FUnDirectSmartSharedReadSumCoverage / (double)FFragmentSize;
	  FUnDirectUniqueReadMeanCoverage = FUnDirectUniqueReadSumCoverage / (double)FFragmentSize;
}

