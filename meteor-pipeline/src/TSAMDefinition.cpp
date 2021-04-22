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

#include "TSAMDefinition.h"
#include <iostream>
#include <cstdlib>

using namespace std;

/// /////////////////////////////////////////////////////////////////////////////

void CIGARToArrayOfByte(const string & aCigarStr, TCigarTypeArray & aCigarArray) {
/*
  this function calculates the alignment length by suming:
  matches/mismatches (M), deletion (D), insertion (I), soft clipping (S)
  in CIGAR string
  example : if CIGAR string = "3M1I3M1D5M" (3 match, 1 insertion, 3 match, 1 deletion, 5 match)
         => aCigarArray : {ctMatch,ctMatch,ctMatch, ctInsertion, ctMatch,ctMatch,ctMatch, ctDeletion, ctMatch,ctMatch,ctMatch,ctMatch,ctMatch}
*/
	char c;

	TCigarItem aCigarItemArray[500];
	int aCigarItemArrayLength;

	// 1- calculate alignment length
	int aLengthAlignment = 0;
	int aCigarStrLength = aCigarStr.size();
	aCigarItemArrayLength = 0;

	int aCigarItemSize = 0;
	// pour chaque caractere de aCigarStrLength
	// code ascii
	for (int i = 0; i<aCigarStrLength; i++)	{

		c = aCigarStr[i];

		//if (c < ':'){ // ascii code comparison (':' comes after '9')
		if (c <= '9'){ // ascii code comparison
		  // x10 to move to the next number rank. e.g.: 19 = 1x10 + 9
		  aCigarItemSize = aCigarItemSize * 10 + (c - '0'); // c - 0 => char to number
		}
		else {
			aLengthAlignment += aCigarItemSize;
			aCigarItemArray[aCigarItemArrayLength].Size = aCigarItemSize;
			switch(c){
				case 'M':
					aCigarItemArray[aCigarItemArrayLength].CigarType = ctMatch;	break;
				case 'S':
					aCigarItemArray[aCigarItemArrayLength].CigarType = ctSoftClipping;	break;
				case 'I':
					aCigarItemArray[aCigarItemArrayLength].CigarType = ctInsertion;	break;
				case 'D':
					aCigarItemArray[aCigarItemArrayLength].CigarType = ctDeletion; break;
				default:
					aCigarItemArray[aCigarItemArrayLength].CigarType = ctUnknown; break;
			}
			aCigarItemSize = 0;
			aCigarItemArrayLength++;
		}
	}
	// 2- transform CIGAR str to CIGAR array
	aCigarArray.clear();
	aCigarArray.resize(aLengthAlignment);

	int k, j = 0;
	for (int i = 0; i < aCigarItemArrayLength; i++) {
		for (k = j; k < j + aCigarItemArray[i].Size; k++) {
			aCigarArray.at(k) = aCigarItemArray[i].CigarType; //// AT
		}
		j = k;
	}
}

// ------------------------------------------------------------------------------

TLocalMatchInfo IsAcceptedLocalMatch(
		const int aReferenceLength, const string & aCigarStr, const int aReadStartLocation,
		const int aAlignmentLengthCutoff, const int aSoftClippingLengthCutoff, const bool aKeepInternalClipping)
{
  TCigarTypeArray aCigarArray;
  TLocalMatchInfo Result;
  size_t i;

  CIGARToArrayOfByte(aCigarStr, aCigarArray);
  int aAlignmentSize = 0;

  for (i = 0; i < aCigarArray.size(); i++) {
    const TCigarType & aCigarType = aCigarArray[i];
    if (aCigarType == ctMatch || aCigarType == ctInsertion || aCigarType == ctDeletion) aAlignmentSize++;
  }
  int aReadEndLocation = aReadStartLocation + aAlignmentSize - 1;

  // initialize result
  Result.IsAccepted = true;
  Result.NewReadStartLocation = aReadStartLocation;
  Result.NewReadEndLocation = aReadEndLocation;
  Result.SupplementaryMismatchesCount = 0;

  bool aIsLocalAlignment = ( aAlignmentSize < (int)aCigarArray.size() );

  if (aIsLocalAlignment)
  {
    bool aIsInternalClipping = false;

    // 1- compute left clipping and check if internal clipping
    if (aCigarArray[0] == ctSoftClipping)
      //aIsInternalClipping = InRange(aReadStartLocation - 1, 1, aReferenceLength); // 1 <= aReadStartLocation -1 <= aReferenceLength
      aIsInternalClipping = (aReadStartLocation > 1 && aReadStartLocation - 1 <= aReferenceLength);
    // compute clipping size
    int aLeftClippingSize = 0;
    i = 0;
    while (aCigarArray[i] == ctSoftClipping)
    {
      aLeftClippingSize++;
      i++;
    }
    //Result.NewReadStartLocation = max(aReadStartLocation - aClippingSize, 1);
    if ( (Result.NewReadStartLocation = aReadStartLocation - aLeftClippingSize) < 1 ) Result.NewReadStartLocation = 1;

    // 2- compute right clipping and check if internal clipping
    if (aCigarArray.back() == ctSoftClipping)
      //aIsInternalClipping = ( aIsInternalClipping || (aReadEndLocation+1 >= 1 && aReadEndLocation+1 <= aReferenceLength) );
      aIsInternalClipping = ( aIsInternalClipping || (aReadEndLocation >= 0 && aReadEndLocation < aReferenceLength) );
    // compute clipping size
    int aRightClippingSize = 0;
    i = aCigarArray.size() - 1;
    while (aCigarArray.at(i) == ctSoftClipping) //// AT vector::at() raises an exception when i is out of bound
    {
      aRightClippingSize++;
      i--;
    }
    //Result.NewReadEndLocation = min(aReadEndLocation + aClippingSize, aReferenceLength);
    if ( (Result.NewReadEndLocation = aReadEndLocation + aRightClippingSize) > aReferenceLength) Result.NewReadEndLocation = aReferenceLength;

    Result.SupplementaryMismatchesCount = (aReadStartLocation - Result.NewReadStartLocation) + (Result.NewReadEndLocation - aReadEndLocation);

    /*
      writeln('alignment size = ', aAlignmentSize);
      writeln('read start location = ', aReadStartLocation);
      writeln('read end location = ', aReadEndLocation);
      writeln('fragment size = ', aReferenceLength);
      writeln('cigar = ', aCigarStr);
      writeln('internal clipping = ', aIsInternalClipping);
      writeln('new read start location = ', Result.NewReadStartLocation);
      writeln('new read end location = ', Result.NewReadEndLocation);
      writeln('supplementary mismatches = ', Result.SupplementaryMismatchesCount);
      Readln;
    */

    // 3- Test acceptable clipping
    //if ( (aAlignmentSize >= aAlignmentLengthCutoff) && ( (! aIsInternalClipping) || aKeepInternalClipping ) )
    if ( (aAlignmentSize >= aAlignmentLengthCutoff) && ( (! aIsInternalClipping) || aKeepInternalClipping ) && aLeftClippingSize+aRightClippingSize <= aSoftClippingLengthCutoff)
        Result.IsAccepted = true;
    else
        Result.IsAccepted = false;
  } // end of if aIsLocalAlignment
  else Result.IsAccepted = true;

  return Result;
}

// ------------------------------------------------------------------------------

