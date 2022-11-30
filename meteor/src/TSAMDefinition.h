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

#ifndef TSAMDEFINITION_H_
#define TSAMDEFINITION_H_

#include <string>
#include <vector>

//using namespace std;

const std::string C_SAM_VERSION = "1.4";
const char C_HEADER_CHAR = '@';
const std::string C_HEADER_TAG = "@HD";
const std::string C_VERSION_TAG = "VN";
const std::string C_SORTING_ORDER_TAG = "SO";
const std::string C_REFERENCE_DICTIONARY_TAG = "@SQ";
const std::string C_REFERENCE_SEQUENCE_NAME_TAG = "SN";
const std::string C_REFERENCE_SEQUENCE_LENGTH_TAG = "LN";

const std::string C_SO_UNKNOWN_STR = "unknown";
const std::string C_SO_UNSORTED_STR = "unsorted";
const std::string C_SO_QUERYNAME_STR = "queryname";
const std::string C_SO_COORDINATE_STR = "coordinate";

// FLAG DESCRIPTION
const int C_NULL_FLAG = 0;
const int C_FLAG_UNMAPPED = 4; // hex = 0x4
const int C_FLAG_STRAND = 16; // hex = 0x10

const int C_SAM_MANDATORY_FIELD_COUNT = 11;

// OPT Fields description
const int C_OPT_FIELD_VALUE_POS = 6; // pos after last ':' (opt:type:value)
const std::string C_MD_OPT_FIELD = "MD:";
const std::string C_NM_OPT_FIELD = "NM:"; // edit distance

//// NO C++ equivalent
//TSAMCharacter = set of 'A' .. 'Z';
//
//TSAMNumber = set of '0' .. '9';

enum TOptMDStatus {mdUnknown, mdMatch, mdMismatch, mdDeletion};

enum TSortingOrderKind {soUnknown, soUnsorted, soQueryName, soCoordinate};

enum TCigarType {ctUnknown, ctMatch, ctInsertion, ctDeletion, ctSoftClipping, chHardClipping};

typedef std::vector<TCigarType> TCigarTypeArray;

typedef struct {
	int Size;
	TCigarType CigarType;
} TCigarItem;

typedef struct {
  bool IsAccepted;
  int NewReadStartLocation;
  int NewReadEndLocation;
  int SupplementaryMismatchesCount;
} TLocalMatchInfo ;

  // Others functions
void CIGARToArrayOfByte(const std::string & aCigarStr, TCigarTypeArray & aCigarArray);

TLocalMatchInfo IsAcceptedLocalMatch(
		const int aReferenceLength, const std::string & aCigarStr, const int aReadStartLocation,
		const int aAlignmentLengthCutoff, const int aSoftClippingLengthCutoff, const bool aKeepInternalClipping); // overload
//TLocalMatchInfo IsAcceptedLocalMatch(
//		const int aReferenceLength, const TCigarTypeArray & aCigarArray, const int aReadStartLocation,
//		const int aAlignmentLengthCutoff, const bool aKeepInternalClipping); // overload //// not used


#endif /* TSAMDEFINITION_H_ */
