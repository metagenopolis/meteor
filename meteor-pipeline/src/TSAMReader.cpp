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

#include "TSAMReader.h"
#include "MeteorConstant.h"
#include "TSAMDefinition.h"
#include <cstring>

using namespace std;

TSAMReader::TSAMReader(): FSAMArray(C_SAM_MANDATORY_FIELD_COUNT+1, ""), FSAMOptField()
{
	//FSAMLine := TStringList.Create; // not used
}
TSAMReader::TSAMReader(const string & aSAMFileName): FSAMFile(aSAMFileName.c_str()), FSAMArray(C_SAM_MANDATORY_FIELD_COUNT+1, ""), FSAMOptField()
{
	////FSAMFile := TStreamReader.Create(aSAMFileName, TEncoding.UTF8, False, 2048);
	//FSAMLine := TStringList.Create; // not used
}

TSAMReader::~TSAMReader() {
	if (FSAMFile.is_open()) FSAMFile.close();
}

// ------------------------------------------------------------------------------
//// TODO, changer le code retour en uint ou uchar:
//// 0 : problem
//// 1 : ok without optional fields
//// 2 : ok
bool TSAMReader::GetMatch_old(const string & s, TSAMMatch & aSAMMatch) {
// return true if retrieve a new match entry

  //// FSAMArray = s.Split([#9, ' ']); // #9 = '\t'
  FSAMArray.clear();
  splitString(FSAMArray, s, '\t');

  //// 1 - essayer s.find_last_of("MD:Z:") et s.find_last_of("NM:i:") au lieu de faire un split
  //// 2 - essayer les regex (c++11) : MD:Z:(.*) et NM:i:(.*) au lieu de faire un split

  int aSAMArrayLength = FSAMArray.size();

  if (aSAMArrayLength < C_SAM_MANDATORY_FIELD_COUNT) return false;

  // extract optional field info
  aSAMMatch.Mismatches = -1; // -1 because unknown
  aSAMMatch.EditDistance = -1; // -1 because unknown

  if (aSAMArrayLength > C_SAM_MANDATORY_FIELD_COUNT) // contains opt
  {
//      extract all different opt fields
//      for an opt entry: 3 lines in the FSAMLine list: field/type/value
//      field, type and value are ':' separated in the SAM format

    // extract MD field
    // get index
    int idx = -1;
    for (int i = aSAMArrayLength - 1; i >= C_SAM_MANDATORY_FIELD_COUNT; i--)
    {
//      string aOptStr = FSAMArray[i].substr(0, 3);
//      if (FSAMArray[i].substr(0, 3) == C_MD_OPT_FIELD) { idx = i; break; }
      if (strncmp(FSAMArray.at(i).c_str(),C_MD_OPT_FIELD.c_str(),3) == 0) { idx = i; break; } //// AT
//      if (FSAMArray[i].substr(0, 3).compare(C_MD_OPT_FIELD) == 0) { idx = i; break; }
    }
    if (idx != -1)                                                                               // 0    5
      aSAMMatch.Mismatches = GetMismatchesCount(FSAMArray.at(idx).substr(C_OPT_FIELD_VALUE_POS-1)); // MD:Z:value
      //// tester surcharge  GetMismatchesCount(const char * aOptMDField)
      //// => appel : GetMismatchesCount( &(FSAMArray[idx][C_OPT_FIELD_VALUE_POS-1]) );

    // extract NM field
    // get index
    idx = -1;
    for (int i = C_SAM_MANDATORY_FIELD_COUNT; i < aSAMArrayLength ; i++)
    {
//      if (FSAMArray[i].substr(0, 3) == C_NM_OPT_FIELD) { idx = i; break; }
    	if (strncmp(FSAMArray.at(i).c_str(),C_NM_OPT_FIELD.c_str(),3) == 0) { idx = i; break; } //// AT
    }
    if (idx != -1)
    {
//      aSAMMatch.EditDistance = atoi(FSAMArray[idx].substr(C_OPT_FIELD_VALUE_POS-1).c_str());
      aSAMMatch.EditDistance = atoi(&(FSAMArray.at(idx).at(C_OPT_FIELD_VALUE_POS-1))); ////
    }
  }
	aSAMMatch.QName = FSAMArray[0];
//	aSAMMatch.Flag = atoi(FSAMArray[1].c_str()); ////
	aSAMMatch.Flag = atoi(&(FSAMArray[1][0]));
	aSAMMatch.RName = FSAMArray[2];
//	aSAMMatch.Pos = atoi(FSAMArray[3].c_str()); ////
	aSAMMatch.Pos = atoi(&(FSAMArray[3][0]));
	aSAMMatch.Cigar = FSAMArray[5];

	//// A que quoi caaa ????
    if ((aSAMMatch.Flag & C_FLAG_STRAND) == C_FLAG_STRAND) aSAMMatch.Strand = 0;
    else aSAMMatch.Strand = 1;
    aSAMMatch.Unmapped = ( (aSAMMatch.Flag & C_FLAG_UNMAPPED) == C_FLAG_UNMAPPED );

    if (aSAMMatch.Unmapped) aSAMMatch.AlignmentLength = 0;
    else aSAMMatch.AlignmentLength = GetAlignmentLength(aSAMMatch.Cigar);

  return true;
}

// ------------------------------------------------------------------------------
//// TODO, changer le code retour en uint ou uchar:
//// 0 : problem
//// 1 : ok without optional fields
//// 2 : ok
bool TSAMReader::GetMatch(const string & s, TSAMMatch & aSAMMatch) {
// return true if retrieve a new match entry

	//// FSAMArray = s.Split([#9, ' ']); // #9 = '\t'

	// last element of FSAMArray is optional fields (unsplitted).
	int aSAMArrayLength = splitString(FSAMArray, s, '\t', C_SAM_MANDATORY_FIELD_COUNT+1);

	//// 1 - essayer s.find_last_of("MD:Z:") et s.find_last_of("NM:i:") au lieu de faire un split
	//// 2 - essayer les regex (c++11) : MD:Z:(.*) et NM:i:(.*) au lieu de faire un split

	if (aSAMArrayLength < C_SAM_MANDATORY_FIELD_COUNT) return false;

	// extract optional fields info
	aSAMMatch.Mismatches = -1; // -1 because unknown
	aSAMMatch.EditDistance = -1; // -1 because unknown

	if (aSAMArrayLength > C_SAM_MANDATORY_FIELD_COUNT) // contains optional fields
	{
	  // one opt entry format: field:type:value

		const string & aOptFields = FSAMArray.back();

		///// extract MD field (for mismatch count)
		size_t begPos = aOptFields.rfind(C_MD_OPT_FIELD);
		size_t aStrSize = 0;
																			 //0000000111
		if (begPos != string::npos) { // MD found                              3456789012
			begPos += C_OPT_FIELD_VALUE_POS - 1; // jump to MD value position (MD:Z:value)
			aStrSize = aOptFields.find_first_of('\t', begPos);
			if (aStrSize == string::npos) aStrSize = aOptFields.size(); // MD is the last field
			aStrSize -= begPos;

			aSAMMatch.Mismatches = GetMismatchesCount(aOptFields.substr(begPos, aStrSize));
			//// tester surcharge  GetMismatchesCount(const char * aOptMDField)
			//// => appel : GetMismatchesCount(&FSAMArray[idx][C_OPT_FIELD_VALUE_POS-1]);
		}

		///// extract NM field (edit distance)
		begPos = aOptFields.find(C_NM_OPT_FIELD);
		if (begPos != string::npos) { // NM found                              3456789012
			begPos += C_OPT_FIELD_VALUE_POS - 1; // jump to NM value position (NM:i:value)
			aStrSize = aOptFields.find_first_of('\t', begPos);
			if (aStrSize == string::npos) aStrSize = aOptFields.size(); // NM is the last field
			aStrSize -= begPos;

			aSAMMatch.EditDistance = atoi(aOptFields.substr(begPos, aStrSize).c_str());
		}
	}
	aSAMMatch.QName = FSAMArray.at(0); //// AT
	//	aSAMMatch.Flag = atoi(FSAMArray[1].c_str()); ////
	aSAMMatch.Flag = atoi(&(FSAMArray.at(1).at(0))); //// AT
	aSAMMatch.RName = FSAMArray.at(2); //// AT
	//	aSAMMatch.Pos = atoi(FSAMArray[3].c_str()); ////
	aSAMMatch.Pos = atoi(&(FSAMArray.at(3).at(0))); //// AT
	aSAMMatch.Cigar = FSAMArray.at(5); //// AT

	// https://broadinstitute.github.io/picard/explain-flags.html
	// 0: read reverse strand
	if ((aSAMMatch.Flag & C_FLAG_STRAND) == C_FLAG_STRAND) aSAMMatch.Strand = 0;
	else aSAMMatch.Strand = 1;
	aSAMMatch.Unmapped = ( (aSAMMatch.Flag & C_FLAG_UNMAPPED) == C_FLAG_UNMAPPED );

	///// get alignment length
	if (aSAMMatch.Unmapped) aSAMMatch.AlignmentLength = 0;
	else aSAMMatch.AlignmentLength = GetAlignmentLength(aSAMMatch.Cigar);

	return true;
}

// ------------------------------------------------------------------------------

//// TODO: changer le code retour en uint ou uchar
bool TSAMReader::FastGetMatch_old(string & s, TSAMMatch & aSAMMatch, vector<StringRef> & SAMArray) const {
// return true if retrieve a new match entry

  //// FSAMArray = s.Split([#9, ' ']); // #9 = '\t'
  SAMArray.clear();
  stringToStringRefVector(s, SAMArray, '\t');

  //// 1 - essayer s.find_last_of("MD:Z:") et s.find_last_of("NM:i:") au lieu de faire un split
  //// 2 - essayer les regex (c++11) : MD:Z:(.*) et NM:i:(.*) au lieu de faire un split

  int aSAMArrayLength = SAMArray.size();

  if (aSAMArrayLength < C_SAM_MANDATORY_FIELD_COUNT) return false;

  // extract optional field info
  aSAMMatch.Mismatches = -1; // -1 because unknown
  aSAMMatch.EditDistance = -1; // -1 because unknown

  if (aSAMArrayLength > C_SAM_MANDATORY_FIELD_COUNT) // contains opt
  {
//      extract all different opt fields
//      for an opt entry: 3 lines in the FSAMLine list: field/type/value
//      field, type and value are ':' separated in the SAM format


    // extract MD field
    // get index
    int idx = -1;
    for (int i = aSAMArrayLength - 1; i >= C_SAM_MANDATORY_FIELD_COUNT; i--)
    {
//      string aOptStr = FSAMArray[i].substr(0, 3);
//      if (FSAMArray[i].substr(0, 3) == C_MD_OPT_FIELD) { idx = i; break; }
      if (strncmp(SAMArray[i].begin(),C_MD_OPT_FIELD.c_str(),3) == 0) { idx = i; break; }
//      if (FSAMArray[i].substr(0, 3).compare(C_MD_OPT_FIELD) == 0) { idx = i; break; }
    }
    if (idx != -1)                                                                               // 0    5
      aSAMMatch.Mismatches = GetMismatchesCount(SAMArray[idx].begin()+C_OPT_FIELD_VALUE_POS-1); // MD:Z:value
      //// tester surcharge  GetMismatchesCount(const char * aOptMDField)
      //// => appel : GetMismatchesCount(&FSAMArray[idx][C_OPT_FIELD_VALUE_POS-1]);

    // extract NM field
    // get index
    idx = -1;
    for (int i = C_SAM_MANDATORY_FIELD_COUNT; i < aSAMArrayLength ; i++)
    {
//      if (FSAMArray[i].substr(0, 3) == C_NM_OPT_FIELD) { idx = i; break; }
    	if (strncmp(SAMArray[i].begin(),C_NM_OPT_FIELD.c_str(),3) == 0) { idx = i; break; }
    }
    if (idx != -1)
    {
      aSAMMatch.EditDistance = atoi(SAMArray[idx].begin()+C_OPT_FIELD_VALUE_POS-1);
    }
  }
	aSAMMatch.QName = SAMArray[0].begin();
	aSAMMatch.Flag = atoi(SAMArray[1].begin());
	aSAMMatch.RName = SAMArray[2].begin();
	aSAMMatch.Pos = atoi(SAMArray[3].begin());
	aSAMMatch.Cigar = SAMArray[5].begin();

	//// A que quoi caaa ????
    if ((aSAMMatch.Flag & C_FLAG_STRAND) == C_FLAG_STRAND) aSAMMatch.Strand = 0;
    else aSAMMatch.Strand = 1;
    aSAMMatch.Unmapped = ( (aSAMMatch.Flag & C_FLAG_UNMAPPED) == C_FLAG_UNMAPPED );

    if (aSAMMatch.Unmapped) aSAMMatch.AlignmentLength = 0;
    else aSAMMatch.AlignmentLength = GetAlignmentLength(aSAMMatch.Cigar);

  return true;
}

// ------------------------------------------------------------------------------

//// TODO: changer le code retour en uint ou uchar
bool TSAMReader::FastGetMatch(string & s, TSAMMatch & aSAMMatch, vector<StringRef> & SAMArray) const {
// return true if retrieve a new match entry

	int aSAMArrayLength = stringToStringRefVector(s, SAMArray, C_SAM_MANDATORY_FIELD_COUNT+1, '\t');

	if (aSAMArrayLength < C_SAM_MANDATORY_FIELD_COUNT) return false;

	// extract optional field info
	aSAMMatch.Mismatches = -1; // -1 because unknown
	aSAMMatch.EditDistance = -1; // -1 because unknown

	if (aSAMArrayLength > C_SAM_MANDATORY_FIELD_COUNT) // contains opt
	{
	// one opt entry format: field:type:value

		const StringRef & aOptFields = SAMArray.back();

		// extract MD field
		const char* begPtr = strstr(aOptFields.begin(), C_MD_OPT_FIELD.c_str());
		const char* endPtr = NULL;
																			 //0000000111
		if (begPtr != NULL) { // MD found                              3456789012
			begPtr += C_OPT_FIELD_VALUE_POS - 1; // jump to MD value position (MD:Z:value)
			endPtr = strchr(begPtr, '\t');
			if (endPtr == NULL) endPtr = begPtr + aOptFields.size(); // MD is the last field

			aSAMMatch.Mismatches = GetMismatchesCount(string(begPtr, endPtr));
			//// tester surcharge  GetMismatchesCount(const char * aOptMDField)
			//// => appel : GetMismatchesCount(&FSAMArray[idx][C_OPT_FIELD_VALUE_POS-1]);
		}
		// extract NM field
		begPtr = strstr(aOptFields.begin(), C_NM_OPT_FIELD.c_str());
		if (begPtr != NULL) { // NM found                              3456789012
			begPtr += C_OPT_FIELD_VALUE_POS - 1; // jump to NM value position (NM:i:value)
			endPtr = strchr(begPtr, '\t');
			if (endPtr == NULL) endPtr = begPtr + aOptFields.size(); // NM is the last field

			aSAMMatch.EditDistance = atoi(begPtr);
			// atoi (ascii to int) converts the string pointed by begPtr to integer.
			// it stops at the first character which is not a valid digit (string end or field delimiter character)
		}
	}
	aSAMMatch.QName = SAMArray[0].begin();
	aSAMMatch.Flag = atoi(SAMArray[1].begin());
	aSAMMatch.RName = SAMArray[2].begin();
	aSAMMatch.Pos = atoi(SAMArray[3].begin());
	aSAMMatch.Cigar = SAMArray[5].begin();

	//// A que quoi caaa ????
    if ((aSAMMatch.Flag & C_FLAG_STRAND) == C_FLAG_STRAND) aSAMMatch.Strand = 0;
    else aSAMMatch.Strand = 1;
    aSAMMatch.Unmapped = ( (aSAMMatch.Flag & C_FLAG_UNMAPPED) == C_FLAG_UNMAPPED );

    if (aSAMMatch.Unmapped) aSAMMatch.AlignmentLength = 0;
    else aSAMMatch.AlignmentLength = GetAlignmentLength(aSAMMatch.Cigar);

  return true;
}

// ------------------------------------------------------------------------------

//// TODO: changer le code retour en uint ou uchar
bool TSAMReader::NoSplitGetMatch(ifstream & infile, TSAMMatch & aSAMMatch) const {
	// return true if retrieve a new match entry

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

	string aOptFields; // optional fields (MD and NM)
	string corbeil;

	try{
		infile >> aSAMMatch.QName; // read 1st column
	}
	catch (const std::exception & e) {
		return false; // reach end of file
	}
	infile >> aSAMMatch.Flag; // read 2nd column
	infile >> aSAMMatch.RName; // read 3rd column
	infile >> aSAMMatch.Pos; // read 4th column
	infile >> corbeil; // read 5th column
	infile >> aSAMMatch.Cigar; // read 6th column

	//// A que quoi caaa ????
	if ((aSAMMatch.Flag & C_FLAG_STRAND) == C_FLAG_STRAND) aSAMMatch.Strand = 0;
	else aSAMMatch.Strand = 1;
	aSAMMatch.Unmapped = (aSAMMatch.Flag & C_FLAG_UNMAPPED) == C_FLAG_UNMAPPED;

	if (aSAMMatch.Unmapped) aSAMMatch.AlignmentLength = 0;
	else aSAMMatch.AlignmentLength = GetAlignmentLength(aSAMMatch.Cigar);

	// flush the 6 unused remaining mandatory fields
	infile >> corbeil >> corbeil >> corbeil >> corbeil >> corbeil >> corbeil;

	// get optional fields in one string then find NM and MD positions
	getline(infile, aOptFields);

  // extract optional fields info

	aSAMMatch.Mismatches = -1; // -1 because unknown
	aSAMMatch.EditDistance = -1; // -1 because unknown

	if (aOptFields.empty()) return false; //// TODO: a changer

	// contains opt
	// extract MD field
	size_t begPos = aOptFields.rfind(C_MD_OPT_FIELD);
	size_t aStrSize = 0;
																		 //0000000111
	if (begPos != string::npos) { // MD found                              3456789012
		begPos += C_OPT_FIELD_VALUE_POS - 1; // jump to MD value position (MD:Z:value)
		aStrSize = aOptFields.find_first_of('\t', begPos);
		if (aStrSize == string::npos) aStrSize = aOptFields.size(); // MD is the last field
		aStrSize -= begPos;

		aSAMMatch.Mismatches = GetMismatchesCount(aOptFields.substr(begPos, aStrSize));
		//// tester surcharge  GetMismatchesCount(const char * aOptMDField)
		//// => appel : GetMismatchesCount(&FSAMArray[idx][C_OPT_FIELD_VALUE_POS-1]);
	}
	// extract NM field
	begPos = aOptFields.find(C_NM_OPT_FIELD);
	if (begPos != string::npos) { // NM found                              3456789012
		begPos += C_OPT_FIELD_VALUE_POS - 1; // jump to NM value position (NM:i:value)
		aStrSize = aOptFields.find_first_of('\t', begPos);
		if (aStrSize == string::npos) aStrSize = aOptFields.size(); // NM is the last field
		aStrSize -= begPos;

		aSAMMatch.EditDistance = atoi(aOptFields.substr(begPos, aStrSize).c_str());
	}
	return true;
}


// ------------------------------------------------------------------------------

int TSAMReader::GetMismatchesCount(const string & aOptMDField) const {

  TOptMDStatus aOptMDStatus;

  int Result = 0;
  aOptMDStatus = mdUnknown;

  // "10A5^AC6" => 10 matches, one mismatch (A in the ref seq), 5 matches, 2 deletions (AC deleted from ref seq), and 6 matches
  // "18^CA19" => 18 matches, 2 deletions (CA deleted from ref seq), and 19 matches => no mismatch
  // pour chaque caractere:
  // MD field always begins and ends with a digit. consecutive mismatches are not valid.
  // We want A|C|T|G preceded by a digit => regex "[0-9][ACTG]"

  //char aMDCharacter[] = "ACTG";
  //char aSAMNumber[] = "0123456789";
  // in ASCII table : '0' < '9' < 'A' < 'Z' < '^'

//// new loop
//  for (size_t i=0; i<aOptMDField.size(); i++)
//   {
// 	 const char & c = aOptMDField[i];
// 	//if (aOptMDStatus in [mdUnknown, mdMatch, mdMismatch])
//     if (aOptMDStatus != mdDeletion)
//     {
//       if (c <= '9') aOptMDStatus = mdMatch;
//       else if (c == '^') aOptMDStatus = mdDeletion;
//       else { // character is in [A-Z]
//    	   aOptMDStatus = mdMismatch;
//    	   Result++;
//       }
//     }
//     else // mdDeletion
//       if (c <= '9') aOptMDStatus = mdMatch;
//       // else [A-Z] => still mdDeletion
//   }

//// equivalent to this old loop :
 for (size_t i=0; i<aOptMDField.size(); i++)
  {
	const char & c = aOptMDField[i];
//	if (aOptMDStatus in [mdUnknown, mdMatch, mdMismatch])
	if (aOptMDStatus == mdUnknown || aOptMDStatus == mdMatch || aOptMDStatus == mdMismatch)
    {
      if (c == '^') aOptMDStatus = mdDeletion;
      //else if (strchr(aMDCharacter, c) != NULL)
      else if (c == 'A' || c == 'C' || c == 'T' || c == 'G')
      {
        aOptMDStatus = mdMismatch;
        Result++;
      }
      // c is a base 10 digit. if (c >= '0' && c <= '9') aOptMDStatus = mdMatch;
      else aOptMDStatus = mdMatch;
    }
    else // mdDeletion
    {
      //if (c in aMDCharacter) aOptMDStatus = mdDeletion;
      if (c == 'A' || c == 'C' || c == 'T' || c == 'G') aOptMDStatus = mdDeletion;
      else aOptMDStatus = mdMatch;
    }
  }
  return Result;
}

// ------------------------------------------------------------------------------
// get alignment length from cigar string
// CIGAR: "3M1I3M1D5M" => 3+1+3+1+5 = 13
int TSAMReader::GetAlignmentLength(const string & aCigarStr) const {

	int Result = 0;
	int aValue = 0;

	for (size_t i=0; i<aCigarStr.size(); i++)
	{
		const char & c = aCigarStr[i];
		if ('0' <= c && c <= '9') aValue = aValue * 10 + c - '0';
//		if (c <= '9') aValue = aValue * 10 + c - '0';
		//else if (c == 'M' || c == 'D' || c == 'I' || c == 'S')
		else if (c == 'M' || c == 'D' || c == 'I') //// LOCAL
		{
			Result += aValue;
			aValue = 0;
		}
		else aValue = 0;
	}
	return Result;
}

void TSAMReader::PrintSAMMatch(const TSAMMatch & aSAMMatch) const {

	cout << "AlignmentLength: " << aSAMMatch.AlignmentLength << endl;
	cout << "Cigar: " << aSAMMatch.Cigar << endl;
	cout << "EditDistance: " << aSAMMatch.EditDistance << endl;
	cout << "Flag: " << aSAMMatch.Flag << endl;
	cout << "MapQ: " << aSAMMatch.MapQ << endl;
	cout << "Mismatches: " << aSAMMatch.Mismatches << endl;
	cout << "Opt: " << aSAMMatch.Opt << endl;
	cout << "OptMD: " << aSAMMatch.OptMD << endl;
	cout << "PNext: " << aSAMMatch.PNext << endl;
	cout << "Pos: " << aSAMMatch.Pos << endl;
	cout << "QName: " << aSAMMatch.QName << endl;
	cout << "Qual: " << aSAMMatch.Qual << endl;
	cout << "RName: " << aSAMMatch.RName << endl;
	cout << "RNext: " << aSAMMatch.RNext << endl;
	cout << "Seq: " << aSAMMatch.Seq << endl;
	cout << "Strand: " << aSAMMatch.Strand << endl;
	cout << "TLen: " << aSAMMatch.TLen << endl;
	cout << "Unmapped: " << aSAMMatch.Unmapped << endl;
	cout << endl;

}

