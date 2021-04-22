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

#ifndef TSAMREADER_H_
#define TSAMREADER_H_
/*
    see SAM/BAM format specification
    http://samtools.sourceforge.net/
    current version = 16ede77 -12 Sep 2014)
*/
#include <string>
#include <fstream>
#include "utils.h"

//using namespace std;

typedef struct {
    std::string QName;
    int Flag;
    std::string RName;
    int Pos; // 1-based coordinate system
    int MapQ;
    std::string Cigar;
    std::string RNext;
    int PNext;
    int TLen;
    std::string Seq;
    std::string Qual;
    std::string Opt;
    std::string OptMD; // String for mismatches positions = read mask
    int Mismatches; // 'OptMD' post-processing (not a SAM field)
    int Strand; // 1 if strand +, 0 else : Flag post-processing
    bool Unmapped; // Flag post-processing
    int AlignmentLength; // Cigar post-processing
    int EditDistance; // 'Opt' post-processing
} TSAMMatch;

class TSAMReader {
	private:
		std::ifstream FSAMFile; // read only stream on a file. Was TStreamReader
		//char line[BUF_SIZE]; istream::getline(char* line,BUF_SIZE,'\n') or non-member function std::getline(istream s, std::string line);

		// std::vector storing a splitted version of a SAM alignment line.
		// Its size is 12: C_SAM_MANDATORY_FIELD_COUNT + 1 (optional fields not splitted)
		std::vector<std::string> FSAMArray; //// was TArray<std::string>
		TStrings FSAMOptField; //// was TStringList

	public:
		TSAMReader();
		TSAMReader(const std::string & aSAMFileName);
		virtual ~TSAMReader();
		bool FastGetMatch_old(std::string & s, TSAMMatch & aSAMMatch, std::vector<StringRef> & SAMArray) const;
		bool FastGetMatch(std::string & s, TSAMMatch & aSAMMatch, std::vector<StringRef> & SAMArray) const;
		bool NoSplitGetMatch(std::ifstream & infile, TSAMMatch & aSAMMatch) const;
		bool GetMatch_old(const std::string & s, TSAMMatch & aSAMMatch);
		bool GetMatch(const std::string & s, TSAMMatch & aSAMMatch);

		int GetMismatchesCount(const std::string & aOptMDField) const;
		int GetAlignmentLength(const std::string & aCigarStr) const;

		void PrintSAMMatch(const TSAMMatch & aSAMMatch) const;
};


#endif /* TSAMREADER_H_ */
