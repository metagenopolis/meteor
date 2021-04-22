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

#ifndef TPROFILETABLE_H_
#define TPROFILETABLE_H_

#include <string>

//using namespace std;

const int INPUT_COLS_NB = 32; // number of column in input files
const int NB_OF_VALUES = 30;  // two first col are fragment id and size

// Store info on a file census.dat
typedef struct TInputFileInfo {
	std::string FileName;
	long Progression;   // progression in the file (number of bytes)
	int WantedColId;    // column index to consider in the file, according to the wanted counted type.
	int OptionalColId;  // column index to process census.dat, where 4 undirect counting types were deduced from total - direct columns.
	                    //   i.e. crtUndirectTotalReads, crtUndirectSharedReads, crtUndirectSmartSharedReads, crtUndirectUniqueReads
	int NbOfColumns;    // number of counting columns in the file (might not be uniform)
	//bool IsOldStyle;  // the file may contain all type of counting. (not usefull ??)
} TInputFileInfo;

class TProfileTable {

	private:

		size_t FSampleCount; // number of columns of F2DArray
		int FLineNumber;
		int FMaxRowsNumber;
		TCountedReadType FCountedReadType;
		std::vector<std::string> FInputFileNames;
		std::vector<long> FProgressionArray;
		std::vector<int> FWantedColumnArray;
		std::vector<int> FNbOfColumnsArray;

		std::vector<TInputFileInfo> FInputFileInfoArray;

		double ** F2DArray;

	public:

		TProfileTable(const std::vector<std::string> & aFileNames, const TCountedReadType aCountedReadType);
		virtual ~TProfileTable();
		void Reset2DArray();
		int MergeColumns(int iWantedCol);
		int MergeColumns();

		void Write2DArray(int aReadLinesNb, std::ofstream & aOutputFile) const;

		int GetMaxRowsNumber() const;
		size_t GetSampleCount() const;
		int GetLineNumber() const;
};

#endif /* TPROFILETABLE_H_ */
