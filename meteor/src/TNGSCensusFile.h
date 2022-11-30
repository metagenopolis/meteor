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

#ifndef TNGSCENSUSFILE_H_
#define TNGSCENSUSFILE_H_

#include "TNGSCensusData.h"
#include "TCountedFragmentList.h"
#include <vector>

const std::string C_CENSUS_TABLE_FILE_EXT = ".dat";
const std::string C_CENSUS_PROFILE_TABLE_FILE_EXT = ".profile.dat";
const int C_INPUT_BUF_ITEM_COUNT = 10000;

const int C_CENSUS_TABLE_FILE_COL_NUMBER = 30;
const std::string C_CENSUS_TABLE_FILE_HEADER =
"FragmentID\t\
FragmentSize\t\
TotalReadCount\t\
SharedReadCount\t\
SmartSharedReadCount\t\
UniqueReadCount\t\
DirectTotalReadCount\t\
DirectSharedReadCount\t\
DirectSmartSharedReadCount\t\
DirectUniqueReadCount\t\
TotalReadCoverage\t\
UniqueReadCoverage\t\
DirectTotalReadCoverage\t\
DirectUniqueReadCoverage\t\
UnDirectTotalReadCoverage\t\
UnDirectUniqueReadCoverage\t\
TotalReadMeanCoverage\t\
SharedReadMeanCoverage\t\
SmartSharedReadMeanCoverage\t\
UniqueReadMeanCoverage\t\
DirectTotalReadMeanCoverage\t\
DirectSharedReadMeanCoverage\t\
DirectSmartSharedReadMeanCoverage\t\
DirectUniqueReadMeanCoverage\t\
UnDirectTotalReadMeanCoverage\t\
UnDirectSharedReadMeanCoverage\t\
UnDirectSmartSharedReadMeanCoverage\t\
UnDirectUniqueReadMeanCoverage\t\
UnDirectTotalReadCount\t\
UnDirectSharedReadCount\t\
UnDirectSmartSharedReadCount\t\
UnDirectUniqueReadCount";

//// every long long was integer (delphi).
//   Size of long long and double is 8
//   This should avoid padding effects when write/read to file in binary mode
//   => sizeof(TNGSCensusItem) should always be 28*8 = 224 Bytes
typedef struct {
    long long FragmentID;
    long long FragmentSize;
    long long TotalReadCount;
    double SharedReadCount;
    double SmartSharedReadCount;
    long long UniqueReadCount;
    long long DirectTotalReadCount;
    double DirectSharedReadCount;
    double DirectSmartSharedReadCount;
    long long DirectUniqueReadCount;
    double TotalReadCoverage;
    double UniqueReadCoverage;
    double DirectTotalReadCoverage;
    double DirectUniqueReadCoverage;
    double UnDirectTotalReadCoverage;
    double UnDirectUniqueReadCoverage;
    double TotalReadMeanCoverage;
    double SharedReadMeanCoverage;
    double SmartSharedReadMeanCoverage;
    double UniqueReadMeanCoverage;
    double DirectTotalReadMeanCoverage;
    double DirectSharedReadMeanCoverage;
    double DirectSmartSharedReadMeanCoverage;
    double DirectUniqueReadMeanCoverage;
    double UnDirectTotalReadMeanCoverage;
    double UnDirectSharedReadMeanCoverage;
    double UnDirectSmartSharedReadMeanCoverage;
    double UnDirectUniqueReadMeanCoverage;

    long long UnDirectTotalReadCount;
    double UnDirectSharedReadCount;
    double UnDirectSmartSharedReadCount;
    long long UnDirectUniqueReadCount;
} TNGSCensusItem;

  //// type binary file of struct TNGSCensusItem
  //TFileOfCensusItem = File of TNGSCensusItem;

class TNGSCensusFile : public TNGSCensusData {

  private:
    std::string FCensusTableFileName;
    std::string FCensusProfileFileName;

  public:

    TNGSCensusFile(const std::string & aNGSCensusTablePathName, const std::string & aNGSCensusTableName);
    virtual ~TNGSCensusFile();

    void WriteCountedFragmentListBinary(TCountedFragmentList & aCountedFragmentList);
    void WriteCountedFragmentListAsciiText(TCountedFragmentList & aCountedFragmentList);
    void WriteCountedFragmentListAsciiText(TCountedFragmentList & aCountedFragmentList, const std::vector<TCountedReadType> & aCountedReadTypeList);

    void WriteCountedFragmentList(TCountedFragmentList & aCountedFragmentList, const std::vector<TCountedReadType> & aCountedReadTypeList, char aWriteMode = 'b'); //override

};

#endif /* TNGSCENSUSFILE_H_ */
