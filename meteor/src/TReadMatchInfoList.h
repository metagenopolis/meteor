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

#ifndef TREADMATCHINFOLIST_H_
#define TREADMATCHINFOLIST_H_

#include <string>
#include <vector>

//using namespace std;

typedef struct {
	int ReadLocation;
	bool ReadStrand;
	int ReadMismatches;
	int ReadEntryID;
	int ReadSAMFlag;
	std::string ReadSAMCigar;
	std::string ReadSAMOpt;
	int ReadMappingLength;
	int ReadEditDistance;
	double ReadDistance;
	bool IsRejected;
} TReadMatchInfoItem;


class TReadMatchInfoList {

	private:
	int FReadID;
	std::vector<TReadMatchInfoItem> FList; // was TObjectList
	size_t FCount;

public:
  TReadMatchInfoList();
  virtual ~TReadMatchInfoList();
  void Initialize(const int & aReadID);
  size_t Count() const; // return FCount
  int ReadID() const; // return FReadID;

  TReadMatchInfoItem & RefItem(size_t i);
  TReadMatchInfoItem & Add();
  void Add(const int aReadLocation, const bool aReadStrand, const int aReadMismatches, const int aReadEntryID);
  void Add(const int aReadLocation, const bool aReadStrand, const int aReadMismatches, const int aReadEntryID, const int aReadMappingLength);
  void Add(const int aReadLocation, const bool aReadStrand, const int aReadMismatches, const int aReadEntryID,
           const int aReadSAMFlag, const std::string & aReadSAMCigar, const std::string & aReadSAMOpt, int aReadMappingLength, int aReadEditDistance);
  void Add(const bool aReadStrand, const int aReadEntryID, const int aReadMappingLength);
  void Add(const bool aReadStrand, const int aReadEntryID, const std::string & aReadSAMCigar);

  TReadMatchInfoItem GetReadMatchInfoItem(int Index) const;


};

#endif /* TREADMATCHINFOLIST_H_ */
