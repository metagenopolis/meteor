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

#include "TReadMatchInfoList.h"
#include <iostream>

using namespace std;

TReadMatchInfoList::TReadMatchInfoList():FList() {
	  FCount = 0;
}

TReadMatchInfoList::~TReadMatchInfoList() {
}

void TReadMatchInfoList::Initialize(const int & aReadID){
  // NB: FCount is only set to zero, but FList is not cleared
  //   => some or all of its elements (TReadMatchInfoItems) might be overwritten by the MatchInfo data of the next read.
  // See TReadMatchInfoList::Add()
  FCount = 0;
  FReadID = aReadID;
}

TReadMatchInfoItem & TReadMatchInfoList::Add(){

	// Warning, Flist size is not necessarly increased.
	//   See TReadMatchInfoList::Initialize(const int & aReadID)
	FCount++;
	if (FCount > FList.size()) FList.push_back(TReadMatchInfoItem());

	TReadMatchInfoItem & aRefReadMatchInfoItem = FList.at(FCount-1);

	aRefReadMatchInfoItem.IsRejected = true;

	return aRefReadMatchInfoItem; // return a reference to the last element.
}

TReadMatchInfoItem & TReadMatchInfoList::RefItem(size_t i) {
	return FList.at(i); // return a reference to the ith element.
}

size_t TReadMatchInfoList::Count() const{
	return FCount;
}

int TReadMatchInfoList::ReadID() const{
	return FReadID;
}

TReadMatchInfoItem TReadMatchInfoList::GetReadMatchInfoItem(int Index) const
{
	return FList.at(Index);
}

