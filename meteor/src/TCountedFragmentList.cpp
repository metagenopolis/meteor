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

#include "TCountedFragmentList.h"

using namespace std;

//TCountedFragmentList::TCountedFragmentList() {
//	// TODO Auto-generated constructor stub
//
//}

TCountedFragmentList::~TCountedFragmentList() {
	// implicitly call base class destructor and FFragmentIndexList destructor.
	// if we had pointers to objects as elements, we would have to delete them first (in a loop)
}

TCountedFragmentList::TCountedFragmentList(const int aFragmentCount):
		vector<TCountedFragment>(),
		FFragmentIndexList(aFragmentCount+1, -1) // +1 because id_fragment=0 does not exist, -1 means not counted.
{
	this->reserve(aFragmentCount/10 +1);
	// As many fragments (~90%) may not be counted at all, we only reserve a capacity of 10% of total fragment count (aFragmentCount/10)
	// Use of method vector::reserve() => further push_back will not induce expensive reallocation

	FFragmentCount = aFragmentCount;
}

const vector<int> & TCountedFragmentList::GetFragmentIndexList() const {
	return FFragmentIndexList;
}

int TCountedFragmentList::GetFragmentCount() const { return FFragmentCount;}
int TCountedFragmentList::GetCountedFragmentCount() const { return this->size(); }

TCountedFragment * TCountedFragmentList::GetItemByFragmentID(int aFragmentID) //// why not return TCountedFragment & (reference) ????
{
  if (FFragmentIndexList.at(aFragmentID) != -1) //// AT
    return &(this->at(FFragmentIndexList.at(aFragmentID))); //// AT

  return NULL;

}

//// TODO : retourne reference ou valeur (copie) ???? a revoir
TCountedFragment & TCountedFragmentList::AddCountedFragment(
//TCountedFragment TCountedFragmentList::AddCountedFragment(
		const int aFragmentID,
		const int aFragmentStartLocation,
		const int aFragmentEndLocation,
		const int aFragmentSize,
		const unsigned char aFragmentStrand)
{
	//TCountedFragment* aCountedFragment;

	if (FFragmentIndexList.at(aFragmentID) == -1) //// AT
	{
	     TCountedFragment aCountedFragment(aFragmentID, aFragmentStartLocation, aFragmentEndLocation, aFragmentSize, aFragmentStrand);

	    // add and update with new index.
	    ////FFragmentIndexList[aFragmentID] = Self.Add(aCountedFragment);
	    this->push_back(aCountedFragment);
	    ////si taille de this peut etre fixee des le depart (dans le constructeur), utiliser TCountedFragment::update()
	    ////this->at[FFragmentIndexList[aFragmentID]].update(aFragmentID, aFragmentStartLocation, aFragmentEndLocation, aFragmentSize, aFragmentStrand)
	    // index of the last inserted element
	    FFragmentIndexList.at(aFragmentID) = this->size() - 1; //// AT

	}
	//else aCountedFragment = this->at(FFragmentIndexList[aFragmentID]); //// pas utile
	//return aCountedFragment;
	//// return a copy of the CountedFragment, is necessary ????
	return this->at(FFragmentIndexList.at(aFragmentID)); //// AT
}

// ------------------------------------------------------------------------------

void TCountedFragmentList::Finalize()
{
	for (size_t i = 0; i < this->size(); i++)
		this->at(i).Finalize();
}

