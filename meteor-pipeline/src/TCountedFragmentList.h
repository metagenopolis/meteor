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

#ifndef TCOUNTEDFRAGMENTLIST_H_
#define TCOUNTEDFRAGMENTLIST_H_

#include "TCountedFragment.h"
#include <vector>

//using namespace std;

//// TODO : pointeur ou pas pointeur ... a revoir
//// TODO : voir si on peut fixer la taille dans le constructeur (appelle le constructeur par défaut de TCountedFragment qui met le flag counted a false)
class TCountedFragmentList : public std::vector<TCountedFragment> {
	private:
    /*
      list of index of counted fragment
      FFragmentIndexList is initialized from 1 to fragmentCount in reference
      -1 => not counted
      else : value correspond to index in self object
    */
    int FFragmentCount;
    std::vector<int> FFragmentIndexList; ////was TIntegerDynArray;
    // FFragmentIndexList[fragment_id] == -1, means that the fragment is not counted

	public:
    TCountedFragmentList(const int aFragmentCount); // set capacity to aFragmentCount with reserve()
    ~TCountedFragmentList();

    const std::vector<int> & GetFragmentIndexList() const; //// const ????
    int GetFragmentCount() const;
    int GetCountedFragmentCount() const;
    TCountedFragment * GetItemByFragmentID(int aFragmentID); //// why not return const TCountedFragment & (const reference) ????

    TCountedFragment & AddCountedFragment(const int aFragmentID, const int aFragmentStartLocation, const int aFragmentEndLocation, const int aFragmentSize, const unsigned char aFragmentStrand);

    void Finalize();
};

#endif /* TCOUNTEDFRAGMENTLIST_H_ */
