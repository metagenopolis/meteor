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

#ifndef TLITEREFERENCEANNOTATIONTABLE_
#define TLITEREFERENCEANNOTATIONTABLE_

#include <vector>
#include "MeteorConstant.h"

class TLiteReferenceAnnotationTable {
	  private:
		std::vector<int>FTable;
	    //TStrings FTempTable; // for reference creation step //// not used
	  public:
	    TLiteReferenceAnnotationTable();
	    ~TLiteReferenceAnnotationTable();

	    int GetFragmentSize(const int aFragmentID) const;
	    int GetFragmentCount() const;
	    void OpenReferenceAnnotationTable(const std::string & aReferenceAnnotationTableFileName);
	    //void AddFragment(const int aFragmentID, const int aFragmentSize); //// not used
	    //void SaveTable(const std::string & aReferenceAnnotationTableFileName); //// not used
};

#endif /* TLITEREFERENCEANNOTATIONTABLE_ */
