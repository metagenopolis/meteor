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

#include "TReferenceIniFile.h"

using namespace std;

TReferenceIniFile::TReferenceIniFile(const string & filename):IniFileManager(filename){
	if (getSectionsNumber() == 0){
		cerr << "Error, no section found in " << filename << endl;
		exit(1);
	}
}

TReferenceIniFile::~TReferenceIniFile() {}

// ------------------------------------------------------------------------------

string TReferenceIniFile::GetReferenceDatabaseDirectory() const
{
	return ReadString(C_REFERENCE_FILE_SECTION, C_REFERENCE_DATABASE_DIR_STR, "");
}

// ------------------------------------------------------------------------------

string TReferenceIniFile::GetReferenceName() const
{
	return ReadString(C_REFERENCE_INFO_SECTION, C_REFERENCE_NAME_STR, "");
}

// ------------------------------------------------------------------------------

TMeteorDBType TReferenceIniFile::GetDatabaseType() const
{
	  string aMeteorDBTypeStr = ReadString(C_REFERENCE_INFO_SECTION, C_DATABASE_TYPE_STR, C_METEOR_DBTYPE_ADT_STR);

	  if (aMeteorDBTypeStr == C_METEOR_DBTYPE_ADT_STR) return mdbADT;
	  if (aMeteorDBTypeStr == C_METEOR_DBTYPE_PARSTREAM_STR) return mdbParstream;
	  if (aMeteorDBTypeStr == C_METEOR_DBTYPE_GIGABASE_STR) return mdbGigabase;
	  if (aMeteorDBTypeStr == C_METEOR_DBTYPE_BINARY_STR) return mdbBinary;
	  else return mdbUnknown;
}

// ------------------------------------------------------------------------------

bool TReferenceIniFile::GetHasLiteInfo() const
{
	return ReadInteger(C_REFERENCE_INFO_SECTION, C_HAS_LITE_INFO, 0) == 1; // default is 0;
}
