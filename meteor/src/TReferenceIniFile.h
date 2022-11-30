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

#ifndef TREFERENCEINIFILE_H_
#define TREFERENCEINIFILE_H_

#include "utils.h"
#include "MeteorConstant.h"
#include "IniFileManager.h"

//using namespace std;

const std::string C_REFERENCE_FILE_SECTION = "reference_file";
const std::string C_REFERENCE_INFO_SECTION = "reference_info";
const std::string C_REFERENCE_DATABASE_DIR_STR = "database_dir";
const std::string C_HAS_LITE_INFO = "has_lite_info";

class TReferenceIniFile : public IniFileManager {

public:
	TReferenceIniFile(const std::string & aReferenceIniFileName);
	virtual ~TReferenceIniFile();

    std::string GetReferenceDatabaseDirectory() const;
    std::string GetReferenceName() const;
    TMeteorDBType GetDatabaseType() const;
    bool GetHasLiteInfo() const;
};

#endif /* TREFERENCEINIFILE_H_ */
