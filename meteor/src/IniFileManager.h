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

#ifndef INIFILEMANAGER_H_
#define INIFILEMANAGER_H_

#include <iostream>
#include <map>
//#include <unordered_map>
//#include <ctime>
#include <string>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include "utils.h"

//using namespace std;

class IniFileManager {

	protected:
		std::string iniFilename;
//		unordered_map< std::string, unordered_map<std::string, std::string> > m_data;
		std::map< std::string, std::map<std::string, std::string> > m_data;
//		std::vector<TSectionKey> m_order;

	public:
		IniFileManager();
		IniFileManager(const std::string & filename);
		virtual ~IniFileManager();

		void loadIniData(const std::string & filename);

		// write to class member iniFilename by default, or stdout if iniFilename is empty.
		void printIniData_other(const std::string & filename = std::string("")) const;
		void printIniData(const std::string & filename = std::string("")) const;

		void setSection(const std::string & section, const std::map<std::string, std::string> & aMapSection);
		bool getSection(const std::string & section, std::map<std::string, std::string> & aMapSection) const;
		size_t getSectionsNumber() const;

		std::string getIniFilename() const;
		//void setIniFilename(const std::string & newFilename);

		std::string ReadString(const std::string & section, const std::string & key, const std::string & default_value) const;
		void   WriteString(const std::string & section, const std::string & key, const std::string & value);

		int  ReadInteger(const std::string & section, const std::string & key, const int & default_value) const;
		void WriteInteger(const std::string & section, const std::string & key, const int & value);

		bool ReadBool(const std::string & section, const std::string & key, const bool & default_value) const;
		void WriteBool(const std::string & section, const std::string & key, const bool & value);

		template <class T> T ReadType(const std::string & section, const std::string & key, const T & default_value) const
		{
			T result;
			std::string str = ReadString(section, key, "_NO_");
			if (str == "_NO_") return default_value;

			std::istringstream iss(str);
			iss >> result;
			if (iss.fail()) { // conversion error
				throw std::logic_error(
					"In file "+iniFilename+": section '["+section+"]', key '"+key+"'\n"
					"Fail to convert '"+str+"' to desired type !");
			}
			return result;
		}
//		time_t getDateTimeValue(const std::string & section, const std::string & key, const time_t & default_value) const;
//		void setDateTimeValue(const std::string & section, const std::string & key, const time_t & value);
};

////////////////////////////////////////////////////////////////////////

#endif /* INIFILEMANAGER_H_ */

