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

#include "IniFileManager.h"

using namespace std;

IniFileManager::IniFileManager() {
}

IniFileManager::IniFileManager(const string & filename) {
	loadIniData(filename);
}

IniFileManager::~IniFileManager() {
}


/**********************************************************************/
void IniFileManager::loadIniData(const string & filename) {

	iniFilename = filename;

	ifstream infile(filename.c_str()); // read only stream on the file
	if (!infile) {
		// ok file does not exist => new file
		return;
	}

	string line, section, key, value;
//	TSectionKey aSectionKey;
	size_t pos, key_size;

	while (getline(infile, line)) { // ending '\n' is removed

		// erase potential ending '\r' (dos formated ini file)
		string::reverse_iterator rit=line.rbegin();
		if (rit != line.rend() && *rit == '\r') line=line.erase(line.size()-1);
		// trim ending blank characters
		line.erase( line.find_last_not_of(" \t")+1 );
		// skip comments and empty lines
		if (line.empty() || line[0]=='#') continue;

		// get section name
		if (line[0] == '['){
			section = line.substr(1, line.size()-2);
		}
		// else search sign '=' position
		else if ( (pos=line.find_first_of('=')) != string::npos ){
			if (pos==0) continue; // nothing before '='

			// get the size of the key (key_size) and the start position of the value (pos).
			key_size = pos;
			// manage blank(s) around sign '='
			while (line[key_size-1]==' ') key_size--;
			while (pos+1<line.size() && line[pos+1]==' ') pos++;
			key   = line.substr(0,key_size);   // substring starts from 0 and spans key_size characters
			value = line.substr(pos+1);          // start from pos until the end. NB: empty value is valid
			//~ cout<<"(key,value): ("<<key<<","<<value<<")"<<endl;
			if(!key.empty()) {
				m_data[section][key] = value;      // update m_data
			}
		}
	}
	infile.close();
}
/**********************************************************************/
string IniFileManager::ReadString(const string & section, const string & key, const string & default_value) const{

	map<string, string>::const_iterator it2;
	map<string, map<string, string> >::const_iterator it1 = m_data.find(section);

	if (it1 != m_data.end()){
		it2 = it1->second.find(key);
		if (it2 != it1->second.end()) return it2->second;
		//throw invalid_argument("In file "+iniFilename+": section '["+section+"]' has no entry '"+key+"'");
	}
	//throw invalid_argument("File "+iniFilename+" has no section '["+section+"]'");
	return default_value;
}
/**********************************************************************/
void IniFileManager::WriteString(const string & section, const string & key, const string & value){
	m_data[section][key] = value;
}
/**********************************************************************/
int IniFileManager::ReadInteger(const string & section, const string & key, const int & default_value) const{
	int result;
	string str = ReadString(section, key, "NOTFOUND");
	if (str == "NOTFOUND") return default_value;

	istringstream iss(str);
	iss >> result;
	if (iss.fail()) { // conversion error
		throw logic_error(
			"In file "+iniFilename+": section '["+section+"]', key '"+key+"'\n"
			"Fail to convert '"+str+"' to integer !");
	}
	return result;
}
/**********************************************************************/
void IniFileManager::WriteInteger(const string & section, const string & key, const int & value){
	m_data[section][key] = numberToString(value);
}
/**********************************************************************/
bool IniFileManager::ReadBool(const string & section, const string & key, const bool & default_value) const {

	string str = ReadString(section, key, "_NO_");
	if (str == "_NO_") return default_value;

	return (str.empty() || str == "0") ? false : true;
}
/**********************************************************************/
void IniFileManager::WriteBool(const string & section, const string & key, const bool & value) {
	m_data[section][key] = value ? "1" : "0";
}
/**********************************************************************/
//time_t IniFileManager::getDateTimeValue(const string & section, const string & key, const time_t & default_value) const {;}
/**********************************************************************/
//void IniFileManager::setDateTimeValue(const string & section, const string & key, const time_t & value){;}
/**********************************************************************/
//// see declaration and definition in header.h
//template <class T>
//T IniFileManager::ReadType(const string & section, const string & key, const T & default_value) const
//{
//	T result;
//	string str = ReadString(section, key, "NOTFOUND");
//	if (str == "NOTFOUND") return default_value;
//
//	istringstream iss(str);
//	iss >> result;
//	if (iss.fail()) { // conversion error
//		throw logic_error(
//			"In file "+iniFilename+": section '["+section+"]', key '"+key+"'\n"
//			"Fail to convert '"+str+"' to integer !");
//	}
//	return result;
//}
/**********************************************************************/
void IniFileManager::printIniData(const string & filename) const {

	// write to class member iniFilename by default
	const string & outputFilename = filename.empty() ? iniFilename : filename;

	std::ofstream ofs (outputFilename.c_str(), std::ofstream::out);

	map<string, string>::const_iterator itr1;
	map<string, map<string, string> >::const_iterator itr2;

	for(itr2 = m_data.begin(); itr2 != m_data.end(); itr2++)
	{
		ofs<<"["<<itr2->first<<"]"<<endl;
		for(itr1 = itr2->second.begin(); itr1 != itr2->second.end(); itr1++)
		{
			ofs<<itr1->first<<" = "<<itr1->second<<endl;
		}
	}
	ofs.close();
}/**********************************************************************/
void IniFileManager::printIniData_other(const string & filename) const {

//	// re-inserting m_data elements one by one into um_tmp will reverse the elements order.
//	unordered_map<string, unordered_map<string, string> > um_tmp;
//
//	// write to class member iniFilename by default
//	const string & outputFilename = filename.empty() ? iniFilename : filename;
//
//	std::ofstream ofs (outputFilename.c_str(), std::ofstream::out);
//
//	unordered_map<string, string>::const_iterator itr1;
//	unordered_map<string, unordered_map<string, string> >::const_iterator itr2;
//
//	// insert element one by one from m_data to um_tmp (=> reverse the order)
//	for(itr2 = m_data.begin(); itr2 != m_data.end(); itr2++)
//	{
//	    for(itr1 = itr2->second.begin(); itr1 != itr2->second.end(); itr1++)
//	    {
//	        um_tmp[itr2->first][itr1->first]=itr1->second;
//	    }
//	}
//	// print um_tmp content
//	//for(itr2 = m_data.begin(); itr2 != m_data.end(); itr2++)
//	for(itr2 = um_tmp.begin(); itr2 != um_tmp.end(); itr2++)
//	{
//	    ofs<<"["<<itr2->first<<"]"<<endl;
//	    for(itr1 = itr2->second.begin(); itr1 != itr2->second.end(); itr1++)
//	    {
//	        ofs<<itr1->first<<" = "<<itr1->second<<endl;
//	    }
//	}
//	ofs.close();
}
/**********************************************************************/
void IniFileManager::setSection(const string & section, const map<string, string> & aMapSection){
	m_data[section] = aMapSection;
}
/**********************************************************************/
bool IniFileManager::getSection(const string & section, map<string, string> & aMapSection) const {
	map<string, map<string, string> >::const_iterator it1 = m_data.find(section);
	if (it1 != m_data.end()){
		aMapSection = it1->second;
		return true;
	}
	//throw invalid_argument("File "+iniFilename+" has no section '["+section+"]'");
	return false; // section not found
}
/**********************************************************************/
size_t IniFileManager::getSectionsNumber() const {
	return m_data.size();
}
string IniFileManager::getIniFilename() const {
	return iniFilename;
}
////////////////////////////////////////////////////////////////////////
