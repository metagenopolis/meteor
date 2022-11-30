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

#include "utils.h"
#include <ctime>
#include <errno.h>

// needed for getFiles
#include <cstdlib>
#include <cstdio>
//#include <errno.h>
//#include <sys/types.h>
#include<dirent.h>
#include<string>
#include<vector>

using namespace std;

// --------------------------------------------------------------------

string ExtractFileDir(const string & filename){
	size_t ilast = filename.size()-1;
	while (filename[ilast]==C_PATH_SEP) ilast--;                // ignore ending separator(s)
	size_t pos = filename.find_last_of(C_PATH_SEP, ilast);      // then find last separtor
	return pos == string::npos ? "." : filename.substr(0, pos); // ending path separator is not included in the resulting substring
}

// --------------------------------------------------------------------

string BaseName(const string & filename){
	size_t ilast = filename.size()-1;
	while (filename[ilast]==C_PATH_SEP) ilast--;             // ignore ending separator(s)
	size_t pos = filename.find_last_of(C_PATH_SEP, ilast);   // then find last separtor
	return pos == string::npos ? filename : filename.substr(pos+1);
}

// --------------------------------------------------------------------
//// TODO : overloaded version with const char* from, to + declare a string copy local to the function
string replaceStr(const std::string & str, const std::string& from, const std::string& to) {
    size_t pos = str.find(from);
    return pos == std::string::npos ? string(str) : string(str).replace(pos, from.length(), to);
}

// --------------------------------------------------------------------

//#include <cstdlib>
//#include <cstdio>
////#include <errno.h>
////#include <sys/types.h>
//#include<dirent.h>
//#include<string>
//#include<vector>

bool getFiles(const string & directory, const string & fileSuffix, vector<string> & fileNameList){

	DIR *dir; struct dirent *ent;
	fileNameList.clear();

	// check if path ends with C_PATH_SEP
	string sep;
	if (directory[directory.size()-1] != C_PATH_SEP) sep+=C_PATH_SEP;

	if ( ( dir = opendir( directory.c_str() ) ) != NULL) {

		// iterate all the files and directories within directory
		while ((ent = readdir(dir)) != NULL) {

			string fileName(ent->d_name);
			size_t pos = fileName.rfind(fileSuffix);

			if (pos != string::npos && pos+fileSuffix.size() == fileName.size()){
				fileNameList.push_back(directory+sep+fileName);
			}
		}
		closedir(dir);
	}
	else {
	  /* could not open directory */
	  perror ("");
	  return false;
	}
	return true;
}

// --------------------------------------------------------------------

// this procedure stores addresse(s) of string chunks in a vector of StringRef
// Chunks are delimited by parameter: char delimiter.
// each delimiter char found in the input string is replaced with '\0'.
// => The input string IS modified and should be destroyed as soon as chunks are exploited.
void stringToStringRefVector(string & str, vector<StringRef> &result, const char & delimiter)
{
	enum State { inSpace, inToken };
	size_t size = str.size();
	State state = inSpace;
	char * pTokenBegin = NULL;	// Init to satisfy compiler.
	//~ for( string::iterator it = str.begin(); it != str.end(); ++it )
	for( size_t it = 0 ; it < size; ++it )
	{
		//~ State const newState = (*it == delimiter? inSpace : inToken);
		State const newState = (str[it] == delimiter? inSpace : inToken);
		if( newState != state )
		{
			switch( newState )
			{
			case inSpace:
				//~ result.push_back( StringRef( pTokenBegin, &*it - pTokenBegin ) );
				result.push_back( StringRef( pTokenBegin, &str[it] - pTokenBegin ) );
				str[it]='\0';
				break;
			case inToken:
				//~ pTokenBegin = &*it;
				pTokenBegin = &str[it];
				//~ if (it-1>0 && str[it-1]==delimiter) printf("c:'%c' ",str[it-1]);//='\0';
				//~ if (it>0) { cout<<"it-1: "<<it-1<<", "; str[it-1]='\0';}
			}
		}
		state = newState;
	}
	// last token
	if( state == inToken )
	{
		//~ result.push_back( StringRef( pTokenBegin, &*str.end() - pTokenBegin ) );
		result.push_back( StringRef( pTokenBegin, &*str.end() - pTokenBegin ) );
	}
}

// --------------------------------------------------------------------

// this procedure stores addresse(s) of string chunks in a vector of StringRef
// Chunks are delimited by parameter: char delimiter.
// each delimiter char found in the input string is replaced with '\0'.
// => The input string IS modified and should be destroyed as soon as chunks are exploited.
size_t stringToStringRefVector( string & str, vector<StringRef> &result, const size_t & chunksNb, const char & delimiter)
{
	enum State { inSpace, inToken };
	size_t size = str.size();
	size_t ivect = 0;
	State state = inSpace;
	char * pTokenBegin = NULL;	// Init to satisfy compiler.
	//~ for( string::iterator it = str.begin(); it != str.end(); ++it )
	for( size_t it = 0 ; it < size; ++it )
	{
		//~ State const newState = (*it == delimiter? inSpace : inToken);
		State const newState = (str[it] == delimiter? inSpace : inToken);
		if( newState != state )
		{
			switch( newState )
			{
				case inSpace:
					// When chunksNb is reached, the remaining part of str is added unchanged (unsplitted) in the vector.
					if(ivect+1 == chunksNb) {
						result[ivect++] = StringRef( pTokenBegin, &*str.end() - pTokenBegin );
						return ivect;
					}
					//result.push_back( StringRef( pTokenBegin, &str[it] - pTokenBegin ) );
					result[ivect++] = StringRef( pTokenBegin, &str[it] - pTokenBegin );
					str[it]='\0';
					break;
				case inToken:
					//~ pTokenBegin = &*it;
					pTokenBegin = &str[it];
			}
		}
		state = newState;
	}
	// last token
	if( state == inToken )
	{
		//~ result.push_back( StringRef( pTokenBegin, &*str.end() - pTokenBegin ) );
		result[ivect++] = StringRef( pTokenBegin, &*str.end() - pTokenBegin );
	}
	return ivect;
}


/**-----------------------------------------------------------------------------
	## FUNCTION:
	int my_mkdir(char *filename)
	-----------------------------------------------------------------------------
	## RETURN: status of the function mkdir. 0 on success, -1 otherwise
	-----------------------------------------------------------------------------
	## PARAMETERS:
	@ char *filename : the name of the file
	-----------------------------------------------------------------------------
	## SPECIFICATION: creation of a directory (for Linux, Windows and Mac)
	-----------------------------------------------------------------------------
*/

int my_mkdir(const char *filename){
	int status = -1;

	#ifdef __WINDOWS__
	status = _mkdir(filename); // needs <direct.h>
//	#elif __linux__
//	status = mkdir(filename, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	#else
	status = mkdir(filename, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); // needs <sys/stat.h>
	#endif

	return status;
}

/**-----------------------------------------------------------------------------
	## FUNCTION:
	int my_mkpath(const string & path)
	-----------------------------------------------------------------------------
	## RETURN: status of the function mkdir. 0 on success, -1 otherwise
	-----------------------------------------------------------------------------
	## PARAMETERS:
	@ const string & path : the path name
	-----------------------------------------------------------------------------
	## SPECIFICATION: create a directory and its parent directories if needed (for Linux, Windows and Mac)
	-----------------------------------------------------------------------------
*/


int ForceDirectories(const string & path){

	int status = -1;
	string tmpPath;
	vector<string> tokens;

	// path starts with C_PATH_SEP
	if (!path.empty() && path[0] == C_PATH_SEP) tmpPath += C_PATH_SEP;

	splitString(tokens, path, C_PATH_SEP);

	for (size_t i = 0; i<tokens.size(); i++){
		tmpPath += tokens[i] + C_PATH_SEP;
		status = my_mkdir(tmpPath.c_str());
		if ( status != 0 && errno != EEXIST ) {
			cerr<<"Error, failed to create directory: "+tmpPath+" while creating path: "+path<<endl;
			perror("");
			return status;
		}
	}
	return 0;
}

/**-----------------------------------------------------------------------------
	## FUNCTION:
	string now()
	-----------------------------------------------------------------------------
	## RETURN: a std::string representing the current local date (yyyy-mm-dd)
	-----------------------------------------------------------------------------
	## PARAMETERS:
	-----------------------------------------------------------------------------
	## SPECIFICATION:
	-----------------------------------------------------------------------------
*/

string str_now(){

	//get current time as time_t (integer)
	time_t t = time(NULL);

	// localtime(&t) convert time_t to tm data type (struct) and return its address.
	// IT IS NOT THREAD SAFE => localtime_r and localtime_s are thread safe but not portable.
	// see https://www.tutorialcup.com/cplusplus/date-time.htm for struct tm description

	// format "%F" get date as yyyy-mm-dd (10 characters)
	// see http://www.cplusplus.com/reference/ctime/strftime/ for other format
	char yyyymmdd[11];
	strftime(yyyymmdd, 11, "%F", localtime(&t));

	return string(yyyymmdd);
}
