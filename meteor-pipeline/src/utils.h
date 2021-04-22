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

#ifndef UTILS_H_
#define UTILS_H_

#include<cstdlib>
#include<string>
#include<sstream>
#include<iostream>
#include<vector>

//using namespace std;

#if defined (_WIN64) || defined(_WIN32) || defined(WIN32) || defined(__WIN32__)
#define __WINDOWS__
#endif

#ifdef __WINDOWS__
	#include <direct.h>    // mkdir()
	#include <io.h>        // define _access()
	#define access _access // _access is used to check if a file exists (windows systems)
	#define F_OK 0	       // access() 2nd argument
#else
	#include <sys/stat.h>  // mkdir()
	#include <unistd.h>    // define the function access to check if a file/directory exists (unix/linux systems)
#endif

#ifdef __WINDOWS__
const char C_PATH_SEP = '\\';
#else
const char C_PATH_SEP = '/';
#endif

typedef std::vector<std::string> TStrings;
typedef long long int64;

// --------------------------------------------------------------------

int my_mkdir(const char *filename); // cross-platform mkdir()
// create a new directory including the parent directories as needed.
//// TODO: test if really cross-platform (=> windows)
int ForceDirectories(const std::string & path);

// --------------------------------------------------------------------

std::string str_now(); // get the current local date as a std::string (yyyy-mm-dd)

// --------------------------------------------------------------------

// Class used to optimize std::string split while reading large files.
// member begin_ is the address of an non-member C std::string.
// StringRef object does NOT own the pointed C std::string.
// => begin_ may be invalidated by further calls to external functions/procedures on this non-member std::string.
class StringRef
{
private:
    const char *     begin_; // address of an external C std::string (that has to end with '\0')
    int             size_;  // size of the pointed C std::string.

public:
    int size() const { return size_; }
    const char * begin() const { return begin_; }
    const char * end() const { return begin_ + size_; }

    StringRef(): begin_(NULL), size_(0){}
    StringRef(char const* const begin, int const size): begin_( begin ), size_( size ){}
};

// --------------------------------------------------------------------

template <typename T> std::string numberToString(T n, int precision = 0) {
     std::ostringstream oss;
     if (precision>0) oss.precision(precision);
     oss << n;
     return oss.str();
}

// --------------------------------------------------------------------

std::string ExtractFileDir(const std::string & filename); // same as dirname
std::string BaseName(const std::string & filename);
std::string replaceStr(const std::string & str, const std::string& from, const std::string& to);
bool getFiles(const std::string & directory, const std::string & fileSuffix, std::vector<std::string> & fileNameList);

// --------------------------------------------------------------------

// these procedures store addresse(s) of std::string chunks in a std::vector of StringRef
// Chunks are delimited by parameter: char delimiter.
// each delimiter char found in the input std::string is replaced with '\0'.
// => The input std::string IS modified and should be destroyed as soon as chunks are exploited.
void stringToStringRefVector(std::string & str, std::vector<StringRef> &result, const char & delimiter = ' ');
size_t stringToStringRefVector(std::string & str, std::vector<StringRef> &result, const size_t & chunksNb, const char & delimiter = ' ');

// --------------------------------------------------------------------

// delimiter type : char (faster) or std::string (=list of delimiters, but is slower)
template <typename T> void splitString(std::vector<std::string> & tokens, const std::string & str, const T & delimiter) {
    // Skip delimiters at beginning
    std::string::size_type lastPos = str.find_first_not_of(delimiter, 0);

    // Find first non-delimiter
    std::string::size_type pos = str.find_first_of(delimiter, lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos) {
        // Found a token, add it to the std::vector
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters
        lastPos = str.find_first_not_of(delimiter, pos);
        // Find next non-delimiter
        pos = str.find_first_of(delimiter, lastPos);
    }
}

// --------------------------------------------------------------------

// overloaded version. Does not push_back strings to empty std::vector tokens => tokens size has to be egal to chunksNb
template <typename T> size_t splitString(std::vector<std::string> & tokens, const std::string & str, const T & delimiter, const size_t & chunksNb) {
    // Skip delimiters at beginning
    std::string::size_type lastPos = str.find_first_not_of(delimiter, 0);

    // Find first non-delimiter
    std::string::size_type pos = str.find_first_of(delimiter, lastPos);

    size_t i = 0;

    while (std::string::npos != pos || std::string::npos != lastPos) {
        // Found a token, add it to the std::vector.
        // When chunksNb is reached, the remaining part of str is added unchanged (unsplitted) in the std::vector.
    	if(i+1 == chunksNb) {
        	tokens[i++] = str.substr(lastPos);
        	return i;
        }
        tokens[i++] = str.substr(lastPos, pos - lastPos);

        // Skip delimiters
        lastPos = str.find_first_not_of(delimiter, pos);
        // Find next non-delimiter
        pos = str.find_first_of(delimiter, lastPos);
    }
    // return the effective number of splitted elements (can be lower than chunksNb)
    return i;
}

#endif /* UTILS_H_ */
