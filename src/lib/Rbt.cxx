/***********************************************************************
 * The rDock program was developed from 1998 - 2006 by the software team
 * at RiboTargets (subsequently Vernalis (R&D) Ltd).
 * In 2006, the software was licensed to the University of York for
 * maintenance and distribution.
 * In 2012, Vernalis and the University of York agreed to release the
 * program as Open Source software.
 * This version is licensed under GNU-LGPL version 3.0 with support from
 * the University of Barcelona.
 * http://rdock.sourceforge.net/
 ***********************************************************************/

// Misc non-member functions in Rbt namespace

#include <dirent.h>  //For directory handling
#include <limits.h>  //For PATH_MAX
#include <stdlib.h>  //For getenv
#include <time.h>    //For time functions
#include <unistd.h>  //For POSIX getcwd

#include <algorithm>  //For sort
#include <fstream>    //For ifstream
//#include <ios>
#include "Rbt.h"
#include "RbtBinaryIO.h"
#include "RbtFileError.h"
#include "RbtResources.h"
#include "RbtVersion.h"
using std::istream;

// GetRbtRoot - returns value of RBT_ROOT env variable
RbtString Rbt::GetRbtRoot() {
    char* szRbtRoot = getenv("RBT_ROOT");
    if (szRbtRoot != (char*)NULL) {
        return RbtString(szRbtRoot);
    } else {
        return GetCurrentDirectory();
    }
}

// DM 02 Aug 2000
// GetRbtHome - returns value of RBT_HOME env variable
// or HOME if RBT_HOME is not defined
// If HOME is undefined, returns current working directory
RbtString Rbt::GetRbtHome() {
    char* szRbtHome = getenv("RBT_HOME");
    if (szRbtHome != (char*)NULL) {
        return RbtString(szRbtHome);
    } else {
        szRbtHome = getenv("HOME");
        if (szRbtHome != (char*)NULL) {
            return RbtString(szRbtHome);
        } else {
            return GetCurrentDirectory();
        }
    }
}

// Rbt::GetCopyright - returns legalese statement
RbtString Rbt::GetCopyright() { return IDS_COPYRIGHT; }
// current version format is 'v<year>.<month>.<build_number>' so we can extract the version number
// by taking the substring from the second character to the 5th character
// Rbt::GetVersion - returns current library version
RbtString Rbt::GetVersion() { return std::string(RBT_VERSION).substr(1, 5); }
// same as above for the build number, but we take the substring from the 7th character to the end
// GetBuildNum - returns current library build number
RbtString Rbt::GetBuild() { return std::string(RBT_VERSION).substr(7); }
// GetProduct - returns library product name
RbtString Rbt::GetProduct() { return IDS_PRODUCT; }
// GetTime - returns current time as an RbtString
RbtString Rbt::GetTime() {
    time_t t = ::time(NULL);                  // Get time in seconds since 1970
    tm* pLocalTime = ::localtime(&t);         // Convert to local time struct
    return RbtString(::asctime(pLocalTime));  // Convert to ascii string
}
// GetCurrentDirectory - returns current working directory
RbtString Rbt::GetCurrentDirectory() {
    RbtString strCwd(".");
    char* szCwd = new char[PATH_MAX + 1];            // Allocate a temp char* array
    if (::getcwd(szCwd, PATH_MAX) != (char*)NULL) {  // Get the cwd
        strCwd = szCwd;
        // strCwd += "/";
    }
    delete[] szCwd;  // Delete the temp array
    return strCwd;
}

// Rbt::GetRbtDirName
// Returns the full path to a subdirectory in the rDock directory structure
//
// For example, if RBT_ROOT environment variable is ~dave/ribodev/molmod/ribodev
// then GetRbtDirName("data") would return ~dave/ribodev/molmod/ribodev/data/
//
RbtString Rbt::GetRbtDirName(const RbtString& strSubDir) {
    RbtString strRbtDir = GetRbtRoot();
    if (strSubDir.size() > 0) {
        strRbtDir += "/";
        strRbtDir += strSubDir;
    }
    return strRbtDir;
}

// Rbt::GetRbtFileName
// DM 17 Dec 1998 - slightly different behaviour
// First check if the file exists in the CWD, if so return this path
// Next check RBT_HOME directory, if so return this path
// Finally, return the path to the file in the rDock directory structure
//(without checking if the file is actually present)
RbtString Rbt::GetRbtFileName(const RbtString& strSubdir, const RbtString& strFile) {
    // First see if the file exists in the current directory
    RbtString strFullPathToFile(strFile);
    // Just open it, don't try and parse it (after all, we don't know what format the file is)
    std::ifstream fileIn(strFullPathToFile.c_str(), std::ios_base::in);
    if (fileIn) {
        fileIn.close();
        return strFullPathToFile;
    } else {
        // DM 02 Aug 200 - check RBT_HOME directory
        strFullPathToFile = GetRbtHome() + "/" + strFile;
        // DM 27 Apr 2005 - under gcc 3.4.3 there are problems reusing the same ifstream object
        // after the first "file open" fails
        // fileIn.open(strFullPathToFile.c_str(),std::ios_base::in);
        std::ifstream fileIn2(strFullPathToFile.c_str(), std::ios_base::in);
        if (fileIn2) {
            fileIn2.close();
            return strFullPathToFile;
        } else {
            return GetRbtDirName(strSubdir) + "/" + strFile;
        }
    }
}

// GetFileType
// Returns the string following the last "." in the file name.
// e.g. GetFileType("receptor.psf") would return "psf"
// If no "." is present, returns the whole file name
RbtString Rbt::GetFileType(const RbtString& strFile) { return strFile.substr(strFile.rfind(".") + 1); }

// Rbt::GetDirList
// Returns a list of files in a directory (strDir) whose names begin with strFilePrefix (optional)
// and whose type is strFileType (optional, as returned by GetFileType)
RbtStringList Rbt::GetDirList(const RbtString& strDir, const RbtString& strFilePrefix, const RbtString& strFileType) {
    RbtStringList dirList;
    RbtBool bMatchPrefix = (strFilePrefix.size() > 0);  // Check if we need to match on file prefix
    RbtBool bMatchType = (strFileType.size() > 0);      // Check if we need to match on file type (suffix)

    DIR* pDir = opendir(strDir.c_str());  // Open directory
    if (pDir != NULL) {
        // DM 6 Dec 1999 - dirent_t appears to be SGI MIPSPro specific
        // At least on Linux g++, it is called dirent
#ifdef __sgi
        dirent_t* pEntry;
#else
        dirent* pEntry;
#endif

        while ((pEntry = readdir(pDir)) != NULL) {  // Read each entry
            RbtString strFile = pEntry->d_name;
            if ((strFile == ".") || (strFile == "..")) continue;  // Don't need to consider . and ..
            RbtBool bMatch = true;
            if (bMatchPrefix && (strFile.find(strFilePrefix) != 0))  // Eliminate those that don't match the prefix
                bMatch = false;
            if (bMatchType && (GetFileType(strFile) != strFileType))  // Eliminate those that don't match the type
                bMatch = false;
            if (bMatch)  // Only store the hits
                dirList.push_back(strFile);
        }
        closedir(pDir);
    }

    // Sort the file list into alphabetical order
    std::sort(dirList.begin(), dirList.end());

    return dirList;
}

// Converts (comma)-delimited string of segment names to segment map
RbtSegmentMap Rbt::ConvertStringToSegmentMap(const RbtString& strSegments, const RbtString& strDelimiter) {
    RbtString::size_type nDelimiterSize = strDelimiter.size();
    RbtSegmentMap segmentMap;

    // Check for null string or null delimiter
    if ((strSegments.size() > 0) && (nDelimiterSize > 0)) {
        RbtString::size_type iBegin = 0;  // indicies into string (sort-of iterators)
        RbtString::size_type iEnd;
        // Not often a do..while loop is used, but in this case it's what we want
        // Even if no delimiter is present, we still need to extract the whole string
        do {
            iEnd = strSegments.find(strDelimiter, iBegin);
            segmentMap[strSegments.substr(iBegin, iEnd - iBegin)] = 0;
            iBegin = iEnd + nDelimiterSize;
        } while (iEnd != RbtString::npos);  // This is the check for delimiter not found
    }

    return segmentMap;
}

// Converts segment map to (comma)-delimited string of segment names
RbtString Rbt::ConvertSegmentMapToString(const RbtSegmentMap& segmentMap, const RbtString& strDelimiter) {
    RbtString strSegments;

    // Check for empty segment map
    if (segmentMap.size() > 0) {
        RbtSegmentMapConstIter iter = segmentMap.begin();
        strSegments += (*iter++).first;  // Add first string
        // Now loop over remaining entries, adding delimiter before each string
        while (iter != segmentMap.end()) {
            strSegments += strDelimiter;
            strSegments += (*iter++).first;
        }
    }
    return strSegments;
}

// Returns a segment map containing the members of map1 which are not in map2
// I know, should really be a template so as to be more universal...one day maybe.
// Or maybe there is already an STL algorithm for doing this.
RbtSegmentMap Rbt::SegmentDiffMap(const RbtSegmentMap& map1, const RbtSegmentMap& map2) {
    RbtSegmentMap map3 = map1;  // Init return value to map1
    for (RbtSegmentMapConstIter iter = map2.begin(); iter != map2.end(); iter++)
        map3.erase((*iter).first);  // Now delete everything in map2
    return map3;
}

// DM 30 Mar 1999
// Converts (comma)-delimited string to string list (similar to ConvertStringToSegmentMap, but returns list not map)
RbtStringList Rbt::ConvertDelimitedStringToList(const RbtString& strValues, const RbtString& strDelimiter) {
    RbtString::size_type nDelimiterSize = strDelimiter.size();
    RbtStringList listOfValues;

    // Check for null string or null delimiter
    if ((strValues.size() > 0) && (nDelimiterSize > 0)) {
        RbtString::size_type iBegin = 0;  // indicies into string (sort-of iterators)
        RbtString::size_type iEnd;
        // Not often a do..while loop is used, but in this case it's what we want
        // Even if no delimiter is present, we still need to extract the whole string
        do {
            iEnd = strValues.find(strDelimiter, iBegin);
            RbtString strValue = strValues.substr(iBegin, iEnd - iBegin);
            listOfValues.push_back(strValue);
            iBegin = iEnd + nDelimiterSize;
        } while (iEnd != RbtString::npos);  // This is the check for delimiter not found
    }
    return listOfValues;
}

// Converts string list to (comma)-delimited string (inverse of above)
RbtString Rbt::ConvertListToDelimitedString(const RbtStringList& listOfValues, const RbtString& strDelimiter) {
    RbtString strValues;

    // Check for empty string list
    if (!listOfValues.empty()) {
        RbtStringListConstIter iter = listOfValues.begin();
        strValues += *iter++;  // Add first string
        // Now loop over remaining entries, adding delimiter before each string
        while (iter != listOfValues.end()) {
            strValues += strDelimiter;
            strValues += *iter++;
        }
    }
    return strValues;
}

////////////////////////////////////////////////////////////////
// I/O ROUTINES
//
// PrintStdHeader
// Print a standard header to an output stream (for log files etc)
// Contains copyright info, library version, date, time etc
// DM 19 Feb 1999 - include executable information
std::ostream& Rbt::PrintStdHeader(std::ostream& s, const RbtString& strExecutable) {
    s << "***********************************************" << std::endl;
    s << GetCopyright() << std::endl;
    if (!strExecutable.empty()) s << "Executable:\t" << strExecutable << std::endl;
    s << "Library:\t" << GetProduct() << "/" << GetVersion() << "/" << GetBuild() << std::endl;
    s << "RBT_ROOT:\t" << GetRbtRoot() << std::endl;
    s << "RBT_HOME:\t" << GetRbtHome() << std::endl;
    s << "Current dir:\t" << GetCurrentDirectory() << std::endl;
    s << "Date:\t\t" << GetTime();
    s << "***********************************************" << std::endl;
    return s;
}

// Used to read RbtCoord. The separator between x y z should be a
// ',', but this takes care of small mistakes, reading any white
// space or commas there is between each variable.
// If necessary, it can be modified to accept the ',' as a
// parameter to be able to use it when other separators are
// needed. Also, it should be possible to implemente it as a
// manipulator of istream.
istream& eatSeps(istream& is) {
    char c;
    while (is.get(c)) {
        if (!isspace(c) && (c != ',')) {
            is.putback(c);
            break;
        }
    }
    return is;
}

void Rbt::ValidateTitle(istream& istr, const std::string& className) {
    std::string title;
    bin_read(istr, title);
    if (title != className) {
        throw RbtFileParseError(
            _WHERE_, "Invalid title string for class " + className + " while deserializing: " + title
        );
    }
}