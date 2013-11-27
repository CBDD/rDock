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

//Miscellaneous non-member functions in Rbt namespace                      
#ifndef _RBT_H_
#define _RBT_H_

#include "RbtTypes.h"
#include "RbtContainers.h"
#include "RbtError.h"

namespace Rbt
{
  ////////////////////////////////////////////////////////////////
  //RESOURCE HANDLING FUNCTIONS
  //
  //GetRbtRoot - returns value of RBT_ROOT env variable
  RbtString GetRbtRoot();
  //GetRbtHome - returns value of RBT_HOME env variable
  //(or HOME if RBT_HOME is undefined)
  RbtString GetRbtHome();
  //GetCopyright - returns legalese statement
  RbtString GetCopyright();
  //GetVersion - returns current library version
  RbtString GetVersion();
  //GetBuildNum - returns current library build number
  RbtString GetBuild();
  //GetProduct - returns library product name
  RbtString GetProduct();
  //GetTime - returns current time and date as an RbtString
  RbtString GetTime();
  //GetCurrentDirectory - returns current working directory
  RbtString GetCurrentDirectory();
  //
  ////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////
  //FILE/DIRECTORY HANDLING FUNCTIONS
  //
  //GetRbtDirName
  //Returns the full path to a subdirectory in the rDock directory structure
  //
  //For example, if RBT_ROOT environment variable is ~dave/ribodev/molmod/ribodev
  //then GetRbtDirName("data") would return ~dave/ribodev/molmod/ribodev/data/
  //
  //If RBT_ROOT is not set, then GetRbtDirName returns ./ irrespective of the
  //subdirectory asked for. Thus the fall-back position is that parameter files
  //are read from the current directory if RBT_ROOT is not defined.
  RbtString GetRbtDirName(const RbtString& strSubdir="");

  //GetRbtFileName
  //As GetRbtDirName but returns the full path to a file in the rDock directory structure
  RbtString GetRbtFileName(const RbtString& strSubdir, const RbtString& strFile);

  //GetFileType
  //Returns the string following the last "." in the file name.
  //e.g. GetFileType("receptor.psf") would return "psf"
  RbtString GetFileType(const RbtString& strFile);

  //GetDirList
  //Returns a list of files in a directory (strDir) whose names begin with strFilePrefix (optional)
  //and whose type is strFileType (optional, as returned by GetFileType)
  RbtStringList GetDirList(const RbtString& strDir, const RbtString& strFilePrefix="", const RbtString& strFileType="");
  //
  ////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////
  //CONVERSION ROUTINES
  //
  //Converts (comma)-delimited string of segment names to segment map
  RbtSegmentMap ConvertStringToSegmentMap(const RbtString& strSegments, const RbtString& strDelimiter=",");
  //Converts segment map to (comma)-delimited string of segment names
  RbtString ConvertSegmentMapToString(const RbtSegmentMap& segmentMap, const RbtString& strDelimiter=",");

  //Returns a segment map containing the members of map1 which are not in map2
  //I know, should really be a template so as to be more universal...one day maybe.
  //Or maybe there is already an STL algorithm for doing this.
  RbtSegmentMap SegmentDiffMap(const RbtSegmentMap& map1, const RbtSegmentMap& map2);

  //DM 30 Mar 1999
  //Converts (comma)-delimited string to string list (similar to ConvertStringToSegmentMap, but returns list not map)
  RbtStringList ConvertDelimitedStringToList(const RbtString& strValues, const RbtString& strDelimiter=",");
  //Converts string list to (comma)-delimited string (inverse of above)
  RbtString ConvertListToDelimitedString(const RbtStringList& listOfValues, const RbtString& strDelimiter=",");

  //
  ////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////
  //I/O ROUTINES
  //
  //PrintStdHeader
  //Print a standard header to an output stream (for log files etc)
  //Contains copyright info, library version, date, time etc
  //DM 19 Feb 1999 - include executable information
  std::ostream& PrintStdHeader(std::ostream& s, const RbtString& strExecutable="");

  //Helper functions to read/write chars from iostreams
  //Throws error if stream state is not Good() before and after the read/write
  //It appears the STL ios_base exception throwing is not yet implemented
  //at least on RedHat 6.1, so this is a temporary workaround (yeah right)
  //
  //DM 25 Sep 2000 - with MIPSPro CC compiler, streamsize is not defined,
  //at least with the iostreams library we are using, so typedef it here
  // it is not in gcc 3.2.1 either SJ
//#ifdef __sgi
  typedef int streamsize;
//#endif

  void WriteWithThrow(std::ostream& ostr, const char* p, streamsize n) throw (RbtError);
  void ReadWithThrow(std::istream& istr, char* p, streamsize n) throw (RbtError);

}

#endif //_RBT_H_
