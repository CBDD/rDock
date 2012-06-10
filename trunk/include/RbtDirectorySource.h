/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtDirectorySource.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

#ifndef _RBTDIRECTORYSOURCE_H_
#define _RBTDIRECTORYSOURCE_H_

#include <sys/stat.h>       // for stat()
#include <dirent.h>         // for scandir()
#include <sys/types.h>
#include <sys/dir.h>
#include <fstream>          // ditto
using std::ifstream;
#include "RbtConfig.h"
#include "RbtFileError.h"	// for exception types


const unsigned int PATH_SIZE  = 1024;	// max path width
const unsigned int LINE_WIDTH = 256;	// max num of chars in a line
const unsigned int ALL        = 0;		// read all files in directory

/**
 * RbtDirectorySource
 * Reads all (or the given number) of files into a <vector>
 * that are in a directory. To extend with other file formasts
 * include other ReadFiles () methods
 */
class RbtDirectorySource {
	protected:
		struct 	stat	fStat;					// stat struct for file access check
#ifdef COMMENT
#ifdef linux
		struct dirent **fNameList;              // file name list in directory
#endif
#ifdef sun										// for gcc 
		struct direct	**fNameList;			// file name list in directory
#endif
#endif //COMMENT
		int             fNum;                   // number of files in the dir
		struct dirent **fNameList;              // file name list in directory
		ifstream        inFile;                 // an actual file to read
		RbtString		thePath;				// leading path 	
		void CheckDirectory(RbtString&);		// to check directory access

		
	public:
		RbtDirectorySource(const RbtString&) throw (RbtError);
		virtual ~RbtDirectorySource();
		
		static RbtString	_CT;
		
		/**
		 * ReadFiles should be re-defined by derived classes to match argument
		 * Ie ReadFile(vector<string> anStrVect) reads file and stores result in anStrVect,
		 * or ReadFile(map< ... > aMap) in aMap, etc.
		 *
		 */
		//virtual void ReadFiles()=0;				
		//virtual void ParseLines()=0;
};

#endif // _RBTDIRECTORYSOURCE_H_
