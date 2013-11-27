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

//#include <stdlib.h>				// for getenv
//#include "RbtFileError.h"
#include <iostream>
#include "RbtDirectorySource.h"
using std::cout;
using std::endl;

RbtString RbtDirectorySource::_CT("RbtDirectorySource");

/**
 * Constructor checks wether the directory exists and gives
 * $RBT_ROOT to the path if required
 */
RbtDirectorySource::RbtDirectorySource(const RbtString& aDirectory) throw (RbtError)
{
//	RbtString theFileStr;								// filename without path after scandir()
//	RbtString theLine;									// line read

	thePath = (RbtString) getenv("RBT_ROOT");			// path string to the dir
	if ( thePath.empty() )								// $RBT_ROOT not defined ...
		throw RbtEnvNotDefined(_WHERE_, "$RBT_ROOT " );
	else if ( aDirectory.size() > PATH_SIZE )			// or string is too long
		throw RbtStringTooLong(_WHERE_," ");
	else if ( stat(aDirectory.c_str(),&fStat) < 0 )		// check in the cwd (or absolute path)
			thePath = aDirectory;						// well, we want to read from cwd/abs path
	else {
		//thePath =+ "/"+aDirectory;					// concatenate
		thePath = aDirectory;							// providing $RBT_ROOT is OK, we do not have to give /
		if (thePath.size() > PATH_SIZE)					// again too long
			throw RbtStringTooLong(_WHERE_, "");
		else if ( stat(thePath.c_str(),&fStat) < 0 )	// if not even in the extended path
			throw RbtDirIsNotAccessible(_WHERE_,"");
	}
	// is it _really_ a directory or just an entry?
	CheckDirectory(thePath);							// throws exceptions if something messy

	
}

void RbtDirectorySource::CheckDirectory(RbtString& aDir)
{
//#ifdef COMENT
	if (!S_ISDIR(fStat.st_mode))
		throw RbtDirIsNotAccessible(_WHERE_,aDir+" : Not a directory. ");
	char *theDir = (char *) aDir.c_str();			// unbelivable annoying buglet in Solaris gcc 2.95.3
	fNum = scandir(theDir,&fNameList,0,alphasort); 	// get dir entries
//	cout << "Not running on Slowaris" << endl;
//	exit(1);
	if (fNum < 1)
		throw RbtNoFileInDir(_WHERE_," No files in directory. ");
	else
		cout << _CT << " Reading files from " << aDir << endl;
//#endif
}

RbtDirectorySource::~RbtDirectorySource()
{
	cout << "RbtDirectorySource destructor" << endl;
}
