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

//Rbt file-based exceptions
         
#ifndef _RBTFILEERROR_H_
#define _RBTFILEERROR_H_

#include "RbtError.h"

const RbtString IDS_FILE_ERROR             = "RBT_FILE_ERROR";
const RbtString IDS_FILE_READ_ERROR        = "RBT_FILE_READ_ERROR";
const RbtString IDS_FILE_WRITE_ERROR       = "RBT_FILE_WRITE_ERROR";
const RbtString IDS_FILE_PARSE_ERROR       = "RBT_FILE_PARSE_ERROR";
const RbtString IDS_FILE_MISSING_PARAMETER = "RBT_FILE_MISSING_PARAMETER";
const RbtString IDS_DIR_NOACCESS           = "RBT_DIR_NOACCESS";
const RbtString IDS_NO_FILEINDIR           = "RBT_NO_FILEINDIR";
const RbtString IDS_NOENV                  = "RBT_NO_ENVIRONMENT";
const RbtString IDS_STRINGTOOLONG          = "RBT_SRINGTOOLONG";

//Unspecified file error
class RbtFileError : public RbtError
{
 public:
  RbtFileError(const RbtString& strFile, RbtInt nLine, const RbtString& strMessage="") :
    RbtError(IDS_FILE_ERROR,strFile,nLine,strMessage) {}
  //Protected constructor to allow derived file error classes to set error name
 protected:
  RbtFileError(const RbtString& strName,const RbtString& strFile, RbtInt nLine, const RbtString& strMessage="") :
    RbtError(strName,strFile,nLine,strMessage) {}
};

//File read error
class RbtFileReadError : public RbtFileError
{
 public:
  RbtFileReadError(const RbtString& strFile, RbtInt nLine, const RbtString& strMessage="") :
    RbtFileError(IDS_FILE_READ_ERROR,strFile,nLine,strMessage) {}
};

//File write error
class RbtFileWriteError : public RbtFileError
{
 public:
  RbtFileWriteError(const RbtString& strFile, RbtInt nLine, const RbtString& strMessage="") :
    RbtFileError(IDS_FILE_WRITE_ERROR,strFile,nLine,strMessage) {}
};

//File parse error
class RbtFileParseError : public RbtFileError
{
 public:
  RbtFileParseError(const RbtString& strFile, RbtInt nLine, const RbtString& strMessage="") :
    RbtFileError(IDS_FILE_PARSE_ERROR,strFile,nLine,strMessage) {}
};

//Missing parameter error
class RbtFileMissingParameter : public RbtFileError
{
 public:
  RbtFileMissingParameter(const RbtString& strFile, RbtInt nLine, const RbtString& strMessage="") :
    RbtFileError(IDS_FILE_MISSING_PARAMETER,strFile,nLine,strMessage) {}
};

// Exceptions for RbtDirectorySource and derived classes:
// 
// Directory access error
class RbtDirIsNotAccessible : public RbtFileError
{
	public:
		RbtDirIsNotAccessible (const RbtString& strFile, RbtInt nLine, const RbtString& strMessage=""):
    		RbtFileError(IDS_DIR_NOACCESS,strFile,nLine,strMessage) {}
};

// Directory contains no file with that extension
class RbtNoFileInDir : public RbtFileError
{
	public:
		RbtNoFileInDir (const RbtString& strFile, RbtInt nLine, const RbtString& strMessage=""):
    		RbtFileError(IDS_NO_FILEINDIR,strFile,nLine,strMessage) {}
};

// Misc exceptions
//
// If path is missing
class RbtEnvNotDefined : public RbtFileError
{
	public:
		RbtEnvNotDefined (const RbtString& strFile, RbtInt nLine, const RbtString& strMessage=""):
			RbtFileError(IDS_NOENV,strFile,nLine,strMessage) {}
};

// If string is too long 
class RbtStringTooLong : public RbtFileError
{
	public:
		RbtStringTooLong (const RbtString& strFile, RbtInt nLine, const RbtString& strMessage=""):
			RbtFileError(IDS_STRINGTOOLONG,strFile,nLine,strMessage) {}
};

#endif //_RBTFILEERROR_H_
