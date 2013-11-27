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

//Rbt exception base class
//RbtErrors can be thrown for unspecified errors, or derived classes can be
//used for particular exception types.
//Has the following methods:
//
// Name      - exception name (set by constructor)
// File      - source file where error occurred
// Line      - source line no.
// Message   - description of the error
// isOK      - if true, status is OK, not an error
// AddMessage- append new message to existing message, e.g. when rethrowing
// Insertion operator defined for dumping exception details to a stream
         
#ifndef _RBTERROR_H_
#define _RBTERROR_H_

#include <iostream>
#include "RbtTypes.h"

//Useful macro for reporting errors
//__FILE__ and __LINE__ are standard macros for source code file name
// and line number respectively
#define _WHERE_ __FILE__, __LINE__

//General errors
const RbtString IDS_ERROR					= "RBT_ERROR";
const RbtString IDS_OK						= "RBT_OK";
const RbtString IDS_INVALID_REQUEST			= "RBT_INVALID_REQUEST";
const RbtString IDS_BAD_ARGUMENT			= "RBT_BAD_ARGUMENT";
const RbtString IDS_ASSERT					= "RBT_ASSERT";
const RbtString IDS_BAD_RECEPTOR_FILE		= "BAD_RECEPTOR_FILE";

//Template for throwing exceptions from assert failures
//(The C++ Programming Language, 3rd Ed, p750)
template <class X, class A> inline void Assert(A assertion)
{
	if (!assertion) {
		throw X();
	}
}
/////////////////////////////////////
//BASE ERROR CLASS
//Note: This is a concrete class so can be instantiated
/////////////////////////////////////

class RbtError
{
  ///////////////////////
  //Public methods
 public:

  //Default constructor
  //Use to create an "error" with status=OK, no line, no file, no message
  //All other constructors set the OK flag to false
  RbtError() : m_strName(IDS_OK),m_strFile(""),m_nLine(0),m_strMessage(""),m_bOK(true) {}

  //Parameterised constructor
  //Use to create unspecified Rbt Errors
  RbtError(const RbtString& strFile, RbtInt nLine, const RbtString& strMessage="") :
    m_strName(IDS_ERROR),m_strFile(strFile),m_nLine(nLine),m_strMessage(strMessage),m_bOK(false) {}

  //Default destructor
  virtual ~RbtError() {}

  //Insertion operator
  friend std::ostream& operator<<(std::ostream& s, const RbtError& error) {
    s << error.Name();
    if (!error.File().empty())
      s << " at " << error.File() << ", line " << error.Line();
    if (!error.Message().empty())
      s << std::endl << error.Message();
    return s;
  }

  //Attribute methods
  RbtString File() const {return m_strFile;} //Get filename
  RbtInt Line() const {return m_nLine;} //Get line number
  RbtString Message() const {return m_strMessage;} //Get message
  RbtString Name() const {return m_strName;} //Get error name
  RbtBool isOK() const {return m_bOK;} //If true, status is OK (not an error)

  //Append new message to existing message
  void AddMessage(const RbtString& strMessage) {m_strMessage+=strMessage;}

  ///////////////////////
  //Protected methods
 protected:
  //Protected constructor to allow derived error classes to set error name
  RbtError(const RbtString& strName, const RbtString& strFile, RbtInt nLine, const RbtString& strMessage="") :
    m_strName(strName),m_strFile(strFile),m_nLine(nLine),m_strMessage(strMessage),m_bOK(false) {}

  ///////////////////////
  //Private data
 private:
  RbtString m_strName;
  RbtString m_strFile;
  RbtString m_strMessage;
  RbtInt m_nLine;
  RbtBool m_bOK;
};

/////////////////////////////////////
//DERIVED ERROR CLASSES - General
/////////////////////////////////////

//Invalid request - object does not supported the request
class RbtInvalidRequest : public RbtError
{
 public:
  RbtInvalidRequest(const RbtString& strFile, RbtInt nLine, const RbtString& strMessage="") :
    RbtError(IDS_INVALID_REQUEST,strFile,nLine,strMessage) {}
};
//Bad argument - e.g. empty atom list
class RbtBadArgument : public RbtError
{
 public:
  RbtBadArgument(const RbtString& strFile, RbtInt nLine, const RbtString& strMessage="") :
    RbtError(IDS_BAD_ARGUMENT,strFile,nLine,strMessage) {}
};
//Assertion failure
//Assert template doesn't pass params to exception class so can't pass file and line number
class RbtAssert : public RbtError
{
 public:
  RbtAssert() :
    RbtError(IDS_ASSERT,"unspecified file",0,"Assertion failed") {}
};

// bad receptor type (ie not PSF/CRD either mol2
class RbtBadReceptorFile : public RbtError
{
 public:
  RbtBadReceptorFile(const RbtString& strFile,RbtInt nLine,const RbtString& strMessage="") : 
    RbtError(IDS_BAD_RECEPTOR_FILE,strFile,nLine,strMessage) {}
};

#endif //_RBTERROR_H_
