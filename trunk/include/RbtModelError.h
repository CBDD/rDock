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

//Rbt model exceptions (e.g. topological errors)
         
#ifndef _RBTMODELERROR_H_
#define _RBTMODELERROR_H_

#include "RbtError.h"

const RbtString IDS_MODEL_ERROR            = "RBT_MODEL_ERROR";

//Unspecified model error
class RbtModelError : public RbtError
{
 public:
  RbtModelError(const RbtString& strFile, RbtInt nLine, const RbtString& strMessage="") :
    RbtError(IDS_MODEL_ERROR,strFile,nLine,strMessage) {}
  //Protected constructor to allow derived model error classes to set error name
 protected:
  RbtModelError(const RbtString& strName,const RbtString& strFile, RbtInt nLine, const RbtString& strMessage="") :
    RbtError(strName,strFile,nLine,strMessage) {}
};

#endif //_RBTMODELERROR_H_
