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

//Rbt docking exceptions, errors during the docking run (e.g. early
//convergence of the GA population)
         
#ifndef _RBTDOCKINGERROR_H_
#define _RBTDOCKINGERROR_H_

#include "RbtError.h"

const RbtString IDS_DOCKING_ERROR            = "RBT_DOCKING_ERROR";

//Unspecified model error
class RbtDockingError : public RbtError
{
 public:
  RbtDockingError(const RbtString& strFile, RbtInt nLine, const RbtString& strMessage="") :
    RbtError(IDS_DOCKING_ERROR,strFile,nLine,strMessage) {}
  //Protected constructor to allow derived docking error classes 
  //to set error name
 protected:
  RbtDockingError(const RbtString& strName,const RbtString& strFile, RbtInt nLine, const RbtString& strMessage="") :
    RbtError(strName,strFile,nLine,strMessage) {}
};

#endif //_RBTDOCKINGERROR_H_
