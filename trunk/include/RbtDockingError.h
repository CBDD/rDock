/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtDockingError.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
