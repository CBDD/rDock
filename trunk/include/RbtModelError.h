/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtModelError.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
