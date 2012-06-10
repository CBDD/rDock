/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtLigandError.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Rbt ligand exceptions, the model is correct, but the ligand is not
//appropiate for this case. (e.g. doesn't comply with all the pharmacophore
//constraints)

#ifndef _RBTLIGANDERROR_H_
#define _RBTLIGANDERROR_H_

#include "RbtError.h"

const RbtString IDS_LIGAND_ERROR            = "RBT_LIGAND_ERROR";

//Unspecified ligand error
class RbtLigandError : public RbtError
{
 public:
  RbtLigandError(const RbtString& strFile, RbtInt nLine, const RbtString& strMessage="") :
    RbtError(IDS_LIGAND_ERROR,strFile,nLine,strMessage) {}
  //Protected constructor to allow derived ligand error classes to set error name
 protected:
  RbtLigandError(const RbtString& strName,const RbtString& strFile, RbtInt nLine, const RbtString& strMessage="") :
    RbtError(strName,strFile,nLine,strMessage) {}
};

#endif //_RBTLIGANDERROR_H_
