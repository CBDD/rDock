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
