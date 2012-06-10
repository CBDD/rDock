/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtDihedralIntraSF.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Ligand intramolecular dihedral scoring function
#ifndef _RBTDIHEDRALINTRASF_H_
#define _RBTDIHEDRALINTRASF_H_

#include "RbtBaseIntraSF.h"
#include "RbtDihedralSF.h"

class RbtDihedralIntraSF : public RbtBaseIntraSF, public RbtDihedralSF
{
 public:
  //Class type string
  static RbtString _CT;
  //Parameter names
  
  RbtDihedralIntraSF(const RbtString& strName = "DIHEDRAL");
  virtual ~RbtDihedralIntraSF();
  
 protected:
  virtual void SetupScore();
  virtual RbtDouble RawScore() const;
  
  //Clear the dihedral list
  //As we are not using smart pointers, there is some memory management to do
  void ClearModel();

 private:
  RbtDihedralList m_dihList;
};

#endif //_RBTDIHEDRALINTRASF_H_
  
