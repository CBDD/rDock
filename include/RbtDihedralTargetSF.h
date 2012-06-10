/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtDihedralTargetSF.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Target intramolecular dihedral scoring function
//for flexible receptors
#ifndef _RBTDIHEDRALTARGETSF_H_
#define _RBTDIHEDRALTARGETSF_H_

#include "RbtBaseInterSF.h"
#include "RbtDihedralSF.h"

class RbtDihedralTargetSF : public RbtBaseInterSF, public RbtDihedralSF
{
 public:
  //Class type string
  static RbtString _CT;
  //Parameter names
  
  RbtDihedralTargetSF(const RbtString& strName = "DIHEDRAL");
  virtual ~RbtDihedralTargetSF();
  
 protected:
  virtual void SetupReceptor();
  virtual void SetupLigand();
  virtual void SetupScore();
  virtual RbtDouble RawScore() const;
  
  //Clear the dihedral list
  //As we are not using smart pointers, there is some memory management to do
  void ClearReceptor();

 private:
  RbtDihedralList m_dihList;
};

#endif //_RBTDIHEDRALTARGETSF_H_
  
