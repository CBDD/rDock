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
  
