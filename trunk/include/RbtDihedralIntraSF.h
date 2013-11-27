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
  
