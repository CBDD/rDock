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

//Simple scoring function which estimates the loss of rotational entropy of the
//ligand from the number of rotable bonds.

#ifndef _RBTROTSF_H_
#define _RBTROTSF_H_

#include "RbtBaseInterSF.h"

class RbtRotSF : public RbtBaseInterSF
{
 public:
  //Static data member for class type
  static RbtString _CT;
  //Parameter names
  //Boolean controlling whether to include bonds to NH3+ in rotable bond count
  static RbtString _INCNH3;
  //Boolean controlling whether to include bonds to OH in  rotable bond count
  static RbtString _INCOH;

  RbtRotSF(const RbtString& strName = "ROT");
  virtual ~RbtRotSF();
  
 protected:
  virtual void SetupReceptor();
  virtual void SetupLigand();
  virtual void SetupScore();
  virtual RbtDouble RawScore() const;
  void ParameterUpdated(const RbtString& strName);

 private:
  RbtInt nRot;
  RbtBool bIncNH3;
  RbtBool bIncOH;
};

#endif //_RBTROTSF_H_
