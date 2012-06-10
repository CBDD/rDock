/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtRotSF.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
