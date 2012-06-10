/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtNmrSF.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Generic nmr restraint scoring function

#ifndef _RBTNMRSF_H_
#define _RBTNMRSF_H_

#include "RbtBaseInterSF.h"
#include "RbtBaseIdxSF.h"
#include "RbtNoeRestraint.h"

class RbtNmrSF : public RbtBaseInterSF, public RbtBaseIdxSF
{
 public:
  //Class type string
  static RbtString _CT;
  //Parameter names
  static RbtString _FILENAME;//Nmr restraint file name
  static RbtString _QUADRATIC;//True = quadratic penalty function; false = linear
  
  RbtNmrSF(const RbtString& strName = "NMR");
  virtual ~RbtNmrSF();
  
 protected:
  virtual void SetupReceptor();
  virtual void SetupLigand();
  virtual void SetupScore();
  virtual RbtDouble RawScore() const;
  void ParameterUpdated(const RbtString& strName);
  
 private:
  RbtDouble NoeDistance(const RbtNoeRestraintAtoms& noe) const;
  RbtDouble StdDistance(const RbtStdRestraintAtoms& std) const;

  RbtBool m_bQuadratic;//synchronised with QUADRATIC named parameter
  RbtNonBondedGridPtr m_spGrid;
  RbtAtomList m_ligAtomList;//All ligand atoms
  RbtNoeRestraintAtomsList m_noeList;//List of all NOE interactions
  RbtStdRestraintAtomsList m_stdList;//List of all STD interactions
};

#endif //_RBTNMRSF_H_
