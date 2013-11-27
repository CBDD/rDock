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
