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

//Post-filtering scoring function for calculating percentage cavity unfilled
//(mode=0) or percentage ligand unbound (mode=1)

#ifndef _RBTCAVITYFILLSF_H_
#define _RBTCAVITYFILLSF_H_

#include "RbtBaseInterSF.h"
#include "RbtFFTGrid.h"

class RbtCavityFillSF : public RbtBaseInterSF
{
 public:
  //Class type string
  static RbtString _CT;
  //Parameter names
  RbtCavityFillSF(const RbtString& strName = "CAVFILL");
  virtual ~RbtCavityFillSF();
  
 protected:
  virtual void SetupReceptor();
  virtual void SetupLigand();
  virtual void SetupScore();
  virtual RbtDouble RawScore() const;
  void ParameterUpdated(const RbtString& strName);
  
 private:
  RbtFFTGridPtr m_spGrid;
  RbtAtomList m_ligAtomList;
};

#endif //_RBTCAVITYFILLSF_H_
