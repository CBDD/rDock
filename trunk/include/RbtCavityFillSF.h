/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtCavityFillSF.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
