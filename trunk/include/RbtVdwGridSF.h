/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtVdwGridSF.h#4 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Precalculated-grid-based intermolecular vdw scoring function

#ifndef _RBTVDWGRIDSF_H_
#define _RBTVDWGRIDSF_H_

#include "RbtBaseInterSF.h"
#include "RbtRealGrid.h"

class RbtVdwGridSF : public RbtBaseInterSF
{
 public:
  //Class type string
  static RbtString _CT;
  //Parameter names
  static RbtString _GRID;//Suffix for grid filename
  static RbtString _SMOOTHED;//Controls whether to smooth the grid values
  
  RbtVdwGridSF(const RbtString& strName = "VDW");
  virtual ~RbtVdwGridSF();
  
 protected:
  virtual void SetupReceptor();
  virtual void SetupLigand();
  virtual void SetupSolvent();
  virtual void SetupScore();
  virtual RbtDouble RawScore() const;
  //DM 25 Oct 2000 - track changes to parameter values in local data members
  //ParameterUpdated is invoked by RbtParamHandler::SetParameter
  void ParameterUpdated(const RbtString& strName);
  
 private:
  //Read grids from input stream
  void ReadGrids(istream& istr) throw (RbtError);
  
  RbtRealGridList m_grids;
  RbtAtomRList m_ligAtomList;
  RbtTriposAtomTypeList m_ligAtomTypes;
  RbtBool m_bSmoothed;
};

#endif //_RBTVDWGRIDSF_H_
