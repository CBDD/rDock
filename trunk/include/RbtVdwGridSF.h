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
