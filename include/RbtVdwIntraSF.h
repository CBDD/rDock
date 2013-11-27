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

//Ligand intramolecular scoring function for vdw potential

#ifndef _RBTVDWINTRASF_H_
#define _RBTVDWINTRASF_H_

#include "RbtBaseIntraSF.h"
#include "RbtVdwSF.h"


class RbtVdwIntraSF : public RbtBaseIntraSF, public RbtVdwSF
{
 public:
  //Class type string
  static RbtString _CT;

  RbtVdwIntraSF(const RbtString& strName = "VDW");
  virtual ~RbtVdwIntraSF();
  
  //Request Handling method
  //Handles the Partition request
  virtual void HandleRequest(RbtRequestPtr spRequest);

 protected:
  virtual void SetupScore();
  virtual RbtDouble RawScore() const;
  
  //DM 25 Oct 2000 - track changes to parameter values in local data members
  //ParameterUpdated is invoked by RbtParamHandler::SetParameter
  void ParameterUpdated(const RbtString& strName);
  
 private:
  RbtAtomRListList m_vdwIntns;//The full list of vdW interactions
  RbtAtomRListList m_prtIntns;//The partitioned interactions (within partition distance)
  RbtAtomRList m_ligAtomList;
};

#endif //_RBTVDWINTRASF_H_
  
