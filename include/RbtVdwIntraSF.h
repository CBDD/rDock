/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtVdwIntraSF.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
  
