/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtPolarIntraSF.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Ligand intramolecular scoring function for all attractive polar
//interactions (HBD,HBA,metal,guanidinium carbon)

#ifndef _RBTPOLARINTRASF_H_
#define _RBTPOLARINTRASF_H_

#include "RbtBaseIntraSF.h"
#include "RbtPolarSF.h"
#include "RbtInteractionGrid.h"

class RbtPolarIntraSF : public RbtBaseIntraSF, public RbtPolarSF
{
 public:
  //Class type string
  static RbtString _CT;
  //Parameter names
  static RbtString _ATTR;
  
  RbtPolarIntraSF(const RbtString& strName = "POLAR");
  virtual ~RbtPolarIntraSF();
  
 protected:
  virtual void SetupScore();
  virtual RbtDouble RawScore() const;
  
  //Request Handling method
  //Handles the Partition request
  virtual void HandleRequest(RbtRequestPtr spRequest);

  //Clear the model lists
  //As we are not using smart pointers, there is some memory management to do
  void ClearModel();
  
  //DM 25 Oct 2000 - track changes to parameter values in local data members
  //ParameterUpdated is invoked by RbtParamHandler::SetParameter
  void ParameterUpdated(const RbtString& strName);
  
 private:
  RbtInteractionCenterList m_posList;
  RbtInteractionCenterList m_negList;
  RbtInteractionListMap m_intns;
  RbtInteractionListMap m_prtIntns;
  RbtBool m_bAttr;
};

#endif //_RBTPOLARINTRASF_H_
  
