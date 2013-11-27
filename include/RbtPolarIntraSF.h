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
  
