/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtPolarIntraSF.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtPolarIntraSF.h"
#include "RbtSFRequest.h"

//Static data members
RbtString RbtPolarIntraSF::_CT("RbtPolarIntraSF");
RbtString RbtPolarIntraSF::_ATTR("ATTR");

//NB - Virtual base class constructor (RbtBaseSF) gets called first,
//implicit constructor for RbtBaseInterSF is called second
RbtPolarIntraSF::RbtPolarIntraSF(const RbtString& strName) : RbtBaseSF(_CT,strName),m_bAttr(true)
{
  //Add parameters
  AddParameter(_ATTR,m_bAttr);
#ifdef _DEBUG
  cout << _CT << " parameterised constructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtPolarIntraSF::~RbtPolarIntraSF() {
  ClearModel();
#ifdef _DEBUG
  cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtPolarIntraSF::SetupScore() {
  ClearModel();
  RbtModelPtr spModel = GetLigand();
  if (spModel.Null())
    return;

  RbtAtomList atomList(spModel->GetAtomList());
  m_posList = CreateDonorInteractionCenters(atomList);
  m_negList = CreateAcceptorInteractionCenters(atomList);

  RbtInt nAtoms = atomList.size();
  m_intns = RbtInteractionListMap(nAtoms,RbtInteractionCenterList());
  m_prtIntns = RbtInteractionListMap(nAtoms,RbtInteractionCenterList());
  if (m_bAttr) {
    BuildIntraMap(m_posList,m_negList,m_intns);
  }
  else {
    BuildIntraMap(m_posList,m_intns);
    BuildIntraMap(m_negList,m_intns);
  }
  //Partition with zero distance is needed to copy all the polar interactions
  //into the partitioned list (this is the list that is scored)
  Partition(m_posList,m_negList,m_intns,m_prtIntns,0.0);
}

RbtDouble RbtPolarIntraSF::RawScore() const {
  return IntraScore(m_posList,m_negList,m_prtIntns,m_bAttr);
}

void RbtPolarIntraSF::ClearModel() {
  //Clear the interaction maps
  for (RbtInteractionListMapIter iter = m_intns.begin(); iter != m_intns.end(); iter++) {
    (*iter).clear();
  }
  m_intns.clear();
  for (RbtInteractionListMapIter iter = m_prtIntns.begin(); iter != m_prtIntns.end(); iter++) {
    (*iter).clear();
  }
  m_prtIntns.clear();

  //Delete the ligand interaction centers
  for (RbtInteractionCenterListIter iter = m_posList.begin(); iter != m_posList.end(); iter++) {
    delete *iter;
  }
  m_posList.clear();
  for (RbtInteractionCenterListIter iter = m_negList.begin(); iter != m_negList.end(); iter++) {
    delete *iter;
  }
  m_negList.clear();
}

//DM 25 Oct 2000 - track changes to parameter values in local data members
//ParameterUpdated is invoked by RbtParamHandler::SetParameter
void RbtPolarIntraSF::ParameterUpdated(const RbtString& strName) {
  //DM 25 Oct 2000 - heavily used params
  if (strName == _ATTR) {
    m_bAttr = GetParameter(_ATTR);
  }
  else {
    RbtPolarSF::OwnParameterUpdated(strName);
    RbtBaseSF::ParameterUpdated(strName);
  }
}

//Request Handling method
//Handles the Partition request
void RbtPolarIntraSF::HandleRequest(RbtRequestPtr spRequest) {
  RbtVariantList params = spRequest->GetParameters();
  RbtInt iTrace = GetTrace();

  switch (spRequest->GetID()) {
    //ID_REQ_SF_PARTITION requests come in two forms:
    //1-param: param[0] = distance => Partition all scoring functions
    //2-param: param[0] = SF Name,
    //         param[1] = distance => Partition a named scoring function
  case ID_REQ_SF_PARTITION:
    if (params.size() == 1) {
      if (iTrace > 2) {
	cout << _CT << "::HandleRequest: Partitioning " << GetFullName() << " at distance=" << params[0] << endl;
      }
      Partition(m_posList,m_negList,m_intns,m_prtIntns,params[0]);
    }
    else if ( (params.size() == 2) && (params[0].String() == GetFullName())) {
      if (iTrace > 2) {
      cout << _CT << "::HandleRequest: Partitioning " << GetFullName() << " at distance=" << params[1] << endl;
      }
      Partition(m_posList,m_negList,m_intns,m_prtIntns,params[1]);
    }
    break;
    
    //Pass all other requests to base handler
  default:
    if (iTrace > 2) {
      cout << _CT << "::HandleRequest: " << GetFullName() << " passing request to base handler" << endl;
    }
    RbtBaseObject::HandleRequest(spRequest);
    break;
  }
}
