/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtVdwIntraSF.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtVdwIntraSF.h"
#include "RbtSFRequest.h"

//Static data members
RbtString RbtVdwIntraSF::_CT("RbtVdwIntraSF");

//NB - Virtual base class constructor (RbtBaseSF) gets called first,
//implicit constructor for RbtBaseInterSF is called second
RbtVdwIntraSF::RbtVdwIntraSF(const RbtString& strName) : RbtBaseSF(_CT,strName)
{
#ifdef _DEBUG
  cout << _CT << " parameterised constructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtVdwIntraSF::~RbtVdwIntraSF() {
#ifdef _DEBUG
  cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

//Request Handling method
//Handles the Partition request
void RbtVdwIntraSF::HandleRequest(RbtRequestPtr spRequest) {
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
      Partition(m_ligAtomList,m_vdwIntns,m_prtIntns,params[0]);
    }
    else if ( (params.size() == 2) && (params[0].String() == GetFullName())) {
      if (iTrace > 2) {
      cout << _CT << "::HandleRequest: Partitioning " << GetFullName() << " at distance=" << params[1] << endl;
      }
      Partition(m_ligAtomList,m_vdwIntns,m_prtIntns,params[1]);
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

void RbtVdwIntraSF::SetupScore() {
  m_ligAtomList.clear();
  for (RbtAtomRListListIter iter = m_vdwIntns.begin(); iter != m_vdwIntns.end() ; iter++)
    (*iter).clear();
  m_vdwIntns.clear();
  for (RbtAtomRListListIter iter = m_prtIntns.begin(); iter != m_prtIntns.end() ; iter++)
    (*iter).clear();
  m_prtIntns.clear();

  RbtModelPtr spModel = GetLigand();
  if (spModel.Null())
    return;

  RbtInt iTrace = GetTrace();
  if (iTrace > 2) {
    cout << _CT << "::SetupScore()" << endl;
  }

  RbtAtomList tmpList = spModel->GetAtomList();
  //Strip off the smart pointers
  std::copy(tmpList.begin(),tmpList.end(),std::back_inserter(m_ligAtomList));

  //Build map of intra-ligand flexible interactions
  m_vdwIntns = RbtAtomRListList(m_ligAtomList.size(),RbtAtomRList());
  m_prtIntns = RbtAtomRListList(m_ligAtomList.size(),RbtAtomRList());
  BuildIntraMap(m_ligAtomList,m_vdwIntns);
  //Partition with zero distance is needed to copy all the vdW interactions
  //into the partitioned list (this is the list that is scored)
  Partition(m_ligAtomList,m_vdwIntns,m_prtIntns,0.0);

}

RbtDouble RbtVdwIntraSF::RawScore() const {
  RbtDouble score = 0.0;//Total score
  //Loop over all ligand atoms
  for (RbtAtomRListConstIter iter = m_ligAtomList.begin(); iter != m_ligAtomList.end(); iter++) {
    RbtInt id = (*iter)->GetAtomId()-1;
    RbtDouble s = VdwScore(*iter,m_prtIntns[id]);
    score += s;
  }
  return score;
}


//DM 25 Oct 2000 - track changes to parameter values in local data members
//ParameterUpdated is invoked by RbtParamHandler::SetParameter
void RbtVdwIntraSF::ParameterUpdated(const RbtString& strName) {
  RbtVdwSF::OwnParameterUpdated(strName);
  RbtBaseSF::ParameterUpdated(strName);
}


