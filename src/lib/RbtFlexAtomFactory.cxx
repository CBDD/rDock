/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtFlexAtomFactory.cxx#2 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtFlexAtomFactory.h"
#include "RbtReceptorFlexData.h"
#include "RbtLigandFlexData.h"
#include "RbtSolventFlexData.h"
#include "RbtChromElement.h"

RbtFlexAtomFactory::RbtFlexAtomFactory(RbtModel* pModel) {
  Visit(pModel);
}

RbtFlexAtomFactory::RbtFlexAtomFactory(RbtModelList modelList) {
  for (RbtModelListConstIter iter = modelList.begin(); iter != modelList.end(); iter++) {
    Visit(*iter);
  }
}

void RbtFlexAtomFactory::Visit(RbtModel* pModel) {
  if (pModel != NULL) {
    RbtFlexData* pFlexData = pModel->GetFlexData();
    if (pFlexData != NULL) {
      pFlexData->Accept(*this);
    }
    else {
      //If flex data is null, then model is rigid by default
      //Therefore all atoms are fixed
      RbtAtomList atomList = pModel->GetAtomList();
      pModel->SetAtomUser2Values(0.0);
      std::copy(atomList.begin(), atomList.end(), std::back_inserter(m_fixedAtomList));
    }
  }
}

void RbtFlexAtomFactory::Clear() {
  m_fixedAtomList.clear();
  m_tetheredAtomList.clear();
  m_freeAtomList.clear();
}

void RbtFlexAtomFactory::VisitReceptorFlexData(RbtReceptorFlexData* pFlexData) {
  RbtModel* pModel = pFlexData->GetModel();
  RbtAtomList atomList = pModel->GetAtomList();
  //Divide the receptor atoms into fixed and tethered (pendant OH/NH3+)
  pModel->SetAtomSelectionFlags(false);
  pModel->SetAtomUser2Values(0.0);
  pModel->SelectFlexAtoms();
  for (RbtAtomListIter iter = atomList.begin(); iter != atomList.end(); ++iter) {
    if ((*iter)->GetSelectionFlag()) {
      m_tetheredAtomList.push_back(*iter);
      //This code assumes that the only flexible atoms in the receptor are
      //terminal OH/NH3+ hydrogens. Max rotational displacement from current
      //coords is equal to twice the bond length.
      RbtAtomList bondedAtomList = Rbt::GetBondedAtomList(*iter);
      if (bondedAtomList.size() == 1) {
	RbtAtom* parent = bondedAtomList.front();
	RbtDouble bondLength = Rbt::Length(parent->GetCoords(),(*iter)->GetCoords());
	(*iter)->SetUser2Value(2.0 * bondLength);
      }
    }
    else {
      m_fixedAtomList.push_back(*iter);
    }
  }
}

void RbtFlexAtomFactory::VisitLigandFlexData(RbtLigandFlexData* pFlexData) {
  RbtModel* pModel = pFlexData->GetModel();
  RbtAtomList atomList = pModel->GetAtomList();
  pModel->SetAtomUser2Values(0.0);
  RbtString transModeStr = pFlexData->GetParameter(RbtLigandFlexData::_TRANS_MODE);
  RbtString rotModeStr = pFlexData->GetParameter(RbtLigandFlexData::_ROT_MODE);
  RbtChromElement::eMode eTransMode = RbtChromElement::StrToMode(transModeStr);
  RbtChromElement::eMode eRotMode = RbtChromElement::StrToMode(rotModeStr);
  //Internally flexible ligand, or free translation
  if (pModel->isFlexible() || (eTransMode == RbtChromElement::FREE) ) {
    std::copy(atomList.begin(), atomList.end(), std::back_inserter(m_freeAtomList));
  }
  else if (eTransMode == RbtChromElement::FIXED) {
    //Completely fixed
    if (eRotMode == RbtChromElement::FIXED) {
      std::copy(atomList.begin(), atomList.end(), std::back_inserter(m_fixedAtomList));
    }
    //Fixed translation; tethered or free rotation
    //Therefore max translation of any atom depends on its distance from the COM
    //For simplicity, we don't distinguish between tethered and free rotation
    //so assume full rotation through +/- 180 deg is possible.
    else {
      RbtPrincipalAxes prAxes = Rbt::GetPrincipalAxes(atomList);
      for (RbtAtomListIter iter = atomList.begin(); iter != atomList.end(); ++iter) {
	m_tetheredAtomList.push_back(*iter);
	RbtDouble radius = Rbt::Length((*iter)->GetCoords(),prAxes.com);
	(*iter)->SetUser2Value(2.0 * radius);
      }
    }
  }
  else if (eTransMode == RbtChromElement::TETHERED) {
    RbtDouble maxTrans = pFlexData->GetParameter(RbtLigandFlexData::_MAX_TRANS);
    //Tethered translation; fixed rotation
    //Therefore max translation of any atom is a constant (+/- maxTrans)
    if (eRotMode == RbtChromElement::FIXED) {
      std::copy(atomList.begin(), atomList.end(), std::back_inserter(m_tetheredAtomList));
      pModel->SetAtomUser2Values(maxTrans);
    }
    //Tethered translation, tethered or free rotation
    //Therefore max translation of any atom depends on its distance from the COM
    //and on maxTrans
    //For simplicity, we don't distinguish between tethered and free rotation
    //so assume full rotation through +/- 180 deg is possible.
    else {
      RbtPrincipalAxes prAxes = Rbt::GetPrincipalAxes(atomList);
      for (RbtAtomListIter iter = atomList.begin(); iter != atomList.end(); ++iter) {
	m_tetheredAtomList.push_back(*iter);
	RbtDouble radius = Rbt::Length((*iter)->GetCoords(),prAxes.com);
	(*iter)->SetUser2Value( (2.0 * radius) + maxTrans);
      }
    }
  }
}

void RbtFlexAtomFactory::VisitSolventFlexData(RbtSolventFlexData* pFlexData) {
  VisitLigandFlexData(pFlexData);
}
