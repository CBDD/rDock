/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtConstSF.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtConstSF.h"

//Static data member for class type
RbtString RbtConstSF::_CT("RbtConstSF");
RbtString RbtConstSF::_SOLVENT_PENALTY("SOLVENT_PENALTY");

RbtConstSF::RbtConstSF(const RbtString& strName) 
  : RbtBaseSF(_CT,strName),m_solventPenalty(0.5) {
  AddParameter(_SOLVENT_PENALTY,m_solventPenalty);
#ifdef _DEBUG
	cout << _CT << " parameterised constructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtConstSF::~RbtConstSF() {
#ifdef _DEBUG
	cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtConstSF::ScoreMap(RbtStringVariantMap& scoreMap) const {
  if (isEnabled()) {
    //Divide the total raw score into "system" and "inter" components.
    RbtDouble rs = InterScore();

    //First deal with the inter score which is stored in its natural location in the map
    RbtString name = GetFullName();
    scoreMap[name] = rs;
    AddToParentMapEntry(scoreMap, rs);

    //Now deal with the system raw score which needs to be stored in SCORE.SYSTEM.CONST
    RbtDouble system_rs = SystemScore();
    RbtString systemName = RbtBaseSF::_SYSTEM_SF + "." + GetName();
    scoreMap[systemName] = system_rs;
    //increment the SCORE.SYSTEM total
    RbtDouble parentScore = scoreMap[RbtBaseSF::_SYSTEM_SF];
    parentScore += system_rs * GetWeight();
    scoreMap[RbtBaseSF::_SYSTEM_SF] = parentScore;
  }
}

RbtDouble RbtConstSF::RawScore() const {
  return InterScore() + SystemScore();
}

RbtDouble RbtConstSF::InterScore() const {
  //With current protocols, ligand should always be enabled, so RawScore is 1
  RbtDouble retVal = GetLigand()->GetEnabled() ? 1.0 : 0.0;
  return retVal;
}
 
RbtDouble RbtConstSF::SystemScore() const {
  RbtDouble retVal = 0.0;
  RbtModelList solventList = GetSolvent();
  for (RbtModelListConstIter iter = solventList.begin(); iter != solventList.end(); ++iter) {
    if ((*iter)->GetEnabled()) {
      retVal += m_solventPenalty;
    }
  }
  return retVal;
}

void RbtConstSF::ParameterUpdated(const RbtString& strName) {
  if (strName == _SOLVENT_PENALTY) {
    m_solventPenalty = GetParameter(_SOLVENT_PENALTY);
  }
  else {
    RbtBaseSF::ParameterUpdated(strName);
  }
}

