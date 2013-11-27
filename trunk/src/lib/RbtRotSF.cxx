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

#include "RbtRotSF.h"

//Static data member for class type
RbtString RbtRotSF::_CT("RbtRotSF");
RbtString RbtRotSF::_INCNH3("INCNH3");
RbtString RbtRotSF::_INCOH("INCOH");

RbtRotSF::RbtRotSF(const RbtString& strName)
  : RbtBaseSF(_CT,strName),nRot(0),bIncNH3(false),bIncOH(false) {
  //Add parameters
  AddParameter(_INCNH3,bIncNH3);
  AddParameter(_INCOH,bIncOH);
#ifdef _DEBUG
  cout << _CT << " parameterised constructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtRotSF::~RbtRotSF() {
#ifdef _DEBUG
  cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtRotSF::SetupReceptor() {}

void RbtRotSF::SetupLigand() {
  if (GetLigand().Null()) {
    nRot = 0;
  }
  else {
    RbtBondList rotatableBondList = Rbt::GetBondList(GetLigand()->GetBondList(),Rbt::isBondRotatable());
    if (!bIncNH3) {
      rotatableBondList = Rbt::GetBondList(rotatableBondList,std::not1(Rbt::isBondToNH3()));
    }
    if (!bIncOH) {
      rotatableBondList = Rbt::GetBondList(rotatableBondList,std::not1(Rbt::isBondToOH()));
    }
    nRot = rotatableBondList.size();
  }
}

void RbtRotSF::SetupScore() {}

RbtDouble RbtRotSF::RawScore() const {
  return RbtDouble(nRot);
}

void RbtRotSF::ParameterUpdated(const RbtString& strName) {
  if (strName == _INCNH3) {
    bIncNH3 = GetParameter(_INCNH3);
  }
  else if (strName == _INCOH) {
    bIncOH = GetParameter(_INCOH);
  }
  else {
    RbtBaseSF::ParameterUpdated(strName);
  }
}
