/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtRotSF.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
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
