/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtDihedralIntraSF.cxx#2 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtDihedralIntraSF.h"

//Static data members
RbtString RbtDihedralIntraSF::_CT("RbtDihedralIntraSF");

RbtDihedralIntraSF::RbtDihedralIntraSF(const RbtString& strName) : RbtBaseSF(_CT,strName)
{
#ifdef _DEBUG
  cout << _CT << " parameterised constructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtDihedralIntraSF::~RbtDihedralIntraSF() {
  ClearModel();
#ifdef _DEBUG
  cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtDihedralIntraSF::SetupScore() {
  ClearModel();
  if (GetLigand().Null())
    return;
  if (GetLigand()->isFlexible()) {
    m_dihList = CreateDihedralList(GetLigand()->GetFlexBonds());
  }
}

RbtDouble RbtDihedralIntraSF::RawScore() const {
  RbtDouble score = 0.0;//Total score
  for (RbtDihedralListConstIter iter = m_dihList.begin(); iter != m_dihList.end(); iter++) {
    score += (**iter)();
  }
  return score;
}

//Clear the dihedral list
//As we are not using smart pointers, there is some memory management to do
void RbtDihedralIntraSF::ClearModel() {
  for (RbtDihedralListIter iter = m_dihList.begin(); iter != m_dihList.end(); iter++) {
    delete *iter;
  }
  m_dihList.clear();
}
