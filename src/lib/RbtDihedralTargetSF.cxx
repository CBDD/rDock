/*This file is part of Rdock.

    Rdock is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Rdock is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Rdock.  If not, see <http://www.gnu.org/licenses/>.*/

#include "RbtDihedralTargetSF.h"

//Static data members
RbtString RbtDihedralTargetSF::_CT("RbtDihedralTargetSF");

RbtDihedralTargetSF::RbtDihedralTargetSF(const RbtString& strName) : RbtBaseSF(_CT,strName)
{
#ifdef _DEBUG
  cout << _CT << " parameterised constructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtDihedralTargetSF::~RbtDihedralTargetSF() {
  ClearReceptor();
#ifdef _DEBUG
  cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtDihedralTargetSF::SetupReceptor() {
  ClearReceptor();
  if (GetReceptor().Null())
    return;
  if (GetReceptor()->isFlexible()) {
    m_dihList = CreateDihedralList(GetReceptor()->GetFlexBonds());
  }
}
void RbtDihedralTargetSF::SetupLigand() {}
void RbtDihedralTargetSF::SetupScore() {}

RbtDouble RbtDihedralTargetSF::RawScore() const {
  RbtDouble score = 0.0;//Total score
  for (RbtDihedralListConstIter iter = m_dihList.begin(); iter != m_dihList.end(); iter++) {
    score += (**iter)();
  }
  return score;
}

//Clear the dihedral list
//As we are not using smart pointers, there is some memory management to do
void RbtDihedralTargetSF::ClearReceptor() {
  for (RbtDihedralListIter iter = m_dihList.begin(); iter != m_dihList.end(); iter++) {
    delete *iter;
  }
  m_dihList.clear();
}
