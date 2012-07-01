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
