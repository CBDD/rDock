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

#include "RbtBaseIntraSF.h"
#include "RbtWorkSpace.h"

//Static data members
RbtString RbtBaseIntraSF::_CT("RbtBaseIntraSF");

RbtBaseIntraSF::RbtBaseIntraSF() : m_zero(0.0)
{
#ifdef _DEBUG
  cout << _CT << " default constructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtBaseIntraSF::~RbtBaseIntraSF()
{
#ifdef _DEBUG
  cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

////////////////////////////////////////
//Public methods
////////////////

RbtModelPtr RbtBaseIntraSF::GetLigand() const {return m_spLigand;}
	
//Override RbtObserver pure virtual
//Notify observer that subject has changed
void RbtBaseIntraSF::Update(RbtSubject* theChangedSubject) {
  RbtWorkSpace* pWorkSpace = GetWorkSpace();
  if (theChangedSubject == pWorkSpace) {
    //Check if ligand has been updated (model #1)
    if (pWorkSpace->GetNumModels() >= 2) {		
      RbtModelPtr spLigand = GetWorkSpace()->GetModel(1);
      if (spLigand != m_spLigand) {
#ifdef _DEBUG
	cout << _CT << "::Update(): Ligand has been updated" << endl;
#endif //_DEBUG
	m_spLigand = spLigand;
	SetupScore();
	//Retain the zero-point offset from the ligand model data if present
	//Otherwise set the zero-point to the current score
	RbtString name = GetFullName() + ".0";
	if (m_spLigand->isDataFieldPresent(name)) {
	  m_zero = m_spLigand->GetDataValue(name);
	}
	else {
	  m_zero = RawScore();
	}
      }
    }
  }
}

//Override RbtBaseSF::ScoreMap to provide additional raw descriptors
void RbtBaseIntraSF::ScoreMap(RbtStringVariantMap& scoreMap) const {
  if (isEnabled()) {
    RbtDouble rs = RawScore() - m_zero;//report the raw score relative to the zero-point offset
    RbtString name = GetFullName();
    scoreMap[name] = rs;
    AddToParentMapEntry(scoreMap, rs);
    scoreMap[name+".0"] = m_zero;
  }
}
