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

#include "RbtBaseUniMolTransform.h"
#include "RbtWorkSpace.h"

//Static data members
RbtString RbtBaseUniMolTransform::_CT("RbtBaseUniMolTransform");

////////////////////////////////////////
//Constructors/destructors
RbtBaseUniMolTransform::RbtBaseUniMolTransform(const RbtString& strClass, const RbtString& strName) :
		RbtBaseTransform(strClass,strName)
{
#ifdef _DEBUG
	cout << _CT << " parameterised constructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtBaseUniMolTransform::~RbtBaseUniMolTransform()
{
#ifdef _DEBUG
	cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

////////////////////////////////////////
//Public methods
////////////////

RbtModelPtr RbtBaseUniMolTransform::GetLigand() const {return m_spLigand;}
	
//Override RbtObserver pure virtual
//Notify observer that subject has changed
void RbtBaseUniMolTransform::Update(RbtSubject* theChangedSubject) {
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
	SetupTransform();
      }
    }
  }
}
