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

#include "RbtSiteMapper.h"
#include "RbtWorkSpace.h"

//Static data members
RbtString RbtSiteMapper::_CT("RbtSiteMapper");

////////////////////////////////////////
//Constructors/destructors
RbtSiteMapper::RbtSiteMapper(const RbtString& strClass, const RbtString& strName) : RbtBaseObject(strClass,strName) {
#ifdef _DEBUG
	cout << _CT << " parameterised constructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtSiteMapper::~RbtSiteMapper() {
#ifdef _DEBUG
	cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

//Override RbtObserver pure virtual
//Notify observer that subject has changed
void RbtSiteMapper::Update(RbtSubject* theChangedSubject) {
  RbtWorkSpace* pWorkSpace = GetWorkSpace();
  RbtInt iTrace = GetTrace();
  if (theChangedSubject == pWorkSpace) {
    //Check if receptor has been updated (model #0)
    if (pWorkSpace->GetNumModels() >= 1) {
      RbtModelPtr spReceptor = GetWorkSpace()->GetModel(0);
      if (spReceptor != m_spReceptor) {
	m_spReceptor = spReceptor;
	if (iTrace > 1) {
	  cout << _CT << "::Update(): Receptor has been updated" << endl;
	}
      }
    }
  }
}
