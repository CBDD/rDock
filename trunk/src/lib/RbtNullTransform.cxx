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

#include "RbtNullTransform.h"
#include "RbtPopulation.h"
#include "RbtWorkSpace.h"

//Static data member for class type
RbtString RbtNullTransform::_CT("RbtNullTransform");

////////////////////////////////////////
//Constructors/destructors
RbtNullTransform::RbtNullTransform(const RbtString& strName) :
										RbtBaseTransform(_CT,strName) {
#ifdef _DEBUG
	cout << _CT << " parameterised constructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtNullTransform::~RbtNullTransform() {
#ifdef _DEBUG
	cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}


////////////////////////////////////////
//Public methods
////////////////
void RbtNullTransform::Update(RbtSubject* theChangedSubject) {}

////////////////////////////////////////
//Private methods
////////////////
void RbtNullTransform::Execute() {
  RbtInt iTrace = GetTrace();
  //Trace level 1 = just output the current conformation (presumed to be the best from the
  //previous transform) to the history file
  if (iTrace == 5) {
    GetWorkSpace()->SaveHistory(true);
    RbtDouble s = GetWorkSpace()->GetSF()->Score();
    cout << "SCORE = " << s << endl;
  }
  //Trace level 2 = output entire GA population to history file, unless pop does not exist,
  //in which case same as level 1
  else if (iTrace > 5) {
    RbtPopulationPtr pop = GetWorkSpace()->GetPopulation();
    if (pop.Ptr() == NULL) {
      GetWorkSpace()->SaveHistory(true);
    }
    else {
	   const RbtGenomeList& genList = pop->GetGenomeList();
	   for (RbtGenomeListConstIter gIter = genList.begin(); gIter != genList.end(); ++gIter) {
	       (*gIter)->GetChrom()->SyncToModel();
	       GetWorkSpace()->SaveHistory(true); 
	   }
    }
  }
}

