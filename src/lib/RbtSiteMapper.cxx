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
