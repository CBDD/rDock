/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtSiteMapper.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
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
