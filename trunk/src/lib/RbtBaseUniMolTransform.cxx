/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtBaseUniMolTransform.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

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
