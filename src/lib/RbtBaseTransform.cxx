/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtBaseTransform.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtBaseTransform.h"
#include "RbtBaseSF.h"
#include "RbtWorkSpace.h"

//Static data members
RbtString RbtBaseTransform::_CT("RbtBaseTransform");

////////////////////////////////////////
//Constructors/destructors
RbtBaseTransform::RbtBaseTransform(const RbtString& strClass, const RbtString& strName) :
					RbtBaseObject(strClass,strName), m_parent(NULL) {
#ifdef _DEBUG
	cout << _CT << " parameterised constructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}


RbtBaseTransform::~RbtBaseTransform() {
	Orphan();//Remove object from parent aggregate
#ifdef _DEBUG
	cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}


////////////////////////////////////////
//Public methods
////////////////

//Fully qualified name, prefixed by all ancestors
RbtString RbtBaseTransform::GetFullName() const {
	//If we have a parent, prefix our short name with our parents full name,
	//else full name is just our short name
	return (m_parent) ? m_parent->GetFullName() + "." + GetName() : GetName();
}

//Actually apply the transform
void RbtBaseTransform::Go()
{
	if (isEnabled()) {
		//Send any stored Scoring Function requests (e.g. to change any params)
		SendSFRequests();
		Execute();
	}
}
		

//Aggregate handling (virtual) methods
//Base class throws an InvalidRequest error

void RbtBaseTransform::Add(RbtBaseTransform*) throw (RbtError) {
	throw RbtInvalidRequest(_WHERE_,"Add() invalid for non-aggregate transforms");
}
void RbtBaseTransform::Remove(RbtBaseTransform*) throw (RbtError) {
	throw RbtInvalidRequest(_WHERE_,"Remove() invalid for non-aggregate transforms");
}

RbtBaseTransform* RbtBaseTransform::GetTransform(RbtUInt i) const throw (RbtError) {
	throw RbtInvalidRequest(_WHERE_,"GetTransform() invalid for non-aggregate transforms");
}
RbtBool RbtBaseTransform::isAgg() const {return false;}
RbtUInt RbtBaseTransform::GetNumTransforms() const {return 0;}

//Aggregate handling (concrete) methods
RbtBaseTransform* RbtBaseTransform::GetParentTransform() const {return m_parent;}

//Force removal from the parent aggregate
void RbtBaseTransform::Orphan() {
	if (m_parent) {
#ifdef _DEBUG
		cout << "RbtBaseTransform::Orphan(): Removing " << GetName() << " from " << m_parent->GetName() << endl;
#endif //_DEBUG
		m_parent->Remove(this);
	}
}

//Scoring function request handling
//Transforms can store up a list of requests to send to the workspace
//scoring function each time Go() is executed
void RbtBaseTransform::AddSFRequest(RbtRequestPtr spRequest) {
	m_SFRequests.push_back(spRequest);
}
void RbtBaseTransform::ClearSFRequests() {
	m_SFRequests.clear();
}
void RbtBaseTransform::SendSFRequests() {
	//Get the current scoring function from the workspace
	RbtWorkSpace* pWorkSpace = GetWorkSpace();
	if (pWorkSpace == NULL) //Return if this transform is not registered
		return;
	RbtBaseSF* pSF = pWorkSpace->GetSF();
	if (pSF == NULL) //Return if workspace does not have a scoring function
		return;
	//Send each request to the scoring function
	for (RbtRequestListConstIter iter = m_SFRequests.begin(); iter != m_SFRequests.end(); iter++) {
		pSF->HandleRequest(*iter);
	}
}
