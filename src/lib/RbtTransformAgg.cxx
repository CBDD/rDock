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

#include "RbtTransformAgg.h"

//Static data member for class type
RbtString RbtTransformAgg::_CT("RbtTransformAgg");

////////////////////////////////////////
//Constructors/destructors
RbtTransformAgg::RbtTransformAgg(const RbtString& strName) : RbtBaseTransform(_CT,strName) {
#ifdef _DEBUG
	cout << _CT << " parameterised constructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtTransformAgg::~RbtTransformAgg() {
	//Delete all our children
#ifdef _DEBUG
	cout << _CT << "::~" << _CT << "(): Deleting child transforms of " << GetName() << endl;
#endif //_DEBUG
	//We need to iterate using a while loop because each deletion will reduce the size of m_sf,
	//hence conventional iterators would become invalid
	while (m_transforms.size() > 0) {
		RbtBaseTransform* pTransform = m_transforms.back();
		//Assertion: parent of child is this object
		Assert<RbtAssert>(!TRANSFORMAGG_CHECK||pTransform->m_parent==this);
#ifdef _DEBUG
		cout << "Deleting " << pTransform->GetName() << " from " << GetName() << endl;
#endif //_DEBUG
		delete pTransform;
	}
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}


////////////////////////////////////////
//Public methods
////////////////

//Aggregate handling methods
void RbtTransformAgg::Add(RbtBaseTransform* pTransform) throw (RbtError) {
	//By first orphaning the transform to be added,
	//we handle attempts to readd existing children automatically,
	pTransform->Orphan();
	pTransform->m_parent = this;
#ifdef _DEBUG
		cout << _CT << "::Add(): Adding " << pTransform->GetName() << " to " << GetName() << endl;
#endif //_DEBUG
	m_transforms.push_back(pTransform);
}

void RbtTransformAgg::Remove(RbtBaseTransform* pTransform) throw (RbtError) {
	RbtBaseTransformListIter iter = std::find(m_transforms.begin(),m_transforms.end(),pTransform);
	if (iter == m_transforms.end()) {
		throw RbtBadArgument(_WHERE_,"Remove(): pTransform is not a member of this aggregate");
	}
	else {
		//Assertion: parent of child is this object
		Assert<RbtAssert>(!TRANSFORMAGG_CHECK||pTransform->m_parent==this);
#ifdef _DEBUG
		cout << _CT << "::Remove(): Removing " << pTransform->GetName() << " from " << GetName() << endl;
#endif //_DEBUG
		m_transforms.erase(iter);
		pTransform->m_parent = NULL;//Nullify the parent pointer of the child that has been removed 
	}
}

RbtBool RbtTransformAgg::isAgg() const {return true;}
RbtUInt RbtTransformAgg::GetNumTransforms() const {return m_transforms.size();}
RbtBaseTransform* RbtTransformAgg::GetTransform(RbtUInt iTransform) const throw (RbtError)
{
	if (iTransform >= m_transforms.size()) {
		throw RbtBadArgument(_WHERE_,"GetTransform(): iTransform out of range");
	}
	else {
		return m_transforms[iTransform];
	}
}

//WorkSpace handling methods
//Register scoring function with a workspace
//Aggregate version registers all children, but NOT itself
//(Aggregates are just containers, and have no need for model information
void RbtTransformAgg::Register(RbtWorkSpace* pWorkSpace) {
#ifdef _DEBUG
		cout << _CT << "::Register(): Registering child transforms of " << GetName() << endl;
#endif //_DEBUG
	for (RbtBaseTransformListIter iter = m_transforms.begin(); iter != m_transforms.end(); iter++) {
		(*iter)->Register(pWorkSpace);
	}
}

//Unregister with a workspace
//Aggregate version unregisters all children, but NOT itself
void RbtTransformAgg::Unregister() {
#ifdef _DEBUG
		cout << _CT << "::Unregister(): Unregistering child transforms of " << GetName() << endl;
#endif //_DEBUG
	for (RbtBaseTransformListIter iter = m_transforms.begin(); iter != m_transforms.end(); iter++) {
		(*iter)->Unregister();
	}
}

//Override RbtObserver pure virtual
//Notify observer that subject has changed
//Does nothing in RbtTransformAgg as aggregates do not require updating
void RbtTransformAgg::Update(RbtSubject* theChangedSubject) {
  //	if (theChangedSubject == GetWorkSpace()) {
  //#ifdef _DEBUG
  //		cout << _CT << "::Update(): " << GetName() << " received update from workspace" << endl;
  //#endif //_DEBUG
  //	}
}

//Request Handling method
//Aggregates handle the request themselves first, then cascade to all children
void RbtTransformAgg::HandleRequest(RbtRequestPtr spRequest) {
	RbtBaseObject::HandleRequest(spRequest);
	for (RbtBaseTransformListIter iter = m_transforms.begin(); iter != m_transforms.end(); iter++) {
		(*iter)->HandleRequest(spRequest);
	}
}

//Virtual function for dumping scoring function details to an output stream
//Called by operator <<
void RbtTransformAgg::Print(ostream& s) const {
	//First print the parameters for this aggregate
	RbtBaseObject::Print(s);
	//Now call Print() for each of the children
	for (RbtBaseTransformListConstIter iter = m_transforms.begin(); iter != m_transforms.end(); iter++) {
		(*iter)->Print(s);
	}
}

//Actually apply the transform
//Aggregate version loops over all child transforms
void RbtTransformAgg::Execute() {
	for (RbtBaseTransformListConstIter iter = m_transforms.begin(); iter != m_transforms.end(); iter++) {
		(*iter)->Go();
	}
}


