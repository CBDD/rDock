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

#include "RbtSubject.h"

//Constructor -does nothing.
RbtSubject::RbtSubject() {
  _RBTOBJECTCOUNTER_CONSTR_("RbtSubject");
}

//Destructor - notify observers of impending destruction
RbtSubject::~RbtSubject() {
#ifdef _DEBUG
	cout << "RbtSubject::~RbtSubject: Notifying observers of impending destruction" << endl;
#endif //_DEBUG
	//We need to iterate using a while loop because call to Deleted will trigger a call back to Detach,
	//reducing the size of m_observers, hence invalidating a conventional iterator
	while (m_observers.size() > 0) {
		RbtObserver* pObserver = m_observers.back();
		pObserver->Deleted(this);
	}
  _RBTOBJECTCOUNTER_DESTR_("RbtSubject");
}

////////////////////////////////////////
//Public methods
////////////////
void RbtSubject::Attach(RbtObserver* pObserver) throw (RbtError) {
	RbtObserverListIter iter = std::find(m_observers.begin(),m_observers.end(),pObserver);
	if (iter != m_observers.end()) {
		throw RbtBadArgument(_WHERE_,"RbtSubject::Attach(): pObserver is already attached to this subject");
	}
	else {
		m_observers.push_back(pObserver);
#ifdef _DEBUG
		cout << "RbtSubject::Attach: Attaching new observer; #observers=" << m_observers.size() << endl;
#endif //_DEBUG
	}
}

void RbtSubject::Detach(RbtObserver* pObserver) throw (RbtError) {
	RbtObserverListIter iter = std::find(m_observers.begin(),m_observers.end(),pObserver);
	if (iter == m_observers.end()) {
		throw RbtBadArgument(_WHERE_,"RbtSubject::Detach(): pObserver not attached to this subject");
	}
	else {
		m_observers.erase(iter);
#ifdef _DEBUG
		cout << "RbtSubject::Detach: Detaching observer; #observers=" << m_observers.size() << endl;
#endif //_DEBUG
	}
}

void RbtSubject::Notify() {
#ifdef _DEBUG
	cout << "RbtSubject::Notify: Notifying observers of change of state" << endl;
#endif //_DEBUG
	for (RbtObserverListIter iter = m_observers.begin(); iter != m_observers.end(); iter++) {
		(*iter)->Update(this);
	}
}

