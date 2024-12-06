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

// Constructor -does nothing.
RbtSubject::RbtSubject() { _RBTOBJECTCOUNTER_CONSTR_("RbtSubject"); }

// Destructor - notify observers of impending destruction
RbtSubject::~RbtSubject() {
    DEBUG("RbtSubject::~RbtSubject: Notifying observers of impending destruction" << endl);
    // We need to iterate using a while loop because call to Deleted will trigger a call back to Detach,
    // reducing the size of m_observers, hence invalidating a conventional iterator
    while (m_observers.size() > 0) {
        RbtObserver* pObserver = m_observers.back();
        pObserver->Deleted(this);
    }
    _RBTOBJECTCOUNTER_DESTR_("RbtSubject");
}

////////////////////////////////////////
// Public methods
////////////////
void RbtSubject::Attach(RbtObserver* pObserver) {
    RbtObserverListIter iter = std::find(m_observers.begin(), m_observers.end(), pObserver);
    if (iter != m_observers.end()) {
        throw RbtBadArgument(_WHERE_, "RbtSubject::Attach(): pObserver is already attached to this subject");
    } else {
        m_observers.push_back(pObserver);
        DEBUG("RbtSubject::Attach: Attaching new observer; #observers=" << m_observers.size() << endl);
    }
}

void RbtSubject::Detach(RbtObserver* pObserver) {
    RbtObserverListIter iter = std::find(m_observers.begin(), m_observers.end(), pObserver);
    if (iter == m_observers.end()) {
        throw RbtBadArgument(_WHERE_, "RbtSubject::Detach(): pObserver not attached to this subject");
    } else {
        m_observers.erase(iter);
        DEBUG("RbtSubject::Detach: Detaching observer; #observers=" << m_observers.size() << endl);
    }
}

void RbtSubject::Notify() {
    DEBUG("RbtSubject::Notify: Notifying observers of change of state" << endl);
    for (RbtObserverListIter iter = m_observers.begin(); iter != m_observers.end(); iter++) {
        (*iter)->Update(this);
    }
}
