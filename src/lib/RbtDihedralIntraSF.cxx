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

#include "RbtDihedralIntraSF.h"

#include "RbtDebug.h"

// Static data members
RbtString RbtDihedralIntraSF::_CT("RbtDihedralIntraSF");

RbtDihedralIntraSF::RbtDihedralIntraSF(const RbtString& strName): RbtBaseSF(_CT, strName) {
    DEBUG_ERR(_CT << " parameterised constructor" << endl);
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtDihedralIntraSF::~RbtDihedralIntraSF() {
    ClearModel();
    DEBUG_ERR(_CT << " destructor" << endl);
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtDihedralIntraSF::SetupScore() {
    ClearModel();
    if (GetLigand().Null()) return;
    if (GetLigand()->isFlexible()) {
        m_dihList = CreateDihedralList(GetLigand()->GetFlexBonds());
    }
}

RbtDouble RbtDihedralIntraSF::RawScore() const {
    RbtDouble score = 0.0;  // Total score
    for (RbtDihedralListConstIter iter = m_dihList.begin(); iter != m_dihList.end(); iter++) {
        score += (**iter)();
    }
    return score;
}

// Clear the dihedral list
// As we are not using smart pointers, there is some memory management to do
void RbtDihedralIntraSF::ClearModel() {
    for (RbtDihedralListIter iter = m_dihList.begin(); iter != m_dihList.end(); iter++) {
        delete *iter;
    }
    m_dihList.clear();
}
