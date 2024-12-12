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

#include "RbtSetupSASF.h"

RbtString RbtSetupSASF::_CT("RbtSetupSASF");

RbtSetupSASF::RbtSetupSASF(const RbtString& strName): RbtBaseSF(_CT, strName) {
    DEBUG_ERR(_CT << "parameterized constructor      <--------" << endl);
    Disable();
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtSetupSASF::~RbtSetupSASF() {
    DEBUG_ERR(_CT << " destructor                    <--------" << endl);
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtSetupSASF::SetupReceptor() {
    DEBUG_ERR(_CT << " SetupReceptor                 <--------" << endl);
    // get receptor atoms
    // theReceptorList = Rbt::GetAtomList(GetReceptor()->GetAtomList(),std::not1(Rbt::isAtomicNo_eq(1)));
    theReceptorList = GetReceptor()->GetAtomList();
    DEBUG_ERR(_CT << "::SetupReceptor(): #ATOMS = " << theReceptorList.size() << endl);
}

void RbtSetupSASF::SetupScore() { DEBUG_ERR(_CT << " SetupScore                    <--------" << endl); }

void RbtSetupSASF::SetupLigand() {
    DEBUG_ERR(_CT << " SetupLigand                   <--------" << endl);
    theLigandList = GetLigand()->GetAtomList();
    DEBUG_ERR(_CT << "::SetupLigand(): #ATOMS = " << theLigandList.size() << endl);
}

RbtDouble RbtSetupSASF::RawScore() const { return 0.0; }
