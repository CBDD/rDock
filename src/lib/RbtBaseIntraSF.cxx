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

#include "RbtBaseIntraSF.h"
#include "RbtDebug.h"
#include "RbtWorkSpace.h"

// Static data members
RbtString RbtBaseIntraSF::_CT("RbtBaseIntraSF");

RbtBaseIntraSF::RbtBaseIntraSF(): m_zero(0.0) {
    DEBUG(_CT << " default constructor" << endl);
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtBaseIntraSF::~RbtBaseIntraSF() {
    DEBUG(_CT << " destructor" << endl);
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

////////////////////////////////////////
// Public methods
////////////////

RbtModelPtr RbtBaseIntraSF::GetLigand() const { return m_spLigand; }

// Override RbtObserver pure virtual
// Notify observer that subject has changed
void RbtBaseIntraSF::Update(RbtSubject* theChangedSubject) {
    RbtWorkSpace* pWorkSpace = GetWorkSpace();
    if (theChangedSubject == pWorkSpace) {
        // Check if ligand has been updated (model #1)
        if (pWorkSpace->GetNumModels() >= 2) {
            RbtModelPtr spLigand = GetWorkSpace()->GetModel(1);
            if (spLigand != m_spLigand) {
                DEBUG(_CT << "::Update(): Ligand has been updated" << endl);
                m_spLigand = spLigand;
                SetupScore();
                // Retain the zero-point offset from the ligand model data if present
                // Otherwise set the zero-point to the current score
                RbtString name = GetFullName() + ".0";
                if (m_spLigand->isDataFieldPresent(name)) {
                    m_zero = m_spLigand->GetDataValue(name);
                } else {
                    m_zero = RawScore();
                }
            }
        }
    }
}

// Override RbtBaseSF::ScoreMap to provide additional raw descriptors
void RbtBaseIntraSF::ScoreMap(RbtStringVariantMap& scoreMap) const {
    if (isEnabled()) {
        RbtDouble rs = RawScore() - m_zero;  // report the raw score relative to the zero-point offset
        RbtString name = GetFullName();
        scoreMap[name] = rs;
        AddToParentMapEntry(scoreMap, rs);
        scoreMap[name + ".0"] = m_zero;
    }
}
