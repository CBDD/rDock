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

#include "RbtBaseBiMolTransform.h"

#include "RbtDebug.h"
#include "RbtWorkSpace.h"

// Static data members
RbtString RbtBaseBiMolTransform::_CT("RbtBaseBiMolTransform");

////////////////////////////////////////
// Constructors/destructors
RbtBaseBiMolTransform::RbtBaseBiMolTransform(const RbtString& strClass, const RbtString& strName):
    RbtBaseTransform(strClass, strName) {
    DEBUG_ERR(_CT << " parameterised constructor" << endl);
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtBaseBiMolTransform::~RbtBaseBiMolTransform() {
    DEBUG_ERR(_CT << " destructor" << endl);
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

////////////////////////////////////////
// Public methods
////////////////

RbtModelPtr RbtBaseBiMolTransform::GetReceptor() const { return m_spReceptor; }
RbtModelPtr RbtBaseBiMolTransform::GetLigand() const { return m_spLigand; }
RbtModelList RbtBaseBiMolTransform::GetSolvent() const { return m_solventList; }

// Override RbtObserver pure virtual
// Notify observer that subject has changed
void RbtBaseBiMolTransform::Update(RbtSubject* theChangedSubject) {
    RbtWorkSpace* pWorkSpace = GetWorkSpace();
    if (theChangedSubject == pWorkSpace) {
        RbtInt numModels = pWorkSpace->GetNumModels();
        // Check if receptor has been updated (model #0)
        if (numModels >= 1) {
            RbtModelPtr spReceptor = pWorkSpace->GetModel(0);
            if (spReceptor != m_spReceptor) {
                DEBUG_ERR(_CT << "::Update(): Receptor has been updated" << endl);
                m_spReceptor = spReceptor;
                SetupReceptor();
            }
        }
        // Check if ligand has been updated (model #1)
        if (numModels >= 2) {
            RbtModelPtr spLigand = pWorkSpace->GetModel(1);
            if (spLigand != m_spLigand) {
                DEBUG_ERR(_CT << "::Update(): Ligand has been updated" << endl);
                m_spLigand = spLigand;
                SetupLigand();
            }
        }
        // check if solvent has been updated (model #2 onwards)
        // At the moment, we only check whether the number of solvent
        // molecules has changed
        RbtModelList solventList;
        if (numModels >= 3) {
            solventList = pWorkSpace->GetModels(2);
        }
        if (solventList.size() != m_solventList.size()) {
            m_solventList = solventList;
            SetupSolvent();
        }
        SetupTransform();
    }
}
