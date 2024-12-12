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

#include "RbtBaseIdxSF.h"

#include "RbtDebug.h"
#include "RbtWorkSpace.h"

// Static data members
RbtString RbtBaseIdxSF::_CT("RbtBaseIdxSF");
RbtString RbtBaseIdxSF::_GRIDSTEP("GRIDSTEP");
RbtString RbtBaseIdxSF::_BORDER("BORDER");

RbtBaseIdxSF::RbtBaseIdxSF(): m_gridStep(0.5), m_border(1.0) {
    DEBUG_ERR(_CT << " default constructor" << endl);
    // Add parameters
    AddParameter(_GRIDSTEP, m_gridStep);
    AddParameter(_BORDER, m_border);
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtBaseIdxSF::~RbtBaseIdxSF() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

RbtDouble RbtBaseIdxSF::GetGridStep() const { return m_gridStep; }

void RbtBaseIdxSF::SetGridStep(RbtDouble step) { SetParameter(_GRIDSTEP, step); }

RbtDouble RbtBaseIdxSF::GetBorder() const { return m_border; }

void RbtBaseIdxSF::SetBorder(RbtDouble border) { SetParameter(_BORDER, border); }

RbtDouble RbtBaseIdxSF::GetMaxError() const {
    // maxError is half a grid diagonal. This is the tolerance we have to allow when indexing the receptor atoms on the
    // grid When retrieving the nearest neighbours to a ligand atom, we do the lookup on the nearest grid point to the
    // ligand atom. However the actual ligand atom - receptor atom distance may be maxError Angstroms closer than the
    // grid point - receptor atom distance used in the indexing.
    return 0.5 * sqrt(3.0) * m_gridStep;
}

// DM 12 Apr 2002
// Returns the maximum range of the scoring function,
// corrected for max grid error, and grid border around docking site
// This should be used by subclasses for selecting the receptor atoms to index
// GetCorrectedRange() = GetRange() + GetMaxError() + GetBorder()
RbtDouble RbtBaseIdxSF::GetCorrectedRange() const { return GetRange() + GetMaxError() + m_border; }

// As this has a virtual base class we need a separate OwnParameterUpdated
// which can be called by concrete subclass ParameterUpdated methods
// See Stroustrup C++ 3rd edition, p395, on programming virtual base classes
void RbtBaseIdxSF::OwnParameterUpdated(const RbtString& strName) {
    if (strName == _GRIDSTEP) {
        m_gridStep = GetParameter(_GRIDSTEP);
    } else if (strName == _BORDER) {
        m_border = GetParameter(_BORDER);
    }
}
