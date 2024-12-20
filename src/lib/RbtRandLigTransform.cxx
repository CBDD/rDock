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

#include "RbtRandLigTransform.h"

// Static data member for class type
RbtString RbtRandLigTransform::_CT("RbtRandLigTransform");
// Parameter names
RbtString RbtRandLigTransform::_TORS_STEP("TORS_STEP");

const RbtRandLigTransform::Config
    RbtRandLigTransform::DEFAULT_CONFIG{};  // Empty initializer to fall back to default values

RbtRandLigTransform::RbtRandLigTransform(const RbtString& strName, const Config& config):
    RbtBaseUniMolTransform(_CT, strName),
    m_rand(Rbt::GetRbtRand()),
    config{config} {
    DEBUG_ERR(_CT << " parameterised constructor" << endl);
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtRandLigTransform::~RbtRandLigTransform() {
    DEBUG_ERR(_CT << " destructor" << endl);
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

////////////////////////////////////////
// Protected methods
///////////////////
void RbtRandLigTransform::SetupTransform() {
    // Clear the rotable bond list from the previous model
    m_rotableBonds.clear();
    // Check for undefined ligand
    if (GetLigand().Null()) return;
    m_rotableBonds = Rbt::GetRotatableBondList(GetLigand()->GetBondList());
}

////////////////////////////////////////
// Private methods
///////////////////
// Pure virtual in RbtBaseTransform
// Actually apply the transform
void RbtRandLigTransform::Execute() {
    RbtModelPtr spLigand = GetLigand();
    if (spLigand.Null()) return;
    for (auto& rotable_bond: m_rotableBonds) {
        RbtDouble thetaDeg = 2.0 * config.torsion_step * m_rand.GetRandom01() - config.torsion_step;
        spLigand->RotateBond(rotable_bond, thetaDeg, false);
    }
}
