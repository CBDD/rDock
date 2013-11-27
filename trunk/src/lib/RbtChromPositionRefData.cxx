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

#include "RbtChromPositionRefData.h"

RbtString RbtChromPositionRefData::_CT = "RbtChromPositionRefData";
const RbtPrincipalAxes RbtChromPositionRefData::CARTESIAN_AXES;

RbtChromPositionRefData::RbtChromPositionRefData(const RbtModel* pModel,
                                const RbtDockingSite* pDockSite,
                                RbtDouble transStepSize,
                                RbtDouble rotStepSize,
                                RbtChromElement::eMode transMode,
                                RbtChromElement::eMode rotMode,
                                RbtDouble maxTrans,
                                RbtDouble maxRot)
        : m_transStepSize(transStepSize),
          m_rotStepSize(rotStepSize),
          m_transMode(transMode),
          m_rotMode(rotMode),
          m_length(6),
          m_xOverLength(2),
          m_maxTrans(maxTrans),
          m_maxRot(maxRot)
{
    RbtAtomList atomList = pModel->GetAtomList();
    //Tethered substructure atom list (may be empty)
    RbtAtomList tetheredAtomList = pModel->GetTetheredAtomList();
    //If we have a tethered substructure, use this as the reference atom list
    //to calculate centre of mass and orientation, else use all atoms
    m_refAtoms = (tetheredAtomList.empty()) ? atomList : tetheredAtomList;
    //All atoms are movable, but use std::copy to strip off the smart pointers
    std::copy(atomList.begin(), atomList.end(), std::back_inserter(m_movableAtoms));
    GetModelValue(m_initialCom, m_initialOrientation);
    m_initialQuat = m_initialOrientation.ToQuat();
    pDockSite->GetCoordList(m_startCoords);
    //Check for zero ranges in TETHERED mode and convert to FIXED
    if ( (m_transMode == RbtChromElement::TETHERED) && (m_maxTrans <= 0.0) ) {
        m_transMode = RbtChromElement::FIXED;
        m_maxTrans = 0.0;
    }
    if ( (m_rotMode == RbtChromElement::TETHERED) && (m_maxRot <= 0.0) ) {
        m_rotMode = RbtChromElement::FIXED;
        m_maxRot = 0.0;
    }
    //Remove degrees of freedom for fixed modes
    if (IsTransFixed()) {
        m_length -= 3;
        m_xOverLength--;
    }
    if (IsRotFixed()) {
        m_length -=3;
        m_xOverLength--;
    }
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtChromPositionRefData::~RbtChromPositionRefData() {
    _RBTOBJECTCOUNTER_DESTR_(_CT);    
}

void RbtChromPositionRefData::GetModelValue(RbtCoord& com,
                                            RbtEuler& orientation) const {
    //Determine the principal axes and centre of mass of the reference atoms
    RbtPrincipalAxes prAxes = Rbt::GetPrincipalAxes(m_refAtoms);
    //Determine the quaternion needed to align Cartesian axes with actual
    //molecule principal axes. This represents the absolute orientation of
    //the molecule.
    RbtQuat q = Rbt::GetQuatFromAlignAxes(CARTESIAN_AXES, prAxes);
    orientation.FromQuat(q);
    com = prAxes.com;
}

void RbtChromPositionRefData::SetModelValue(const RbtCoord& com,
                                            const RbtEuler& orientation) {
    //Determine the principal axes and centre of mass of the reference atoms
    RbtPrincipalAxes prAxes = Rbt::GetPrincipalAxes(m_refAtoms);
    //Determine the overall rotation required.
    //1) Go back to realign with Cartesian axes
    RbtQuat qBack = Rbt::GetQuatFromAlignAxes(prAxes, CARTESIAN_AXES);
    //2) Go forward to the desired orientation
    RbtQuat qForward = orientation.ToQuat();
    //3 Combine the two rotations
    RbtQuat q = qForward * qBack;
    for (RbtAtomRListIter iter = m_movableAtoms.begin();
                                 iter != m_movableAtoms.end();
                                 ++iter) {
        (*iter)->Translate(-prAxes.com);//Move to origin
        (*iter)->RotateUsingQuat(q);//Rotate
        (*iter)->Translate(com);//Move to new centre of mass
    }
}
