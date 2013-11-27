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

#include "RbtChromPositionElement.h"

RbtString RbtChromPositionElement::_CT = "RbtChromPositionElement";

RbtChromPositionElement::RbtChromPositionElement(
                                const RbtModel* pModel,
                                const RbtDockingSite* pDockSite,
                                RbtDouble transStepSize,
                                RbtDouble rotStepSize,
                                RbtChromElement::eMode transMode,
                                RbtChromElement::eMode rotMode,
                                RbtDouble maxTrans,
                                RbtDouble maxRot)
{
    m_spRefData = new RbtChromPositionRefData(
                                pModel,
                                pDockSite,
                                transStepSize,
                                rotStepSize,
                                transMode,
                                rotMode,
                                maxTrans,
                                maxRot);
    SyncFromModel();
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtChromPositionElement::RbtChromPositionElement(
                                RbtChromPositionRefDataPtr spRefData,
                                const RbtCoord& com,
                                const RbtEuler& orientation)
        : m_spRefData(spRefData), m_com(com), m_orientation(orientation)
{
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtChromPositionElement::~RbtChromPositionElement() {
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtChromPositionElement::Reset() {
    m_com = m_spRefData->GetInitialCOM();
    m_orientation = m_spRefData->GetInitialOrientation();    
}

void RbtChromPositionElement::Randomise() {
    RandomiseCOM();
    RandomiseOrientation();
}

void RbtChromPositionElement::RandomiseCOM() {
    RbtDouble dist;
    RbtVector axis;
    RbtInt numStartCoords;
    RbtInt iCoord;
    switch (m_spRefData->GetTransMode()) {
        //TETHERED: Perform a single mutation from the initial COM
        //up to the maximum permitted
         case RbtChromElement::TETHERED:
            dist = m_spRefData->GetMaxTrans() * GetRand().GetRandom01();
            axis = GetRand().GetRandomUnitVector();
            m_com = m_spRefData->GetInitialCOM() + dist * axis;
            break;
        //FREE: Randomise across docking volume coords
        case RbtChromElement::FREE:
            numStartCoords = m_spRefData->GetNumStartCoords();
            iCoord = GetRand().GetRandomInt(numStartCoords);
            m_com = m_spRefData->GetStartCoord(iCoord);
            break;
        //FIXED: Revert to initial COM
        default:
            m_com = m_spRefData->GetInitialCOM();
            break;
    }
}

void RbtChromPositionElement::RandomiseOrientation() {
    RbtDouble theta;
    RbtVector axis;
    RbtDouble heading, attitude, bank;
    switch (m_spRefData->GetRotMode()) {
        //TETHERED: Perform a single mutation from the initial orientation
        //up to the maximum permitted
         case RbtChromElement::TETHERED:
            m_orientation = m_spRefData->GetInitialOrientation();
            theta = m_spRefData->GetMaxRot() * GetRand().GetRandom01();
            axis = GetRand().GetRandomUnitVector();
            m_orientation.Rotate(axis, theta);
            break;
        //FREE: completely scramble the initial orientation
        case RbtChromElement::FREE:
            heading = 2.0*M_PI*GetRand().GetRandom01() - M_PI;
            attitude =    M_PI*GetRand().GetRandom01() - 0.5*M_PI;
            bank =    2.0*M_PI*GetRand().GetRandom01() - M_PI;
            m_orientation = RbtEuler(heading, attitude, bank);
            break;
        //FIXED: Revert to initial orientation
        default:
            m_orientation = m_spRefData->GetInitialOrientation();
            break;
    }
}

void RbtChromPositionElement::Mutate(RbtDouble relStepSize) {
    MutateCOM(relStepSize);
    MutateOrientation(relStepSize);
}

void RbtChromPositionElement::MutateCOM(RbtDouble relStepSize) {
    RbtDouble absTransStepSize;
    RbtDouble dist;
    RbtVector axis;
    switch (m_spRefData->GetTransMode()) {
         case RbtChromElement::TETHERED:
            absTransStepSize = relStepSize * m_spRefData->GetTransStepSize();
            if (absTransStepSize > 0) {
               dist = absTransStepSize * GetRand().GetRandom01();
               axis = GetRand().GetRandomUnitVector();
               m_com += dist * axis;
               CorrectTetheredCOM();
            }
            break;
         case RbtChromElement::FREE:
            absTransStepSize = relStepSize * m_spRefData->GetTransStepSize();
            if (absTransStepSize > 0) {
               dist = absTransStepSize * GetRand().GetRandom01();
               axis = GetRand().GetRandomUnitVector();
               m_com += dist * axis;
            }
            break;
         //FIXED: Do nothing
         default:
            break;
    }
}

void RbtChromPositionElement::MutateOrientation(RbtDouble relStepSize) {
    RbtDouble absRotStepSize;
    RbtDouble theta;
    RbtVector axis;
    switch (m_spRefData->GetRotMode()) {
         case RbtChromElement::TETHERED:
            absRotStepSize = relStepSize * m_spRefData->GetRotStepSize();
            if (absRotStepSize > 0) {
                theta = absRotStepSize * GetRand().GetRandom01();
                axis = GetRand().GetRandomUnitVector();
                m_orientation.Rotate(axis, theta);
                CorrectTetheredOrientation();
             }
            break;
         case RbtChromElement::FREE:
            absRotStepSize = relStepSize * m_spRefData->GetRotStepSize();
            if (absRotStepSize > 0) {
                theta = absRotStepSize * GetRand().GetRandom01();
                axis = GetRand().GetRandomUnitVector();
                m_orientation.Rotate(axis, theta);
            }
            break;
         //FIXED: Do nothing
         default:
            break;
    }
}

void RbtChromPositionElement::SyncFromModel() {
    m_spRefData->GetModelValue(m_com, m_orientation);
}

void RbtChromPositionElement::SyncToModel() {
    m_spRefData->SetModelValue(m_com, m_orientation);
}

RbtChromElement* RbtChromPositionElement::clone() const {
    return new RbtChromPositionElement(m_spRefData, m_com, m_orientation);
}

void RbtChromPositionElement::GetVector(RbtDoubleList& v) const {
    if (!m_spRefData->IsTransFixed()) {
        v.push_back(m_com.x);
        v.push_back(m_com.y);
        v.push_back(m_com.z);
    }
    if (!m_spRefData->IsRotFixed()) {
        v.push_back(m_orientation.GetHeading());
        v.push_back(m_orientation.GetAttitude());
        v.push_back(m_orientation.GetBank());
    }
}

void RbtChromPositionElement::GetVector(RbtXOverList& v) const {
    if (!m_spRefData->IsTransFixed()) {
        RbtXOverElement comElement;
        comElement.push_back(m_com.x);
        comElement.push_back(m_com.y);
        comElement.push_back(m_com.z);
        v.push_back(comElement);
    }
    if (!m_spRefData->IsRotFixed()) {
        RbtXOverElement orientationElement;
        orientationElement.push_back(m_orientation.GetHeading());
        orientationElement.push_back(m_orientation.GetAttitude());
        orientationElement.push_back(m_orientation.GetBank());
        v.push_back(orientationElement);
    }
}

void RbtChromPositionElement::SetVector(const RbtDoubleList& v, RbtInt& i) throw (RbtError) {
    if (VectorOK(v,i)) {
        if (!m_spRefData->IsTransFixed()) {
            //2013Nov26 (DM) Bug fix to unsafe code.
            //We cannot assume that the multiple increments will be processed left->right.
            //This assumption is broken in g++ 4. Safer to break into separate statements.
            //m_com = RbtCoord(v[i++], v[i++], v[i++]);
            RbtDouble x(v[i++]);
            RbtDouble y(v[i++]);
            RbtDouble z(v[i++]);
            m_com = RbtCoord(x, y, z);
        }
        if (!m_spRefData->IsRotFixed()) {
            //2013Nov26 (DM) Bug fix to unsafe code.
            //We cannot assume that the multiple increments will be processed left->right.
            //This assumption is broken in g++ 4. Safer to break into separate statements.
            //m_orientation = RbtEuler(v[i++], v[i++], v[i++]);
            RbtDouble heading(v[i++]);
            RbtDouble attitude(v[i++]);
            RbtDouble bank(v[i++]);
            m_orientation = RbtEuler(heading, attitude, bank);
            m_orientation.Standardise();
        }
    }
    else {
        throw RbtBadArgument(_WHERE_,"Index out of range or insufficient elements remaining");
    }
 }
 
 void RbtChromPositionElement::SetVector(const RbtXOverList& v, RbtInt& i) throw (RbtError) {
    if (VectorOK(v,i)) {
        if (!m_spRefData->IsTransFixed()) {
            RbtXOverElement comElement(v[i++]);
            if (comElement.size() == 3) {
                //As we crossover an intact COM vector there should be no need to check
                //for tethered bounds
                m_com = RbtCoord(comElement[0], comElement[1], comElement[2]);
            }
            else {
                throw RbtBadArgument(_WHERE_, "comElement vector is of incorrect length");
            }
        }
        if (!m_spRefData->IsRotFixed()) {
            RbtXOverElement orientationElement(v[i++]);
            if (orientationElement.size() == 3) {
                //As we crossover an intact orientation vector there should be no need to check
                //for tethered bounds
                m_orientation = RbtEuler(orientationElement[0],
                                            orientationElement[1],
                                            orientationElement[2]);
            }
            else {
                throw RbtBadArgument(_WHERE_, "orientationElement vector is of incorrect length");
            }
        }
    }
    else {
        throw RbtBadArgument(_WHERE_, "Index out of range or insufficient elements remaining");
    }
}

void RbtChromPositionElement::GetStepVector(RbtDoubleList& v) const {
    if (!m_spRefData->IsTransFixed()) {
        RbtDouble transStepSize = m_spRefData->GetTransStepSize();
        for (RbtInt i = 0; i < 3; ++i) {
            v.push_back(transStepSize);
        }
    }
    if (!m_spRefData->IsRotFixed()) {
        RbtDouble rotStepSize = m_spRefData->GetRotStepSize();
        for (RbtInt i = 0; i < 3; ++i) {
           v.push_back(rotStepSize);
        }
    }
}

RbtDouble RbtChromPositionElement::CompareVector(const RbtDoubleList& v, RbtInt& i) const {
    RbtDouble retVal(0.0);
    if (VectorOK(v,i)) {
        if (!m_spRefData->IsTransFixed()) {
            //2013Nov26 (DM) Bug fix to unsafe code.
            //We cannot assume that the multiple increments will be processed left->right.
            //This assumption is broken in g++ 4. Safer to break into separate statements.
            //RbtCoord otherCom = RbtCoord(v[i++], v[i++], v[i++]);
            RbtDouble x(v[i++]);
            RbtDouble y(v[i++]);
            RbtDouble z(v[i++]);
            RbtCoord otherCom = RbtCoord(x, y, z);
            RbtDouble transStepSize = m_spRefData->GetTransStepSize();
            //Compare distance between centres of mass
            if (transStepSize > 0.0) {
                RbtDouble absDiff = Rbt::Length(m_com, otherCom);
                RbtDouble relDiff = absDiff / transStepSize;
                retVal = std::max(retVal, relDiff);
            }
        }
        if (!m_spRefData->IsRotFixed()) {
            //2013Nov26 (DM) Bug fix to unsafe code.
            //We cannot assume that the multiple increments will be processed left->right.
            //This assumption is broken in g++ 4. Safer to break into separate statements.
            //RbtEuler otherOrientation = RbtEuler(v[i++], v[i++], v[i++]);
            RbtDouble heading(v[i++]);
            RbtDouble attitude(v[i++]);
            RbtDouble bank(v[i++]);
            RbtEuler otherOrientation = RbtEuler(heading, attitude, bank);
            RbtDouble rotStepSize = m_spRefData->GetRotStepSize();
            //Compare orientations
            if (rotStepSize > 0.0) {
                //Determine the difference between the two orientations
                //in terms of the axis/angle needed to align them
                //q.s = cos(phi / 2)
                RbtQuat qAlign = otherOrientation.ToQuat()
                                    * m_orientation.ToQuat().Conj();
                RbtDouble cosHalfTheta = qAlign.s;
                if (cosHalfTheta < -1.0) {
                    cosHalfTheta = -1.0;
                }
                else if (cosHalfTheta > 1.0) {
                    cosHalfTheta = 1.0;
                }
                RbtDouble absDiff = fabs(StandardisedValue(2.0 * acos(cosHalfTheta)));
                RbtDouble relDiff = absDiff / rotStepSize;
                retVal = std::max(retVal, relDiff);
            }
        }      
    }
    else {
        retVal = -1.0;
    }
    return retVal;
}

void RbtChromPositionElement::Print(ostream& s) const {
    s << "COM " << m_com << endl;
    s << "EULER " << m_orientation << endl;
}

RbtDouble RbtChromPositionElement::StandardisedValue(RbtDouble rotationAngle) {
   while (rotationAngle >= M_PI) {
      rotationAngle -= 2.0 * M_PI;
   }
   while (rotationAngle < -M_PI) {
       rotationAngle += 2.0 * M_PI;
   }
   return rotationAngle;
}

void RbtChromPositionElement::CorrectTetheredCOM() {
   //If we are out of bounds
   //revert to the initial COM and translate as far as we can
   //towards the new point until we hit the boundary
   RbtDouble maxTrans = m_spRefData->GetMaxTrans();
   //vector from initial COM to new point
   RbtVector axis = m_com - m_spRefData->GetInitialCOM();
   //avoid square roots until necessary
   if (axis.Length2() > (maxTrans*maxTrans)) {
        //Stay just inside the boundary sphere
        m_com = m_spRefData->GetInitialCOM()
                + (0.999 * maxTrans * axis.Unit()); 
    }
}

void RbtChromPositionElement::CorrectTetheredOrientation() {
    //Check for orientation out of bounds
    RbtDouble maxRot = m_spRefData->GetMaxRot();
    RbtQuat qAlign = m_spRefData->GetInitialQuat() * m_orientation.ToQuat().Conj();
    RbtDouble cosHalfTheta = qAlign.s;
    if (cosHalfTheta < -1.0) {
        cosHalfTheta = -1.0;
    }
    else if (cosHalfTheta > 1.0) {
        cosHalfTheta = 1.0;
    }
    //Theta is the rotation angle required to realign with reference 
    RbtDouble theta = StandardisedValue(2.0 * acos(cosHalfTheta));
    //Deal with pos and neg angles independently as we have to
    //invert the rotation axis for negative angles
    if (theta < -maxRot) {
        RbtVector axis = -qAlign.v / sin(theta / 2.0);
        //Adjust theta to bring the orientation just inside the tethered bound
        theta += 0.999 * maxRot;
        m_orientation.Rotate(axis, theta);                    
    }
    else if (theta > maxRot) {
        RbtVector axis = qAlign.v / sin(theta / 2.0);
        //Adjust theta to bring the orientation just inside the tethered bound
        theta -= 0.999 * maxRot;
        m_orientation.Rotate(axis, theta);                    
    }
}


