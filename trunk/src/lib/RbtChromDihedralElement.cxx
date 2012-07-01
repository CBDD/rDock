/*This file is part of Rdock.

    Rdock is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Rdock is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Rdock.  If not, see <http://www.gnu.org/licenses/>.*/

#include "RbtChromDihedralElement.h"

RbtString RbtChromDihedralElement::_CT = "RbtChromDihedralElement";

RbtChromDihedralElement::RbtChromDihedralElement(RbtBondPtr spBond,
                                RbtAtomList tetheredAtoms,
                                RbtDouble stepSize,
                                RbtChromElement::eMode mode,
                                RbtDouble maxDihedral)
        : m_value(0.0) {
	m_spRefData = new RbtChromDihedralRefData(spBond,
                                                tetheredAtoms,
                                                stepSize,
                                                mode,
                                                maxDihedral);
    //Set the initial genotype to match the current phenotype
	SyncFromModel();
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtChromDihedralElement::RbtChromDihedralElement(
        RbtChromDihedralRefDataPtr spRefData, RbtDouble value)
        : m_spRefData(spRefData), m_value(value) {
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtChromDihedralElement::~RbtChromDihedralElement() {
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtChromDihedralElement::Reset() {
    m_value = m_spRefData->GetInitialValue();    
}

void RbtChromDihedralElement::Randomise() {
    RbtDouble maxDelta = m_spRefData->GetMaxDihedral();
    RbtDouble delta;
    switch (m_spRefData->GetMode()) {
        case RbtChromElement::TETHERED:
            delta = 2.0 * maxDelta * GetRand().GetRandom01() - maxDelta;
            m_value = StandardisedValue(m_spRefData->GetInitialValue() + delta);
            break;
        case RbtChromElement::FREE:
            m_value = 360.0 * GetRand().GetRandom01() - 180.0;
            break;
        default:
            m_value = m_spRefData->GetInitialValue();
            break;
    }
}

void RbtChromDihedralElement::Mutate(RbtDouble relStepSize) {
	RbtDouble absStepSize = relStepSize * m_spRefData->GetStepSize();
    RbtDouble delta;
    if (absStepSize > 0) {
        switch (m_spRefData->GetMode()) {
            case RbtChromElement::TETHERED:
        	   delta = 2.0 * absStepSize * GetRand().GetRandom01() - absStepSize;
        	   m_value = StandardisedValue(m_value + delta);
               CorrectTetheredDihedral();
               break;
            case RbtChromElement::FREE:
               delta = 2.0 * absStepSize * GetRand().GetRandom01() - absStepSize;
               m_value = StandardisedValue(m_value + delta);
               break;
            default:
                break;
        }
    }
}

void RbtChromDihedralElement::SyncFromModel() {
	m_value = m_spRefData->GetModelValue();
}

void RbtChromDihedralElement::SyncToModel() {
	m_spRefData->SetModelValue(m_value);
}

RbtChromElement* RbtChromDihedralElement::clone() const {
	return new RbtChromDihedralElement(m_spRefData, m_value);
}

void RbtChromDihedralElement::GetVector(RbtDoubleList& v) const {
    v.push_back(m_value);
}

void RbtChromDihedralElement::GetVector(RbtXOverList& v) const {
    RbtXOverElement dihedralElement;
    dihedralElement.push_back(m_value);
    v.push_back(dihedralElement);
}

void RbtChromDihedralElement::SetVector(const RbtDoubleList& v, RbtInt& i) throw (RbtError) {
    if (VectorOK(v,i)) {
        m_value = StandardisedValue(v[i++]);
    }
    else {
        throw RbtBadArgument(_WHERE_, "Index out of range or insufficient elements remaining");
    }
}

void RbtChromDihedralElement::SetVector(const RbtXOverList& v, RbtInt& i) throw (RbtError) {
    if (VectorOK(v,i)) {
        RbtXOverElement dihedralElement(v[i++]);
        if (dihedralElement.size() == 1) {
            m_value = dihedralElement[0];
        }
        else {
            throw RbtBadArgument(_WHERE_, "dihedralElement vector is of incorrect length");
        }
    }
    else {
        throw RbtBadArgument(_WHERE_, "Index out of range or insufficient elements remaining");
    }
}

void RbtChromDihedralElement::GetStepVector(RbtDoubleList& v) const {
    v.push_back(m_spRefData->GetStepSize());
}

RbtDouble RbtChromDihedralElement::CompareVector(const RbtDoubleList& v, RbtInt& i) const {
    RbtDouble retVal(0.0);
    if (!VectorOK(v,i)) {
        retVal = -1.0;
    }
    else {
        RbtDouble otherAngle = v[i++];
        RbtDouble stepSize = m_spRefData->GetStepSize();
        if (stepSize > 0.0) {
            RbtDouble absDiff = fabs(StandardisedValue(m_value - otherAngle));
            retVal = absDiff / stepSize;
        }
    }
    return retVal;
}

void RbtChromDihedralElement::Print(ostream& s) const {
    s << "DIHEDRAL " << m_value << endl;
}

RbtDouble RbtChromDihedralElement::StandardisedValue(RbtDouble dihedralAngle) {
   while (dihedralAngle >= 180) {
      dihedralAngle -= 360.0;
   }
   while (dihedralAngle < -180.0) {
       dihedralAngle += 360.0;
   }
   return dihedralAngle;
}

void RbtChromDihedralElement::CorrectTetheredDihedral() {
    RbtDouble maxDelta = m_spRefData->GetMaxDihedral();
    RbtDouble initialValue = m_spRefData->GetInitialValue();
    RbtDouble delta = StandardisedValue(m_value - initialValue);
    if (delta > maxDelta) {
        m_value = StandardisedValue(initialValue + maxDelta);
    }
    else if (delta < -maxDelta) {
        m_value = StandardisedValue(initialValue - maxDelta);
    }
}
