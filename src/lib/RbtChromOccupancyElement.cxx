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

#include "RbtChromOccupancyElement.h"

RbtString RbtChromOccupancyElement::_CT = "RbtChromOccupancyElement";

RbtChromOccupancyElement::RbtChromOccupancyElement(RbtModel* pModel,
						   RbtDouble stepSize,
                                                   RbtDouble threshold)
        : m_value(1.0) {
  m_spRefData = new RbtChromOccupancyRefData(pModel, stepSize, threshold);
  //Set the initial genotype to match the current phenotype
  SyncFromModel();
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtChromOccupancyElement::RbtChromOccupancyElement(
  RbtChromOccupancyRefDataPtr spRefData, RbtDouble value)
  : m_spRefData(spRefData), m_value(value) {
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtChromOccupancyElement::~RbtChromOccupancyElement() {
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtChromOccupancyElement::Reset() {
  m_value = m_spRefData->GetInitialValue();    
}

void RbtChromOccupancyElement::Randomise() {
  m_value = GetRand().GetRandom01();
}

void RbtChromOccupancyElement::Mutate(RbtDouble relStepSize) {
  RbtDouble absStepSize = relStepSize * m_spRefData->GetStepSize();
  RbtDouble delta;
  if (absStepSize > 0) {
    delta = 2.0 * absStepSize * GetRand().GetRandom01() - absStepSize;
    m_value = StandardisedValue(m_value + delta);
  }
}

void RbtChromOccupancyElement::SyncFromModel() {
  m_value = m_spRefData->GetModelValue();
}

void RbtChromOccupancyElement::SyncToModel() {
  m_spRefData->SetModelValue(m_value);
}

RbtChromElement* RbtChromOccupancyElement::clone() const {
  return new RbtChromOccupancyElement(m_spRefData, m_value);
}

void RbtChromOccupancyElement::GetVector(RbtDoubleList& v) const {
  v.push_back(m_value);
}

void RbtChromOccupancyElement::GetVector(RbtXOverList& v) const {
    RbtXOverElement occupancyElement;
    occupancyElement.push_back(m_value);
    v.push_back(occupancyElement);
}

void RbtChromOccupancyElement::SetVector(const RbtDoubleList& v, RbtInt& i) throw (RbtError) {
    if (VectorOK(v,i)) {
        m_value = StandardisedValue(v[i++]);
    }
    else {
        throw RbtBadArgument(_WHERE_, "Index out of range or insufficient elements remaining");
    }
}

void RbtChromOccupancyElement::SetVector(const RbtXOverList& v, RbtInt& i) throw (RbtError) {
    if (VectorOK(v,i)) {
        RbtXOverElement occupancyElement(v[i++]);
        if (occupancyElement.size() == 1) {
            m_value = occupancyElement[0];
        }
        else {
            throw RbtBadArgument(_WHERE_, "occupancyElement vector is of incorrect length");
        }
    }
    else {
        throw RbtBadArgument(_WHERE_, "Index out of range or insufficient elements remaining");
    }
}

void RbtChromOccupancyElement::GetStepVector(RbtDoubleList& v) const {
    v.push_back(m_spRefData->GetStepSize());
}

RbtDouble RbtChromOccupancyElement::CompareVector(const RbtDoubleList& v, RbtInt& i) const {
    RbtDouble retVal(0.0);
    if (!VectorOK(v,i)) {
        retVal = -1.0;
    }
    else {
        RbtDouble otherOccupancy = v[i++];
        RbtDouble stepSize = m_spRefData->GetStepSize();
        if (stepSize > 0.0) {
            RbtDouble absDiff = fabs(m_value - otherOccupancy);
            retVal = absDiff / stepSize;
        }
    }
    return retVal;
}

void RbtChromOccupancyElement::Print(ostream& s) const {
    s << "OCCUPANCY " << m_value << endl;
}

RbtDouble RbtChromOccupancyElement::StandardisedValue(RbtDouble occupancy) {
  if (occupancy < 0.0) {
    occupancy = 0.0;
  }
  else if (occupancy > 1.0) {
    occupancy = 1.0;
  }
  return occupancy;
}
