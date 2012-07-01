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

#include "RbtChromOccupancyRefData.h"
#include "RbtModel.h"

RbtString RbtChromOccupancyRefData::_CT = "RbtChromOccupancyRefData";

RbtChromOccupancyRefData::RbtChromOccupancyRefData(RbtModel* pModel,
						   RbtDouble stepSize,
                                                   RbtDouble threshold)
  : m_pModel(pModel),m_stepSize(stepSize),m_threshold(threshold)
{
    m_initialValue = GetModelValue();
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtChromOccupancyRefData::~RbtChromOccupancyRefData() {
    _RBTOBJECTCOUNTER_DESTR_(_CT);    
}

RbtDouble RbtChromOccupancyRefData::GetModelValue() const {
  return (m_pModel) ? m_pModel->GetOccupancy() : 0.0;
}

void RbtChromOccupancyRefData::SetModelValue(RbtDouble occupancy) {
  if (m_pModel) {
    m_pModel->SetOccupancy(occupancy,m_threshold);
  }
}
