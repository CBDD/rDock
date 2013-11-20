/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtChromOccupancyRefData.cxx#2 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

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
