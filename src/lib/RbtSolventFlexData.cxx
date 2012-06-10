/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtSolventFlexData.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtSolventFlexData.h"

const RbtString& RbtSolventFlexData::_OCCUPANCY = "OCCUPANCY"; 
const RbtString& RbtSolventFlexData::_OCCUPANCY_STEP = "OCCUPANCY_STEP"; 

RbtSolventFlexData::RbtSolventFlexData(RbtDockingSite* pDockSite)
        : RbtLigandFlexData(pDockSite) {
    AddParameter(_OCCUPANCY, 1.0);
    AddParameter(_OCCUPANCY_STEP, 1.0);
}
