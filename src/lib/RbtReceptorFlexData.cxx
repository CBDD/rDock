/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtReceptorFlexData.cxx#2 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtReceptorFlexData.h"

const RbtString& RbtReceptorFlexData::_FLEX_DISTANCE = "FLEX_DISTANCE"; 
const RbtString& RbtReceptorFlexData::_DIHEDRAL_STEP = "DIHEDRAL_STEP"; 

RbtReceptorFlexData::RbtReceptorFlexData(RbtDockingSite* pDockSite)
        : RbtFlexData(pDockSite) {
    AddParameter(_FLEX_DISTANCE, 3.0);
    AddParameter(_DIHEDRAL_STEP, 30.0);
}
