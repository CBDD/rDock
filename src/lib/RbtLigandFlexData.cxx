/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtLigandFlexData.cxx#2 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtLigandFlexData.h"
#include "RbtChromElement.h"

const RbtString& RbtLigandFlexData::_TRANS_STEP = "TRANS_STEP"; 
const RbtString& RbtLigandFlexData::_ROT_STEP = "ROT_STEP"; 
const RbtString& RbtLigandFlexData::_DIHEDRAL_STEP = "DIHEDRAL_STEP"; 
const RbtString& RbtLigandFlexData::_TRANS_MODE = "TRANS_MODE"; 
const RbtString& RbtLigandFlexData::_ROT_MODE = "ROT_MODE"; 
const RbtString& RbtLigandFlexData::_DIHEDRAL_MODE = "DIHEDRAL_MODE"; 
const RbtString& RbtLigandFlexData::_MAX_TRANS = "MAX_TRANS"; 
const RbtString& RbtLigandFlexData::_MAX_ROT = "MAX_ROT"; 
const RbtString& RbtLigandFlexData::_MAX_DIHEDRAL = "MAX_DIHEDRAL"; 

RbtLigandFlexData::RbtLigandFlexData(RbtDockingSite* pDockSite)
        : RbtFlexData(pDockSite) {
    AddParameter(_TRANS_STEP, 2.0);
    AddParameter(_ROT_STEP, 30.0);
    AddParameter(_DIHEDRAL_STEP, 30.0);
    AddParameter(_TRANS_MODE, RbtChromElement::ModeToStr(RbtChromElement::FREE));
    AddParameter(_ROT_MODE, RbtChromElement::ModeToStr(RbtChromElement::FREE));
    AddParameter(_DIHEDRAL_MODE, RbtChromElement::ModeToStr(RbtChromElement::FREE));
    AddParameter(_MAX_TRANS, 1.0);
    AddParameter(_MAX_ROT, 30.0);
    AddParameter(_MAX_DIHEDRAL, 30.0);
}
