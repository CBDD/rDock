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
