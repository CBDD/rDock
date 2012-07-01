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

#include "RbtSolventFlexData.h"

const RbtString& RbtSolventFlexData::_OCCUPANCY = "OCCUPANCY"; 
const RbtString& RbtSolventFlexData::_OCCUPANCY_STEP = "OCCUPANCY_STEP"; 

RbtSolventFlexData::RbtSolventFlexData(RbtDockingSite* pDockSite)
        : RbtLigandFlexData(pDockSite) {
    AddParameter(_OCCUPANCY, 1.0);
    AddParameter(_OCCUPANCY_STEP, 1.0);
}
