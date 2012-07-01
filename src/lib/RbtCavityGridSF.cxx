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

#include "RbtCavityGridSF.h"
#include "RbtWorkSpace.h"
#include "RbtReceptorFlexData.h"
#include "RbtLigandFlexData.h"
#include "RbtSolventFlexData.h"
#include "RbtChromPositionRefData.h"

//Static data members
RbtString RbtCavityGridSF::_CT("RbtCavityGridSF");
RbtString RbtCavityGridSF::_RMAX("RMAX");
RbtString RbtCavityGridSF::_QUADRATIC("QUADRATIC");
	
//NB - Virtual base class constructor (RbtBaseSF) gets called first,
//implicit constructor for RbtBaseInterSF is called second
RbtCavityGridSF::RbtCavityGridSF(const RbtString& strName)
  : RbtBaseSF(_CT,strName),m_maxDist(0.0),m_rMax(0.1),m_bQuadratic(false) {
  //Add parameters
  AddParameter(_RMAX,m_rMax);
  AddParameter(_QUADRATIC,m_bQuadratic);
#ifdef _DEBUG
  cout << _CT << " parameterised constructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtCavityGridSF::~RbtCavityGridSF() {
#ifdef _DEBUG
  cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}
	

void RbtCavityGridSF::SetupReceptor() {
  //DM 09 Apr 2002 - we can reuse the precalculated distance grid directly from the docking site
  m_spGrid = RbtRealGridPtr();
  RbtDockingSitePtr spDS = GetWorkSpace()->GetDockingSite();
  if (spDS.Null())
    return;
  m_spGrid = spDS->GetGrid();
  if (m_spGrid.Null())
    return;
  m_maxDist = m_spGrid->MaxValue();
}

void RbtCavityGridSF::SetupLigand() {}

void RbtCavityGridSF::SetupSolvent() {}

void RbtCavityGridSF::SetupScore() {
    m_atomList.clear();
    RbtWorkSpace* pWorkSpace = GetWorkSpace();
    if (pWorkSpace) {
	   HeavyAtomFactory atomFactory(pWorkSpace->GetModels());
       m_atomList = atomFactory.GetAtomList();
    }
}

RbtDouble RbtCavityGridSF::RawScore() const {
  RbtDouble score(0.0);

  //Check grid is defined
  if (m_spGrid.Null())
    return score;

  //Quadratic penalty function
  if (m_bQuadratic) {
    for (RbtAtomRListConstIter iter = m_atomList.begin(); iter != m_atomList.end(); iter++) {
      const RbtCoord& c = (*iter)->GetCoords();
      //Check if atom is off grid, if so default to max grid value
      RbtDouble r = m_spGrid->isValid(c) ? m_spGrid->GetSmoothedValue(c) : m_maxDist;
      RbtDouble dr = r - m_rMax;
      if (dr > 0.0) {
        score += dr*dr;
      }
    }
  }
  //Linear penalty function
  else {
    for (RbtAtomRListConstIter iter = m_atomList.begin(); iter != m_atomList.end(); iter++) {
      const RbtCoord& c = (*iter)->GetCoords();
      //Check if atom is off grid, if so default to max grid value
      RbtDouble r = m_spGrid->isValid(c) ? m_spGrid->GetSmoothedValue(c) : m_maxDist;
      RbtDouble dr = r - m_rMax;
      if (dr > 0.0) {
        score += dr;
      }
    }
  }
  return score;
}

//DM 25 Oct 2000 - track changes to parameter values in local data members
//ParameterUpdated is invoked by RbtParamHandler::SetParameter
void RbtCavityGridSF::ParameterUpdated(const RbtString& strName) {
  if (strName == _RMAX) {
    m_rMax = GetParameter(_RMAX);
  }
  else if (strName == _QUADRATIC) {
    m_bQuadratic = GetParameter(_QUADRATIC);
  }
  else {
    RbtBaseSF::ParameterUpdated(strName);
  }
}

RbtCavityGridSF::HeavyAtomFactory::HeavyAtomFactory(RbtModelList modelList) {
    for (RbtModelListConstIter iter = modelList.begin(); iter != modelList.end(); iter++) {
        if ( (*iter).Ptr() != NULL) {
            RbtFlexData* pFlexData = (*iter)->GetFlexData();
            if (pFlexData != NULL) {
                pFlexData->Accept(*this);
            }
        }
    }    
}


void RbtCavityGridSF::HeavyAtomFactory::VisitReceptorFlexData(RbtReceptorFlexData* pFlexData) {
    //Receptor atoms: do nothing (the only movable atoms are hydrogens at present)    
}

void RbtCavityGridSF::HeavyAtomFactory::VisitLigandFlexData(RbtLigandFlexData* pFlexData) {
    //Extract all non-hydrogen atoms from ligand
    RbtAtomList atomList = pFlexData->GetModel()->GetAtomList();
    std::copy_if(atomList.begin(), atomList.end(),
                std::back_inserter(m_atomList),
                std::not1(Rbt::isAtomicNo_eq(1)));
}

void RbtCavityGridSF::HeavyAtomFactory::VisitSolventFlexData(RbtSolventFlexData* pFlexData) {
    //Extract all non-hydrogen atoms from solvent but only if the solvent can translate.
    //As there are no dihedrals to worry about, the water oxygens are fixed in position
    //for orientation-only flexibility.
    RbtString transModeStr = pFlexData->GetParameter(RbtLigandFlexData::_TRANS_MODE);
    RbtDouble maxTrans = pFlexData->GetParameter(RbtLigandFlexData::_MAX_TRANS);
    RbtChromElement::eMode eTransMode = RbtChromElement::StrToMode(transModeStr);
    if ( (eTransMode == RbtChromElement::FREE)
        || ((eTransMode == RbtChromElement::TETHERED) && (maxTrans > 0.0)) ) {
        RbtAtomList atomList = pFlexData->GetModel()->GetAtomList();
        std::copy_if(atomList.begin(), atomList.end(),
                    std::back_inserter(m_atomList),
                    std::not1(Rbt::isAtomicNo_eq(1)));
    }
}
 
