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

#include "RbtPolarIdxSF.h"
#include "RbtWorkSpace.h"
#include "RbtPlane.h"

//Static data members
RbtString RbtPolarIdxSF::_CT("RbtPolarIdxSF");
RbtString RbtPolarIdxSF::_INCR("INCR");
RbtString RbtPolarIdxSF::_ATTR("ATTR");
RbtString RbtPolarIdxSF::_THRESHOLD_POS("THRESHOLD_POS");
RbtString RbtPolarIdxSF::_THRESHOLD_NEG("THRESHOLD_NEG");

//NB - Virtual base class constructor (RbtBaseSF) gets called first,
//implicit constructor for RbtBaseInterSF is called second
RbtPolarIdxSF::RbtPolarIdxSF(const RbtString& strName)
  : RbtBaseSF(_CT,strName),m_bAttr(true),m_bFlexRec(false),m_bSolvent(false),
  m_nPos(0),m_nNeg(0),m_posThreshold(0.25),m_negThreshold(0.25)
{
  //Add parameters
  AddParameter(_INCR,2.4);
  AddParameter(_ATTR,m_bAttr);
  AddParameter(_THRESHOLD_POS,m_posThreshold);
  AddParameter(_THRESHOLD_NEG,m_negThreshold);
#ifdef _DEBUG
  cout << _CT << " parameterised constructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtPolarIdxSF::~RbtPolarIdxSF() {
  ClearReceptor();
  ClearLigand();
  ClearSolvent();
#ifdef _DEBUG
  cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}
	
//Override RbtBaseSF::ScoreMap to provide additional raw descriptors
void RbtPolarIdxSF::ScoreMap(RbtStringVariantMap& scoreMap) const {
  if (isEnabled()) {
    //EnableAnnotations(true);
    //ClearAnnotationList();
    //We can only annotate the ligand-receptor interactions
    //as the rDock Viewer annotation format is hardwired to expect
    //ligand-receptor atom indices
    RbtDouble rs = InterScore();
    //EnableAnnotations(false);
    rs += LigandSolventScore();
    
    //First deal with the inter score which is stored in its natural location in the map
    RbtString name = GetFullName();
    scoreMap[name] = rs;
    AddToParentMapEntry(scoreMap, rs);
    
    //Now deal with the system raw scores which need to be stored in SCORE.INTER.POLAR
    RbtDouble system_rs = ReceptorScore() + SolventScore() + ReceptorSolventScore();
    if (system_rs != 0.0) {
        RbtString systemName = RbtBaseSF::_SYSTEM_SF + "." + GetName();
        scoreMap[systemName] = system_rs;
        //increment the parent SCORE.SYSTEM total with the weighted score
        RbtDouble parentScore = scoreMap[RbtBaseSF::_SYSTEM_SF];
        parentScore += system_rs * GetWeight();
        scoreMap[RbtBaseSF::_SYSTEM_SF] = parentScore;
    }

    //scoreMap[name+".nhbd"] = m_nPos;
    //scoreMap[name+".nhba"] = m_nNeg;
    //RbtStringList annList;
    //RenderAnnotationList(GetName(),annList);
    //scoreMap[RbtAnnotationHandler::_ANNOTATION_FIELD] += annList;
    //ClearAnnotationList();
  }
}

void RbtPolarIdxSF::SetupReceptor() {
  ClearReceptor();
  if (GetReceptor().Null())
    return;
  m_bFlexRec = GetReceptor()->isFlexible();
  RbtDouble idxIncr = GetParameter(_INCR);
  RbtDouble maxError = GetMaxError();
  idxIncr += maxError;
  RbtDouble flexDist = 2.0;//Additional indexing increment for flexibile OH and NH3
  RbtDockingSitePtr spDS = GetWorkSpace()->GetDockingSite();
  RbtInt iTrace = GetTrace();

  RbtInt nCoords = GetReceptor()->GetNumSavedCoords()-1;
  if (nCoords > 0) {
    RbtAtomList atomList = GetReceptor()->GetAtomList();
    m_spPosGrid = CreateInteractionGrid();
    m_recepPosList = CreateDonorInteractionCenters(atomList);
    m_spNegGrid = CreateInteractionGrid();
    m_recepNegList = CreateAcceptorInteractionCenters(atomList);
    for (RbtInt i = 1; i <= nCoords; i++) {
      if (iTrace > 0) {
	cout << _CT << ": Indexing receptor coords # " << i << endl;
      }
      GetReceptor()->RevertCoords(i);
      for (RbtInteractionCenterListConstIter iter = m_recepPosList.begin(); iter != m_recepPosList.end(); iter++) {
	RbtDouble rvdw = (*iter)->GetAtom1Ptr()->GetVdwRadius();
	m_spPosGrid->SetInteractionLists(*iter,rvdw+idxIncr);
      }
      m_spPosGrid->UniqueInteractionLists();
      for (RbtInteractionCenterListConstIter iter = m_recepNegList.begin(); iter != m_recepNegList.end(); iter++) {
	RbtDouble rvdw = (*iter)->GetAtom1Ptr()->GetVdwRadius();
	m_spNegGrid->SetInteractionLists(*iter,rvdw+idxIncr);
      }
      m_spNegGrid->UniqueInteractionLists();
    }
  }
  else {
    RbtAtomList atomList = spDS->GetAtomList(GetReceptor()->GetAtomList(),0.0,GetCorrectedRange());
    m_spPosGrid = CreateInteractionGrid();
    m_spNegGrid = CreateInteractionGrid();
    m_recepPosList = CreateDonorInteractionCenters(atomList);
    m_recepNegList = CreateAcceptorInteractionCenters(atomList);

    //For flexible receptors, separate the interaction centers into rigid and flexible
    if (m_bFlexRec) {
      GetReceptor()->SetAtomSelectionFlags(false);
      GetReceptor()->SelectFlexAtoms();//This leaves all moveable atoms selected
      Rbt::isInteractionCenterSelected isSel;
      RbtInteractionCenterListIter fIter = std::stable_partition(m_recepPosList.begin(),m_recepPosList.end(),std::not1(isSel));
      std::copy(fIter,m_recepPosList.end(),std::back_inserter(m_flexRecPosList));
      m_recepPosList.erase(fIter,m_recepPosList.end());

      fIter = std::stable_partition(m_recepNegList.begin(),m_recepNegList.end(),std::not1(isSel));
      std::copy(fIter,m_recepNegList.end(),std::back_inserter(m_flexRecNegList));
      m_recepNegList.erase(fIter,m_recepNegList.end());

      //Build map of intra-protein interactions (similar to intra-ligand)
      //Include flexible-flexible and flexible-rigid
      RbtInt nAtoms = GetReceptor()->GetNumAtoms();
      m_flexRecIntns = RbtInteractionListMap(nAtoms,RbtInteractionCenterList());
      m_flexRecPrtIntns = RbtInteractionListMap(nAtoms,RbtInteractionCenterList());
      if (m_bAttr) {
	BuildIntraMap(m_flexRecPosList,m_flexRecNegList,m_flexRecIntns);//flexible donor - flexible acceptor
	BuildIntraMap(m_flexRecPosList,m_recepNegList,m_flexRecIntns);//flexible donor - rigid acceptor
	BuildIntraMap(m_flexRecNegList,m_recepPosList,m_flexRecIntns);//flexible acceptor - rigid donor
      }
      else {
	BuildIntraMap(m_flexRecPosList,m_flexRecIntns);//flexible donor - flexible donor
	BuildIntraMap(m_flexRecPosList,m_recepPosList,m_flexRecIntns);//flexible donor - rigid donor
	BuildIntraMap(m_flexRecNegList,m_flexRecIntns);//flexible acceptor - flexible acceptor
	BuildIntraMap(m_flexRecNegList,m_recepNegList,m_flexRecIntns);//flexible acceptor - rigid acceptor
      }
      //We can get away with partitioning the variable interactions just at the beginning.
      //For grosser receptor flexibility we would have to partition periodically during docking
      RbtDouble partitionDist = GetRange() + flexDist;
      Partition(m_flexRecPosList,m_flexRecNegList,m_flexRecIntns,m_flexRecPrtIntns,partitionDist);
      //Index the flexible interaction centers over a larger radius
      //NOTE: WE ASSUME ONLY -OH and -NH3 rotation here (protons can't move more than 2.0A at most)
      //Grosser rotations will require a different approach
      for (RbtInteractionCenterListConstIter iter = m_flexRecPosList.begin(); iter != m_flexRecPosList.end(); iter++) {
	RbtDouble rvdw = (*iter)->GetAtom1Ptr()->GetVdwRadius();
	m_spPosGrid->SetInteractionLists(*iter,rvdw+idxIncr+flexDist);
      }
      for (RbtInteractionCenterListConstIter iter = m_flexRecNegList.begin(); iter != m_flexRecNegList.end(); iter++) {
	RbtDouble rvdw = (*iter)->GetAtom1Ptr()->GetVdwRadius();
	m_spNegGrid->SetInteractionLists(*iter,rvdw+idxIncr+flexDist);
      }
      if (iTrace > 0) {
	RbtDouble score = ReceptorScore();
	cout << GetWorkSpace()->GetName() << " " << GetFullName() << ": Intra-receptor score = " << score << endl;
      }
    }
    
    //Index the rigid interaction centers as usual
    for (RbtInteractionCenterListConstIter iter = m_recepPosList.begin(); iter != m_recepPosList.end(); iter++) {
      RbtDouble rvdw = (*iter)->GetAtom1Ptr()->GetVdwRadius();
      m_spPosGrid->SetInteractionLists(*iter,rvdw+idxIncr);
    }
    for (RbtInteractionCenterListConstIter iter = m_recepNegList.begin(); iter != m_recepNegList.end(); iter++) {
      RbtDouble rvdw = (*iter)->GetAtom1Ptr()->GetVdwRadius();
      m_spNegGrid->SetInteractionLists(*iter,rvdw+idxIncr);
    }
  }
}

void RbtPolarIdxSF::SetupLigand() {
  ClearLigand();
  if (GetLigand().Null())
    return;

  RbtAtomList atomList(GetLigand()->GetAtomList());
  m_ligPosList = CreateDonorInteractionCenters(atomList);
  m_ligNegList = CreateAcceptorInteractionCenters(atomList);
}

void RbtPolarIdxSF::SetupSolvent() {
  ClearSolvent();
  RbtModelList solventList = GetSolvent();
  if (solventList.empty())
  	return;
  
  m_bSolvent = true;
  
  //Concatenate all solvent atoms into a single list
  RbtAtomList solventAtomList;
  for (RbtModelListConstIter iter = solventList.begin(); iter != solventList.end(); ++iter) {
  	RbtAtomList atomList = (*iter)->GetAtomList();
   	std::copy(atomList.begin(),atomList.end(),std::back_inserter(solventAtomList));
  }
  
  //Create the interaction centers
  m_solventPosList = CreateDonorInteractionCenters(solventAtomList);
  m_solventNegList = CreateAcceptorInteractionCenters(solventAtomList);
  m_solventIntns = RbtInteractionListMap(solventAtomList.size(),RbtInteractionCenterList());

  //Build the indexed interaction map
  if (m_bAttr) {
    BuildIntraMap(m_solventPosList,m_solventNegList,m_solventIntns);
  }
  else {
    BuildIntraMap(m_solventPosList,m_solventIntns);
    BuildIntraMap(m_solventNegList,m_solventIntns);
  }
}

void RbtPolarIdxSF::SetupScore() {
  //No further setup required
}


//Overall score is the sum of:
//ligand-receptor
//ligand-solvent
//intra-receptor
//intra-solvent
//receptor-solvent
RbtDouble RbtPolarIdxSF::RawScore() const {
  return InterScore() + LigandSolventScore() + ReceptorScore() + SolventScore() +
  		ReceptorSolventScore();
}


//Clear the receptor and ligand grids and lists respectively
//As we are not using smart pointers, there is some memory management to do
void RbtPolarIdxSF::ClearReceptor() {
  m_spPosGrid = RbtInteractionGridPtr();
  m_spNegGrid = RbtInteractionGridPtr();
  m_flexRecIntns.clear();
  m_flexRecPrtIntns.clear();
  m_bFlexRec = false;
  DeleteList(m_recepPosList);
  DeleteList(m_flexRecPosList);
  DeleteList(m_recepNegList);
  DeleteList(m_flexRecNegList);
}

void RbtPolarIdxSF::ClearLigand() {
	DeleteList(m_ligPosList);
	DeleteList(m_ligNegList);
}

void RbtPolarIdxSF::ClearSolvent() {
	m_solventIntns.clear();
	m_bSolvent = false;
	DeleteList(m_solventPosList);
	DeleteList(m_solventNegList);
}

//Deletes the interaction center objects contained in the list, then clears the list
void RbtPolarIdxSF::DeleteList(RbtInteractionCenterList& icList) {
 for (RbtInteractionCenterListIter iter = icList.begin(); iter != icList.end(); ++iter) {
    delete *iter;
  }
  icList.clear();
}

//DM 25 Oct 2000 - track changes to parameter values in local data members
//ParameterUpdated is invoked by RbtParamHandler::SetParameter
void RbtPolarIdxSF::ParameterUpdated(const RbtString& strName) {
  //DM 25 Oct 2000 - heavily used params
  if (strName == _ATTR) {
    m_bAttr = GetParameter(_ATTR);
  }
  else if (strName == _THRESHOLD_POS) {
    m_posThreshold = GetParameter(_THRESHOLD_POS);
  }
  else if (strName == _THRESHOLD_NEG) {
    m_negThreshold = GetParameter(_THRESHOLD_NEG);
  }
  else {
    RbtPolarSF::OwnParameterUpdated(strName);
    RbtBaseIdxSF::OwnParameterUpdated(strName);
    RbtBaseSF::ParameterUpdated(strName);
  }
}

//Intra-receptor
RbtDouble RbtPolarIdxSF::ReceptorScore() const {
  return (m_bFlexRec) ? IntraScore(m_flexRecPosList,m_flexRecNegList,m_flexRecPrtIntns,m_bAttr) : 0.0;
}

//Intra-solvent
RbtDouble RbtPolarIdxSF::SolventScore() const {
  return (m_bSolvent) ? IntraScore(m_solventPosList,m_solventNegList,m_solventIntns,m_bAttr) : 0.0;
}

//Ligand-receptor
RbtDouble RbtPolarIdxSF::InterScore() const {
	return InterScore(m_ligPosList, m_ligNegList, true);
}

//Receptor-solvent
RbtDouble RbtPolarIdxSF::ReceptorSolventScore() const {
  return (m_bSolvent) ? InterScore(m_solventPosList,m_solventNegList, false) : 0.0;
}

//Ligand-solvent score: very inefficient
RbtDouble RbtPolarIdxSF::LigandSolventScore() const {
  RbtDouble score = 0.0;
  if (!m_bSolvent)
    return score;
  RbtPolarSF::f1prms Rprms = GetRprms();//Distance params
  RbtPolarSF::f1prms A1prms = GetA1prms();//Donor angle params
  RbtPolarSF::f1prms A2prms = GetA2prms();//Acceptor angle params
  if (m_bAttr) {
    //Pos-neg
    for (RbtInteractionCenterListConstIter iter = m_solventPosList.begin(); iter != m_solventPosList.end(); iter++) {
       score += (*iter)->GetAtom1Ptr()->GetUser1Value() * 
       						PolarScore(*iter,m_ligNegList,Rprms,A1prms,A2prms);
    }
    //Neg-pos
    for (RbtInteractionCenterListConstIter iter = m_solventNegList.begin(); iter != m_solventNegList.end(); iter++) {
       score += (*iter)->GetAtom1Ptr()->GetUser1Value() * 
       						PolarScore(*iter,m_ligPosList,Rprms,A2prms,A1prms);
    }
  }
  else {
    //pos-pos
    for (RbtInteractionCenterListConstIter iter = m_solventPosList.begin(); iter != m_solventPosList.end(); iter++) {
       score += (*iter)->GetAtom1Ptr()->GetUser1Value() * 
       						PolarScore(*iter,m_ligPosList,Rprms,A1prms,A1prms);
    }
    //neg-neg
    for (RbtInteractionCenterListConstIter iter = m_solventNegList.begin(); iter != m_solventNegList.end(); iter++) {
       score += (*iter)->GetAtom1Ptr()->GetUser1Value() * 
       						PolarScore(*iter,m_ligNegList,Rprms,A2prms,A2prms);
    }
  }
  return score;
}

//Reusable method for receptor-ligand and receptor-solvent scores
//bCount controls whether to count the positive and negative interaction scores
RbtDouble RbtPolarIdxSF::InterScore(const RbtInteractionCenterList& posList,
  						const RbtInteractionCenterList& negList,
  						RbtBool bCount) const {
  RbtDouble score = 0.0;//Total score
  if (bCount) {
  	m_nPos = 0;
  	m_nNeg = 0;
  }
  
  //Check grid is defined
  if (m_spPosGrid.Null() || m_spNegGrid.Null())
    return score;

  RbtPolarSF::f1prms Rprms = GetRprms();//Distance params
  RbtPolarSF::f1prms A1prms = GetA1prms();//Donor angle params
  RbtPolarSF::f1prms A2prms = GetA2prms();//Acceptor angle params

  RbtDouble s(0.0);//Partial scores
  //Ligand HBA
  for (RbtInteractionCenterListConstIter lIter = negList.begin(); lIter != negList.end(); lIter++) {
    RbtAtom* pLig1 = (*lIter)->GetAtom1Ptr();
    const RbtCoord& cLig1 = pLig1->GetCoords();
    //If this is an attractive potential we calculate the score with all adjacent +ve centres (HBD/M+/guan)
    if (m_bAttr) {
      const RbtInteractionCenterList& rList = m_spPosGrid->GetInteractionList(cLig1);
      s = PolarScore(*lIter,rList,Rprms,A2prms,A1prms);
    }
    else {
      //If this is an repulsive potential we calculate the score with all adjacent HBA
      const RbtInteractionCenterList& rList = m_spNegGrid->GetInteractionList(cLig1);
      s = PolarScore(*lIter,rList,Rprms,A2prms,A2prms);
    }
    s *= pLig1->GetUser1Value();
    if (bCount && (fabs(s) > m_negThreshold)) {
      m_nNeg++;
    }
    score += s;
  }

  //Ligand HBD
  for (RbtInteractionCenterListConstIter lIter = posList.begin(); lIter != posList.end(); lIter++) {
    RbtAtom* pLig1 = (*lIter)->GetAtom1Ptr();
    const RbtCoord& cLig1 = pLig1->GetCoords();
    //If this is an attractive potential we calculate the score with all adjacent HBA
    if (m_bAttr) {
      const RbtInteractionCenterList& rList = m_spNegGrid->GetInteractionList(cLig1);
      s = PolarScore(*lIter,rList,Rprms,A1prms,A2prms);
    }
    else {
      //If this is an repulsive potential we calculate the score with all adjacent +ve centres (HBD/M+/guan)
      const RbtInteractionCenterList& rList = m_spPosGrid->GetInteractionList(cLig1);
      s = PolarScore(*lIter,rList,Rprms,A1prms,A1prms);
    }
    s *= pLig1->GetUser1Value();
    if (bCount && (fabs(s) > m_posThreshold)) {
      m_nPos++;
    }
    score += s;
  }
  return score;
}
