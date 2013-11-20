/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtAromIdxSF.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtAromIdxSF.h"
#include "RbtWorkSpace.h"

//Static data members
RbtString RbtAromIdxSF::_CT("RbtAromIdxSF");
RbtString RbtAromIdxSF::_INCR("INCR");
RbtString RbtAromIdxSF::_R12("R12");
RbtString RbtAromIdxSF::_DR12MIN("DR12MIN");
RbtString RbtAromIdxSF::_DR12MAX("DR12MAX");
RbtString RbtAromIdxSF::_DAMIN("DAMIN");
RbtString RbtAromIdxSF::_DAMAX("DAMAX");
RbtString RbtAromIdxSF::_THRESHOLD("THRESHOLD");
	
RbtAromIdxSF::RbtAromIdxSF(const RbtString& strName)
  : RbtBaseSF(_CT,strName),
    m_R12(3.5),m_DR12Min(0.25),m_DR12Max(0.6),m_DAMin(20.0),m_DAMax(30.0),
    m_nArom(0),m_nGuan(0),m_threshold(0.25)
{
  //Add parameters
  AddParameter(_INCR,5.0);
  AddParameter(_R12,m_R12);
  AddParameter(_DR12MIN,m_DR12Min);
  AddParameter(_DR12MAX,m_DR12Max);
  AddParameter(_DAMIN,m_DAMin);
  AddParameter(_DAMAX,m_DAMax);
  AddParameter(_THRESHOLD,m_threshold);  
#ifdef _DEBUG
  cout << _CT << " parameterised constructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtAromIdxSF::~RbtAromIdxSF() {
  ClearReceptor();
  ClearLigand();
#ifdef _DEBUG
  cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}
	
//Override RbtBaseSF::ScoreMap to provide additional raw descriptors
/*void RbtAromIdxSF::ScoreMap(RbtStringVariantMap& scoreMap) const {
  if (isEnabled()) {
    EnableAnnotations(true);
    ClearAnnotationList();
    RbtDouble rs = RawScore();
    EnableAnnotations(false);

    RbtDouble w = GetWeight();
    RbtString name = GetFullName();
    RbtDouble s = w * rs;
    //DM 31 Mar 2003
    //New convention - report the raw scores and the weight in the SD file
    //rather than the weighted score and the raw score
    //Raw scores are less likely to change, whereas weights are quite volatile
    scoreMap[name] = rs;
    scoreMap[name+".w"] = w;
    scoreMap[name+".narom"] = m_nArom;
    scoreMap[name+".nguan"] = m_nGuan;
    //scoreMap[name+".ss"] = m_ss;
    //scoreMap[name+".sp"] = m_sp;
    //scoreMap[name+".pp"] = m_pp;
    RbtStringList annList;
    RenderAnnotationList(GetName(),annList);
    scoreMap[RbtAnnotationHandler::_ANNOTATION_FIELD] += annList;
    ClearAnnotationList();
    return s;
  }
  else {
    return 0.0;
  }
}*/

void RbtAromIdxSF::SetupReceptor() {
  ClearReceptor();
  if (GetReceptor().Null())
    return;

  m_spAromGrid = CreateInteractionGrid();
  m_spGuanGrid = CreateInteractionGrid();

  RbtDouble idxIncr = GetParameter(_INCR);//vdw Radius increment for indexing
  RbtDouble maxError = GetMaxError();
  idxIncr += maxError;

  //First, index the receptor aromatic rings
  //Limit to those within RANGE of the docking site
  //Have to filter manually as GetReceptorAtomsInRange not appropriate for ring centers
  RbtAtomListList recepRingLists = GetReceptor()->GetRingAtomLists();
  RbtDockingSitePtr spDS = GetWorkSpace()->GetDockingSite();

  RbtInt nCoords = GetReceptor()->GetNumSavedCoords()-1;
  if (nCoords > 0) {
    for (RbtAtomListListConstIter rIter = recepRingLists.begin(); rIter != recepRingLists.end(); rIter++) {
      if (Rbt::GetNumAtoms(*rIter,Rbt::isPiAtom()) == (*rIter).size()) {
	RbtPseudoAtomPtr spPseudoAtom = GetReceptor()->AddPseudoAtom(*rIter);
	RbtAtomList recepPiAtoms = spPseudoAtom->GetAtomList();
	RbtAtomPtr spRecepPiAtom1 = recepPiAtoms[0];
	RbtAtomPtr spRecepPiAtom2 = recepPiAtoms[1];
	RbtInteractionCenter* pIntnCenter(new RbtInteractionCenter(spPseudoAtom,spRecepPiAtom1,spRecepPiAtom2));
	m_recepAromList.push_back(pIntnCenter);
      }
    }
    RbtAtomList recepGuanList = Rbt::GetAtomList(GetReceptor()->GetAtomList(),Rbt::isAtomGuanidiniumCarbon());
    for (RbtAtomListConstIter iter = recepGuanList.begin(); iter != recepGuanList.end(); iter++) {
      RbtAtomList recepPiAtoms = Rbt::GetBondedAtomList(*iter);
      RbtAtomPtr spRecepPiAtom1 = recepPiAtoms[0];
      RbtAtomPtr spRecepPiAtom2 = recepPiAtoms[1];
      RbtInteractionCenter* pIntnCenter(new RbtInteractionCenter(*iter,spRecepPiAtom1,spRecepPiAtom2));
      m_recepGuanList.push_back(pIntnCenter);//Store the interaction center
    }

    for (RbtInt i = 1; i <= nCoords; i++) {
      cout << _CT << ": Indexing receptor coords # " << i << endl;
      GetReceptor()->RevertCoords(i);
      for (RbtInteractionCenterListConstIter iter = m_recepAromList.begin(); iter != m_recepAromList.end(); iter++) {
	m_spAromGrid->SetInteractionLists(*iter,idxIncr);
      }
      m_spAromGrid->UniqueInteractionLists();
      for (RbtInteractionCenterListConstIter iter = m_recepGuanList.begin(); iter != m_recepGuanList.end(); iter++) {
	m_spGuanGrid->SetInteractionLists(*iter,idxIncr);
      }
      m_spGuanGrid->UniqueInteractionLists();
    }
  }
  else {
    RbtDockingSite::isAtomInRange bIsInRange(spDS->GetGrid(),0.0,GetCorrectedRange());
    for (RbtAtomListListConstIter rIter = recepRingLists.begin(); rIter != recepRingLists.end(); rIter++) {
      //Check that all ring atoms are pi-atoms (crude test for aromaticity)
      if (Rbt::GetNumAtoms(*rIter,Rbt::isPiAtom()) == (*rIter).size()) {
	RbtPseudoAtomPtr spPseudoAtom = GetReceptor()->AddPseudoAtom(*rIter);
	if (bIsInRange(spPseudoAtom)) {
	  RbtAtomList recepPiAtoms = spPseudoAtom->GetAtomList();
	  RbtAtomPtr spRecepPiAtom1 = recepPiAtoms[0];
	  RbtAtomPtr spRecepPiAtom2 = recepPiAtoms[1];
	  RbtInteractionCenter* pIntnCenter(new RbtInteractionCenter(spPseudoAtom,spRecepPiAtom1,spRecepPiAtom2));
	  m_recepAromList.push_back(pIntnCenter);
	  m_spAromGrid->SetInteractionLists(pIntnCenter,idxIncr);
	}
      }
    }

    //Now the guanidinium groups (positively charged carbons)
    //Limit to those within RANGE of the docking site
    RbtAtomList atomList = spDS->GetAtomList(GetReceptor()->GetAtomList(),0.0,GetCorrectedRange());
    RbtAtomList recepGuanList = Rbt::GetAtomList(atomList,Rbt::isAtomGuanidiniumCarbon());
    
    for (RbtAtomListConstIter iter = recepGuanList.begin(); iter != recepGuanList.end(); iter++) {
      RbtAtomList recepPiAtoms = Rbt::GetBondedAtomList(*iter);
      RbtAtomPtr spRecepPiAtom1 = recepPiAtoms[0];
      RbtAtomPtr spRecepPiAtom2 = recepPiAtoms[1];
      RbtInteractionCenter* pIntnCenter(new RbtInteractionCenter(*iter,spRecepPiAtom1,spRecepPiAtom2));
      m_recepGuanList.push_back(pIntnCenter);//Store the interaction center
      m_spGuanGrid->SetInteractionLists(pIntnCenter,idxIncr);
    }
  }
}

void RbtAromIdxSF::SetupLigand() {
  ClearLigand();
  if (GetLigand().Null())
    return;

  //First, create the ligand aromatic interaction centers
  RbtAtomListList ligRingLists = GetLigand()->GetRingAtomLists();
  for (RbtAtomListListConstIter rIter = ligRingLists.begin(); rIter != ligRingLists.end(); rIter++) {
    //Check that all ring atoms are pi-atoms (crude test for aromaticity)
    if (Rbt::GetNumAtoms(*rIter,Rbt::isPiAtom()) == (*rIter).size()) {
      RbtPseudoAtomPtr spPseudoAtom = GetLigand()->AddPseudoAtom(*rIter);
      //RbtAtomList ligPiAtoms = Rbt::GetBondedAtomList(spPseudoAtom);
      RbtAtomList ligPiAtoms = spPseudoAtom->GetAtomList();
      RbtAtomPtr spLigPiAtom1 = ligPiAtoms[0];
      RbtAtomPtr spLigPiAtom2 = ligPiAtoms[1];
      RbtInteractionCenter* pIntnCenter(new RbtInteractionCenter(spPseudoAtom,spLigPiAtom1,spLigPiAtom2));
      m_ligAromList.push_back(pIntnCenter);
    }      
  }

  //Now the guanidinium groups (positively charged carbons)
  RbtAtomList ligGuanList = Rbt::GetAtomList(GetLigand()->GetAtomList(),Rbt::isAtomGuanidiniumCarbon());
  for (RbtAtomListConstIter iter = ligGuanList.begin(); iter != ligGuanList.end(); iter++) {
    RbtAtomList ligPiAtoms = Rbt::GetBondedAtomList(*iter);
    RbtAtomPtr spLigPiAtom1 = ligPiAtoms[0];
    RbtAtomPtr spLigPiAtom2 = ligPiAtoms[1];
    RbtInteractionCenter* pIntnCenter(new RbtInteractionCenter(*iter,spLigPiAtom1,spLigPiAtom2));
    m_ligGuanList.push_back(pIntnCenter);//Store the interaction center
  }
}

void RbtAromIdxSF::SetupScore() {
  //No further setup required
}

RbtDouble RbtAromIdxSF::RawScore() const {
  RbtDouble score = 0.0;
  m_nArom = 0;
  m_nGuan = 0;

  //Check grids are defined
  if (m_spAromGrid.Null() || m_spGuanGrid.Null())
    return score;

  f1prms Rprms = GetRprms();//Distance params
  f1prms Aprms = GetAprms();//Donor angle params  

  //Loop over ligand aromatic interaction centers
  for (RbtInteractionCenterListConstIter ligIter = m_ligAromList.begin(); ligIter != m_ligAromList.end(); ligIter++) {
    const RbtCoord& cLig1 = (*ligIter)->GetAtom1Ptr()->GetCoords();
    //Get the list of nearby receptor aromatic centers
    const RbtInteractionCenterList& recepAromList = m_spAromGrid->GetInteractionList(cLig1);
    //Get the list of nearby receptor guanidinium centers
    const RbtInteractionCenterList& recepGuanList = m_spGuanGrid->GetInteractionList(cLig1);
    
    RbtDouble s = AromScore(*ligIter,recepAromList,Rprms,Aprms);
    s += AromScore(*ligIter,recepGuanList,Rprms,Aprms);
    //RbtDouble s = PiScore(*ligIter,recepAromList);
    if (s > m_threshold) {
      m_nArom++;
    }
    score += s;
  }
  
  //Loop over ligand guanidinium interaction centers
  for (RbtInteractionCenterListConstIter ligIter = m_ligGuanList.begin(); ligIter != m_ligGuanList.end(); ligIter++) {
    const RbtCoord& cLig1 = (*ligIter)->GetAtom1Ptr()->GetCoords();
    //Get the list of nearby receptor aromatic centers
    const RbtInteractionCenterList& recepAromList = m_spAromGrid->GetInteractionList(cLig1);
    
    RbtDouble s = AromScore(*ligIter,recepAromList,Rprms,Aprms);
    //RbtDouble s = 0.0;
    if (s > m_threshold) {
      m_nGuan++;
    }
    score += s;
  }
  return score;
}


//Clear the receptor and ligand grids and lists respectively
//As we are not using smart pointers, there is some memory management to do
void RbtAromIdxSF::ClearReceptor() {
  //Wipe the grids
  m_spAromGrid = RbtInteractionGridPtr();
  m_spGuanGrid = RbtInteractionGridPtr();
  //Delete the receptor interaction centers
  for (RbtInteractionCenterListIter iter = m_recepAromList.begin(); iter != m_recepAromList.end(); iter++) {
    delete *iter;
  }
  m_recepAromList.clear();
  for (RbtInteractionCenterListIter iter = m_recepGuanList.begin(); iter != m_recepGuanList.end(); iter++) {
    delete *iter;
  }
  m_recepGuanList.clear();
}
void RbtAromIdxSF::ClearLigand() {
  //Delete the ligand interaction centers
  for (RbtInteractionCenterListIter iter = m_ligAromList.begin(); iter != m_ligAromList.end(); iter++) {
    delete *iter;
  }
  m_ligAromList.clear();
  for (RbtInteractionCenterListIter iter = m_ligGuanList.begin(); iter != m_ligGuanList.end(); iter++) {
    delete *iter;
  }
  m_ligGuanList.clear();
}

//DM 25 Oct 2000 - track changes to parameter values in local data members
//ParameterUpdated is invoked by RbtParamHandler::SetParameter
void RbtAromIdxSF::ParameterUpdated(const RbtString& strName) {
  if (strName == _R12) {
    m_R12 = GetParameter(_R12);
  }
  else if (strName == _DR12MIN) {
    m_DR12Min = GetParameter(_DR12MIN);
  }
  else if (strName == _DR12MAX) {
    m_DR12Max = GetParameter(_DR12MAX);
  }
  else if (strName == _DAMIN) {
    m_DAMin = GetParameter(_DAMIN);
  }
  else if (strName == _DAMAX) {
    m_DAMax = GetParameter(_DAMAX);
  }
  else if (strName == _THRESHOLD) {
    m_threshold = GetParameter(_THRESHOLD);
  }
  else {
    RbtBaseIdxSF::OwnParameterUpdated(strName);
    RbtBaseSF::ParameterUpdated(strName);
  }
}

//The actual aromatic score, between a given interaction center and a list of near neighbour centers
RbtDouble RbtAromIdxSF::AromScore(const RbtInteractionCenter* pIC1, const RbtInteractionCenterList& IC2List,
				  const f1prms& Rprms, const f1prms& Aprms) const {
  RbtDouble s(0.0);
  if (IC2List.empty()) {
    return s;
  }

  RbtBool bAnnotate = isAnnotationEnabled();

  RbtAtom* pAtom1_1 = pIC1->GetAtom1Ptr();
  const RbtCoord& cAtom1_1 = pAtom1_1->GetCoords();
  RbtAtom* pAtom1_2 = pIC1->GetAtom2Ptr();
  RbtAtom* pAtom1_3 = pIC1->GetAtom3Ptr();
  RbtPlane pl1 = RbtPlane(cAtom1_1,pAtom1_2->GetCoords(),pAtom1_3->GetCoords());

  for (RbtInteractionCenterListConstIter IC2Iter = IC2List.begin(); IC2Iter != IC2List.end(); IC2Iter++) {
    RbtAtom* pAtom2_1 = (*IC2Iter)->GetAtom1Ptr();
    const RbtCoord& cAtom2_1 = pAtom2_1->GetCoords();
    RbtAtom* pAtom2_2 = (*IC2Iter)->GetAtom2Ptr();
    RbtAtom* pAtom2_3 = (*IC2Iter)->GetAtom3Ptr();
    RbtPlane pl2 = RbtPlane(cAtom2_1,pAtom2_2->GetCoords(),pAtom2_3->GetCoords());

    //Calculate average perp. distance from each pi center to the other ring plane
    RbtDouble R = (fabs(Rbt::DistanceFromPointToPlane(cAtom1_1, pl2)) +
                       fabs(Rbt::DistanceFromPointToPlane(cAtom2_1, pl1)))/2.0;
    RbtDouble f = f1(R-Rprms.R0,Rprms);
    //Only calculate average slip angle if f  > 0
    if (f > 0.0) {
      RbtVector v = cAtom2_1-cAtom1_1;
      RbtDouble sa = (acos(fabs(Rbt::Dot(v.Unit(),pl1.VNorm())))+acos(fabs(Rbt::Dot(v.Unit(),pl2.VNorm()))))*90.0/M_PI;
      f *= f1(sa,Aprms);
      if (f > 0.0) {
	s += f;
	if (bAnnotate) {
	  //Convention is that favourable annotations are negative, so invert the raw score here
	  RbtAnnotationPtr spAnnotation(new RbtAnnotation(pAtom1_1,pAtom2_1,R,-f));
	  AddAnnotation(spAnnotation);
	}
      }
    }
  }
  return s;
}


//The actual aromatic score, between a given interaction center and a list of near neighbour centers
RbtDouble RbtAromIdxSF::PiScore(const RbtInteractionCenter* pIC1, const RbtInteractionCenterList& IC2List) const {
  m_ss = 0.0;//sigma-sigma
  m_pp = 0.0;//pi-pi
  m_sp = 0.0;//sigma-pi

  if (IC2List.empty()) {
    return 0.0;
  }

  RbtDouble delta = 0.47;

  RbtAtom* pAtom1_1 = pIC1->GetAtom1Ptr();
  RbtAtom* pAtom1_2 = pIC1->GetAtom2Ptr();
  RbtAtom* pAtom1_3 = pIC1->GetAtom3Ptr();
  RbtPlane pl1 = RbtPlane(pAtom1_1->GetCoords(),pAtom1_2->GetCoords(),pAtom1_3->GetCoords());
  RbtVector v1 = delta*pl1.VNorm();
  //First calculate and store the sigma and pi-system coords (above and below plane)
  RbtPseudoAtom* pseudo1 = dynamic_cast<RbtPseudoAtom*>(pAtom1_1);
  RbtCoordList sigma1;
  RbtCoordList pi1;
  if (pseudo1) {
    RbtAtomList al1 = pseudo1->GetAtomList();
    for (RbtAtomListConstIter iter = al1.begin(); iter != al1.end(); iter++) {
      sigma1.push_back((*iter)->GetCoords());
      pi1.push_back((*iter)->GetCoords() + v1);
      pi1.push_back((*iter)->GetCoords() - v1);
    }
  }

  for (RbtInteractionCenterListConstIter IC2Iter = IC2List.begin(); IC2Iter != IC2List.end(); IC2Iter++) {
    RbtAtom* pAtom2_1 = (*IC2Iter)->GetAtom1Ptr();
    RbtAtom* pAtom2_2 = (*IC2Iter)->GetAtom2Ptr();
    RbtAtom* pAtom2_3 = (*IC2Iter)->GetAtom3Ptr();
    RbtPlane pl2 = RbtPlane(pAtom2_1->GetCoords(),pAtom2_2->GetCoords(),pAtom2_3->GetCoords());
    RbtVector v2 = delta*pl2.VNorm();
    //First calculate and store the sigma and pi-system coords (above and below plane)
    RbtPseudoAtom* pseudo2 = dynamic_cast<RbtPseudoAtom*>(pAtom2_1);
    RbtCoordList sigma2;
    RbtCoordList pi2;
    if (pseudo2) {
      RbtAtomList al2 = pseudo2->GetAtomList();
      for (RbtAtomListConstIter iter = al2.begin(); iter != al2.end(); iter++) {
	sigma2.push_back((*iter)->GetCoords());
	pi2.push_back((*iter)->GetCoords() + v2);
	pi2.push_back((*iter)->GetCoords() - v2);
      }
    }


    for (RbtCoordListConstIter cIter1 = sigma1.begin(); cIter1 != sigma1.end(); cIter1++) {
      for (RbtCoordListConstIter cIter2 = sigma2.begin(); cIter2 != sigma2.end(); cIter2++) {
	RbtDouble r = Rbt::Length(*cIter1,*cIter2) + 0.0001;//prevent div by zero
	m_ss += (1.0/r);
      }
      for (RbtCoordListConstIter cIter2 = pi2.begin(); cIter2 != pi2.end(); cIter2++) {
	RbtDouble r = Rbt::Length(*cIter1,*cIter2) + 0.0001;//prevent div by zero
	m_sp += (1.0/r);
      }
    }
    for (RbtCoordListConstIter cIter1 = pi1.begin(); cIter1 != pi1.end(); cIter1++) {
      for (RbtCoordListConstIter cIter2 = pi2.begin(); cIter2 != pi2.end(); cIter2++) {
	RbtDouble r = Rbt::Length(*cIter1,*cIter2) + 0.0001;//prevent div by zero
	m_pp += (1.0/r);
      }
      for (RbtCoordListConstIter cIter2 = sigma2.begin(); cIter2 != sigma2.end(); cIter2++) {
	RbtDouble r = Rbt::Length(*cIter1,*cIter2) + 0.0001;//prevent div by zero
	m_sp += (1.0/r);
      }
    }
  }
  //Correct the totals for the absolute charges (sigma = +1; pi = -0.5) and convert to kcal
  m_ss *= 332.17;
  m_sp *= (-0.5*332.17);
  m_pp *= (0.25*332.17);
  return m_ss + m_sp + m_pp;
}
