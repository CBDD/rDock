/***********************************************************************
* The rDock program was developed from 1998 - 2006 by the software team 
* at RiboTargets (subsequently Vernalis (R&D) Ltd).
* In 2006, the software was licensed to the University of York for 
* maintenance and distribution.
* In 2012, Vernalis and the University of York agreed to release the 
* program as Open Source software.
* This version is licensed under GNU-LGPL version 3.0 with support from
* the University of Barcelona.
* http://rdock.sourceforge.net/
***********************************************************************/

#include "RbtVdwIdxSF.h"
#include "RbtWorkSpace.h"
#include "RbtFlexAtomFactory.h"

//Static data members
RbtString RbtVdwIdxSF::_CT("RbtVdwIdxSF");
RbtString RbtVdwIdxSF::_THRESHOLD_ATTR("THRESHOLD_ATTR");
RbtString RbtVdwIdxSF::_THRESHOLD_REP("THRESHOLD_REP");
RbtString RbtVdwIdxSF::_ANNOTATION_LIPO("ANNOTATION_LIPO");
RbtString RbtVdwIdxSF::_ANNOTATE("ANNOTATE");
RbtString RbtVdwIdxSF::_FAST_SOLVENT("FAST_SOLVENT");

//NB - Virtual base class constructor (RbtBaseSF) gets called first,
//implicit constructor for RbtBaseInterSF is called second
RbtVdwIdxSF::RbtVdwIdxSF(const RbtString& strName)
  : RbtBaseSF(_CT,strName),m_nAttr(0),m_nRep(0),m_attrThreshold(-0.5),m_repThreshold(0.5),
    m_lipoAnnot(-0.1),m_bAnnotate(true),m_bFlexRec(false),m_bFastSolvent(true)
{
  AddParameter(_THRESHOLD_ATTR,m_attrThreshold);
  AddParameter(_THRESHOLD_REP,m_repThreshold);
  AddParameter(_ANNOTATION_LIPO,m_lipoAnnot);//Threshold for outputting lipo vdW annotations
  AddParameter(_ANNOTATE,m_bAnnotate);//Threshold for outputting lipo vdW annotations
  AddParameter(_FAST_SOLVENT,m_bFastSolvent);//Controls solvent performance enhancements
#ifdef _DEBUG
  cout << _CT << " parameterised constructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtVdwIdxSF::~RbtVdwIdxSF() {
#ifdef _DEBUG
  cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}
	
void RbtVdwIdxSF::ScoreMap(RbtStringVariantMap& scoreMap) const {
  if (isEnabled()) {
//XB uncommented next line
//    EnableAnnotations(m_bAnnotate);//DM 10 Apr 2003 - only annotate if required
//    EnableAnnotations(false);//DM 17 Jan 2006 - disable for now (terse mode)
//XB uncommented next line
//    ClearAnnotationList();
    //We can only annotate the ligand-receptor interactions
    //as the rDock Viewer annotation format is hardwired to expect
    //ligand-receptor atom indices
    //Divide the total raw score into "system" and "inter" components.
    RbtDouble rs = InterScore();
//    EnableAnnotations(false);
//XB uncommented next line
    rs += LigandSolventScore();//lig-solvent belongs with the receptor-ligand inter component
    
    //First deal with the inter score which is stored in its natural location in the map
    RbtString name = GetFullName();
    scoreMap[name] = rs;
    AddToParentMapEntry(scoreMap, rs);
    
    //Now deal with the system raw scores which need to be stored in SCORE.INTER.VDW
    RbtDouble system_rs = ReceptorScore() + SolventScore() + ReceptorSolventScore();
    if (system_rs != 0.0) {
        RbtString systemName = RbtBaseSF::_SYSTEM_SF + "." + GetName();
        scoreMap[systemName] = system_rs;
        //increment the SCORE.SYSTEM total
        RbtDouble parentScore = scoreMap[RbtBaseSF::_SYSTEM_SF];
        parentScore += system_rs * GetWeight();
        scoreMap[RbtBaseSF::_SYSTEM_SF] = parentScore;
    }

//XB uncomented next 8 lines
//    scoreMap[name+".nattr"] = m_nAttr;
//    scoreMap[name+".nrep"] = m_nRep;
//    if (m_bAnnotate) {
//      RbtStringList annList;
//      RenderAnnotationsByResidue(annList);//Summarised by residue
//      scoreMap[RbtAnnotationHandler::_ANNOTATION_FIELD] += annList;
//      ClearAnnotationList();
//    }
  }
}

void RbtVdwIdxSF::SetupReceptor() {
  m_spGrid = RbtNonBondedGridPtr();
  m_recAtomList.clear();
  m_recRigidAtomList.clear();
  m_recFlexAtomList.clear();
  m_bFlexRec = false;
  m_recFlexIntns.clear();
  m_recFlexPrtIntns.clear();
  if (GetReceptor().Null())
    return;
  m_bFlexRec = GetReceptor()->isFlexible();

  m_recAtomList = GetReceptor()->GetAtomList();
  m_spGrid = CreateNonBondedGrid();
  RbtDouble maxError = GetMaxError();
  RbtDouble flexDist = 2.0;
  RbtDockingSitePtr spDS = GetWorkSpace()->GetDockingSite();
  RbtInt iTrace = GetTrace();

  RbtInt nCoords = GetReceptor()->GetNumSavedCoords()-1;
  if (nCoords > 0) {
    for (RbtInt i = 1; i <= nCoords; i++) {
      if (iTrace > 0) {
	cout << _CT << ": Indexing receptor coords # " << i << endl;
      }
      GetReceptor()->RevertCoords(i);
      RbtAtomList atomList = spDS->GetAtomList(m_recAtomList,0.0,GetCorrectedRange());
      for (RbtAtomListConstIter iter = atomList.begin(); iter != atomList.end(); iter++) {
	RbtDouble range = MaxVdwRange(*iter);
	m_spGrid->SetAtomLists(*iter,range+maxError);
      }
      m_spGrid->UniqueAtomLists();
    }
  }
  else {
    RbtAtomList atomList = spDS->GetAtomList(m_recAtomList,0.0,GetCorrectedRange());
    std::copy(atomList.begin(),atomList.end(),std::back_inserter(m_recRigidAtomList));
    //For flexible receptors, separate the site atoms into rigid and flexible
    if (m_bFlexRec) {
      GetReceptor()->SetAtomSelectionFlags(false);
      GetReceptor()->SelectFlexAtoms();//This leaves all moveable atoms selected
      Rbt::isAtomSelected isSel;
      RbtAtomRListIter fIter = std::stable_partition(m_recRigidAtomList.begin(),m_recRigidAtomList.end(),std::not1(isSel));
      std::copy(fIter,m_recRigidAtomList.end(),std::back_inserter(m_recFlexAtomList));
      m_recRigidAtomList.erase(fIter,m_recRigidAtomList.end());

      //Build map of intra-protein interactions (similar to intra-ligand)
      //Include flexible-flexible and flexible-rigid
      RbtInt nAtoms = GetReceptor()->GetNumAtoms();
      m_recFlexIntns = RbtAtomRListList(nAtoms,RbtAtomRList());
      m_recFlexPrtIntns = RbtAtomRListList(nAtoms,RbtAtomRList());
      BuildIntraMap(m_recFlexAtomList,m_recFlexIntns);//Flexible-flexible
      BuildIntraMap(m_recFlexAtomList,m_recRigidAtomList,m_recFlexIntns);//flexible-rigid
      //We can get away with partitioning the variable interactions just at the beginning.
      //For grosser receptor flexibility we would have to partition periodically during docking
      RbtDouble partitionDist = MaxVdwRange(RbtTriposAtomType::H_P) + (2.0 * flexDist);
      Partition(m_recFlexAtomList, m_recFlexIntns, m_recFlexPrtIntns, partitionDist);

      //Index the flexible atoms over a larger radius
      //NOTE: WE ASSUME ONLY -OH and -NH3 rotation here (protons can't move more than 2.0A at most)
      //Grosser rotations will require a different approach
      for (RbtAtomRListConstIter iter = m_recFlexAtomList.begin(); iter != m_recFlexAtomList.end(); iter++) {
	RbtDouble range = MaxVdwRange(*iter);
	m_spGrid->SetAtomLists(*iter,range+maxError+flexDist);
      }
      if (iTrace > 0) {
	RbtDouble score = ReceptorScore();
	cout << GetWorkSpace()->GetName() << " " << GetFullName() << ": Intra-receptor score = " << score << endl;
      }
    }
    //Index the rigid atoms as usual
    for (RbtAtomRListConstIter iter = m_recRigidAtomList.begin(); iter != m_recRigidAtomList.end(); iter++) {
      RbtDouble range = MaxVdwRange(*iter);
      m_spGrid->SetAtomLists(*iter,range+maxError);
    }
  }
}

void RbtVdwIdxSF::SetupLigand() {
  m_ligAtomList.clear();
  if (GetLigand().Null())
    return;
  
  RbtAtomList tmpList = GetLigand()->GetAtomList();
  //Strip off the smart pointers
  std::copy(tmpList.begin(),tmpList.end(),std::back_inserter(m_ligAtomList));
}

// DM 13 June 2006 - performance enhancements to take account of fixed/tethered/free solvent
// Approach:
// 1) Divide the solvent atoms into two lists
//    a) fixed (no displacement) or tethered (predictable (small) displacements)
//    b) free (unpredicatable, large displacements)
// 2) Index the fixed/tethered atoms on an indexing grid, as is done for the receptor
// 3) Build an interaction map for (fixed/tethered - fixed/tethered) solvent intns.
//    This can be partitioned in advance based on the maximum displacement of any of the atoms
// 4) Build an interaction map for (free - free) solvent intns. This can not be partitioned
//    due to the large allowed displacements.
// 5) Score calculations are as follows:
//    a) Receptor - solvent. Use the receptor indexing grid (as used for receptor - ligand).
//    b) Solvent (fix/teth) - solvent (fix/teth). Partitioned intn map.
//       We can not use the solvent indexing grid here, due to self interactions.
//    c) Solvent (free) - solvent (fix/teth). Solvent indexing grid.
//    d) Solvent (free) - solvent (free). Unpartitioned intn map.
//    e) Ligand - solvent(fix/teth). Solvent indexing grid.
//    f) Ligand - solvent(free). Brute force loop.
// Note: 5d and 5f are still inefficient, but free solvent is unlikely to be used very often
// in practice.
void RbtVdwIdxSF::SetupSolvent() {
  m_spSolventGrid = RbtNonBondedGridPtr();
  m_solventAtomList.clear();
  m_solventFixTethAtomList.clear();
  m_solventFreeAtomList.clear();
  m_solventFixTethIntns.clear();
  m_solventFixTethPrtIntns.clear();
  m_solventFreeIntns.clear();
  RbtModelList solventList = GetSolvent();
  if (solventList.empty()) {
    return;
  }
  for (RbtModelListConstIter iter = solventList.begin(); iter != solventList.end(); ++iter) {
    RbtAtomList atomList = (*iter)->GetAtomList();
    std::copy(atomList.begin(),atomList.end(),std::back_inserter(m_solventAtomList));
    //For faster score calculations, divide the solvent atoms into fixed/tethered and free atom lists
    if (m_bFastSolvent) {
      RbtFlexAtomFactory flexAtomFactory(*iter);
      RbtAtomRList fixedAtomList = flexAtomFactory.GetFixedAtomList();
      RbtAtomRList tetheredAtomList = flexAtomFactory.GetTetheredAtomList();
      RbtAtomRList freeAtomList = flexAtomFactory.GetFreeAtomList();
      std::copy(fixedAtomList.begin(), fixedAtomList.end(), std::back_inserter(m_solventFixTethAtomList));
      std::copy(tetheredAtomList.begin(), tetheredAtomList.end(), std::back_inserter(m_solventFixTethAtomList));
      std::copy(freeAtomList.begin(), freeAtomList.end(), std::back_inserter(m_solventFreeAtomList));
    }
    //The slower alternative (for testing comparisons) is to treat all solvent atoms as free
    //i.e. make no assumptions as to their coordinates
    else {
      std::copy(atomList.begin(),atomList.end(),std::back_inserter(m_solventFreeAtomList));
    }
  }

  //Index the fixed/tethered atoms on the solvent indexing grid.
  //The flexAtomFactory has cached the maximum translations possible for each atom in the User2Value attribute
  //Also build the interaction map for fixed/tethered - fixed/tethered intns
  //The interaction map can be partitioned in advance, based on the maximum displacement
  //of any of the tethered atoms
  if (!m_solventFixTethAtomList.empty()) {
    m_spSolventGrid = CreateNonBondedGrid();
    RbtDouble maxError = GetMaxError();
    RbtDouble maxFlexDist = 0.0;
    for (RbtAtomRListConstIter iter = m_solventFixTethAtomList.begin(); iter != m_solventFixTethAtomList.end(); iter++) {
      RbtDouble range = MaxVdwRange(*iter);
      RbtDouble flexDist = (*iter)->GetUser2Value();
      maxFlexDist = std::max(flexDist, maxFlexDist);
      m_spSolventGrid->SetAtomLists(*iter,range+maxError+flexDist);
    }
    m_solventFixTethIntns = RbtAtomRListList(m_solventAtomList.size(),RbtAtomRList());
    m_solventFixTethPrtIntns = RbtAtomRListList(m_solventAtomList.size(),RbtAtomRList());
    BuildIntraMap(m_solventFixTethAtomList,m_solventFixTethIntns);
    RbtDouble partitionDist = MaxVdwRange(RbtTriposAtomType::O_3) + (2.0 * maxFlexDist);
    Partition(m_solventFixTethAtomList, m_solventFixTethIntns, m_solventFixTethPrtIntns, partitionDist);
    if (GetTrace() > 0) {
      cout << endl << "Faster calculation of vdW scores involving fixed/tethered solvent is enabled..." << endl;
      cout << "#Fixed/tethered solvent atoms = " << m_solventFixTethAtomList.size() << endl;
      cout << "Max translation of any fixed/tethered solvent atom = " << maxFlexDist << "A" << endl;
    }
  }
  //Build the interaction map for the free solvent - free solvent interactions only
  if (!m_solventFreeAtomList.empty()) {
    m_solventFreeIntns = RbtAtomRListList(m_solventAtomList.size(),RbtAtomRList());
    BuildIntraMap(m_solventFreeAtomList,m_solventFreeIntns);
    if (GetTrace() > 0) {
      if (m_bFastSolvent) {
	cout << endl << "Calculation of vdW scores involving freely translating solvent can not be optimised..." << endl;
	cout << "#Free solvent atoms = " << m_solventFreeAtomList.size() << endl;
      }
      else {
	cout << endl << "Faster calculation of vdW scores involving fixed/tethered solvent is disabled..." << endl;
      }
    }
  }
}

void RbtVdwIdxSF::SetupScore() {
  //No further setup required
}



RbtDouble RbtVdwIdxSF::RawScore() const {
  return InterScore() + LigandSolventScore() + ReceptorScore() + SolventScore() +
  		ReceptorSolventScore();
}

//DM 25 Oct 2000 - track changes to parameter values in local data members
//ParameterUpdated is invoked by RbtParamHandler::SetParameter
void RbtVdwIdxSF::ParameterUpdated(const RbtString& strName) {
  if (strName == _THRESHOLD_ATTR) {
    m_attrThreshold = GetParameter(_THRESHOLD_ATTR);
  }
  else if (strName == _THRESHOLD_REP) {
    m_repThreshold = GetParameter(_THRESHOLD_REP);
  }
  else if (strName == _ANNOTATION_LIPO) {
    m_lipoAnnot = GetParameter(_ANNOTATION_LIPO);
  }
  else if (strName == _ANNOTATE) {
    m_bAnnotate = GetParameter(_ANNOTATE);
  }
  else if (strName == _FAST_SOLVENT) {
    m_bFastSolvent = GetParameter(_FAST_SOLVENT);
  }
  else {
    RbtVdwSF::OwnParameterUpdated(strName);
    RbtBaseIdxSF::OwnParameterUpdated(strName);
    RbtBaseSF::ParameterUpdated(strName);
  }
}

//DM 06 Feb 2003
//This method processes the raw vdw annotation list (all non-zero pair-wise scores)
//and outputs three subsets
//VDW_SUM - sum of all vdW scores by residue (vector is between amino acid CA and nearest ligand atom)
//VDW_REP - all repulsive (+ve) pair-wise vdW interactions
//VDW_LIPO - all attractive (-ve) pair-wise vdW interactions between apolar C/H atoms
//           with scores better (more -ve) than ANNOTATION_LIPO parameter
void RbtVdwIdxSF::RenderAnnotationsByResidue(RbtStringList& retVal) const {
  RbtString strSum = GetName()+"_SUM";
  RbtString strLipo = GetName()+"_LIPO";
  RbtString strRepul = GetName()+"_REP";

  RbtAnnotationList annList = GetAnnotationList();
  std::sort(annList.begin(),annList.end(),Rbt::RbtAnn_Cmp_AtomId2());

  RbtString oldResName("UNLIKELY_TO_MATCH_BY_CHANCE");
  RbtAnnotationPtr spAnn;
 
  for (RbtAnnotationListConstIter aIter = annList.begin(); aIter != annList.end(); ++aIter) {
    RbtString resName = (*aIter)->GetFQResName();
    //New residue encountered
    if (resName != oldResName) {
      //Render the previous residue annotation (unless this is the first time through)
      if (spAnn.Ptr() != NULL) {
	//For residue summaries, we copy the score into the distance attribute. The viewer 2.0 only displays distance labels
	//so we fool it here into displaying the vdw summary scores (more useful)
	spAnn->SetDistance(spAnn->GetScore());
	retVal.push_back(strSum + "," + spAnn->Render());
      }
      
      oldResName = resName;     
      //Find all the atoms in this residue, so that we can locate the central atom for annotation purposes
      //Amino acid = CA; nucleic acid = C1'; other residues = first non-hydrogen atom (or first atom)
      RbtAtomList resAtoms = Rbt::GetMatchingAtomList(m_recAtomList,resName);
      //If assertion fails, implies that atom 2 is not in the receptor
      Assert<RbtAssert>(!resAtoms.empty());
      RbtAtomList centralAtoms = Rbt::GetMatchingAtomList(resAtoms,"CA");
      RbtAtom* centralAtom;
      if (!centralAtoms.empty()) {
	centralAtom = centralAtoms.front();//Amino acid
      }
      else {
	centralAtoms = Rbt::GetMatchingAtomList(resAtoms,"C1'");
	if (!centralAtoms.empty()) {
	  centralAtom = centralAtoms.front();//Nucleic acid
	}
	else {
	  centralAtoms = Rbt::GetAtomList(resAtoms,std::not1(Rbt::isAtomicNo_eq(1)));
	  centralAtom = (!centralAtoms.empty()) ? centralAtoms.front() : resAtoms.front();//other residue (e.g. solvent)
	}
      }
      //Create a new annotation, with atom 2 = central atom
      spAnn = new RbtAnnotation((*aIter)->GetAtom1Ptr(),centralAtom,(*aIter)->GetDistance(),(*aIter)->GetScore());
    }
    
    //Previously encountered residue
    else {
      //Accumulate the annotation
      *(spAnn) += (**aIter);
    }

    //Output the raw atom-atom annotation if
    //a) it is repulsive (score > 0)
    //b) it is attractive and between two lipo atoms (C,H)
    RbtDouble s((*aIter)->GetScore());
    if (s > 0.0) {
      retVal.push_back(strRepul + "," + (*aIter)->Render());
    }
    else if ((s < m_lipoAnnot) && (*aIter)->GetAtom1Ptr()->GetUser1Flag() &&
	     (*aIter)->GetAtom2Ptr()->GetUser1Flag()) {
      retVal.push_back(strLipo + "," + (*aIter)->Render());
    }
  }
  
  //Render the final annotation
  if (spAnn.Ptr() != NULL) {
    //For residue summaries, we copy the score into the distance attribute. The viewer 2.0 only displays distance labels
    //so we fool it here into displaying the vdw summary scores (more useful)
    spAnn->SetDistance(spAnn->GetScore());
    retVal.push_back(strSum + "," + spAnn->Render());
  }
}

RbtDouble RbtVdwIdxSF::InterScore() const {
  RbtDouble score = 0.0;
  m_nAttr = 0;
  m_nRep = 0;

  //Check grid is defined
  if (m_spGrid.Null())
    return score;

  //Loop over all ligand atoms
  for (RbtAtomRListConstIter iter = m_ligAtomList.begin(); iter != m_ligAtomList.end(); iter++) {
    const RbtCoord& c = (*iter)->GetCoords();
    const RbtAtomRList& recepAtomList = m_spGrid->GetAtomList(c);
    RbtDouble s = VdwScore(*iter,recepAtomList);
    score += s;
    if (s > m_repThreshold) {
      m_nRep++;
    }
    else if (s < m_attrThreshold) {
      m_nAttr++;
    }
  }
  return score;
}

//Intra-receptor
RbtDouble RbtVdwIdxSF::ReceptorScore() const {
  if (!m_bFlexRec) return 0.0;
  RbtDouble score = 0.0;//Total score
  //Loop over all flexible site atoms
  for (RbtAtomRListConstIter iter = m_recFlexAtomList.begin(); iter != m_recFlexAtomList.end(); iter++) {
    RbtInt id = (*iter)->GetAtomId()-1;
//XB changed call from "VdwScore" to "VdwScoreIntra" and created new function
// in "RbtVdwSF.cxx" to avoid using reweighting terms for intra
    //RbtDouble s = VdwScoreIntra(*iter,m_recFlexPrtIntns[id]);
    RbtDouble s = VdwScore(*iter,m_recFlexPrtIntns[id]);
    score += s;
  }
  return score;
}

//Intra-solvent
RbtDouble RbtVdwIdxSF::SolventScore() const {
  RbtDouble score = 0.0;
  //Use the partitioned intn map for fixed/tethered - fixed/tethered intns
  for (RbtAtomRListConstIter iter = m_solventFixTethAtomList.begin(); iter != m_solventFixTethAtomList.end(); iter++) {
    RbtInt id = (*iter)->GetAtomId()-1;
    score += VdwScoreEnabledOnly(*iter,m_solventFixTethPrtIntns[id]);
  }
  if (!m_spSolventGrid.Null()) {
    //Use the indexing grid for free - fixed/tethered intns
    for (RbtAtomRListConstIter iter = m_solventFreeAtomList.begin(); iter != m_solventFreeAtomList.end(); iter++) {
      const RbtCoord& c = (*iter)->GetCoords();
      const RbtAtomRList& atomList = m_spSolventGrid->GetAtomList(c);
      score += VdwScoreEnabledOnly(*iter,atomList);
    }
  }
  //Use the intn map for free - free intns (this is still inefficient, but unlikely to be used in practice)
  for (RbtAtomRListConstIter iter = m_solventFreeAtomList.begin(); iter != m_solventFreeAtomList.end(); iter++) {
    RbtInt id = (*iter)->GetAtomId()-1;
    score += VdwScoreEnabledOnly(*iter,m_solventFreeIntns[id]);
  }
  return score;
}

//Receptor-solvent
RbtDouble RbtVdwIdxSF::ReceptorSolventScore() const {
  RbtDouble score = 0.0;
  if (m_spGrid.Null())
    return score;
  for (RbtAtomRListConstIter iter = m_solventAtomList.begin(); iter != m_solventAtomList.end(); iter++) {
    //DM 7 June 2006 - take into account the enabled state of each solvent atom
    if ((*iter)->GetEnabled()) {
      const RbtCoord& c = (*iter)->GetCoords();
      const RbtAtomRList& recepAtomList = m_spGrid->GetAtomList(c);
//XB changed call from "VdwScore" to "VdwScoreIntra" and created new function
// in "RbtVdwSF.cxx" to avoid using reweighting terms for intra
      //score += VdwScoreIntra(*iter,recepAtomList);
      score += VdwScore(*iter,recepAtomList);
    }
  }
  return score;
}

//Ligand-solvent
RbtDouble RbtVdwIdxSF::LigandSolventScore() const {
  RbtDouble score = 0.0;
  //Use the solvent indexing grid for ligand - fixed/tethered solvent intns
  if (!m_spSolventGrid.Null()) {
    for (RbtAtomRListConstIter iter = m_ligAtomList.begin(); iter != m_ligAtomList.end(); iter++) {
      const RbtCoord& c = (*iter)->GetCoords();
      const RbtAtomRList& atomList = m_spSolventGrid->GetAtomList(c);
      score += VdwScoreEnabledOnly(*iter,atomList);
    }
  }
  //Use inefficient brute force for ligand - free solvent intns
  for (RbtAtomRListConstIter iter = m_solventFreeAtomList.begin(); iter != m_solventFreeAtomList.end(); iter++) {
    //DM 7 June 2006 - take into account the enabled state of each solvent atom
    if ((*iter)->GetEnabled()) {
      score += VdwScore(*iter,m_ligAtomList);
    }
  }
  return score;
}
