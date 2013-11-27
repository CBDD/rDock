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

#include "RbtSAIdxSF.h"
#include "RbtWorkSpace.h"
#include "RbtSFRequest.h"

RbtString RbtSAIdxSF::_CT ("RbtSAIdxSF");
RbtString RbtSAIdxSF::_INCR ("INCR");

RbtSAIdxSF::RbtSAIdxSF(const RbtString& aName) :
  RbtBaseSF(_CT,aName), m_maxR(2.0), m_bFlexRec(false), m_lig_0(0.0), m_lig_free(0.0), m_lig_bound(0.0),
  m_site_0(0.0), m_site_free(0.0), m_site_bound(0.0),
  m_solvent_0(0.0), m_solvent_free(0.0), m_solvent_bound(0.0)
{
  //INCR = increment to be added to radius of each atom for indexing on the near-neighbour grid
  //Used to calculate maximum range of scoring function for each atom
  //Will be adjusted dynamically in Setup, based on max radius of any atom type
  //r_s = solvent probe radius (constant 0.6)
  AddParameter(_INCR,m_maxR + 2*HHS_Solvation::r_s);
  m_spSolvSource = RbtParameterFileSourcePtr(new RbtParameterFileSource(Rbt::GetRbtFileName("data/sf","solvation_asp.prm")));
  Setup();
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtSAIdxSF::~RbtSAIdxSF()
{
  ClearReceptor();
  ClearLigand();
  ClearSolvent();
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtSAIdxSF::SetupReceptor() { 
  ClearReceptor();
  if (GetReceptor().Null())
    return;
  RbtInt iTrace = GetTrace();

  //Trap multiple receptor conformations here: this SF does not support them yet
  RbtBool bEnsemble = (GetReceptor()->GetNumSavedCoords() > 1);
  if (bEnsemble) {
    throw RbtInvalidRequest(_WHERE_,"Solvation scoring function does not support multiple receptor conformations yet");
  }

  //MAKE THE ASSUMPTION that the only flexible receptor atoms are terminal OH/NH3
  //and that they don't move very far (up to 2A)
  m_bFlexRec = GetReceptor()->isFlexible();
  RbtDouble flexDist = 2.0;
  theIdxGrid = CreateNonBondedHHSGrid();
  RbtDouble idxIncr = GetParameter(_INCR).Double() + GetMaxError();
  RbtDouble flexIncr = idxIncr + flexDist;

  // At this stage we deal with all receptor atoms, to avoid edge effects when determining
  // surface areas for atoms on the periphery of the docking site
  // At the end we filter the interaction center list to just those receptor atoms
  // in range of the docking site
  RbtAtomList theReceptorList = GetReceptor()->GetAtomList();
  theRSPList = CreateInteractionCenters(theReceptorList);

  //For flexible receptors, separate the interaction centers into rigid and flexible
  //When we build up the intra-protein variable distances (BuildIntraMap), ensure that
  //the variable interactions are stored on the (relatively small number) of flexible interaction centers
  //rather than being scattered around the entire receptor list
  //This is more efficient when it comes to RawScore, as fewer interaction centers
  //have to have their variable distances updated
  if (m_bFlexRec) {
    GetReceptor()->SetAtomSelectionFlags(false);
    GetReceptor()->SelectFlexAtoms();//This leaves all moveable atoms selected
    Rbt::isHHSSelected isSel;
    HHS_SolvationRListIter fIter = std::stable_partition(theRSPList.begin(),theRSPList.end(),std::not1(isSel));
    std::copy(fIter,theRSPList.end(),std::back_inserter(theFlexList));
    theRSPList.erase(fIter,theRSPList.end());
    BuildIntraMap(theFlexList);//flexible-flexible
    BuildIntraMap(theFlexList,theRSPList);//flexible-rigid
    //Store the per-atom invariant free areas for later retrieval
    Rbt::SaveHHS saveInvariantArea;
    std::for_each(theFlexList.begin(),theFlexList.end(),saveInvariantArea);
    // Index the flexible interaction centers within range of the docking site
    // Use a larger increment
    for(HHS_SolvationRListConstIter iter = theFlexList.begin(); iter != theFlexList.end(); iter++) {
      RbtAtom* pAtom = (*iter)->GetAtom();
      theIdxGrid->SetHHSLists(*iter,(*iter)->GetR_i()+flexIncr);
    }
  }

  BuildIntraMap(theRSPList);//rigid-rigid
  //Store the per-atom invariant free areas for later retrieval
  Rbt::SaveHHS saveInvariantArea;
  std::for_each(theRSPList.begin(),theRSPList.end(),saveInvariantArea);

  if (m_bFlexRec) {
    //Do a one-shot partitioning of the variable distances
    //For grosser flexibility than OH/NH3 rotation, would have to partition more frequently
    //during docking
    RbtDouble dist = GetR_i(RbtHHSType::HNp) + GetParameter(_INCR).Double() + flexDist;
    Partition(theFlexList,dist);
    Rbt::OverlapVariableHHS updateVariableArea;
    std::for_each(theFlexList.begin(),theFlexList.end(),updateVariableArea);
    //This nasty piece of code sets the selection flag to true for all atoms
    //that are on the receiving end of a stored variable interaction
    //Use this flag to detect the peripheral rigid atoms (beyond the range of the docking site)
    //that are in range of one of the flexible atoms
    GetReceptor()->SetAtomSelectionFlags(false);
    for (HHS_SolvationRListConstIter iIter = theFlexList.begin(); iIter != theFlexList.end(); ++iIter) {
      const HHS_SolvationRList& varList = (*iIter)->GetVariable();
      for (HHS_SolvationRListConstIter jIter = varList.begin(); jIter != varList.end(); ++jIter) {
	(*jIter)->GetAtom()->SetSelectionFlag(true);
      }
    }
  }
      
  // Identify the rigid atoms within range of the docking site and manage as a separate list
  // These should be the only atoms that can be buried by docked ligands
  RbtDockingSite::isAtomInRange isNearCavity(GetWorkSpace()->GetDockingSite()->GetGrid(),0.0,GetCorrectedRange());
  Rbt::isAtomSelected isSelected;
  for(HHS_SolvationRListConstIter iter = theRSPList.begin(); iter != theRSPList.end(); iter++) {
    RbtAtom* pAtom = (*iter)->GetAtom();
    if (isNearCavity(pAtom)) {
      theCavList.push_back(*iter);
    }
    else if (m_bFlexRec && isSelected(pAtom)) {
      thePeriphList.push_back(*iter);
      if (iTrace > 1) {
	cout << _CT << ": Peripheral rigid atom within range of flexible atoms: " << pAtom->GetFullAtomName() << endl;
      }
    }
  }

  // Index the rigid interaction centers within range of the docking site
  for(HHS_SolvationRListConstIter iter = theCavList.begin(); iter != theCavList.end(); iter++) {
    RbtAtom* pAtom = (*iter)->GetAtom();
    theIdxGrid->SetHHSLists(*iter,(*iter)->GetR_i()+idxIncr);
  }

  //Initial solvation free energy (rigid and flexible atom contributions)
  m_site_0 = TotalEnergy(theCavList);
  if (m_bFlexRec) {
    m_site_0 += TotalEnergy(theFlexList);
  }

  if (iTrace > 1) {
    cout << endl << _CT << ": Rigid receptor atoms within range of docking site" << endl << endl;
    cout << _CT << ": ATOM, TYPE, ASP, SOLV" << endl;
    RbtHHSType hhsType;
    for(HHS_SolvationRListConstIter iter = theCavList.begin(); iter != theCavList.end(); iter++) {
      RbtHHSType::eType t = (*iter)->GetHHSType();
      RbtAtom* pAtom = (*iter)->GetAtom();
      RbtDouble sasa = (*iter)->GetArea();
      RbtDouble asp = (*iter)->GetSigma();
      RbtDouble energy =(*iter)->GetEnergy();
      cout << _CT << ": " << pAtom->GetFullAtomName() << ", " << hhsType.Type2Str(t) << ", "
	   << sasa << ", " << asp << ", " << energy << endl;
    }
    if (m_bFlexRec) {
      cout << endl << _CT << ": Flexible receptor atoms within range of docking site" << endl << endl;
      cout << _CT << ": ATOM, TYPE, ASP, SOLV" << endl;
      RbtHHSType hhsType;
      for(HHS_SolvationRListConstIter iter = theFlexList.begin(); iter != theFlexList.end(); iter++) {
	RbtHHSType::eType t = (*iter)->GetHHSType();
	RbtAtom* pAtom = (*iter)->GetAtom();
	RbtDouble sasa = (*iter)->GetArea();
	RbtDouble asp = (*iter)->GetSigma();
	RbtDouble energy =(*iter)->GetEnergy();
	cout << _CT << ": " << pAtom->GetFullAtomName() << ", " << hhsType.Type2Str(t) << ", "
	     << sasa << ", " << asp << ", " << energy << endl;
      }
    }
  }
  
  if (iTrace > 0) {
    cout << endl;
    cout << "Initial site dG(solv)        : " << m_site_0 << " kcal/mol" << endl;
    //Count the number of undefined types
    Rbt::isHHSType_eq isUndefined(RbtHHSType::UNDEFINED);
    RbtInt nUndefRec = std::count_if(theRSPList.begin(),theRSPList.end(),isUndefined);
    RbtInt nUndefSite = std::count_if(theCavList.begin(),theCavList.end(),isUndefined);

    cout << "#UNDEFINED TYPES (rigid recep)  : " << nUndefRec << endl;
    cout << "#UNDEFINED TYPES (rigid site)   : " << nUndefSite << endl;
    if (m_bFlexRec) {
      RbtInt nUndefFlex = std::count_if(theFlexList.begin(),theFlexList.end(),isUndefined);
      cout << "#UNDEFINED TYPES (flex site)    : " << nUndefFlex << endl;
    }
  }
}
/////////////////////////////////////////////////////////////////
//
// SetupLigand:
//
/////////////////////////////////////////////////////////////////
void RbtSAIdxSF::SetupLigand()
{
  ClearLigand();
  RbtModelPtr spModel = GetLigand();
  if (spModel.Null())
    return;
  RbtAtomList theLigandList = spModel->GetAtomList();
  theLSPList = CreateInteractionCenters(theLigandList);
  BuildIntraMap(theLSPList);
  //Store the per-atom invariant free areas for later retrieval
  Rbt::SaveHHS saveInvariantArea;
  std::for_each(theLSPList.begin(),theLSPList.end(),saveInvariantArea);
  RbtInt iTrace = GetTrace();
  //The solvation term also describes the intramolecular solvation energy of the ligand
  //so, as with the intramolecular terms, we want to take the "zero" point energy from
  //the initial (Corina) conformation of the ligand. Read the zero (free ligand) value
  //from the input model if it exists, otherwise calculate here
  RbtString name = RbtBaseSF::_INTRA_SF + "." + GetName() + ".lig_0";
  if (spModel->isDataFieldPresent(name)) {
    if (iTrace > 0) {
      cout << "Restoring initial ligand solvation energy from file..." << endl;
    }
    m_lig_0 = spModel->GetDataValue(name);
  }
  else {
    if (iTrace > 0) {
      cout << "Calculating initial ligand solvation energy..." << endl;
    }
    //Variable distances
    Rbt::OverlapVariableHHS updateVariableArea;
    std::for_each(theLSPList.begin(),theLSPList.end(),updateVariableArea);
    //Total solv energy for initial ligand conf
    m_lig_0 = TotalEnergy(theLSPList);
    if (iTrace > 2) {
      PrintWeightMatrix();
    }
  }
  
  if (iTrace > 0) {
    cout << endl;
    cout << "Initial ligand dG(solv):     " << m_lig_0 << " kcal/mol" << endl;
    //Count the number of undefined types
    Rbt::isHHSType_eq isUndefined(RbtHHSType::UNDEFINED);
    RbtInt nUndef = std::count_if(theLSPList.begin(),theLSPList.end(),isUndefined);
    cout << "#UNDEFINED TYPES (lig)    : " << nUndef << endl;
  }
}

void RbtSAIdxSF::SetupSolvent()
{
  ClearSolvent();
  RbtModelList solventModelList = GetSolvent();
  RbtInt iTrace = GetTrace();
  //Process the solvent models individually in order to build up each intra-solvent interaction map
  //correctly.
  for (RbtModelListConstIter iter = solventModelList.begin(); iter != solventModelList.end(); ++iter) {
    HHS_SolvationRList intnList = CreateInteractionCenters((*iter)->GetAtomList());
    BuildIntraMap(intnList);
    //Store the per-atom invariant free areas for later retrieval
    std::for_each(intnList.begin(),intnList.end(),Rbt::SaveHHS());
    //Variable distances to previous solvent models already processed
    //These are strictly intermolecular distances, but as we consider the total set of solvent models
    //as a single entity for the purposes of partitioning the total score, it makes to store them here.
    for (HHS_SolvationRListConstIter iIter = theSolventList.begin(); iIter != theSolventList.end(); ++iIter) {
        for (HHS_SolvationRListConstIter jIter = intnList.begin(); jIter != intnList.end(); ++jIter) {
            (*iIter)->AddVariable(*jIter);
        }
    }
    //Concatenate all the interaction centers from each solvent model into a single list
    std::copy(intnList.begin(),intnList.end(),std::back_inserter(theSolventList));
  }
  //Calculate the initial solvation score for the entire set of solvent models
  //DM 9 June 2006 - don't include the intra-solvent interactions in the zero-point calculation
  //as we don't know whether the individual solvent models will be enabled or not.
  //This gives a natural penalty desolvation score to all enabled solvent models in the
  //docking site, which may reduce the need for the additional solvent penalty term in RbtConstSF.
  //The zero-point score now just represents the desolvation score for isolated solvent models.
  //std::for_each(theSolventList.begin(),theSolventList.end(),Rbt::OverlapVariableHHS());
  m_solvent_0 = TotalEnergy(theSolventList);
  if (iTrace > 0) {
    cout << endl;
    cout << "Initial solvent dG(solv):     " << m_solvent_0 << " kcal/mol" << endl;
    //Count the number of undefined types
    Rbt::isHHSType_eq isUndefined(RbtHHSType::UNDEFINED);
    RbtInt nUndef = std::count_if(theSolventList.begin(),theSolventList.end(),isUndefined);
    cout << "#UNDEFINED TYPES (sol)    : " << nUndef << endl;
  }
}

void RbtSAIdxSF::SetupScore()
{
}

RbtDouble RbtSAIdxSF::RawScore(void) const
{
  // restore invariant surface areas for ligand, solvent and rigid receptor
  for (HHS_SolvationRListConstIter iter = theLSPList.begin(); iter != theLSPList.end(); ++iter)
    (*iter)->Restore();
  for (HHS_SolvationRListConstIter iter = theCavList.begin(); iter != theCavList.end(); ++iter)
    (*iter)->Restore();
  for (HHS_SolvationRListConstIter iter = theSolventList.begin(); iter != theSolventList.end(); ++iter)
    (*iter)->Restore();

  //INTRA-LIGAND interactions
  for (HHS_SolvationRListConstIter iter = theLSPList.begin(); iter != theLSPList.end(); ++iter) {
    (*iter)->OverlapVariable();
  }

  //INTRA-SOLVENT interactions (take account of solvent enabled state)
  for (HHS_SolvationRListConstIter iter = theSolventList.begin(); iter != theSolventList.end(); ++iter) {
    (*iter)->OverlapVariableEnabledOnly();
  }
  
  //INTRA-SITE interactions (if receptor is flexible)
  if (m_bFlexRec) {
    for (HHS_SolvationRListConstIter iter = thePeriphList.begin(); iter != thePeriphList.end(); ++iter)
      (*iter)->Restore();
    for (HHS_SolvationRListConstIter iter = theFlexList.begin(); iter != theFlexList.end(); ++iter)
      (*iter)->Restore();
    for (HHS_SolvationRListConstIter iter = theFlexList.begin(); iter != theFlexList.end(); ++iter)
      (*iter)->OverlapVariable();
  }

  //SITE-SOLVENT interactions (take account of solvent enabled state)
  for(HHS_SolvationRListConstIter iIter = theSolventList.begin(); iIter != theSolventList.end(); iIter++) {
    RbtAtom* pSolventAtom = (*iIter)->GetAtom();
    if (pSolventAtom->GetEnabled()) {
      const RbtCoord& rAtomCoords     = pSolventAtom->GetCoords();
      const HHS_SolvationRList& rList = theIdxGrid->GetHHSList(rAtomCoords);
      for (HHS_SolvationRListConstIter jIter = rList.begin(); jIter != rList.end(); ++jIter)
	(*iIter)->Overlap(*jIter,HHS_Solvation::Pij_14);
    }
  }

  //Optional - record the "free" desolvation energies here, for use by ScoreMap
  //These include Intra-ligand, Intra-solvent, Intra-site and Site-solvent interactions
  //i.e. everything except ligand-site and ligand-solvent
  //We use the AnnotationHandler just to have a boolean state we can enable
  //No annotations are recorded.
  if (isAnnotationEnabled()) {
    m_lig_free = TotalEnergy(theLSPList);  
    m_solvent_free = TotalEnergy(theSolventList);
    m_site_free = TotalEnergy(theCavList)+TotalEnergy(theFlexList);
  }
  else {
    //Clear the intermediate scores just to avoid keeping stale values around 
    m_lig_free = 0.0;
    m_solvent_free = 0.0;
    m_site_free = 0.0;
  }

  //Now the intermolecular binding interactions themselves
  //
  //LIGAND-SITE
  //Retrieve the receptor near-neighbours from the indexing grid
  for(HHS_SolvationRListConstIter iIter = theLSPList.begin(); iIter != theLSPList.end(); iIter++) {
    const RbtCoord& rAtomCoords		= (*iIter)->GetAtom()->GetCoords();
    const HHS_SolvationRList& rList	= theIdxGrid->GetHHSList(rAtomCoords);
    for (HHS_SolvationRListConstIter jIter = rList.begin(); jIter != rList.end(); ++jIter)
      (*iIter)->Overlap(*jIter,HHS_Solvation::Pij_14);
  }
  
  
  //LIGAND-SOLVENT (VERY SLOW) (take account of solvent enabled state)
  //TODO: index the solvent intn centers on a grid, providing solvent position is tethered
  for(HHS_SolvationRListConstIter iIter = theLSPList.begin(); iIter != theLSPList.end(); iIter++) {
    for (HHS_SolvationRListConstIter jIter = theSolventList.begin(); jIter != theSolventList.end(); ++jIter) {
      RbtAtom* pSolventAtom = (*jIter)->GetAtom();
      if (pSolventAtom->GetEnabled()) {
	(*iIter)->Overlap(*jIter,HHS_Solvation::Pij_14);
      }
    }
  }
  
  //DM 8 June 2006 - record the absolute scores here to make ScoreMap analysis easier
  //(previously were relative to _0 scores)

  //Total solv energy for current ligand conf (bound)
  m_lig_bound = TotalEnergy(theLSPList);

  //Total solv energy for site (bound
  m_site_bound = TotalEnergy(theCavList)+TotalEnergy(theFlexList);
  
  //Total solv energy for solvent (bound)
  m_solvent_bound = TotalEnergy(theSolventList);

  //Total score is overall change in solvation energy for ligand, site and solvent
  //relative to initial unbound scores
  RbtDouble theScore = (m_lig_bound - m_lig_0) + (m_site_bound - m_site_0)
    + (m_solvent_bound - m_solvent_0);
  return theScore;
}

RbtDouble RbtSAIdxSF::GetP_i(RbtHHSType::eType theType) const {
  return m_solvTable[theType].p;
}

RbtDouble RbtSAIdxSF::GetR_i(RbtHHSType::eType theType) const {
  return m_solvTable[theType].r;
}

RbtDouble RbtSAIdxSF::GetASP(RbtHHSType::eType theType, RbtDouble chg) const {
  RbtDouble asp = m_solvTable[theType].asp;
  return m_solvTable[theType].chg_scaling ? asp*fabs(chg) : asp;
}

// helper function for parameter estimation
void RbtSAIdxSF::PrintWeightMatrix(void)  const
{
  RbtIntList count;		// count of atoms contributing the area
  RbtDoubleList	area;		// total SASA by a single atom type
  
  RbtHHSType hhsType;
  cout << "TYPE";
  // go through on each atom type
  for(int i = 0 ; i < RbtHHSType::MAXTYPES; ++i) {
    count.push_back(0);
    area.push_back(0.0);
    cout << "," << hhsType.Type2Str(RbtHHSType::eType(i));
  }
  cout << endl;
  for(HHS_SolvationRListConstIter iter = theLSPList.begin(); iter != theLSPList.end(); iter++) {
    RbtHHSType::eType hhsType = (*iter)->GetHHSType();
    count[hhsType]++;
    RbtDouble sa = (*iter)->GetArea();
    if (m_solvTable[hhsType].chg_scaling) {
      RbtDouble chg = fabs((*iter)->GetAtom()->GetGroupCharge());
      if (chg > 0.0)
	sa *= chg;
    }
    area[hhsType] += sa;
  }
  cout << "CT"; 	// stands for Weight Matrix (and helps grepping)
  for(int i = 0 ; i < RbtHHSType::MAXTYPES; ++i) {
    cout << "," << count[i];
  }
  cout << endl;
  cout << "WM"; 	// stands for Weight Matrix (and helps grepping)
  for(int i = 0 ; i < RbtHHSType::MAXTYPES; ++i) {
    cout << "," << area[i];
  }
  cout << endl;
}

void RbtSAIdxSF::ClearReceptor(void) {
  theIdxGrid = RbtNonBondedGridPtr();
  for (HHS_SolvationRListIter iter = theRSPList.begin(); iter != theRSPList.end(); iter++) {
    delete *iter;
  }
  for (HHS_SolvationRListIter iter = theFlexList.begin(); iter != theFlexList.end(); iter++) {
    delete *iter;
  }
  theRSPList.clear();
  theFlexList.clear();
  theCavList.clear();
  thePeriphList.clear();
  m_bFlexRec = false;
  m_site_0 = 0.0;
  m_site_free = 0.0;
  m_site_bound = 0.0;
}

void RbtSAIdxSF::ClearLigand(void)
{
  for (HHS_SolvationRListIter iter = theLSPList.begin(); iter != theLSPList.end(); iter++) {
    delete *iter;
  }
  theLSPList.clear();
  m_lig_0 = 0.0;
  m_lig_free = 0.0;
  m_lig_bound = 0.0;
}

void RbtSAIdxSF::ClearSolvent(void)
{
  for (HHS_SolvationRListIter iter = theSolventList.begin(); iter != theSolventList.end(); iter++) {
    delete *iter;
  }
  theSolventList.clear();
  m_solvent_0 = 0.0;
  m_solvent_free = 0.0;
  m_solvent_bound = 0.0;
}

void RbtSAIdxSF::ScoreMap(RbtStringVariantMap& scoreMap) const {
  if (isEnabled()) {
    //DM 8 June 2006. Divide the total solvation score into changes in INTER, INTRA and SYSTEM
    //THIS IS A POTENTIALLY SIGNIFICANT CHANGE AND REQUIRES CAREFUL VALIDATION
    //The total of the INTER, INTRA and SYSTEM partitions should equal the overall desolvation
    //raw score, but the reporting of the contributions has changed.

    //First we need to force the calculation of all the raw score components
    EnableAnnotations(true);
    RbtDouble rs = RawScore();
    EnableAnnotations(false);
    RbtString name = GetFullName();
    
    //INTER - the change in desolvation score for all components (site, ligand, solvent)
    //for the current internal conformations when intermolecular desolvation interactions
    //are taken into account. i.e. the difference between bound and free states. 
    RbtDouble inter_rs = (m_site_bound - m_site_free) + (m_lig_bound - m_lig_free)
      + (m_solvent_bound - m_solvent_free); 
    scoreMap[name] = inter_rs;
    AddToParentMapEntry(scoreMap, inter_rs);

    //INTRA - the change in the internal desolvation score for the ligand between
    //the initial ligand conformation and the current ligand conformation.
    //This is nothing to do with the binding event, and so belongs with the other
    //ligand intramolecular scores.
    RbtDouble intra_rs = m_lig_free - m_lig_0;
    RbtString intraName = RbtBaseSF::_INTRA_SF + "." + GetName();
    scoreMap[intraName] = intra_rs;
    scoreMap[intraName+".lig_0"] = m_lig_0;
    //increment the SCORE.INTRA total
    RbtDouble parentScore = scoreMap[RbtBaseSF::_INTRA_SF];
    parentScore += intra_rs * GetWeight();
    scoreMap[RbtBaseSF::_INTRA_SF] = parentScore;

    //SYSTEM - the change in the internal desolvation score for the site and solvent between
    //the initial conformations and the current conformations.
    //This is nothing to do with the binding event, and so belongs with the other
    //system scores.
    RbtDouble system_rs = (m_site_free - m_site_0) + (m_solvent_free - m_solvent_0);
    RbtString systemName = RbtBaseSF::_SYSTEM_SF + "." + GetName();
    scoreMap[systemName] = system_rs;
    //increment the SCORE.SYSTEM total
    parentScore = scoreMap[RbtBaseSF::_SYSTEM_SF];
    parentScore += system_rs * GetWeight();
    scoreMap[RbtBaseSF::_SYSTEM_SF] = parentScore;
  }
}

void RbtSAIdxSF::Setup() {
  RbtHHSType hhsType;
  RbtInt iTrace = GetTrace();
  m_maxR = 0.0;//keep track of maximum radius for any atom type
  if (iTrace > 1) {
    cout << endl << _CT << "::Setup()" << endl;
    cout << "TYPE,R,P,ASP,CHG_SCALING" << endl;
  }
  //Dummy read to force parsing of file, otherwise the first SetSection is overridden
  RbtStringList secList = m_spSolvSource->GetSectionList();
  RbtString _R("R");
  RbtString _P("P");
  RbtString _ASP("ASP");
  RbtString _CHG_SCALING("CHG_SCALING");
  m_solvTable.clear();

  for (RbtInt i = RbtHHSType::UNDEFINED; i < RbtHHSType::MAXTYPES; i++) {
    RbtString stri = hhsType.Type2Str(RbtHHSType::eType(i));
    m_spSolvSource->SetSection(stri);
    RbtDouble r = m_spSolvSource->GetParameterValue(_R);
    RbtDouble p = m_spSolvSource->GetParameterValue(_P);
    RbtDouble asp = m_spSolvSource->GetParameterValue(_ASP);
    RbtBool chg_scaling = m_spSolvSource->isParameterPresent(_CHG_SCALING);
    m_solvTable.push_back(solvprms(r,p,asp,chg_scaling));
    if (iTrace > 1) {
      cout << stri << "," << r << "," << p << "," << asp << "," << chg_scaling << endl;
    }
    m_maxR = std::max(m_maxR,r);
  }
  //Overall maximum range of the scoring function
  SetRange(2 * (m_maxR + HHS_Solvation::r_s));
  //Similar, but this represents the increment to any particular atom radius
  //to give the maximum range for that atom (used for indexing)
  SetParameter(_INCR,m_maxR + 2*HHS_Solvation::r_s);
  if (iTrace > 1) {
    cout << _CT << ": Maximum radius of any atom type = " << m_maxR << endl;
    cout << _CT << "::RANGE = " << GetRange() << endl;
    cout << _CT << "::INCR = " << GetParameter(_INCR) << endl;
  }
}

//Sum the exposed areas and energies (ASP*area) for the list of solvation interaction centers
//energy and area are the return values (by reference)
RbtDouble RbtSAIdxSF::TotalEnergy(const HHS_SolvationRList& intnCenters) const {
  //return std::accumulate(intnCenters.begin(),intnCenters.end(),0.0,Rbt::AccumHHSEnergy);
  RbtDouble energy(0.0);
  for (HHS_SolvationRListConstIter iter = intnCenters.begin(); iter != intnCenters.end(); ++iter)
    energy += (*iter)->GetEnergy();
  return energy;
}

//Request Handling method
//Handles the Partition request
void RbtSAIdxSF::HandleRequest(RbtRequestPtr spRequest) {
  RbtVariantList params = spRequest->GetParameters();
  RbtInt iTrace = GetTrace();

  switch (spRequest->GetID()) {
    
    //ID_REQ_SF_PARTITION requests come in two forms:
    //1-param: param[0] = distance => Partition all scoring functions
    //2-param: param[0] = SF Name,
    //         param[1] = distance => Partition a named scoring function
  case ID_REQ_SF_PARTITION:
    if (params.size() == 1) {
      if (iTrace > 2) {
	cout << _CT << "::HandleRequest: Partitioning " << GetFullName() << " at distance=" << params[0] << endl;
      }
      Partition(theLSPList,params[0]);
    }
    else if ( (params.size() == 2) && (params[0].String() == GetFullName())) {
      if (iTrace > 2) {
      cout << _CT << "::HandleRequest: Partitioning " << GetFullName() << " at distance=" << params[1] << endl;
      }
      Partition(theLSPList,params[1]);
    }
    break;
    
    //Pass all other requests to base handler
  default:
    if (iTrace > 2) {
      cout << _CT << "::HandleRequest: " << GetFullName() << " passing request to base handler" << endl;
    }
    RbtBaseObject::HandleRequest(spRequest);
    break;
  }
}
void RbtSAIdxSF::Partition(HHS_SolvationRList& intnCenters, RbtDouble dist) {
  //Rbt::PartitionHHS partition(dist);
  //std::for_each(intnCenters.begin(),intnCenters.end(),partition);
  for (HHS_SolvationRListIter iter = intnCenters.begin(); iter != intnCenters.end(); ++iter)
    (*iter)->Partition(dist);
}

HHS_SolvationRList RbtSAIdxSF::CreateInteractionCenters(const RbtAtomList& atomList) const {
  HHS_SolvationRList retList;
  // assign atom types and ASP parameters
  RbtHHSType hhsType;
  RbtInt iTrace = GetTrace();
  for (RbtAtomListConstIter iter = atomList.begin(); iter != atomList.end(); iter++) {
    RbtHHSType::eType t = hhsType(*iter);
    //Solvation model is calibrated for implicit non-polar H's, so ignore atoms of type H
    //Carbons are typed as though with implicit hydrogens, whether the H's are implicit or not,
    //so all should work well even for all-atom models
    if (t != RbtHHSType::H) {
      RbtDouble p_i = GetP_i(t);
      RbtDouble r_i = GetR_i(t);
      RbtDouble sigma = GetASP(t,(*iter)->GetGroupCharge());
      retList.push_back(new HHS_Solvation(t,*iter,p_i,r_i,sigma));
      if (t == RbtHHSType::UNDEFINED)
	cout << _CT << ": WARNING - CANNOT ASSIGN SOLVATION ATOM TYPE FOR " << (*iter)->GetFullAtomName() << endl;
    }
    else if (iTrace > 1)
      cout << _CT << ": INFO - ignoring non-polar H " << (*iter)->GetFullAtomName() << endl;
  }

  return retList;
}

void RbtSAIdxSF::BuildIntraMap(HHS_SolvationRList& ICList) const {
  //Calculate the invariant partial surface areas for all atoms in the list
  //Have to distinguish 1-2, 1-3, and 1-4+ invariant connections, and 1-4+ variable connections
  //If we come across a variable distance atom pair then store in the Variable list of the first atom
  for(HHS_SolvationRListIter iIter = ICList.begin(); iIter != ICList.end(); iIter++) {
    RbtAtom* pAtom1 = (*iIter)->GetAtom();
    Rbt::isAtom_12Connected is12(pAtom1);
    Rbt::isAtom_13Connected is13(pAtom1);
    Rbt::isAtomSelected is14variable;
    //Set the selection flags for all flexible interactions to this atom (1-4 variable)
    RbtModel* pModel = pAtom1->GetModelPtr();
    pModel->SetAtomSelectionFlags(false);
    pModel->SelectFlexAtoms(pAtom1);
    for(HHS_SolvationRListIter jIter = iIter+1; jIter != ICList.end(); jIter++) {
      RbtAtom* pAtom2 = (*jIter)->GetAtom();
      if (is12(pAtom2))
	(*iIter)->Overlap(*jIter,HHS_Solvation::Pij_12);
      else if (is13(pAtom2))
	(*iIter)->Overlap(*jIter,HHS_Solvation::Pij_13);
      else if (is14variable(pAtom2)) {
	(*iIter)->AddVariable(*jIter);
      }
      else //must be 1-4+ invariant
	(*iIter)->Overlap(*jIter,HHS_Solvation::Pij_14);
    }
  }
}

void RbtSAIdxSF::BuildIntraMap(HHS_SolvationRList& ICList1, HHS_SolvationRList& ICList2) const {
  //As above, but between two lists of interaction centers
  //Variable distances are stored in members of the first list
  for(HHS_SolvationRListIter iIter = ICList1.begin(); iIter != ICList1.end(); iIter++) {
    RbtAtom* pAtom1 = (*iIter)->GetAtom();
    Rbt::isAtom_12Connected is12(pAtom1);
    Rbt::isAtom_13Connected is13(pAtom1);
    Rbt::isAtomSelected is14variable;
    //Set the selection flags for all flexible interactions to this atom (1-4 variable)
    RbtModel* pModel = pAtom1->GetModelPtr();
    pModel->SetAtomSelectionFlags(false);
    pModel->SelectFlexAtoms(pAtom1);
    for(HHS_SolvationRListIter jIter = ICList2.begin(); jIter != ICList2.end(); jIter++) {
      RbtAtom* pAtom2 = (*jIter)->GetAtom();
      if (is12(pAtom2))
	(*iIter)->Overlap(*jIter,HHS_Solvation::Pij_12);
      else if (is13(pAtom2))
	(*iIter)->Overlap(*jIter,HHS_Solvation::Pij_13);
      else if (is14variable(pAtom2)) {
	(*iIter)->AddVariable(*jIter);
      }
      else //must be 1-4+ invariant
	(*iIter)->Overlap(*jIter,HHS_Solvation::Pij_14);
    }
  }
}
