/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtNmrSF.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtNmrSF.h"
#include "RbtNmrRestraintFileSource.h"
#include "RbtWorkSpace.h"

//Static data members
RbtString RbtNmrSF::_CT("RbtNmrSF");
RbtString RbtNmrSF::_FILENAME("FILENAME");
RbtString RbtNmrSF::_QUADRATIC("QUADRATIC");
	
//NB - Virtual base class constructor (RbtBaseSF) gets called first,
//implicit constructors for RbtBaseInterSF is called second
RbtNmrSF::RbtNmrSF(const RbtString& strName)
  : RbtBaseSF(_CT,strName), m_bQuadratic(true)
{
  //Add parameters
  AddParameter(_FILENAME,RbtString("default.noe"));
  AddParameter(_QUADRATIC,m_bQuadratic);
  SetRange(6.0);
  SetTrace(1);
#ifdef _DEBUG
  cout << _CT << " parameterised constructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtNmrSF::~RbtNmrSF() {
#ifdef _DEBUG
  cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}
	

void RbtNmrSF::SetupReceptor() {
  m_spGrid = RbtNonBondedGridPtr();
  if (GetReceptor().Null())
    return;

  RbtInt iTrace = GetTrace();

  //Create indexed grid for STD restraints penalty function
  m_spGrid = CreateNonBondedGrid();
  //Fixed distance range used for indexing the relevant receptor atoms
  RbtDouble range = GetRange() + GetMaxError();
  RbtDockingSitePtr spDS = GetWorkSpace()->GetDockingSite();

  //Find the non-polar hydrogens (nph) in the receptor
  //If we are using an implicit hydrogen model, need to include extended carbons (xc) also
  RbtAtomList recepAtomList = GetReceptor()->GetAtomList();
  RbtAtomList nphList = Rbt::GetAtomList(recepAtomList,Rbt::isAtomicNo_eq(1));
  nphList = Rbt::GetAtomList(nphList,std::not1(Rbt::isAtomHBondDonor()));
  RbtAtomList xcList = Rbt::GetAtomList(recepAtomList,Rbt::isAtomExtended());
  std::copy(xcList.begin(), xcList.end(), std::back_inserter(nphList));

  //Now limit the list to those atoms in the vicinity of the docking site
  nphList = spDS->GetAtomList(nphList,0.0,GetCorrectedRange());

  if (iTrace > 0) {
    cout << _CT << "::SetupReceptor() - " << nphList.size() << " nonpolar Hs and extended carbons near docking site" << endl;
    for (RbtAtomListConstIter iter = nphList.begin(); iter != nphList.end(); iter++) {
      cout << (**iter) << endl;
    }
  }
  
  for (RbtAtomListConstIter iter = nphList.begin(); iter != nphList.end(); iter++) {
    m_spGrid->SetAtomLists(*iter,range);
  }  
}

void RbtNmrSF::SetupLigand() {
  m_ligAtomList.clear();
  if (GetLigand().Null())
    return;
  
  m_ligAtomList = GetLigand()->GetAtomList();
}

void RbtNmrSF::SetupScore() {
  m_noeList.clear();
  m_stdList.clear();
  if (GetLigand().Null() || GetReceptor().Null())
    return;
  
  RbtInt iTrace = GetTrace();
  
  //Create a restraint filesource with the appropriate filename
  RbtString strRestrFile = Rbt::GetRbtFileName("data/receptors",GetParameter(_FILENAME));
  RbtNmrRestraintFileSourcePtr spRestrSource(new RbtNmrRestraintFileSource(strRestrFile));
  RbtString strName = GetName();
  
  //Now convert from restraint names to an actual restraint atoms list by matching against
  //the receptor and ligand atom lists, and setup the noe restraint interactions
  RbtAtomList atomList = GetReceptor()->GetAtomList();
  std::copy(m_ligAtomList.begin(),m_ligAtomList.end(),std::back_inserter(atomList));

  // (1) NOE Restraints
  //Get the noe restraint names list from the restraint file source
  RbtNoeRestraintNamesList noeNamesList = spRestrSource->GetNoeRestraintList();
  for (RbtNoeRestraintNamesListIter rIter = noeNamesList.begin(); rIter != noeNamesList.end(); rIter++) {
    //Create the NOE restraint - constructor matches the atoms in atomList with the names in *rIter
    RbtNoeRestraintAtoms noe(*rIter,atomList);
    if (iTrace > 0) {
      cout << endl << *rIter << endl << noe << endl;
    }
    //Only store if NOE is OK - i.e. has found the restraint names in the atom list
    if (noe.isOK()) {
      m_noeList.push_back(noe);
    }
    else {
      cout << "** WARNING - unable to match NOE restraint names to atom list" << endl;
      cout << noe << " not added" << endl;
    }
  }

  // (2) STD Restraints
  //Get the std restraint names list from the restraint file source
  RbtStdRestraintNamesList stdNamesList = spRestrSource->GetStdRestraintList();
  for (RbtStdRestraintNamesListIter rIter = stdNamesList.begin(); rIter != stdNamesList.end(); rIter++) {
    //Create the STD restraint - constructor matches the atoms in m_ligAtomList with the names in *rIter
    //Note for STD restraints it only makes senses to check the names against the ligand atoms
    RbtStdRestraintAtoms std(*rIter,m_ligAtomList);
    if (iTrace > 0) {
      cout << endl << *rIter << endl << std << endl;
    }
    //Only store if STD is OK - i.e. has found the restraint names in the atom list
    //and maxDist is within the range of the scoring function
    if (std.maxDist > GetRange()) {
      cout << "** WARNING - maxDist is beyond the range of the scoring function index" << endl;
      cout << "** Increase the RANGE parameter in " << GetFullName() << endl;
      cout << std << " not added" << endl;
    }
    else if (std.isOK()) {
      m_stdList.push_back(std);
    }
    else {
      cout << "** WARNING - unable to match STD restraint names to atom list" << endl;
      cout << std << " not added" << endl;
    }
  }
}

RbtDouble RbtNmrSF::RawScore() const {
  RbtDouble score(0.0);
  //QUADRATIC POTENTIAL
  if (m_bQuadratic) {
    for (RbtNoeRestraintAtomsListConstIter iter = m_noeList.begin(); iter != m_noeList.end(); iter++) {
      RbtDouble r = NoeDistance(*iter);
      RbtDouble dr = r - (*iter).maxDist;//delta(R)
      RbtDouble s = (dr > 0.0) ? dr*dr : 0.0;
      //cout << "(NOE) R,RMAX,DR,S=" << r << "," << (*iter).maxDist << "," << dr << "," << s << endl;
      score += s;
    }
    for (RbtStdRestraintAtomsListConstIter iter = m_stdList.begin(); iter != m_stdList.end(); iter++) {
      RbtDouble r = StdDistance(*iter);
      RbtDouble dr = r - (*iter).maxDist;//delta(R)
      RbtDouble s = (dr > 0.0) ? dr*dr : 0.0;
      //cout << "(STD) R,RMAX,DR,S=" << r << "," << (*iter).maxDist << "," << dr << "," << s << endl;
      score += s;
    }
  }
  //LINEAR POTENTIAL
  else {
    for (RbtNoeRestraintAtomsListConstIter iter = m_noeList.begin(); iter != m_noeList.end(); iter++) {
      RbtDouble r = NoeDistance(*iter);
      RbtDouble dr = r - (*iter).maxDist;//delta(R)
      RbtDouble s = (dr > 0.0) ? dr : 0.0;
      //cout << "(NOE) R,RMAX,DR,S=" << r << "," << (*iter).maxDist << "," << dr << "," << s << endl;
      score += s;
    }
    for (RbtStdRestraintAtomsListConstIter iter = m_stdList.begin(); iter != m_stdList.end(); iter++) {
      RbtDouble r = StdDistance(*iter);
      RbtDouble dr = r - (*iter).maxDist;//delta(R)
      RbtDouble s = (dr > 0.0) ? dr : 0.0;
      //cout << "(STD) R,RMAX,DR,S=" << r << "," << (*iter).maxDist << "," << dr << "," << s << endl;
      score += s;
    }
  }
  return score;
}

RbtDouble RbtNmrSF::NoeDistance(const RbtNoeRestraintAtoms& noe) const {
  //Simple, unambiguous restraint - just return the interatomic distance
  if (noe.isSimple()) {
    return Rbt::BondLength(noe.from.atoms.front(),noe.to.atoms.front());
  }
  //More complex cases
  //Compile the list of coords at each end
  //If the restraint end is defined of type MEAN, just store the center of mass of
  //the atoms in the list
  //If the restraint end is defined of type OR or AND, store all coords
  else {
    RbtCoordList fromCoords;
    if (noe.from.type==Rbt::NOE_MEAN)
      fromCoords.push_back(Rbt::GetCenterOfMass(noe.from.atoms));
    else
      fromCoords = Rbt::GetCoordList(noe.from.atoms);
    
    RbtCoordList toCoords;
    if (noe.to.type==Rbt::NOE_MEAN)
      toCoords.push_back(Rbt::GetCenterOfMass(noe.to.atoms));
    else
      toCoords = Rbt::GetCoordList(noe.to.atoms);
    
    RbtDouble dist_sq(0.0);//Keep track of minimum distance**2
    //Iterate over coords in each list and return the appropriate distance**2 between any of them
    for (RbtCoordListConstIter fIter = fromCoords.begin(); fIter != fromCoords.end(); fIter++) {
      //dist1_sq is the appropriate distance**2 between the current "from" coord
      //and all the coords in the "to" list.
      //i.e. if to.type==NOE_AND, dist1_sq is the max distance**2 to any atom in the "to" list
      //else dist1_sq is the min distance**2 to any atom in the "to" list
      RbtDouble dist1_sq(0.0);
      for (RbtCoordListConstIter tIter = toCoords.begin(); tIter != toCoords.end(); tIter++) {
	RbtDouble r12_sq = Rbt::Length2(*fIter,*tIter);
	dist1_sq = (tIter==toCoords.begin()) ? r12_sq
	  : (noe.to.type==Rbt::NOE_AND) ? std::max(dist1_sq,r12_sq)
	  : std::min(dist1_sq,r12_sq);
      }
      //dist_sq is the appropriate overall distance**2 for all the calculated dist1_sq's
      //i.e. if from.type==NOE_AND, dist_sq is the max of all the dist1_sq's
      //else dist_sq is the min of all the dist1_sq's
      dist_sq = (fIter==fromCoords.begin()) ? dist1_sq
	: (noe.from.type==Rbt::NOE_AND) ? std::max(dist_sq,dist1_sq)
	: std::min(dist_sq,dist1_sq);
    }
    return sqrt(dist_sq);
  }
}

RbtDouble RbtNmrSF::StdDistance(const RbtStdRestraintAtoms& std) const {
  //Compile the list of coords at each end
  //If the restraint end is defined of type MEAN, just store the center of mass of
  //the atoms in the list
  //If the restraint end is defined of type OR or AND, store all coords
  RbtCoordList fromCoords;
  if (std.from.type==Rbt::NOE_MEAN)
    fromCoords.push_back(Rbt::GetCenterOfMass(std.from.atoms));
  else
    fromCoords = Rbt::GetCoordList(std.from.atoms);
  
  RbtDouble dist_sq(999.9);//Keep track of minimum distance**2
  //Iterate over coords in each list and return the appropriate distance**2 between any of them
  for (RbtCoordListConstIter fIter = fromCoords.begin(); fIter != fromCoords.end(); fIter++) {
    //dist1_sq is the minimum distance**2 between the current "from" coord
    //and all the coords in the "to" list.
    RbtDouble dist1_sq(999.9);
    //For STD restraints, the list of "to" coords comes from the indexing grid
    const RbtAtomRList& toAtoms = m_spGrid->GetAtomList(*fIter);
    for (RbtAtomRListConstIter tIter = toAtoms.begin(); tIter != toAtoms.end(); tIter++) {
      RbtDouble r12_sq = Rbt::Length2(*fIter,(*tIter)->GetCoords());
      dist1_sq = (tIter==toAtoms.begin()) ? r12_sq : std::min(dist1_sq,r12_sq);
    }
    //dist_sq is the appropriate overall distance**2 for all the calculated dist1_sq's
    //i.e. if from.type==NOE_AND, dist_sq is the max of all the dist1_sq's
    //else dist_sq is the min of all the dist1_sq's
    dist_sq = (fIter==fromCoords.begin()) ? dist1_sq
      : (std.from.type==Rbt::NOE_AND) ? std::max(dist_sq,dist1_sq)
      : std::min(dist_sq,dist1_sq);
  }
  return sqrt(dist_sq);
}

void RbtNmrSF::ParameterUpdated(const RbtString& strName) {
  if (strName == _QUADRATIC) {
    m_bQuadratic = GetParameter(_QUADRATIC);
  }
  else {
    RbtBaseSF::ParameterUpdated(strName);
  }
}
