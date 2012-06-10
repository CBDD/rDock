/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtModelMutator.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtModelMutator.h"
#include "RbtModel.h"
#include "RbtAtomFuncs.h"

#ifdef _DEBUG
namespace Tmp {
  void PrintAtomList(RbtAtomRListConstIter i1, RbtAtomRListConstIter i2) {
    for (RbtAtomRListConstIter iter = i1; iter!= i2; iter++) {
      cout << "\t" << (*iter)->GetFullAtomName() << " ("  << (*iter)->GetAtomId() << ")";
    }
    cout << endl;
  }
}
#endif //_DEBUG

////////////////////////////////////////
//Constructors/destructors
RbtModelMutator::RbtModelMutator(RbtModel* pModel, const RbtBondList& rotBonds, const RbtAtomList& tetheredAtoms)
  : m_pModel(pModel), m_rotBonds(rotBonds),m_tetheredAtoms(tetheredAtoms) {
  if (m_pModel != NULL) {
    Setup();
  }
  _RBTOBJECTCOUNTER_CONSTR_("RbtModelMutator");
}

RbtModelMutator::~RbtModelMutator() {
  _RBTOBJECTCOUNTER_DESTR_("RbtModelMutator");
}

////////////////////////////////////////
//Public methods
////////////////
  
const RbtAtomRListList& RbtModelMutator::GetFlexIntns() const {
  return m_flexIntns;
}

const RbtAtomRListList& RbtModelMutator::GetFlexAtoms() const {
  return m_rotAtoms;
}

RbtBondList RbtModelMutator::GetFlexBonds() const {
  return m_rotBonds;
}

////////////////////////////////////////
//Private methods
////////////////
void RbtModelMutator::Setup() {
  m_dih1Atoms.clear();
  m_dih2Atoms.clear();
  m_dih3Atoms.clear();
  m_dih4Atoms.clear();
  m_rotAtoms.clear();
  m_flexIntns.clear();

  if (m_pModel == NULL)
    return;

  RbtInt nAtoms = m_pModel->GetNumAtoms();
  RbtInt nTethered = m_tetheredAtoms.size();
  m_flexIntns = RbtAtomRListList(nAtoms,RbtAtomRList());

  //Loop over each rotable bond
  for (RbtBondListConstIter iter = m_rotBonds.begin(); iter != m_rotBonds.end(); iter++) {
    RbtAtom* pAtom2 = (*iter)->GetAtom1Ptr();
    RbtAtom* pAtom3 = (*iter)->GetAtom2Ptr();    
    //The following lines get the bonded atom lists on each end of the rotable bond, 
    //taking care not to include the atoms actually in the rotable bond
    RbtAtomList bondedAtoms2 = Rbt::GetAtomList(Rbt::GetBondedAtomList(pAtom2),
                               std::not1(std::bind2nd(Rbt::isAtomPtr_eq(),pAtom3)));
    RbtAtomList bondedAtoms3 = Rbt::GetAtomList(Rbt::GetBondedAtomList(pAtom3),
                               std::not1(std::bind2nd(Rbt::isAtomPtr_eq(),pAtom2)));
    //Assertion - check bonded atom lists are not empty
    Assert<RbtAssert>(!MUT_CHECK || !(bondedAtoms2.empty() || bondedAtoms3.empty()));

    Rbt::ToSpin(*iter,m_pModel->m_atomList,m_pModel->m_bondList);
    pAtom2->SetSelectionFlag(false);
    pAtom3->SetSelectionFlag(false);
    //If we have selected over half the molecule to rotate then invert the selection
    //to minimise the number of atoms to rotate
    //DM 2 Jul 2002 - in tethered mode we want to ensure that the minimum number of *tethered* atoms are rotated
    //(preferably none). i.e. we rotate the free end of the bond, even if this is over half the molecule
    RbtInt nSelected = (nTethered == 0) ? Rbt::GetNumSelectedAtoms(m_pModel->m_atomList) : Rbt::GetNumSelectedAtoms(m_tetheredAtoms);
    RbtInt nHalf = (nTethered == 0) ? (nAtoms-2)/2 : (nTethered-2)/2;
    if (nSelected > nHalf) {
      //cout << "Over half the molecule selected: " << nSelected << " atoms" << endl;
      Rbt::InvertAtomSelectionFlags(m_pModel->m_atomList);
      pAtom2->SetSelectionFlag(false);
      pAtom3->SetSelectionFlag(false);
      m_dih1Atoms.push_back(bondedAtoms3.front());
      m_dih2Atoms.push_back(pAtom3);
      m_dih3Atoms.push_back(pAtom2);
      m_dih4Atoms.push_back(bondedAtoms2.front());
      //cout << "Inverted: " << Rbt::GetNumSelectedAtoms(m_pModel->m_atomList) << " atoms now selected" << endl;
      //cout << "Dihedral spec: " << bondedAtoms3.front()->GetAtomName() << "\t" << pAtom3->GetAtomName() << "\t"
      //	   << pAtom2->GetAtomName() << "\t" << bondedAtoms2.front()->GetAtomName() << endl;
    }
    else {
      m_dih1Atoms.push_back(bondedAtoms2.front());
      m_dih2Atoms.push_back(pAtom2);
      m_dih3Atoms.push_back(pAtom3);
      m_dih4Atoms.push_back(bondedAtoms3.front());
      //cout << "Dihedral spec: " << bondedAtoms2.front()->GetAtomName() << "\t" << pAtom2->GetAtomName() << "\t"
      //	   << pAtom3->GetAtomName() << "\t" << bondedAtoms3.front()->GetAtomName() << endl;
    }

    //Store the smaller atom list (or free atom list in tethered mode) for this rotable bond
    RbtAtomRList sList;
    std::copy_if(m_pModel->m_atomList.begin(),m_pModel->m_atomList.end(),std::back_inserter(sList),Rbt::isAtomSelected());
    m_rotAtoms.push_back(sList);
    
    //DM 25 Apr 2002
    //Work out the flexible interactions across the bond
    //Method is as follows.
    //sList = selected atoms on one side of the bond
    //uList = unselected atoms on other side of the bond
    //By definition, all atom pairs from {sList}-{uList} are flexible interactions
    //but we need to avoid counting interaction pairs more than once
    //foreach rotable bond {
    //  Determine sList, uList {
    //  foreach sAtom in sList {
    //    foreach uAtom in uList {
    //      append uAtom to sAtom's interaction list
    //      append sAtom to uAtom's interaction list
    //    }
    //  }
    //}
    //foreach atom in model {
    //  sort interaction list for atom
    //  remove duplicates from interaction list
    //}
    Rbt::InvertAtomSelectionFlags(m_pModel->m_atomList);
    pAtom2->SetSelectionFlag(false);
    pAtom3->SetSelectionFlag(false);
    RbtAtomRList uList;
    std::copy_if(m_pModel->m_atomList.begin(),m_pModel->m_atomList.end(),std::back_inserter(uList),Rbt::isAtomSelected());
    //Sort the lists of selected and unselected atoms by atom ID
    std::sort(sList.begin(),sList.end(),Rbt::RbtAtomPtrCmp_AtomId());
    std::sort(uList.begin(),uList.end(),Rbt::RbtAtomPtrCmp_AtomId());

    //Loop over selected atom list
    //and append the unselected atom list to the interaction list for each atom
    for (RbtAtomRListConstIter sIter = sList.begin(); sIter != sList.end(); sIter++) {
      RbtInt id = (*sIter)->GetAtomId()-1;
      //Assertion - check id is within range
      Assert<RbtAssert>(!MUT_CHECK || (id >= 0 && id < m_flexIntns.size()));
      RbtAtomRListListIter lIter = m_flexIntns.begin() + id;
      std::copy(uList.begin(),uList.end(),std::back_inserter(*lIter));
    }

    //Ditto for unselected atom list
    for (RbtAtomRListConstIter uIter = uList.begin(); uIter != uList.end(); uIter++) {
      RbtInt id = (*uIter)->GetAtomId()-1;
      //Assertion - check id is within range
      Assert<RbtAssert>(!MUT_CHECK || (id >= 0 && id < m_flexIntns.size()));
      RbtAtomRListListIter lIter = m_flexIntns.begin() + id;
      std::copy(sList.begin(),sList.end(),std::back_inserter(*lIter));
    }
  }
  //Now clean up the interaction lists for each atom to remove duplicates
  for (RbtAtomRListListIter lIter = m_flexIntns.begin(); lIter != m_flexIntns.end(); lIter++) {
    if (!(*lIter).empty()) {
      std::sort((*lIter).begin(),(*lIter).end(),Rbt::RbtAtomPtrCmp_AtomId());
      //std::unique will shuffle all the duplicate elements to the end of the vector
      //and return an iterator to the end of the unique elements
      RbtAtomRListIter uniqIter = std::unique((*lIter).begin(),(*lIter).end());
#ifdef _DEBUG
      RbtInt id = lIter - m_flexIntns.begin();
      RbtAtom* pAtom = m_pModel->m_atomList[id];
      cout << "Atom " << pAtom->GetFullAtomName() << " (" << id+1 << "): with dups=" << (*lIter).size();
#endif //_DEBUG
      //Remove duplicates
      (*lIter).erase(uniqIter,(*lIter).end());
#ifdef _DEBUG
      cout << "; unique=" << (*lIter).size() << endl;
      Tmp::PrintAtomList((*lIter).begin(),(*lIter).end());
#endif //_DEBUG
    }
  }
}
