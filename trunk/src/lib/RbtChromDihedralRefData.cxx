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

#include "RbtChromDihedralRefData.h"
#include "RbtModel.h"
#include "RbtAtomFuncs.h"

RbtString RbtChromDihedralRefData::_CT = "RbtChromDihedralRefData";

RbtChromDihedralRefData::RbtChromDihedralRefData(RbtBondPtr spBond,
                                                     RbtAtomList tetheredAtoms,
                                                     RbtDouble stepSize,
                                                     RbtChromElement::eMode mode,
                                                     RbtDouble maxDihedral)
    : m_stepSize(stepSize),
      m_mode(mode),
      m_maxDihedral(maxDihedral)
{
    Setup(spBond, tetheredAtoms);
    m_initialValue = GetModelValue();
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtChromDihedralRefData::~RbtChromDihedralRefData() {
    _RBTOBJECTCOUNTER_DESTR_(_CT);    
}

RbtDouble RbtChromDihedralRefData::GetModelValue() const {
	return Rbt::BondDihedral(m_atom1, m_atom2, m_atom3, m_atom4);
}

void RbtChromDihedralRefData::SetModelValue(RbtDouble dihedralAngle) {
	RbtDouble delta = dihedralAngle - GetModelValue();
	//Only rotate if delta is non-zero
	if (fabs(delta)>0.001) {
	    //Coords of atom 1
	    RbtCoord coord1(m_atom2->GetCoords());
	    //Vector along the bond between atom 1 and atom 2 (rotation axis)
	    RbtVector bondVector(m_atom3->GetCoords() - coord1);
	    RbtVector toOrigin = -coord1;
	    RbtQuat quat(bondVector, delta*M_PI/180.0);
	    for (RbtAtomRListIter iter = m_rotAtoms.begin();
	    		                         iter != m_rotAtoms.end();
	    		                         ++iter) {
	      (*iter)->Translate(toOrigin);
	      (*iter)->RotateUsingQuat(quat);
	      (*iter)->Translate(coord1);
	    }
	}
}

void RbtChromDihedralRefData::Setup(RbtBondPtr spBond,
										const RbtAtomList& tetheredAtoms)
{
    RbtAtom* pAtom2 = spBond->GetAtom1Ptr();
    RbtAtom* pAtom3 = spBond->GetAtom2Ptr();
    //WARNING - we assume that the bond atoms are registered with a model
    RbtModel* pModel = pAtom2->GetModelPtr();
	RbtAtomList atomList = pModel->GetAtomList();
	RbtBondList bondList = pModel->GetBondList();
	RbtInt nAtoms = atomList.size();
	RbtInt nTethered = tetheredAtoms.size();
    //The following lines get the bonded atom lists on each end of the rotable bond, 
    //taking care not to include the atoms actually in the rotable bond
    RbtAtomList bondedAtoms2 = Rbt::GetAtomList(Rbt::GetBondedAtomList(pAtom2),
                               std::not1(std::bind2nd(Rbt::isAtomPtr_eq(),pAtom3)));
    RbtAtomList bondedAtoms3 = Rbt::GetAtomList(Rbt::GetBondedAtomList(pAtom3),
                               std::not1(std::bind2nd(Rbt::isAtomPtr_eq(),pAtom2)));
    //Assertion - check bonded atom lists are not empty
    Assert<RbtAssert>(!MUT_CHECK || !(bondedAtoms2.empty() || bondedAtoms3.empty()));

    Rbt::ToSpin(spBond,atomList,bondList);
    pAtom2->SetSelectionFlag(false);
    pAtom3->SetSelectionFlag(false);
    //If we have selected over half the molecule to rotate then invert the selection
    //to minimise the number of atoms to rotate
    //DM 2 Jul 2002 - in tethered mode we want to ensure that the minimum number of *tethered* atoms are rotated
    //(preferably none). i.e. we rotate the free end of the bond, even if this is over half the molecule
    RbtInt nSelected = (nTethered == 0)
    		? Rbt::GetNumSelectedAtoms(atomList)
    		: Rbt::GetNumSelectedAtoms(tetheredAtoms);
    RbtInt nHalf = (nTethered == 0)
    		? (nAtoms-2)/2
    		: (nTethered-2)/2;
    if (nSelected > nHalf) {
      //cout << "Over half the molecule selected: " << nSelected << " atoms" << endl;
      Rbt::InvertAtomSelectionFlags(atomList);
      pAtom2->SetSelectionFlag(false);
      pAtom3->SetSelectionFlag(false);
      m_atom1 = bondedAtoms3.front();
      m_atom2 = pAtom3;
      m_atom3 = pAtom2;
      m_atom4 = bondedAtoms2.front();
      //cout << "Inverted: " << Rbt::GetNumSelectedAtoms(m_pModel->m_atomList) << " atoms now selected" << endl;
      //cout << "Dihedral spec: " << bondedAtoms3.front()->GetAtomName() << "\t" << pAtom3->GetAtomName() << "\t"
      //	   << pAtom2->GetAtomName() << "\t" << bondedAtoms2.front()->GetAtomName() << endl;
    }
    else {
      m_atom1 = bondedAtoms2.front();
      m_atom2 = pAtom2;
      m_atom3 = pAtom3;
      m_atom4 = bondedAtoms3.front();
      //cout << "Dihedral spec: " << bondedAtoms2.front()->GetAtomName() << "\t" << pAtom2->GetAtomName() << "\t"
      //	   << pAtom3->GetAtomName() << "\t" << bondedAtoms3.front()->GetAtomName() << endl;
    }

    //Store the smaller atom list (or free atom list in tethered mode) for this rotable bond
    m_rotAtoms.clear();
    std::copy_if(atomList.begin(), atomList.end(), std::back_inserter(m_rotAtoms),
    				Rbt::isAtomSelected());
}
