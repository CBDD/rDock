/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtModelMutator.h#4 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Legacy class, largely replaced by RbtChromElement subclasses
//The only remaining purpose of RbtModelMutator is to generate lists
//of flexible intramolecular interactions.
//i.e. intramolecular atom pairs across rotatable bonds
#ifndef _RBTMODELMUTATOR_H_
#define _RBTMODELMUTATOR_H_

#include "RbtAtom.h"
#include "RbtBond.h"
#include "RbtPrincipalAxes.h"

class RbtModel;//forward declaration

//Only check request assertions in debug build
#ifdef _NDEBUG
const RbtBool MUT_CHECK = true;
#else
const RbtBool MUT_CHECK = true;
#endif //_NDEBUG

class RbtModelMutator
{
 public:
  ////////////////////////////////////////
  //Constructors/destructors
  RbtModelMutator(RbtModel* pModel, const RbtBondList& rotBonds, const RbtAtomList& tetheredAtoms);
  virtual ~RbtModelMutator(); //Default destructor
  
  ////////////////////////////////////////
  //Public methods
  ////////////////
  const RbtAtomRListList& GetFlexIntns() const;
  const RbtAtomRListList& GetFlexAtoms() const;
  RbtBondList GetFlexBonds() const;
  
 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////
  
  
 private:
  ////////////////////////////////////////
  //Private methods
  /////////////////
  
  RbtModelMutator(const RbtModelMutator&);//Copy constructor disabled by default
  RbtModelMutator& operator=(const RbtModelMutator&);//Copy assignment disabled by default
  
  void Setup();

 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////
  
  
 private:
  ////////////////////////////////////////
  //Private data
  //////////////
  RbtModel* m_pModel;//Parent model
  RbtBondList m_rotBonds;//A record of the defined rotatable bonds passed to the constructor
  RbtAtomRList m_dih1Atoms;//Atom 1 in each rotable bond dihedral definition
  RbtAtomRList m_dih2Atoms;//Atom 2 in each rotable bond dihedral definition
  RbtAtomRList m_dih3Atoms;//Atom 3 in each rotable bond dihedral definition
  RbtAtomRList m_dih4Atoms;//Atom 4 in each rotable bond dihedral definition
  RbtAtomRListList m_rotAtoms;//List of atom lists on smallest side of each rotable bond
  //DM 25 Apr 2002
  //m_flexIntns is a vector of size N(atoms)
  //Each element is an RbtAtomList (i.e. vector<RbtAtomPtr>)
  //
  //*** The vector is indexed by AtomId
  //*** We assume all atom IDs in the molecule are assigned sequentially from 1 ***
  //*** For molecules read from RbtMdlFileSource and RbtPsfFileSource this is guaranteed
  //
  //m_flexIntns[i-1] is a vector of all the atoms which can interact with atom i across at least one rotable bond:
  //(these interaction distances are conformation dependent)
  //This is a symmetric matrix -
  //For a flexible interaction i,j, atom j will appear in the RbtAtomList for atom i, and vice versa.
  //Atoms which are NOT in this list are therefore in the same rigid fragment as atom i:
  //(these interaction distances are fixed and can therefore be ignored by any scoring function)
  RbtAtomRListList m_flexIntns;

  //DM 2 Jul 2002 - if model has tethered atoms we need to modify the behaviour of RbtModelMutator
  //Could be subclassed if necessary
  RbtAtomList m_tetheredAtoms;
};

//Useful typedefs
typedef SmartPtr<RbtModelMutator> RbtModelMutatorPtr;//Smart pointer

#endif //_RBTMODELMUTATOR_H_
