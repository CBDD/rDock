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

//Common functions involving atoms and bonds

#ifndef _RBTATOMFUNCS_H_
#define _RBTATOMFUNCS_H_

#include "RbtAtom.h"
#include "RbtBond.h"


namespace Rbt
{
  //DM 31 Oct 2000
  //Given a bond, determines if it is in a ring (cutdown version of ToSpin)
  RbtBool FindCyclic(RbtBondPtr spBond, RbtAtomList& atomList, RbtBondList& bondList);

  //DM 4 Dec 1998
  //Given a bond, set the selection flags for all atoms which are connected to atom 2 of the bond
  //Returns true if bond is in a ring (i.e. if atom 1's flag gets set also)
  //DM 8 Feb 2000 - standalone version (formerly only available as RbtModel method)
  //WARNING - no check that spBond is actually present in bondList, or that atom and bond lists
  //are consistent
  RbtBool ToSpin(RbtBondPtr spBond, RbtAtomList& atomList, RbtBondList& bondList);
  
  //DM 7 Dec 1998
  //Set the atom and bond cyclic flags for all atoms and bonds in the model
  //DM 8 Feb 2000 - standalone version (formerly only available as RbtModel method)
  void SetAtomAndBondCyclicFlags(RbtAtomList& atomList, RbtBondList& bondList);

	//Find the smallest ring containing the given atom
	//Assumes SetAtomAndBondCyclicFlags has already been called
	//Note: if the atom is a member of two equally sized rings, only one will be returned
	//30 Oct 2000 - Standalone version
	RbtAtomList FindRing(RbtAtomPtr spAtom, RbtBondList& bondList);
	
	//30 Oct 2000 (DM) - Find all rings, standalone version
	void FindRings(RbtAtomList& atomList, RbtBondList& bondList, RbtAtomListList& ringList);
}


#endif //_RBTATOMFUNCS_H_
