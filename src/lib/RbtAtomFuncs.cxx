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

#include "RbtAtomFuncs.h"

//DM 31 Oct 2000
//Given a bond, determines if it is in a ring (cutdown version of ToSpin)
RbtBool Rbt::FindCyclic(RbtBondPtr spBond, RbtAtomList& atomList, RbtBondList& bondList)
{
  //Max no. of atoms to process before giving up and return bCyclic=false
  //DM 8 Nov 2000 - effectively remove the cap
  const RbtInt maxProc = 999999;
	
  RbtAtomPtr spAtom1 = spBond->GetAtom1Ptr();
  RbtAtomPtr spAtom2 = spBond->GetAtom2Ptr();

  Rbt::SetAtomSelectionFlags(atomList,false);//Clear the atom selection flags
  Rbt::SetBondSelectionFlags(bondList,false);//Clear the bond selection flags

  //Temporary atom list containing atoms to be processed
  //Note: this is a true list (not a vector) as we will be making numerous insertions and deletions
  RbtAtomTrueList pendingAtomList;

  //Seed the pending list with atom 2
  spAtom2->SetSelectionFlag(true);
  spBond->SetSelectionFlag(true);
  pendingAtomList.push_back(spAtom2);

	RbtInt nProc=0;
	
  //While we still have atoms to process
  while (!pendingAtomList.empty() && (nProc < maxProc)) {
    //Take the last atom from the list and remove it
    RbtAtomPtr spA1 = pendingAtomList.back();
    pendingAtomList.pop_back();
    nProc++;

    //Get the bonds this atom is in
    const RbtBondMap& bondMap = spA1->GetBondMap();

    //Loop round each bond
    for (RbtBondMapConstIter bIter = bondMap.begin(); bIter != bondMap.end(); bIter++) {
      RbtBond* pBnd = (*bIter).first;//Remember not to use smart pointers!
      //If we haven't considered this bond before
      if ( !pBnd->GetSelectionFlag() ) {
				//If this atom is first in the bond (bIter->second == true), get atom 2 ptr, else get atom 1 ptr 
				RbtAtomPtr spA2 = ((*bIter).second) ? pBnd->GetAtom2Ptr() : pBnd->GetAtom1Ptr();
				if (spA2 == spAtom1) {
#ifdef _DEBUG
					cout << spAtom1->GetFullAtomName() << " is cyclic" << endl;
#endif //_DEBUG
	  			return true;
				}
				else if (!spA2->GetSelectionFlag()) {
	  			//Mark the atom and bond as selected and add the atom to the pending list
	  			spA2->SetSelectionFlag(true);
	  			pBnd->SetSelectionFlag(true);
	  			pendingAtomList.push_back(spA2);
				}
      }
    }
  }
  return false;
}

//DM 4 Dec 1998
//Given a bond, set the selection flags for all atoms which are connected to atom 2 of the bond
//Returns true if bond is in a ring (i.e. if atom 1's flag gets set also)
//DM 8 Feb 2000 - standalone version (formerly only available as RbtModel method)
RbtBool Rbt::ToSpin(RbtBondPtr spBond, RbtAtomList& atomList, RbtBondList& bondList)
{
  RbtAtomPtr spAtom1 = spBond->GetAtom1Ptr();
  RbtAtomPtr spAtom2 = spBond->GetAtom2Ptr();

  //Check if both atoms are actually in the model
  //if ( (spAtom1->GetModelPtr() != this) || (spAtom2->GetModelPtr() != this) )
  //  return false;

  Rbt::SetAtomSelectionFlags(atomList,false);//Clear the atom selection flags
  Rbt::SetBondSelectionFlags(bondList,false);//Clear the bond selection flags
  RbtBool bIsCyclic(false);

  //Temporary atom list containing atoms to be processed
  //Note: this is a true list (not a vector) as we will be making numerous insertions and deletions
  RbtAtomTrueList pendingAtomList;

  //Seed the pending list with atom 2
  spAtom2->SetSelectionFlag(true);
  spBond->SetSelectionFlag(true);
  pendingAtomList.push_back(spAtom2);

#ifdef _DEBUG
  //cout << "ToSpin: Bond ID=" << spBond->GetBondId() << endl;
#endif //_DEBUG

  //While we still have atoms to process
  while (!pendingAtomList.empty()) {
    //Take the last atom from the list and remove it
    RbtAtomPtr spA1 = pendingAtomList.back();
    pendingAtomList.pop_back();
#ifdef _DEBUG
    //cout << "ToSpin: Checking atom " << spA1->GetAtomId() << endl;
#endif //_DEBUG

    //Get the bonds this atom is in
    const RbtBondMap& bondMap = spA1->GetBondMap();

    //Loop round each bond
    for (RbtBondMapConstIter bIter = bondMap.begin(); bIter != bondMap.end(); bIter++) {
      RbtBond* pBnd = (*bIter).first;//Remember not to use smart pointers!
      //If we haven't considered this bond before
      if ( !pBnd->GetSelectionFlag() ) {
	//If this atom is first in the bond (bIter->second == true), get atom 2 ptr, else get atom 1 ptr 
	RbtAtomPtr spA2 = ((*bIter).second) ? pBnd->GetAtom2Ptr() : pBnd->GetAtom1Ptr();
#ifdef _DEBUG
	//	cout << "ToSpin: Checking bond " << spBnd->GetBondId() << " to atom " << spA2->GetAtomId() << endl;
#endif //_DEBUG
	//If we've got back to atom 1 of the bond passed into ToSpin, then the bond must be in a ring
	//so set the cyclic flag but don't consider this atom any further (we don't want to set the selection
	//flag otherwise we'd end up selecting the entire molecule).
	if (spA2 == spAtom1) {
	  bIsCyclic = true;
#ifdef _DEBUG
	  //cout << "ToSpin: We've hit " << spAtom1->GetAtomId() << " so must be cyclic!!" << endl;
#endif //_DEBUG
	}
	else if (!spA2->GetSelectionFlag()) {
	  //Mark the atom and bond as selected and add the atom to the pending list
	  spA2->SetSelectionFlag(true);
	  pBnd->SetSelectionFlag(true);
	  pendingAtomList.push_back(spA2);
#ifdef _DEBUG
	  //	  cout << "ToSpin: Adding atom " << spA2->GetAtomId() << " to the pending list" << endl;
#endif //_DEBUG
	}
      }
    }
  }
  return bIsCyclic;
}


//DM 7 Dec 1998
//Set the atom and bond cyclic flags for all atoms and bonds in the model
//DM 8 Feb 2000 - standalone version (formerly only available as RbtModel method)
void Rbt::SetAtomAndBondCyclicFlags(RbtAtomList& atomList, RbtBondList& bondList)
{
  //Clear the flags
  Rbt::SetAtomCyclicFlags(atomList,false);
  Rbt::SetBondCyclicFlags(bondList,false);

  //A bond is cyclic if Rbt::FindCyclic returns true
  //An atom is cyclic if it is involved in a cyclic bond
  //So loop around each bond, call FindCyclic and set the flags for the bond, atom 1 and atom 2 accordingly
  for (RbtBondListIter iter = bondList.begin(); iter != bondList.end(); iter++) {
    RbtBondPtr spBond(*iter);
    if (Rbt::FindCyclic(spBond,atomList,bondList)) {
      spBond->SetCyclicFlag(true);
      spBond->GetAtom1Ptr()->SetCyclicFlag(true);
      spBond->GetAtom2Ptr()->SetCyclicFlag(true);
    }
  }
}

//Find the smallest ring containing the given atom
//Assumes SetAtomAndBondCyclicFlags has already been called
//Note: if the atom is a member of two equally sized rings, only one will be returned
//30 Oct 2000 - Standalone version
RbtAtomList Rbt::FindRing(RbtAtomPtr spAtom, RbtBondList& bondList)
{
  //Check that atom is cyclic
  if (!spAtom->GetCyclicFlag())
    return RbtAtomList();

  Rbt::SetBondSelectionFlags(bondList,false);//Clear the bond selection flags

  //Get the cyclic bonds this atom is in
  RbtBondMap cyclicBondMap = spAtom->GetCyclicBondMap();
  RbtBondMapIter bIter = cyclicBondMap.begin();
  (*bIter).first->SetSelectionFlag(true);
  RbtAtomPtr spAtom2 = Rbt::GetBondedAtomPtr(*bIter);
  //Seed the partial ring lists with the first two atoms
  RbtAtomList firstRing;
  firstRing.push_back(spAtom);
  firstRing.push_back(spAtom2);

  //This is a list of atom lists (actually a vector of atom vectors)
  RbtAtomListList partialRings;
  partialRings.push_back(firstRing);

  //A bit risky, but loop forever and trust that we will find the seed atom eventually
  //After all, it is in a ring!
  //Whichever partial ring makes it first must be the smallest ring
  for (;;) {
    //First replicate any partial rings whose head atom is at a ring branch point
    //(defined by having >1 unselected cyclic bond)
    RbtUInt nRings = partialRings.size();
#ifdef _DEBUG
    //cout << "FindRing: " << nRings << " rings prior to forking" << endl;
#endif //_DEBUG

    //Use indices rather than iterators, as we may be increasing the number of partial rings in the ring list
    //(Iterators are unstable if the container increases in size inside the loop)
    for (RbtUInt ringIter = 0; ringIter < nRings; ringIter++) {
      //Get the head atom
      RbtAtomPtr spHeadAtom = (partialRings[ringIter]).back();
      //Get the cyclic bonds this atom is in
      RbtBondMap cyclicBondMap = spHeadAtom->GetCyclicBondMap();
      //Count how many are unselected
      RbtInt nUnSelected(0);
      //Loop round each bond
      for (RbtBondMapIter bIter = cyclicBondMap.begin(); bIter != cyclicBondMap.end(); bIter++) {
				if (!(*bIter).first->GetSelectionFlag())
	  			nUnSelected++;
      }
#ifdef _DEBUG
      //cout << "FindRing: #Unselected cyclic bonds = " << nUnSelected << " for " << ringIter << " at atom " << spHeadAtom->GetAtomId() << endl;
#endif //_DEBUG

      //Add a new copy of the partial ring for each unselected cyclic bond from this atom above the normal one
      for (RbtInt forkIter = 0; forkIter < (nUnSelected - 1); forkIter++) {
#ifdef _DEBUG
	//cout << "FindRing: Forking ring path " << ringIter << " at atom " << spHeadAtom->GetAtomId() << endl;
#endif //_DEBUG
				partialRings.push_back(partialRings[ringIter]);
      }
    }
    
    //Now grow each partial ring by one atom from the head atom
    for (RbtAtomListListIter rIter = partialRings.begin(); rIter != partialRings.end(); rIter++) {
      RbtAtomPtr spHeadAtom = (*rIter).back();
      //Get the cyclic bonds this atom is in
      RbtBondMap cyclicBondMap = spHeadAtom->GetCyclicBondMap();
      //Loop round each bond
      for (RbtBondMapIter bIter = cyclicBondMap.begin(); bIter != cyclicBondMap.end(); bIter++) {
				RbtBond* pBnd = (*bIter).first;//Remember not to use smart pointers!
				//Check we haven't considered this bond already
				if (!pBnd->GetSelectionFlag()) {
	  			//If this atom is first in the bond (bIter->second == true), get atom 2 ptr, else get atom 1 ptr 
  				RbtAtomPtr spA2 = Rbt::GetBondedAtomPtr(*bIter);
	  			//We've hit atom 1 so this must be the smallest ring
	  			if (spA2 == spAtom) {
#ifdef _DEBUG
	    			//cout << "FindRing: Found seed atom " << spAtom->GetAtomId() << endl;
#endif //_DEBUG
	    			return *rIter;
	  			}
	  			else {
#ifdef _DEBUG
	    			//cout << "FindRing: Adding atom " << spA2->GetAtomId() << " to " << spHeadAtom->GetAtomId() << endl;
#endif //_DEBUG
	    			(*rIter).push_back(spA2);
	    			pBnd->SetSelectionFlag(true);
	   			  break;//Break because we only want to add one atom to each head atom
	  			}
				}
      }
    }
  }
}


//30 Oct 2000 (DM) - Find all rings, standalone version
void Rbt::FindRings(RbtAtomList& atomList, RbtBondList& bondList, RbtAtomListList& ringList)
{
  ringList.clear();
  //Determine the cyclic atoms and bonds
  Rbt::SetAtomAndBondCyclicFlags(atomList,bondList);

  //First set the selection flags to the inverse of the cyclic flag
  //so we initially have all the cyclic atoms unselected
  for (RbtAtomListIter iter = atomList.begin(); iter != atomList.end(); iter++)
    (*iter)->SetSelectionFlag(!(*iter)->GetCyclicFlag());
  
  //Loop over all atoms and use each unselected atom as the seed for FindRing
  for (RbtAtomListIter iter = atomList.begin(); iter != atomList.end(); iter++) {
    if ( !(*iter)->GetSelectionFlag() ) {
      RbtAtomList ringAtoms = Rbt::FindRing(*iter,bondList);
      if (!ringAtoms.empty()) {
				//Now we've found a ring, set the selection flags to true for each ring member
				Rbt::SetAtomSelectionFlags(ringAtoms,true);
				//Store in our list of rings
				ringList.push_back(ringAtoms);
      }
    }
  }     
}

