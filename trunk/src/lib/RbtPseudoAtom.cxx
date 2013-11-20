/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtPseudoAtom.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtPseudoAtom.h"

///////////////////////////////////////////////
//Constructors / destructors

//Constructor supplying atom list of constituent atoms
//Note: the base RbtAtom constructor sets up default attributes for the pseudoatom only
//(the constituent atoms are already constructed so we don't need to worry about them)
//Use Atomic No = -1 to signify a pseudoatom
RbtPseudoAtom::RbtPseudoAtom(const RbtAtomList& atomList, RbtInt nAtomId, RbtString strAtomName) throw (RbtError)
  : RbtAtom(nAtomId,-1,strAtomName)
{
  //Check atom list isn't empty (should be an assert)
  if (atomList.empty())
    throw RbtBadArgument(_WHERE_,"Null atom list passed to RbtPseudoAtom constructor");

  //Store the atom list
  m_atomList = atomList;
  
  //Base the pseudoatom attributes on the first atom in the list
  //Note: this is a write-once action; if the attributes of the first
  //atom change, the pseudoatom attributes will not.
  //Assumes that all constituent atoms are in the same subunit and segment
  RbtAtomPtr spAtom1(atomList.front());//1st atom in list
  RbtAtomPtr spAtomN(atomList.back());//Last atom in list
  //Pseudo-atom name is i the format [C1..N4]
  strAtomName += "["+spAtom1->GetAtomName()+".."+spAtomN->GetAtomName()+"]";
  SetAtomName(strAtomName);
  SetSubunitId(spAtom1->GetSubunitId());
  SetSubunitName(spAtom1->GetSubunitName());
  SetSegmentName(spAtom1->GetSegmentName());
  UpdateCoords();
  _RBTOBJECTCOUNTER_CONSTR_("RbtPseudoAtom");
}

//Default destructor
RbtPseudoAtom::~RbtPseudoAtom()
{
  _RBTOBJECTCOUNTER_DESTR_("RbtPseudoAtom");
}

//Virtual function for dumping atom details to an output stream
//Derived classes (e.g. pseudoatom) can override if required
ostream& RbtPseudoAtom::Print(ostream& s) const
{
  RbtAtom::Print(s);//First call base class Print on the pseudo-atom itself
  s << endl << "Pseudo-atom with " << m_atomList.size() << " constituent atoms:" << endl;
  //Now dump the atom details for all the constituent atoms
  for (RbtAtomListConstIter iter = m_atomList.begin(); iter != m_atomList.end(); iter++)
    s << (**iter) << endl;
  return s;
}

///////////////////////////////////////////////
//Public accessor functions
//3-D attributes

//Coords
//DM 11 Jul 2000 - base GetCoords method is now non-virtual again
//Pseudoatoms must update their mean coords by calling UpdateCoords
//each time the underlying atoms are moved
void RbtPseudoAtom::UpdateCoords() {
  RbtCoord coord;
  RbtUInt nSize = m_atomList.size();
  if (nSize > 0) {
    for (RbtAtomListConstIter iter = m_atomList.begin(); iter != m_atomList.end(); iter++) {
      coord += (*iter)->GetCoords();
    }
    coord /= nSize;
  }
  SetCoords(coord);
}

