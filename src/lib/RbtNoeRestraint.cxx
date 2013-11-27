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

#include "RbtNoeRestraint.h"

//Returns a string representing the NOE restraint type
RbtString Rbt::ConvertNoeRestraintTypeToString(const eNoeType type) {
  switch (type) {
  case Rbt::NOE_OR:
    return "OR";
  case Rbt::NOE_MEAN:
    return "MEAN";
  case Rbt::NOE_AND:
    return "AND";
  default:
    return "UNDEFINED";
  }
}

///////////////////////////////////////////////////////////////////////
//Simple struct for holding one end of an NOE restraint
//Has an atom name string list and a restraint type
RbtBool RbtNoeEndNames::isOK() const {
  return (!names.empty() && (type!=Rbt::NOE_UNDEFINED));
}

//Insertion operator for the above
ostream& operator<<(ostream& s, const RbtNoeEndNames& n) {
  if (n.names.size() == 1)
    s << n.names.front();
  else {
    s << Rbt::ConvertNoeRestraintTypeToString(n.type);
    s << "( ";
    for (RbtStringListConstIter iter = n.names.begin(); iter != n.names.end(); iter++)
      s << *iter << " ";
    s << ")";
  }
  return s;
}

///////////////////////////////////////////////////////////////////////
//As above, but with a real atom list in place of the string list

//Constructor accepting an RbtNoeEndNames and an atom list to match against
RbtNoeEndAtoms::RbtNoeEndAtoms(const RbtNoeEndNames& n, const RbtAtomList& atomList) {
  atoms = Rbt::GetMatchingAtomList(atomList,n.names);
  type = n.type;
}
RbtBool RbtNoeEndAtoms::isOK() const {
  return (!atoms.empty() && (type!=Rbt::NOE_UNDEFINED));
}

//Insertion operator for the above
ostream& operator<<(ostream& s, const RbtNoeEndAtoms& n) {
  if (n.atoms.size() == 1)
    s << n.atoms.front()->GetFullAtomName();
  else {
    s << Rbt::ConvertNoeRestraintTypeToString(n.type);
    s << "( ";
    if (n.atoms.size() <= 4) {
      for (RbtAtomListConstIter iter = n.atoms.begin(); iter != n.atoms.end(); iter++)
	s << (*iter)->GetFullAtomName() << " ";
    }
    else {
      s << n.atoms.front()->GetFullAtomName() << " ... "
	<< n.atoms.back()->GetFullAtomName() << " <"
	<< n.atoms.size() << " atoms> ";
    }
    s << ")";
  }
  return s;
}


///////////////////////////////////////////////////////////////////////
//Simple struct for holding a complete NOE restraint definition
//Has two RbtNoeEndNames, and max distance

RbtBool RbtNoeRestraintNames::isOK() const {
  return from.isOK() && to.isOK();
}

//Insertion operator for the above
ostream& operator<<(ostream& s, const RbtNoeRestraintNames& noe) {
  s << noe.from << " - " << noe.to << ": dist<" << noe.maxDist;
  return s;
}

///////////////////////////////////////////////////////////////////////
//As above, but with real RbtNoeEndAtoms in place of RbtNoeEndNames

//Constructor accepting an RbtNoeRestraintNames and an atom list to match against
RbtNoeRestraintAtoms::RbtNoeRestraintAtoms(const RbtNoeRestraintNames& n, const RbtAtomList& atomList)
  : from(n.from,atomList),to(n.to,atomList),maxDist(n.maxDist) {}

RbtBool RbtNoeRestraintAtoms::isOK() const {
  return from.isOK() && to.isOK();
}

//true if exactly one atom at each end of the NOE (unambiguous)
RbtBool RbtNoeRestraintAtoms::isSimple() const {
  return (from.atoms.size()==1) && (to.atoms.size()==1);
}

//Insertion operator for the above
ostream& operator<<(ostream& s, const RbtNoeRestraintAtoms& noe) {
  s << noe.from << " - " << noe.to << ": dist<" << noe.maxDist;
  return s;
}

///////////////////////////////////////////////////////////////////////
//Simple struct for holding a complete STD restraint definition
//Has one RbtNoeEndNames, and max distance from receptor

RbtBool RbtStdRestraintNames::isOK() const {
  return from.isOK();
}

//Insertion operator for the above
ostream& operator<<(ostream& s, const RbtStdRestraintNames& std) {
  s << std.from << " - receptor: dist<" << std.maxDist;
  return s;
}

///////////////////////////////////////////////////////////////////////
//As above, but with real RbtNoeEndAtoms in place of RbtNoeEndNames

//Constructor accepting an RbtStdRestraintNames and an atom list to match against
RbtStdRestraintAtoms::RbtStdRestraintAtoms(const RbtStdRestraintNames& s, const RbtAtomList& atomList)
  : from(s.from,atomList),maxDist(s.maxDist) {}

RbtBool RbtStdRestraintAtoms::isOK() const {
  return from.isOK();
}

//true if exactly one atom at in the STD definition (unambiguous)
RbtBool RbtStdRestraintAtoms::isSimple() const {
  return (from.atoms.size()==1);
}

//Insertion operator for the above
ostream& operator<<(ostream& s, const RbtStdRestraintAtoms& std) {
  s << std.from << " - receptor: dist<" << std.maxDist;
  return s;
}
