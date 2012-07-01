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
