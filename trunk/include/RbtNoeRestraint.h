/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtNoeRestraint.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Miscellaneous helper classes/structs for handling NOE interactions

#ifndef _RBTNOERESTRAINT_H_
#define _RBTNOERESTRAINT_H_

#include "RbtAtom.h"

namespace Rbt
{
  //NOE Restraint type
  //UNDEFINED: undefined restraint type
  //OR: The distance to each of the atoms in the list is calculated separately,
  //    and any of these can satisfy the restraint.
  //    In practice, this means only the shortest distance is scored
  //MEAN: The distance to the mean coords of the atoms in the list is calculated,
  //      and used to evaluate the restraint
  //OR: The distance to each of the atoms in the list is calculated separately,
  //    and all of these should satisfy the restraint.
  //    In practice, this means only the longest distance is scored
  enum eNoeType {
    NOE_UNDEFINED,
    NOE_OR,      
    NOE_MEAN,
    NOE_AND
  };

  //Returns a string representing the NOE restraint type
  RbtString ConvertNoeRestraintTypeToString(const eNoeType type);
}

///////////////////////////////////////////////////////////////////////
//Simple struct for holding one end of an NOE restraint
//Has an atom name string list and a restraint type
class RbtNoeEndNames {
 public:
  RbtNoeEndNames() {}
  RbtBool isOK() const;
  //Public data
  RbtStringList names;
  Rbt::eNoeType type;
};

//Insertion operator for the above
ostream& operator<<(ostream& s, const RbtNoeEndNames& n);

///////////////////////////////////////////////////////////////////////
//As above, but with a real atom list in place of the string list
class RbtNoeEndAtoms {
 public:
  //Default constructor
  RbtNoeEndAtoms() {}
  //Constructor accepting an RbtNoeEndNames and an atom list to match against
  RbtNoeEndAtoms(const RbtNoeEndNames& n, const RbtAtomList& atomList);
  RbtBool isOK() const;
  //Public data
  RbtAtomList atoms;
  Rbt::eNoeType type;
};

//Insertion operator for the above
ostream& operator<<(ostream& s, const RbtNoeEndAtoms& n);


///////////////////////////////////////////////////////////////////////
//Simple struct for holding a complete NOE restraint definition
//Has two RbtNoeEndNames, and max distance
class RbtNoeRestraintNames {
 public:
  RbtNoeRestraintNames() : maxDist(5.0) {}
  RbtBool isOK() const;
  //Public data
  RbtNoeEndNames from;
  RbtNoeEndNames to;
  RbtDouble maxDist;
};

//Insertion operator for the above
ostream& operator<<(ostream& s, const RbtNoeRestraintNames& noe);

///////////////////////////////////////////////////////////////////////
//As above, but with real RbtNoeEndAtoms in place of RbtNoeEndNames
class RbtNoeRestraintAtoms {
 public:
  //Default constructor
  RbtNoeRestraintAtoms() : maxDist(5.0) {}
  //Constructor accepting an RbtNoeRestraintNames and an atom list to match against
  RbtNoeRestraintAtoms(const RbtNoeRestraintNames& n, const RbtAtomList& atomList);
  RbtBool isOK() const;//true if at least one atom at each end of the NOE and NOE type is not undefined
  RbtBool isSimple() const;//true if exactly one atom at each end of the NOE (unambiguous)
  //Public data
  RbtNoeEndAtoms from;
  RbtNoeEndAtoms to;
  //RbtDouble minDist;
  RbtDouble maxDist;
};

//Insertion operator for the above
ostream& operator<<(ostream& s, const RbtNoeRestraintAtoms& noe);

///////////////////////////////////////////////////////////////////////
//Simple struct for holding a complete STD restraint definition
//Has one RbtNoeEndNames, and max distance (from receptor)
class RbtStdRestraintNames {
 public:
  RbtStdRestraintNames() : maxDist(5.0) {}
  RbtBool isOK() const;
  //Public data
  RbtNoeEndNames from;
  RbtDouble maxDist;
};

//Insertion operator for the above
ostream& operator<<(ostream& s, const RbtStdRestraintNames& std);

///////////////////////////////////////////////////////////////////////
//As above, but with real RbtNoeEndAtoms in place of RbtNoeEndNames
class RbtStdRestraintAtoms {
 public:
  //Default constructor
  RbtStdRestraintAtoms() : maxDist(5.0) {}
  //Constructor accepting an RbtStdRestraintNames and an atom list to match against
  RbtStdRestraintAtoms(const RbtStdRestraintNames& n, const RbtAtomList& atomList);
  RbtBool isOK() const;//true if at least one atom in the definition and type is not undefined
  RbtBool isSimple() const;//true if exactly one atom at in the STD definition (unambiguous)
  //Public data
  RbtNoeEndAtoms from;
  RbtDouble maxDist;
};

//Insertion operator for the above
ostream& operator<<(ostream& s, const RbtStdRestraintAtoms& std);

//Useful typedefs
typedef vector<RbtNoeRestraintNames> RbtNoeRestraintNamesList;
typedef RbtNoeRestraintNamesList::iterator RbtNoeRestraintNamesListIter;
typedef RbtNoeRestraintNamesList::const_iterator RbtNoeRestraintNamesListConstIter;

typedef vector<RbtNoeRestraintAtoms> RbtNoeRestraintAtomsList;
typedef RbtNoeRestraintAtomsList::iterator RbtNoeRestraintAtomsListIter;
typedef RbtNoeRestraintAtomsList::const_iterator RbtNoeRestraintAtomsListConstIter;

typedef vector<RbtStdRestraintNames> RbtStdRestraintNamesList;
typedef RbtStdRestraintNamesList::iterator RbtStdRestraintNamesListIter;
typedef RbtStdRestraintNamesList::const_iterator RbtStdRestraintNamesListConstIter;

typedef vector<RbtStdRestraintAtoms> RbtStdRestraintAtomsList;
typedef RbtStdRestraintAtomsList::iterator RbtStdRestraintAtomsListIter;
typedef RbtStdRestraintAtomsList::const_iterator RbtStdRestraintAtomsListConstIter;

#endif //_RBTNOERESTRAINT_H_
