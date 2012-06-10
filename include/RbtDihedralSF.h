/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtDihedralSF.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Base implementation class for all dihedral scoring functions. Provides
//RbtDihedral struct for dihedral definitions and score, plus a method for
//generating lists of dihedral definitions from a bond list.
//This version hardwired to use Tripos 5.2 atom types and params.

#ifndef _RBTDIHEDRALSF_H_
#define _RBTDIHEDRALSF_H_

#include "RbtBaseSF.h"
#include "RbtAtom.h"
#include "RbtBond.h"
#include "RbtTriposAtomType.h"
#include "RbtParameterFileSource.h"

//Class for holding dihedral atom specifiers and force field params
//Main method is operator() - calculates single dihedral score for current conformation
//DM (9 Dec 2003): Extended to allow additional terms to be added to the dihedral potential (AddTerm)
//Two main uses for AddTerm:
//  i) to define a Fourier series expansion, rather than a simple cosine potential
// ii) implicit hydrogen correction (ghost terms are added for the missing hydrogen dihedrals
//                                   defined as offsets from the real heavy-atom dihedral specifiers)

class RbtDihedral {
 public:

  //Simple struct for storing the force field params for a single term of the dihedral potential
  class prms {
  public:
    prms(RbtDouble ss=3.0, RbtDouble kk=0.2, RbtDouble o=0.0) :
      s(fabs(ss)), k(kk), sign( (ss > 0.0) ? +1.0 : -1.0), offset(o) {}
    RbtDouble s;//rotational degeneracy
    RbtDouble k;//barrier height
    RbtDouble sign;//+1 or -1, to invert potential
    RbtDouble offset;//dihedral offset in degrees (useful for ghost implicit hydrogens)
  };

  //Constructor takes the real atom specifiers, plus the first term of the potential
  RbtDihedral(RbtAtom* pAtom1, RbtAtom* pAtom2, RbtAtom* pAtom3, RbtAtom* pAtom4, const prms& dihprms);
  RbtDouble operator() () const;//Calculate dihedral score for this interaction
  RbtAtom* GetAtom1Ptr() const {return m_pAtom1;}
  RbtAtom* GetAtom2Ptr() const {return m_pAtom2;}
  RbtAtom* GetAtom3Ptr() const {return m_pAtom3;}
  RbtAtom* GetAtom4Ptr() const {return m_pAtom4;}
  //Add additional terms to the potential
  void AddTerm(const prms& dihprms);

 private:
  RbtDihedral();//Forbid default constructor
  RbtAtom* m_pAtom1;
  RbtAtom* m_pAtom2;
  RbtAtom* m_pAtom3;
  RbtAtom* m_pAtom4;
  vector<prms> m_prms;
};

//Useful typedefs
typedef vector<RbtDihedral*> RbtDihedralList;//Vector of smart pointers
typedef RbtDihedralList::iterator RbtDihedralListIter;
typedef RbtDihedralList::const_iterator RbtDihedralListConstIter;

class RbtDihedralSF : public virtual RbtBaseSF
{
 public:
  //Class type string
  static RbtString _CT;
  //Parameter names
  static RbtString _IMPL_H_CORR;//DM 6 Aug 2002 - option to correct barrier heights for implicit H's

  virtual ~RbtDihedralSF();
  
 protected:
  RbtDihedralSF();
  //Creates a vector of pointers to dihedral objects from a supplied bond list
  //Dihedral params are set using Tripos 5.2 atom types + params
  //NOTE: It is the responsibility of the subclass to delete the dihedral objects which are created
  RbtDihedralList CreateDihedralList(const RbtBondList& bondList);

 private:
  //Lookup the dihedral params for a given set of types
  RbtDihedral::prms FindDihedralParams(RbtTriposAtomType::eType t1,
				       RbtTriposAtomType::eType t2,
				       RbtTriposAtomType::eType t3,
				       RbtTriposAtomType::eType t4);

  //Determines the list of bonded atoms for pAtom1, and the list of dihedral offsets for any implicit hydrogens
  //pAtom2 is the other central atom in the bond, is excluded from the returned list of bonded atoms
  void CalcBondedAtoms(RbtAtom* pAtom1, RbtAtom* pAtom2, RbtAtomList&  bondedAtoms, RbtDoubleList& offsets);

  RbtParameterFileSourcePtr m_spDihedralSource;
  RbtStringList m_centralPairs;
  RbtTriposAtomType m_triposType;
};

#endif //_RBTDIHEDRALSF_H_
