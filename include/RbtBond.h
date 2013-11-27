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

//Handles all the properties of a bond

#ifndef _RBTBOND_H_
#define _RBTBOND_H_

#include "RbtConfig.h"
#include "RbtAtom.h"

class RbtBond
{
 public:
  ////////////////////////////////////////
  //Constructors/destructors
  ////////////////////////////////////////

  //Default constructor
  RbtBond();

  //Cconstructor supplying all parameters
  RbtBond(RbtInt nBondId,
	  RbtAtomPtr& spAtom1,
	  RbtAtomPtr& spAtom2,
	  RbtInt nFormalBondOrder=1);

  virtual ~RbtBond(); //Default destructor

  RbtBond(const RbtBond& bond);//Copy constructor

  RbtBond& operator=(const RbtBond& bond);//Copy assignment


  ///////////////////////////////////////////////
  //Stream functions
  ///////////////////////////////////////////////

  //Insertion operator (primarily for debugging)
  //Note: needs to be friend so can access bond private data
  //without using the accessor functions
  friend ostream& operator<<(ostream& s, const RbtBond& bond);


  ////////////////////////////////////////
  //Public accessor functions
  ////////////////

  //Bond ID
  RbtInt GetBondId() const {return m_nBondId;}
  void SetBondId(const RbtInt nBondId) {m_nBondId = nBondId;}

  //Atom pointers
  RbtAtomPtr GetAtom1Ptr() const {return m_spAtom1;}
  RbtAtomPtr GetAtom2Ptr() const {return m_spAtom2;}
  void SetAtom1Ptr(RbtAtomPtr& spAtom1) {m_spAtom1 = spAtom1;}
  void SetAtom2Ptr(RbtAtomPtr& spAtom2) {m_spAtom2 = spAtom2;}

  //Formal Bond order
  RbtInt GetFormalBondOrder() const {return m_nFormalBondOrder;}
  void SetFormalBondOrder(const RbtInt nFormalBondOrder) {m_nFormalBondOrder = nFormalBondOrder;}

  //Partial Bond order
  RbtDouble GetPartialBondOrder() const {return m_dPartialBondOrder;}
  void SetPartialBondOrder(const RbtDouble dPartialBondOrder) {m_dPartialBondOrder = dPartialBondOrder;}

  //CyclicFlag - flag to indicate bond is in a ring (set by RbtModel::FindRing)
  RbtBool GetCyclicFlag() const {return m_bCyclic;}
  void SetCyclicFlag(RbtBool bCyclic=true) {m_bCyclic = bCyclic;}

  //SelectionFlag - general purpose flag can be set/cleared by various search algorithms (e.g. FindRings)
  RbtBool GetSelectionFlag() const {return m_bSelected;}
  void SetSelectionFlag(RbtBool bSelected=true) {m_bSelected = bSelected;}


  ////////////////////////////////////////
  //Public methods
  ////////////////
  //Returns bond length
  RbtDouble Length() const;

 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////


 private:
  ////////////////////////////////////////
  //Private methods
  /////////////////

 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////


 private:
  ////////////////////////////////////////
  //Private data
  //////////////
  RbtInt m_nBondId;//Original bond ID in PSF file
  RbtAtomPtr m_spAtom1;//Smart pointer to atom 1
  RbtAtomPtr m_spAtom2;//Smart pointer to atom 2
  RbtInt m_nFormalBondOrder;//Formal bond order (1,2,3, no aromatic bond orders of 1.5)
  RbtDouble m_dPartialBondOrder;//Partial bond order (1.0, 1.5, 2.0, 3.0 etc)
  RbtBool m_bCyclic;//Is the bond in a ring ?
  RbtBool m_bSelected; //Can be set/cleared by various search algorithms (e.g. FindRings)
};

//Useful typedefs
typedef SmartPtr<RbtBond> RbtBondPtr;//Smart pointer
typedef vector<RbtBondPtr> RbtBondList;//Vector of smart pointers
typedef RbtBondList::iterator RbtBondListIter;
typedef RbtBondList::const_iterator RbtBondListConstIter;


///////////////////////////////////////////////
// Non-member functions (in Rbt namespace)
//////////////////////////////////////////
namespace Rbt
{
  ////////////////////////////////////////////
  //Predicates (for use by STL algorithms)
  ////////////////////////////////////////////
  //DM 08 May 2002 - convert to regular RbtBond* functions
  //More universal as RbtBondPtr parameters will be automatically degraded to
  //regular pointers anyway
  typedef std::unary_function<RbtBond*,RbtBool> RbtBondUnaryPredicate;
  typedef std::binary_function<RbtBond*,RbtBond*,RbtBool> RbtBondBinaryPredicate;

  //Is bond selected ?
  class isBondSelected : public RbtBondUnaryPredicate {
  public:
    explicit isBondSelected() {}
    RbtBool operator() (RbtBond* pBond) const {return pBond->GetSelectionFlag();}
  };

  //Is bond cyclic ?
  class isBondCyclic : public RbtBondUnaryPredicate {
  public:
    explicit isBondCyclic() {}
    RbtBool operator() (RbtBond* pBond) const {return pBond->GetCyclicFlag();}
  };

  //Is bond rotatable ?
  class isBondRotatable : public RbtBondUnaryPredicate {
  public:
    explicit isBondRotatable() {}
    RbtBool operator() (RbtBond*) const;
  };

  //Is bond to a terminal NH3+ group?
  class isBondToNH3 : public RbtBondUnaryPredicate {
  public:
    explicit isBondToNH3() {}
    RbtBool operator() (RbtBond*) const;
  };

  //Is bond to a terminal OH group?
  class isBondToOH : public RbtBondUnaryPredicate {
  public:
    explicit isBondToOH() {}
    RbtBool operator() (RbtBond* pBond) const;
  };
  
  //DM 1 April 1999
  //Is bond2 equal to bond1 (checks if underlying regular pointers match)
  //Note: this is a binary rather than unary predicate
  //class isBondPtr_eq : public RbtBondBinaryPredicate {
  //public:
  //  explicit isBondPtr_eq() {}
  //  RbtBool operator() (RbtBondPtr spBond1, RbtBondPtr spBond2) const {return spBond1.Ptr() == spBond2.Ptr();}
  //};

  //Is bond.Ptr() equal to RbtBond* (checks if underlying regular pointers match)
  class isBondPtr_eq : public RbtBondUnaryPredicate {
    RbtBond* p;
  public:
    explicit isBondPtr_eq(RbtBond* pp) : p(pp) {}
    RbtBool operator() (RbtBond* pBond) const {return pBond == p;}
  };

  //DM 2 Aug 1999
  //Is bond2 equal to bond1 (checks if bond ID, atom1 and atom2 match)
  class isBond_eq : public RbtBondUnaryPredicate {
    RbtBond* p;
    Rbt::isAtom_eq bIsAtomEqual;
  public:
    explicit isBond_eq(RbtBond* pp) : p(pp) {}
    RbtBool operator() (RbtBond* pBond) const {
      return ( (pBond->GetBondId() == p->GetBondId()) &&
	       (bIsAtomEqual(pBond->GetAtom1Ptr(),p->GetAtom1Ptr())) &&
	       (bIsAtomEqual(pBond->GetAtom2Ptr(),p->GetAtom2Ptr())) );
    }
  };

  //DM 7 June 1999
  //Is bond an amide bond?
  class isBondAmide : public RbtBondUnaryPredicate {
  public:
    explicit isBondAmide() {}
    RbtBool operator() (RbtBond*) const;
  };
  
  
  ////////////////////////////////////////////
  //Bond list functions (implemented as STL algorithms)
  ////////////////////////////////////////////

  //Generic template version of GetNumBonds, passing in your own predicate
  template<class Predicate> RbtUInt GetNumBonds(const RbtBondList& bondList, const Predicate& pred) {
    return std::count_if(bondList.begin(),bondList.end(),pred);
  }

  //Generic template version of GetBondList, passing in your own predicate
  template<class Predicate> RbtBondList GetBondList(const RbtBondList& bondList, const Predicate& pred) {
    RbtBondList newBondList;
    std::copy_if(bondList.begin(),bondList.end(),std::back_inserter(newBondList),pred);
    return newBondList;
  }

  //Generic template version of FindBond, passing in your own predicate
  template<class Predicate> RbtBondListIter FindBond(RbtBondList& bondList, const Predicate& pred) {
    return std::find_if(bondList.begin(),bondList.end(),pred);
  }

  //Selected bonds
  void SetBondSelectionFlags(RbtBondList& bondList,RbtBool bSelected=true);
  inline RbtUInt GetNumSelectedBonds(const RbtBondList& bondList) {
    return GetNumBonds(bondList,Rbt::isBondSelected());
  }

  inline RbtBondList GetSelectedBondList(const RbtBondList& bondList) {
    return GetBondList(bondList,Rbt::isBondSelected());
  }

  //Cyclic bonds
  void SetBondCyclicFlags(RbtBondList& bondList,RbtBool bCyclic=true);
  inline RbtUInt GetNumCyclicBonds(const RbtBondList& bondList) {
    return GetNumBonds(bondList,Rbt::isBondCyclic());
  }
  inline RbtBondList GetCyclicBondList(const RbtBondList& bondList) {
    return GetBondList(bondList,Rbt::isBondCyclic());
  }

  //Rotatable bonds
  inline RbtUInt GetNumRotatableBonds(const RbtBondList& bondList) {
    return GetNumBonds(bondList,Rbt::isBondRotatable());
  }
  inline RbtBondList GetRotatableBondList(const RbtBondList& bondList) {
    return GetBondList(bondList,Rbt::isBondRotatable());
  }
}

#endif //_RBTBOND_H_
