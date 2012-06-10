/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtInteractionGrid.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Grid for indexing interactions involving more than 2 atoms

#ifndef _RBTINTERACTIONGRID_H_
#define _RBTINTERACTIONGRID_H_

#include "RbtBaseGrid.h"
#include "RbtAtom.h"

//simple container for up to 3 atoms, to hold one half of an interaction
//i.e. receptor atoms or ligand atoms.
class RbtInteractionCenter {
 public:
  enum eLP {NONE,PLANE,LONEPAIR};
  RbtInteractionCenter(RbtAtom* pAtom1=NULL, RbtAtom* pAtom2=NULL, RbtAtom* pAtom3=NULL, eLP LP=NONE) :
    m_pAtom1(pAtom1),m_pAtom2(pAtom2),m_pAtom3(pAtom3),m_LP(LP)
  {
    _RBTOBJECTCOUNTER_CONSTR_("RbtInteractionCenter");
  }
  ~RbtInteractionCenter()
  {
    _RBTOBJECTCOUNTER_DESTR_("RbtInteractionCenter");
  }
  RbtAtom* GetAtom1Ptr() const {return m_pAtom1;}
  RbtAtom* GetAtom2Ptr() const {return m_pAtom2;}
  RbtAtom* GetAtom3Ptr() const {return m_pAtom3;}
  eLP LP() const {return m_LP;}
  //Returns list of constituent atoms
  //(deconvolutes pseudoatoms into their constituent RbtAtom* lists)
  RbtAtomRList GetAtomList() const;
  RbtBool isSelected() const;

 private:
  //Could be a useful general function
  //If pAtom is a pseudo atom, then pushes all the constituent atoms onto list
  //else, push pAtom onto the list
  void AccumulateAtomList(const RbtAtom* pAtom, RbtAtomRList& atomList) const;
  RbtAtom* m_pAtom1;
  RbtAtom* m_pAtom2;
  RbtAtom* m_pAtom3;
  eLP m_LP;
};

typedef vector<RbtInteractionCenter*> RbtInteractionCenterList;//Vector of regular pointers
typedef RbtInteractionCenterList::iterator RbtInteractionCenterListIter;
typedef RbtInteractionCenterList::const_iterator RbtInteractionCenterListConstIter;

namespace Rbt {
  //Less than operator for sorting RbtInteractionCenter* by pointer value
  class InteractionCenterCmp {
  public:
    RbtBool operator()(const RbtInteractionCenter* pIC1, const RbtInteractionCenter* pIC2) const {
      return pIC1 < pIC2;
    }
  };
  //Is interaction center selected ?
  class isInteractionCenterSelected : public std::unary_function<RbtInteractionCenter*,RbtBool> {
  public:
    explicit isInteractionCenterSelected() {}
    RbtBool operator() (const RbtInteractionCenter* pIC) const {return pIC->isSelected();}
  };

  //Is the distance between interaction centers less than a given value ?
  //Function checks d**2 to save performing a sqrt
  class isInteractionD_lt : public std::unary_function<RbtInteractionCenter*,RbtBool> {
    RbtDouble d_sq;
    RbtAtom* a;
  public:
    explicit isInteractionD_lt(const RbtInteractionCenter* pIC1, RbtDouble dd) : d_sq(dd*dd),a(pIC1->GetAtom1Ptr()) {}
    RbtBool operator() (const RbtInteractionCenter* pIC2) const {return Rbt::Length2(pIC2->GetAtom1Ptr()->GetCoords(),a->GetCoords()) < d_sq;}
  };

  //Select/deselect the interaction center (selects all constituent atoms)
  class SelectInteractionCenter {
    RbtBool b;
  public:
    explicit SelectInteractionCenter(RbtBool bb) : b(bb) {}
    void operator() (RbtInteractionCenter* pIC);
  };
}

//A map of interaction centers indexed by unsigned int
//Used to store the receptor atom lists at each grid point
//DM 11 Jul 2000 - use map of regular RbtAtom* list (not RbtAtomPtr smart pointer list)

//DM 3 Nov 2000 - replace map by vector for faster lookup
//typedef map<RbtUInt,RbtInteractionCenterList> RbtInteractionListMap;
typedef vector<RbtInteractionCenterList> RbtInteractionListMap;

typedef RbtInteractionListMap::iterator RbtInteractionListMapIter;
typedef RbtInteractionListMap::const_iterator RbtInteractionListMapConstIter;

class RbtInteractionGrid : public RbtBaseGrid
{
 public:
	//Class type string
	static RbtString _CT;
  ////////////////////////////////////////
  //Constructors/destructors
  //Construct a NXxNYxNZ grid running from gridMin at gridStep resolution
  RbtInteractionGrid(const RbtCoord& gridMin, const RbtCoord& gridStep,
	     RbtUInt NX, RbtUInt NY, RbtUInt NZ, RbtUInt NPad=0);

  //Constructor reading all params from binary stream
  RbtInteractionGrid(istream& istr);

  ~RbtInteractionGrid(); //Default destructor

  //Copy constructor
  RbtInteractionGrid(const RbtInteractionGrid&);
  //Copy constructor taking base class argument
  RbtInteractionGrid(const RbtBaseGrid&);
  //Copy assignment
  RbtInteractionGrid& operator=(const RbtInteractionGrid&);
  //Copy assignment taking base class argument
  RbtInteractionGrid& operator=(const RbtBaseGrid&);

  ////////////////////////////////////////
  //Virtual functions for reading/writing grid data to streams in
  //text and binary format
  //Subclasses should provide their own private OwnPrint,OwnWrite, OwnRead
  //methods to handle subclass data members, and override the public
  //Print,Write and Read methods
	virtual void Print(ostream& ostr) const;//Text output
	virtual void Write(ostream& ostr) const;//Binary output (serialisation)
	virtual void Read(istream& istr);//Binary input, replaces existing grid

    ////////////////////////////////////////
  //Public methods
  ////////////////

  /////////////////////////
  //Get attribute functions
  /////////////////////////
  const RbtInteractionCenterList& GetInteractionList(RbtUInt iXYZ) const;
  const RbtInteractionCenterList& GetInteractionList(const RbtCoord& c) const;

  /////////////////////////
  //Set attribute functions
  /////////////////////////
  void SetInteractionLists(RbtInteractionCenter* pIntn, RbtDouble radius);
  void ClearInteractionLists();
  void UniqueInteractionLists();


 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////
  //Protected method for writing data members for this class to text stream
	void OwnPrint(ostream& ostr) const;
  //Protected method for writing data members for this class to binary stream
	void OwnWrite(ostream& ostr) const;
	//Protected method for reading data members for this class from binary stream
	void OwnRead(istream& istr) throw (RbtError);


 private:
  ////////////////////////////////////////
  //Private methods
  /////////////////
  RbtInteractionGrid(); //Disable default constructor

  //Helper function called by copy constructor and assignment operator
  void CopyGrid(const RbtInteractionGrid&);
  //DM 3 Nov 2000 - create InteractionListMap of the appropriate size
  void CreateMap();

 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////


 private:
  ////////////////////////////////////////
  //Private data
  //////////////
  RbtInteractionListMap m_intnMap;//Used to store the interaction center lists at each grid point
  const RbtInteractionCenterList m_emptyList;//Dummy list used by GetAtomList
};

//Useful typedefs
typedef SmartPtr<RbtInteractionGrid> RbtInteractionGridPtr;//Smart pointer

#endif //_RBTINTERACTIONGRID_H_
