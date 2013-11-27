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

//Grid for indexing lipophilic interactions

#ifndef _RBTNONBONDEDGRID_H_
#define _RBTNONBONDEDGRID_H_

#include "RbtBaseGrid.h"
#include "RbtAtom.h"

//A map of atom vectors indexed by unsigned int
//Used to store the receptor atom lists at each grid point
//DM 11 Jul 2000 - use map of regular RbtAtom* list (not RbtAtomPtr smart pointer list)

//DM 3 Nov 2000 - replace map by vector for faster lookup
//typedef map<RbtUInt,RbtAtomRList> RbtAtomListMap;
typedef vector<RbtAtomRList> RbtAtomListMap;

typedef RbtAtomListMap::iterator RbtAtomListMapIter;
typedef RbtAtomListMap::const_iterator RbtAtomListMapConstIter;

class RbtNonBondedGrid : public RbtBaseGrid
{
 public:
	//Class type string
	static RbtString _CT;
  ////////////////////////////////////////
  //Constructors/destructors
  //Construct a NXxNYxNZ grid running from gridMin at gridStep resolution
  RbtNonBondedGrid(const RbtCoord& gridMin, const RbtCoord& gridStep,
	     RbtUInt NX, RbtUInt NY, RbtUInt NZ, RbtUInt NPad=0);

  //Constructor reading all params from binary stream
  RbtNonBondedGrid(istream& istr);

  ~RbtNonBondedGrid(); //Default destructor

  //Copy constructor
  RbtNonBondedGrid(const RbtNonBondedGrid&);
  //Copy constructor taking base class argument
  RbtNonBondedGrid(const RbtBaseGrid&);
  //Copy assignment
  RbtNonBondedGrid& operator=(const RbtNonBondedGrid&);
  //Copy assignment taking base class argument
  RbtNonBondedGrid& operator=(const RbtBaseGrid&);

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
  //RbtAtomList GetAtomList(RbtUInt iXYZ) const;
  //RbtAtomList GetAtomList(const RbtCoord& c) const;
  const RbtAtomRList& GetAtomList(RbtUInt iXYZ) const;
  const RbtAtomRList& GetAtomList(const RbtCoord& c) const;

  /////////////////////////
  //Set attribute functions
  /////////////////////////
  void SetAtomLists(RbtAtom* pAtom, RbtDouble radius);
  void ClearAtomLists();
  void UniqueAtomLists();


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
  RbtNonBondedGrid(); //Disable default constructor

  //Helper function called by copy constructor and assignment operator
  void CopyGrid(const RbtNonBondedGrid&);
  //DM 6 Nov 2000 - create AtomListMap of the appropriate size
  void CreateMap();

 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////


 private:
  ////////////////////////////////////////
  //Private data
  //////////////
  RbtAtomListMap m_atomMap;//Used to store the receptor atom lists at each grid point
  const RbtAtomRList m_emptyList;//Dummy list used by GetAtomList
};

//Useful typedefs
typedef SmartPtr<RbtNonBondedGrid> RbtNonBondedGridPtr;//Smart pointer

#endif //_RBTNONBONDEDGRID_H_
