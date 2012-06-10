/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtPseudoAtom.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Derived class for pseudo-atoms.

#ifndef _RBTPSEUDOATOM_H_
#define _RBTPSEUDOATOM_H_

#include "RbtAtom.h"

class RbtPseudoAtom : public RbtAtom
{
 public:
  ///////////////////////////////////////////////
  //Constructors / destructors

  //Constructor supplying atom list of constituent atoms
  RbtPseudoAtom(const RbtAtomList& atomList, RbtInt nAtomId = -1, RbtString strAtomName = "PSEUDO") throw (RbtError);

  //Default destructor
  virtual ~RbtPseudoAtom();

  ///////////////////////////////////////////////
  //Public accessor functions
  //Virtual function for dumping pseudoatom details to an output stream
  virtual ostream& Print(ostream& s) const;

  ///////////////////////////////////////////////
  //Public accessor functions
  //2-D attributes

  //Derived class methods for returning the constituent atom list
  RbtUInt GetNumAtoms() const {return m_atomList.size();}
  RbtAtomList GetAtomList() const {return m_atomList;}

  ///////////////////////////////////////////////
  //Public accessor functions
  //3-D attributes

  //Coords
  //DM 8 Dec 1998 - override base RbtAtom methods
  //virtual RbtCoord GetCoords() const;
  //virtual RbtDouble GetX() const;
  //virtual RbtDouble GetY() const;
  //virtual RbtDouble GetZ() const;
  //
  //DM 11 Jul 2000 - base GetCoords method is now non-virtual again
  //Pseudoatoms must update their mean coords by calling UpdateCoords
  //each time the underlying atoms are moved
  void UpdateCoords();

  ///////////////////////////////////////////////

  //Other public methods


 protected:

 private:
  //Private methods
  RbtPseudoAtom();//Disable default constructor
  RbtPseudoAtom(const RbtPseudoAtom& pseudoAtom);//Disable copy constructor
  RbtPseudoAtom& operator=(const RbtPseudoAtom& pseudoAtom);//Disable copy assignment


 private:
  //Private data
  RbtAtomList m_atomList; //List of smart pointers to constituent atoms
};

//Useful typedefs
typedef SmartPtr<RbtPseudoAtom> RbtPseudoAtomPtr;//Smart pointer
typedef vector<RbtPseudoAtomPtr> RbtPseudoAtomList;//Vector of smart pointers
typedef RbtPseudoAtomList::iterator RbtPseudoAtomListIter;
typedef RbtPseudoAtomList::const_iterator RbtPseudoAtomListConstIter;

#endif //_RBTPSEUDOATOM_H_
