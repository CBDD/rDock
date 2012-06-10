/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtNmrRestraintFileSource.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//File reader for external NMR restraints. Supports NOE distance restraints
//and STD ligand-receptor contact surface restraints.
//NOE Syntax:
//<ATOM_LIST_1> <ATOM_LIST_2> <MAX DISTANCE>
//
//STD Syntax:
//STD <ATOM_LIST_1> <MAX DISTANCE> (this is the literal string "STD")
//
//Comment lines are supported (# first character)
//Atom names in list should be comma-separated, with NO SPACES
//
//For both types of restraint, brackets can be used around each atom list to
//denote various levels of ambiguity:
//
//(ATOM_LIST) = apply restraint to mean coords of atoms in list (MEAN)
//[ATOM_LIST] = apply restraint to all atoms in list simultaneously (AND)
//ATOM_LIST = restraint can be satisfied by any atom in list  (OR)
//
//e.g.
//MOL_1:C3 (_22:C5,_23:N1) 5.0 #MOL_1:C3 within 5A of the MEAN of _22:C5,_23:N1
//MOL_1:C3 [_22:C5,_23:N1] 5.0 #MOL_1:C3 within 5A of _22:C5 AND 5A of _23:N1
//MOL_1:C3 _22:C5,_23:N1   5.0 #MOL_1:C3 within 5A of _22:C5 OR 5A of _23:N1
//STD MOL_1:C3 5.0             #MOL_1:C3 within 5A of any part of the receptor

#ifndef _RBTNMRRESTRAINTFILESOURCE_H_
#define _RBTNMRRESTRAINTFILESOURCE_H_

#include "RbtBaseFileSource.h"
#include "RbtNoeRestraint.h"

class RbtNmrRestraintFileSource : public RbtBaseFileSource
{
 public:
  ////////////////////////////////////////
  //Constructors/destructors
  RbtNmrRestraintFileSource(const RbtString& fileName);

  //Destructor
  virtual ~RbtNmrRestraintFileSource();


  ////////////////////////////////////////
  //Public methods
  ////////////////
  RbtUInt GetNumNoeRestraints();//Number of NOE restraints in file
  RbtNoeRestraintNamesList GetNoeRestraintList();//List of NOE restraints
  RbtUInt GetNumStdRestraints();//Number of STD restraints in file
  RbtStdRestraintNamesList GetStdRestraintList();//List of STD restraints

 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////


 private:
  ////////////////////////////////////////
  //Private methods
  /////////////////

  RbtNmrRestraintFileSource(); //Disable default constructor
  RbtNmrRestraintFileSource(const RbtNmrRestraintFileSource&);//Copy constructor disabled by default
  RbtNmrRestraintFileSource& operator=(const RbtNmrRestraintFileSource&);//Copy assignment disabled by default

  //Pure virtual in RbtBaseFileSource - needs to be defined here
  virtual void Parse() throw (RbtError);
  void ClearRestraintCache();

  //Returns NOE restraint type and modifies the atom name string accordingly
  //Returns UNDEFINED if the atom name string has bad syntax
  Rbt::eNoeType NoeRestraintType(RbtString& strAtomNames);


 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////


 private:
  ////////////////////////////////////////
  //Private data
  //////////////
  RbtNoeRestraintNamesList m_noeRestraintList;
  RbtStdRestraintNamesList m_stdRestraintList;
};

//Useful typedefs
typedef SmartPtr<RbtNmrRestraintFileSource> RbtNmrRestraintFileSourcePtr;//Smart pointer
typedef vector<RbtNmrRestraintFileSourcePtr> RbtNmrRestraintFileSourceList;//Vector of smart pointers
typedef RbtNmrRestraintFileSourceList::iterator RbtNmrRestraintFileSourceListIter;
typedef RbtNmrRestraintFileSourceList::const_iterator RbtNmrRestraintFileSourceListConstIter;

#endif //_RBTNMRRESTRAINTFILESOURCE_H_
