/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtPdbFileSource.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Handles retrieval of atomic (coordinate) info
//from XPLOR PDB files (with segment name in last column).

#ifndef _RBTPDBFILESOURCE_H_
#define _RBTPDBFILESOURCE_H_

#include "RbtBaseMolecularFileSource.h"

class RbtPdbFileSource : public RbtBaseMolecularFileSource
{
 public:
  //Constructors
  RbtPdbFileSource(const char* fileName);
  RbtPdbFileSource(const RbtString fileName);

  //Default destructor
  virtual ~RbtPdbFileSource();

  ////////////////////////////////////////
  //Override public methods from RbtBaseMolecularDataSource
  virtual RbtBool isTitleListSupported() {return true;};
  virtual RbtBool isAtomListSupported() {return true;};
  virtual RbtBool isCoordinatesSupported() {return true;};
  virtual RbtBool isBondListSupported() {return false;};

 protected:
  //Pure virtual in RbtBaseFileSource - needs to be defined here
  virtual void Parse() throw (RbtError);

 private:
  //Private methods
  RbtPdbFileSource();//Disable default constructor
  RbtPdbFileSource(const RbtPdbFileSource&);//Copy constructor disabled by default
  RbtPdbFileSource& operator=(const RbtPdbFileSource&);//Copy assignment disabled by default
};

//useful typedefs
typedef SmartPtr<RbtPdbFileSource> RbtPdbFileSourcePtr;//Smart pointer

#endif //_RBTPDBFILESOURCE_H_
