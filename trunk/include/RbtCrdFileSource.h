/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtCrdFileSource.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Handles retrieval of atomic (coordinate) info from Charmm CRD files

#ifndef _RBTCRDFILESOURCE_H_
#define _RBTCRDFILESOURCE_H_

#include "RbtBaseMolecularFileSource.h"

class RbtCrdFileSource : public RbtBaseMolecularFileSource
{
 public:
  //Constructors
  RbtCrdFileSource(const char* fileName);
  RbtCrdFileSource(const RbtString fileName);

  //Default destructor
  virtual ~RbtCrdFileSource();

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
  RbtCrdFileSource();//Disable default constructor
  RbtCrdFileSource(const RbtCrdFileSource&);//Copy constructor disabled by default
  RbtCrdFileSource& operator=(const RbtCrdFileSource&);//Copy assignment disabled by default
};

//useful typedefs
typedef SmartPtr<RbtCrdFileSource> RbtCrdFileSourcePtr;//Smart pointer

#endif //_RBTCRDFILESOURCE_H_
