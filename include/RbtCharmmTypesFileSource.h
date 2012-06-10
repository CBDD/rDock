/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtCharmmTypesFileSource.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//File reader for Charmm MASSES.RTF atom type file
//Provides the atom type list as a vector of structs

#ifndef _RBTCHARMMTYPESFILESOURCE_H_
#define _RBTCHARMMTYPESFILESOURCE_H_

#include "RbtBaseFileSource.h"

//Simple struct for holding the type info in masses.rtf
class CharmmType {
public:
	CharmmType() : nAtomType(0),strAtomType(""),mass(0.0),element(""),comment("") {
	}
	RbtInt nAtomType;
	RbtString strAtomType;
	RbtDouble mass;
	RbtString element;
	RbtString comment;
};

typedef vector<CharmmType> RbtCharmmTypeList;
typedef RbtCharmmTypeList::iterator RbtCharmmTypeListIter;
  
class RbtCharmmTypesFileSource : public RbtBaseFileSource
{
 public:
  //Constructors
  RbtCharmmTypesFileSource(const char* fileName);
  RbtCharmmTypesFileSource(const RbtString fileName);

  //Destructor
  virtual ~RbtCharmmTypesFileSource();
 
  ////////////////////////////////////////
  //Public methods
  RbtInt GetNumTypes();
  RbtCharmmTypeList GetTypeList();

 protected:
  //Protected methods

 private:
  //Private methods
  RbtCharmmTypesFileSource();//Disable default constructor
  RbtCharmmTypesFileSource(const RbtCharmmTypesFileSource&);//Copy constructor disabled by default
  RbtCharmmTypesFileSource& operator=(const RbtCharmmTypesFileSource&);//Copy assignment disabled by default

  //Pure virtual in RbtBaseFileSource - needs to be defined here
  virtual void Parse() throw (RbtError);
  void ClearTypesCache();

 protected:
  //Protected data

 private:
  //Private data
  RbtCharmmTypeList m_typesList;
};

//useful typedefs
typedef SmartPtr<RbtCharmmTypesFileSource> RbtCharmmTypesFileSourcePtr;//Smart pointer

#endif //_RBTCHARMMTYPESFILESOURCE_H_
