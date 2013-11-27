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
