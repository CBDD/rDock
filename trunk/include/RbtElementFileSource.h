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

//File reader for Rbt elemental data file (atomic no,element name,vdw radii)
//Provides the data as a vector of structs

#ifndef _RBTELEMENTFILESOURCE_H_
#define _RBTELEMENTFILESOURCE_H_

#include "RbtBaseFileSource.h"

//Simple struct for holding the element data
class RbtElementData {
public:
	RbtElementData() : atomicNo(0),element(""),minVal(0),maxVal(0),commonVal(0),
		mass(0.0),vdwRadius(0.0) {}
	RbtInt atomicNo;
	RbtString element;
	RbtInt minVal;
	RbtInt maxVal;
	RbtInt commonVal;
	RbtDouble mass;
	RbtDouble vdwRadius;//Regular vdW radius
};

//Map with element data indexed by element name
typedef map<RbtString,RbtElementData> RbtStringElementDataMap;
typedef RbtStringElementDataMap::iterator RbtStringElementDataMapIter;
typedef RbtStringElementDataMap::const_iterator RbtStringElementDataMapConstIter;
//Map with element data indexed by atomic number
typedef map<RbtInt,RbtElementData> RbtIntElementDataMap;
typedef RbtIntElementDataMap::iterator RbtIntElementDataMapIter;
typedef RbtIntElementDataMap::const_iterator RbtIntElementDataMapConstIter;

class RbtElementFileSource : public RbtBaseFileSource
{
 public:
  //Constructors
  //RbtElementFileSource(const char* fileName);
  RbtElementFileSource(const RbtString& fileName);

  //Destructor
  virtual ~RbtElementFileSource();
 
  ////////////////////////////////////////
  //Public methods
  RbtString GetTitle();
  RbtString GetVersion();
  RbtUInt GetNumElements();
  RbtStringList GetElementNameList();//List of element names
  RbtIntList GetAtomicNumberList();//List of atomic numbers
  //Get element data for a given element name, throws error if not found
  RbtElementData GetElementData(const RbtString& strElementName) throw (RbtError);
  //Get element data for a given atomic number, throws error if not found
  RbtElementData GetElementData(RbtInt nAtomicNumber) throw (RbtError);
  //Check if given element name is present
  RbtBool isElementNamePresent(const RbtString& strElementName);
  //Check if given atomic number is present
  RbtBool isAtomicNumberPresent(RbtInt nAtomicNumber);

  //Get vdw radius increment for hydrogen-bonding donors
  RbtDouble GetHBondRadiusIncr();
  //Get vdw radius increment for atoms with implicit hydrogens
  RbtDouble GetImplicitRadiusIncr();

 protected:
  //Protected methods

 private:
  //Private methods
  RbtElementFileSource();//Disable default constructor
  RbtElementFileSource(const RbtElementFileSource&);//Copy constructor disabled by default
  RbtElementFileSource& operator=(const RbtElementFileSource&);//Copy assignment disabled by default

  //Pure virtual in RbtBaseFileSource - needs to be defined here
  virtual void Parse() throw (RbtError);
  void ClearElementDataCache();

 protected:
  //Protected data

 private:
  //Private data
  RbtString m_strTitle;
  RbtString m_strVersion;
  RbtDouble m_dHBondRadiusIncr;//Increment to add to vdW radius for H-bonding hydrogens
  RbtDouble m_dImplicitRadiusIncr;//Increment to add to vdW radius for atoms with implicit hydrogens
  RbtStringElementDataMap m_elementNameMap;
  RbtIntElementDataMap m_atomicNumberMap;
};

//useful typedefs
typedef SmartPtr<RbtElementFileSource> RbtElementFileSourcePtr;//Smart pointer

#endif //_RBTELEMENTFILESOURCE_H_
