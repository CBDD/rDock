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

//Abstract base class (interface) for molecular data sources

#ifndef _RBTBASEMOLECULARDATASOURCE_H_
#define _RBTBASEMOLECULARDATASOURCE_H_

#include "RbtConfig.h"
#include "RbtAtom.h"
#include "RbtBond.h"
#include "RbtVariant.h"

class RbtBaseMolecularDataSource
{
 public:
  ////////////////////////////////////////
  //Constructors/destructors

  //Parameterised constructor
  //strName is a descriptive text label for the type of data source
  //Derived classes should use this constructor to set a suitable name
  RbtBaseMolecularDataSource(RbtString strName) : m_strName(strName) {}

  virtual ~RbtBaseMolecularDataSource() {} //Default destructor

  ////////////////////////////////////////
  //Public methods
  ////////////////
  RbtString GetName() {return m_strName;} //Returns source identifier string

  virtual RbtBool isTitleListSupported() {return false;};
  virtual RbtBool isAtomListSupported() {return false;};
  virtual RbtBool isCoordinatesSupported() {return false;};
  virtual RbtBool isBondListSupported() {return false;};

  //Pure virtual - derived classes must override
  
  //Reset source so that next time it is used, a new set of atom and bond objects
  //are created
  virtual void Reset() = 0;

  virtual RbtInt GetNumTitles() throw (RbtError) = 0;
  virtual RbtInt GetNumAtoms() throw (RbtError) = 0;
  virtual RbtInt GetNumBonds() throw (RbtError) = 0;
  virtual RbtInt GetNumSegments() throw (RbtError) = 0;

  virtual RbtStringList GetTitleList() throw (RbtError) = 0;
  virtual RbtAtomList GetAtomList() throw (RbtError) = 0;
  virtual RbtBondList GetBondList() throw (RbtError) = 0;
  virtual RbtSegmentMap GetSegmentMap() throw (RbtError) = 0;

  //DM 12 May 1999 - support for data records (e.g. SD file)
  //Does source support data records (default=false)
  virtual RbtBool isDataSupported() {return false;};
  //Get number of data fields
  virtual RbtInt GetNumData() throw (RbtError) = 0;
  //Get list of field names as string list
  virtual RbtStringList GetDataFieldList() throw (RbtError) = 0;
  //Get all data as map of key=field name, value=variant (double,string or string list)
  virtual RbtStringVariantMap GetDataMap() throw (RbtError) = 0;
  //Query as to whether a particular data field name is present
  virtual RbtBool isDataFieldPresent(const RbtString& strDataField) = 0;
  //Get a particular data value
  virtual RbtVariant GetDataValue(const RbtString& strDataField) throw (RbtError) = 0;


 private:
  ////////////////////////////////////////
  //Private methods
  /////////////////
  RbtBaseMolecularDataSource(); //Disable default constructor
  RbtBaseMolecularDataSource(const RbtBaseMolecularDataSource&);//Copy constructor disabled by default
  RbtBaseMolecularDataSource& operator=(const RbtBaseMolecularDataSource&);//Copy assignment disabled by default

 private:
  ////////////////////////////////////////
  //Private data
  //////////////
  RbtString m_strName;//Source identifier string
};

#endif //_RBTBASEMOLECULARDATASOURCE_H_
