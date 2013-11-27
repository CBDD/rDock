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

//Factory class for creating docking site mapper objects.

#ifndef _RBTSITEMAPPERFACTORY_H_
#define _RBTSITEMAPPERFACTORY_H_

#include "RbtConfig.h"
#include "RbtSiteMapper.h"
#include "RbtParameterFileSource.h"

class RbtSiteMapperFactory
{
  //Parameter name which identifies a site mapper definition 
  static RbtString _MAPPER;
 public:
  ////////////////////////////////////////
  //Constructors/destructors
  
  RbtSiteMapperFactory(); //Default constructor  
  virtual ~RbtSiteMapperFactory(); //Default destructor
  
  
  ////////////////////////////////////////
  //Public methods
  ////////////////
  
  //Creates a single site mapper object of named class
  virtual RbtSiteMapper* Create(const RbtString& strMapperClass, const RbtString& strName) throw (RbtError);

  //Creates a single site mapper object whose class is defined
  //as the value of the SITE_MAPPER parameter in the strName section of the
  //parameter file source argument. Also sets the site mapper parameters from the
  //remaining parameter values in the current section
  //Note: the current section is restored to its original value upon exit
  virtual RbtSiteMapper* CreateFromFile(RbtParameterFileSourcePtr spPrmSource, const RbtString& strName) throw (RbtError);
  
  
 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////
  
  
 private:
  ////////////////////////////////////////
  //Private methods
  /////////////////
  
  RbtSiteMapperFactory(const RbtSiteMapperFactory&);//Copy constructor disabled by default  
  RbtSiteMapperFactory& operator=(const RbtSiteMapperFactory&);//Copy assignment disabled by default
  
  
 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////
  
  
 private:
  ////////////////////////////////////////
  //Private data
  //////////////
};

//Useful typedefs
typedef SmartPtr<RbtSiteMapperFactory> RbtSiteMapperFactoryPtr;//Smart pointer

#endif //_RBTSITEMAPPERFACTORY_H_
