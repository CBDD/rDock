/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtSiteMapperFactory.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
