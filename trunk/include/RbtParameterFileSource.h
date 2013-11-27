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

//File reader for Rbt parameter files

#ifndef _RBTPARAMETERFILESOURCE_H_
#define _RBTPARAMETERFILESOURCE_H_

#include "RbtBaseFileSource.h"
#include "RbtVariant.h"

class RbtParameterFileSource : public RbtBaseFileSource
{
 public:
  //Constructors
  RbtParameterFileSource(const char* fileName);
  RbtParameterFileSource(const RbtString& fileName);

  //Destructor
  virtual ~RbtParameterFileSource();
 
  ////////////////////////////////////////
  //Public methods
  RbtString GetTitle();
  RbtString GetVersion();
	//DM 06 June 2000 - limit #parameters to those in current section
  RbtUInt GetNumParameters();
  //DM 4 Feb 1999 Renamed from GetParameters()
  //RbtStringDoubleMap GetParameterMap();
  //DM 12 Feb 1999 - only return the list of parameter names, not their values
  //DM 06 Jun 2000 - limits params to those in named section
  RbtStringList GetParameterList();
  //DM 4 Feb 1999 Get a particular named parameter value as a double
  RbtDouble GetParameterValue(const RbtString& strParamName) throw (RbtError);
  //DM 12 Feb 1999 Get a particular named parameter value as a string
  RbtString GetParameterValueAsString(const RbtString& strParamName) throw (RbtError);
  //DM 11 Feb 1999 Check if parameter is present
  RbtBool isParameterPresent(const RbtString& strParamName);

  //DM 11 Feb 1999 - section handling
  //Parameters can be grouped into named sections
  //such that the same parameter name can appear in multiple sections
  //Default namespace is the unnamed section.
  //Main use is for simulation protocols which may need to define a variable number
  //of phases - e.g. high temperature sampling, cooling phase, minimisation phase
  //and need the same parameters to appear in each
  RbtInt GetNumSections();//Number of named sections
  RbtStringList GetSectionList();//List of section names
  RbtString GetSection() const;//Get current section name
  void SetSection(const RbtString& strSection="");//Set current section name
  

 protected:
  //Protected methods

 private:
  //Private methods
  RbtParameterFileSource();//Disable default constructor
  RbtParameterFileSource(const RbtParameterFileSource&);//Copy constructor disabled by default
  RbtParameterFileSource& operator=(const RbtParameterFileSource&);//Copy assignment disabled by default

  //Pure virtual in RbtBaseFileSource - needs to be defined here
  virtual void Parse() throw (RbtError);
  void ClearParamsCache();
  //Add a new section name
  //Throws an error if section name is empty, or is a duplicate
  void AddSection(const RbtString& strSection) throw (RbtError);
  //Returns the fully qualified parameter name (<section>::<parameter name>)
  //Checks if name already contains a section name, if so just returns the name unchanged
  //If not, prefixes the name with the current section name
  RbtString GetFullParameterName(const RbtString& strParamName);

 protected:
  //Protected data

 private:
  //Private data
  RbtString m_strTitle;
  RbtString m_strVersion;
  RbtStringVariantMap m_paramsMap;
  RbtStringList m_sectionNames;//List(vector) of section names
  RbtString m_strSection;//Current section
};

//useful typedefs
typedef SmartPtr<RbtParameterFileSource> RbtParameterFileSourcePtr;//Smart pointer

#endif //_RBTPARAMETERFILESOURCE_H_
