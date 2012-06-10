/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtBaseMolecularFileSource.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Abstract base class for molecular file-based data sources

#ifndef _RBTBASEMOLECULARFILESOURCE_H_
#define _RBTBASEMOLECULARFILESOURCE_H_

#include "RbtBaseFileSource.h"
#include "RbtBaseMolecularDataSource.h"
#include "RbtParameterFileSource.h"

class RbtBaseMolecularFileSource : public RbtBaseFileSource,
				   public RbtBaseMolecularDataSource
{
 public:
  //Constructors
  //RbtBaseMolecularFileSource(const char* fileName, const char* sourceName);
  //Single-record version
  RbtBaseMolecularFileSource(const RbtString& fileName, const RbtString& sourceName);
//Multi-record version, with record delimiter passed as an argument
  RbtBaseMolecularFileSource(const RbtString& fileName, const RbtString& strRecDelim, const RbtString& sourceName);

  //Destructor
  virtual ~RbtBaseMolecularFileSource();

  //Pure virtual methods from RbtBaseMolecularDataSource
  void Reset();
  RbtInt GetNumTitles() throw (RbtError);
  RbtInt GetNumAtoms() throw (RbtError);
  RbtInt GetNumBonds() throw (RbtError);
  RbtInt GetNumSegments() throw (RbtError);

  RbtStringList GetTitleList() throw (RbtError);
  RbtAtomList GetAtomList() throw (RbtError);
  RbtBondList GetBondList() throw (RbtError);
  RbtSegmentMap GetSegmentMap() throw (RbtError);

  //DM 12 May 1999 - support for data records (e.g. SD file)
  //Get number of data fields
  virtual RbtInt GetNumData() throw (RbtError);
  //Get list of field names as string list
  virtual RbtStringList GetDataFieldList() throw (RbtError);
  //Get all data as map of key=field name, value=variant (double,string or string list)
  virtual RbtStringVariantMap GetDataMap() throw (RbtError);
  //Query as to whether a particular data field name is present
  virtual RbtBool isDataFieldPresent(const RbtString& strDataField);
  //Get a particular data value
  virtual RbtVariant GetDataValue(const RbtString& strDataField) throw (RbtError);

  //These methods allow filtering of the data source by segment name
  //So for example, we could just read the segment named TAR from the source
  //The filter modifies the behaviour of GetAtomList, GetBondList etc
  RbtSegmentMap GetSegmentFilterMap();
  void SetSegmentFilterMap(const RbtSegmentMap& segmentFilterMap);
  void ClearSegmentFilterMap();
  RbtBool isSegmentFilterMapDefined();
 
 protected:
  //Protected methods
  void ClearMolCache();//Clear the atom, bond and title lists
  void RemoveAtom(RbtAtomPtr spAtom);//DM 16 Feb 2000 - moved from RbtMdlFileSource
  void RenumberAtomsAndBonds();//DM 30 Oct 2000 - now independent from RemoveAtom
  //Moved from RbtMOL2FileSource - sets partial charges on a per-residue basis
  //Automatically determines protonation states of ASP/GLU/HIS using Mandatory and Forbidden
  //atoms in spParamSource
  void SetupPartialIonicGroups(RbtAtomList& atoms, RbtParameterFileSourcePtr spParamSource);

 private:
  //Private methods
  RbtBaseMolecularFileSource();//Disable default constructor
  RbtBaseMolecularFileSource(const RbtBaseMolecularFileSource&);//Copy constructor disabled by default
  RbtBaseMolecularFileSource& operator=(const RbtBaseMolecularFileSource&);//Copy assignment disabled by default

  //Get methods with segment filter applied
  RbtInt GetNumAtomsWithFilter();
  RbtInt GetNumBondsWithFilter();
  RbtAtomList GetAtomListWithFilter();
  RbtBondList GetBondListWithFilter();

 protected:
  //Protected data
  RbtStringList m_titleList;
  RbtAtomList m_atomList;
  RbtBondList m_bondList;
  RbtSegmentMap m_segmentMap;
  RbtStringVariantMap m_dataMap;//DM 12 May 1999 - for storing data records

 private:
  //Private data
  RbtSegmentMap m_segmentFilterMap;
};

//useful typedefs
typedef SmartPtr<RbtBaseMolecularFileSource> RbtMolecularFileSourcePtr;//Smart pointer

#endif //_RBTBASEMOLECULARFILESOURCE_H_
