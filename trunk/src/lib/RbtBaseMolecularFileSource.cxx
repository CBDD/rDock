/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtBaseMolecularFileSource.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtBaseMolecularFileSource.h"
#include "RbtFileError.h"

//Constructors
//RbtBaseMolecularFileSource::RbtBaseMolecularFileSource(const char* fileName, const char* sourceName) :
//  RbtBaseFileSource(fileName), RbtBaseMolecularDataSource(sourceName)
//{
//  _RBTOBJECTCOUNTER_CONSTR_("RbtBaseMolecularFileSource");
//}

//Single-record version
RbtBaseMolecularFileSource::RbtBaseMolecularFileSource(const RbtString& fileName, const RbtString& sourceName) :
  RbtBaseFileSource(fileName), RbtBaseMolecularDataSource(sourceName)
{
  _RBTOBJECTCOUNTER_CONSTR_("RbtBaseMolecularFileSource");
}

//Multi-record version, with record delimiter passed as an argument
RbtBaseMolecularFileSource::RbtBaseMolecularFileSource(const RbtString& fileName, const RbtString& strRecDelim,
						       const RbtString& sourceName) :
  RbtBaseFileSource(fileName,strRecDelim), RbtBaseMolecularDataSource(sourceName)
{
  _RBTOBJECTCOUNTER_CONSTR_("RbtBaseMolecularFileSource");
}

//Destructors
RbtBaseMolecularFileSource::~RbtBaseMolecularFileSource()
{
  ClearMolCache();//Tidy up the atom, bond, title and segment lists
  _RBTOBJECTCOUNTER_DESTR_("RbtBaseMolecularFileSource");
}

//Public methods from RbtBaseMolecularDataSource

//Reset source so that next time it is used, a new set of atom and bond objects
//are created
void RbtBaseMolecularFileSource::Reset()
{
  ClearMolCache();
}

RbtInt RbtBaseMolecularFileSource::GetNumTitles() throw (RbtError)
{
  if (isTitleListSupported()) {
    Parse();
    return m_titleList.size();
  }
  else
    throw RbtInvalidRequest(_WHERE_,"Title list not supported by " + GetName());
}

RbtInt RbtBaseMolecularFileSource::GetNumAtoms() throw (RbtError)
{
  if (isAtomListSupported()) {
    Parse();
    if (isSegmentFilterMapDefined())
      return GetNumAtomsWithFilter();
    else
      return m_atomList.size();
  }
  else
    throw RbtInvalidRequest(_WHERE_,"Atom list not supported by " + GetName());
}

RbtInt RbtBaseMolecularFileSource::GetNumBonds() throw (RbtError)
{
  if (isBondListSupported()) {
    Parse();
    if (isSegmentFilterMapDefined())
      return GetNumBondsWithFilter();
    else
      return m_bondList.size();
  }
  else
    throw RbtInvalidRequest(_WHERE_,"Bond list not supported by " + GetName());
}


RbtInt RbtBaseMolecularFileSource::GetNumSegments() throw (RbtError)
{
  if (isAtomListSupported()) {
    Parse();
    return m_segmentMap.size();
  }
  else
    throw RbtInvalidRequest(_WHERE_,"Atom list not supported by " + GetName());
}

RbtStringList RbtBaseMolecularFileSource::GetTitleList() throw (RbtError)
{
  if (isTitleListSupported()) {
    Parse();
    return m_titleList;
  }
  else
    throw RbtInvalidRequest(_WHERE_,"Title list not supported by " + GetName());
}

RbtAtomList RbtBaseMolecularFileSource::GetAtomList() throw (RbtError)
{
  if (isAtomListSupported()) {
    Parse();
    if (isSegmentFilterMapDefined())
      return GetAtomListWithFilter();
    else
      return m_atomList;
  }
  else
    throw RbtInvalidRequest(_WHERE_,"Atom list not supported by " + GetName());
}


RbtBondList RbtBaseMolecularFileSource::GetBondList() throw (RbtError)
{
  if (isBondListSupported()) {
    Parse();
    if (isSegmentFilterMapDefined())
      return GetBondListWithFilter();
    else
      return m_bondList;
  }
  else
    throw RbtInvalidRequest(_WHERE_,"Bond list not supported by " + GetName());
}


RbtSegmentMap RbtBaseMolecularFileSource::GetSegmentMap() throw (RbtError)
{
  if (isAtomListSupported()) {
    Parse();
    return m_segmentMap;
  }
  else
    throw RbtInvalidRequest(_WHERE_,"Atom list not supported by " + GetName());
}
 
//DM 12 May 1999 - support for data records (e.g. SD file)
//Get number of data fields
RbtInt RbtBaseMolecularFileSource::GetNumData() throw (RbtError)
{
  if (isDataSupported()) {
    Parse();
    return m_dataMap.size();
  }
  else
    throw RbtInvalidRequest(_WHERE_,"Data records not supported by " + GetName());
}

//Get list of field names as string list
RbtStringList RbtBaseMolecularFileSource::GetDataFieldList() throw (RbtError)
{
  if (isDataSupported()) {
    Parse();
    RbtStringList dataFieldList;
    dataFieldList.reserve(m_dataMap.size());
    for (RbtStringVariantMapConstIter iter = m_dataMap.begin(); iter != m_dataMap.end(); iter++)
      dataFieldList.push_back((*iter).first);
    return dataFieldList;
  }
  else
    throw RbtInvalidRequest(_WHERE_,"Data records not supported by " + GetName());
}

//Get all data as map of key=field name, value=variant (double,string or string list)
RbtStringVariantMap RbtBaseMolecularFileSource::GetDataMap() throw (RbtError)
{
  if (isDataSupported()) {
    Parse();
    return m_dataMap;
  }
  else
    throw RbtInvalidRequest(_WHERE_,"Data records not supported by " + GetName());
}

//Query as to whether a particular data field name is present
RbtBool RbtBaseMolecularFileSource::isDataFieldPresent(const RbtString& strDataField)
{
  if (isDataSupported()) {
    Parse();
    return m_dataMap.find(strDataField) != m_dataMap.end();
  }
  else
    throw RbtInvalidRequest(_WHERE_,"Data records not supported by " + GetName());
}

//Get a particular data value
RbtVariant RbtBaseMolecularFileSource::GetDataValue(const RbtString& strDataField) throw (RbtError)
{
  if (isDataSupported()) {
    Parse();
    RbtStringVariantMapConstIter iter = m_dataMap.find(strDataField);
    if (iter != m_dataMap.end())
      return (*iter).second;
    else
      throw RbtFileMissingParameter(_WHERE_,strDataField + " data field not found in " + GetFileName());
  }
  else
    throw RbtInvalidRequest(_WHERE_,"Data records not supported by " + GetName());
}


//These methods allow filtering of the data source by segment name
//So for example, we could just read the segment named TAR from the source
//The filter modifies the behaviour of getNumAtoms(), GetAtomList,
//GetNumBonds(), GetBondList()
RbtSegmentMap RbtBaseMolecularFileSource::GetSegmentFilterMap()
{
  return m_segmentFilterMap;
}

void RbtBaseMolecularFileSource::SetSegmentFilterMap(const RbtSegmentMap& segmentFilterMap)
{
  m_segmentFilterMap = segmentFilterMap;
}

void RbtBaseMolecularFileSource::ClearSegmentFilterMap()
{
  m_segmentFilterMap.clear();
}

RbtBool RbtBaseMolecularFileSource::isSegmentFilterMapDefined()
{
  return !m_segmentFilterMap.empty();
}

////////////////////////////////////////////////////////////////////
//Protected methods

//Clear the atom, bond, title and segment lists
void RbtBaseMolecularFileSource::ClearMolCache()
{
  m_titleList.clear();
  m_atomList.clear();
  m_bondList.clear();
  m_segmentMap.clear();
  m_dataMap.clear();
  m_bParsedOK = false;
}

//DM 16 Feb 2000 - formerly in RbtMdlFileSource
//Removes an atom and all bonds to the atom, from the atom and bond lists
//Does not adjust the attributes of the atoms whose bonds have been broken
//DM 30 Oct 2000 - no longer renumbers the atoms and bonds
//RenumberAtomsAndBonds should be called after all atoms have been removed
void RbtBaseMolecularFileSource::RemoveAtom(RbtAtomPtr spAtom)
{
	//DM 31 Oct 2000 - return by const ref
  const RbtBondMap& bondMap = spAtom->GetBondMap();

  //First remove all bonds from the atom
  for (RbtBondMapConstIter mapIter = bondMap.begin(); mapIter != bondMap.end(); mapIter++) {
    RbtBondListIter bIter = Rbt::FindBond(m_bondList,Rbt::isBond_eq((*mapIter).first));
    if (bIter != m_bondList.end()) {
#ifdef _DEBUG
      cout << "Removing bond #" << (*bIter)->GetBondId()
	   << " (" << (*bIter)->GetAtom1Ptr()->GetAtomName()
	   << "-" << (*bIter)->GetAtom2Ptr()->GetAtomName()
	   << ")" << endl;
#endif //_DEBUG
      m_bondList.erase(bIter);//Erase the bond
    }
  }

  //Now we have an isolated atom we can remove it
  //Find the atom in the FileSource atom list
  //DM 2 Aug 1999 - search by atom  attributes, not by memory location (v risky)
  //RbtAtomListIter aIter = Rbt::FindAtom(m_atomList,std::bind2nd(Rbt::isAtomPtr_eq(),spAtom));
  RbtAtomListIter aIter = Rbt::FindAtom(m_atomList,std::bind2nd(Rbt::isAtom_eq(),spAtom));
  if (aIter != m_atomList.end()) {
#ifdef _DEBUG
    cout << "Removing atom #" << (*aIter)->GetAtomId() << ", " << (*aIter)->GetAtomName() << endl;
#endif //_DEBUG
    m_atomList.erase(aIter);//Erase the atom
  }
}

//DM 30 Oct 2000 - now independent from RemoveAtom
void RbtBaseMolecularFileSource::RenumberAtomsAndBonds() {
  //Renumber the bond and atom IDs from 1
  RbtInt nAtomId(1);
  for (RbtAtomListIter aIter = m_atomList.begin(); aIter != m_atomList.end(); aIter++,nAtomId++) {
    (*aIter)->SetAtomId(nAtomId);
  }
  RbtInt nBondId(1);
  for (RbtBondListIter bIter = m_bondList.begin(); bIter != m_bondList.end(); bIter++,nBondId++) {
    (*bIter)->SetBondId(nBondId);
  }
}

////////////////////////////////////////////////////////////////////
//Private methods

//This version of GetNumAtoms applies the segment filter
//Called by the public GetNumAtoms
//Assume that Parse has already been called
RbtInt RbtBaseMolecularFileSource::GetNumAtomsWithFilter()
{
  RbtInt nAtoms(0);
  //The segment map data member (m_segmentMap) already contains the atom counts for each segment present
  //but we need to iterate over the segment filter map to see how many of those requested are
  //actually present
  for (RbtSegmentMapIter iter = m_segmentFilterMap.begin(); iter != m_segmentFilterMap.end(); iter++) {
    RbtSegmentMapIter match = m_segmentMap.find((*iter).first);//try to find a match
    if (match != m_segmentMap.end())
      nAtoms += (*match).second; //The value (.second) of the map element is the number of atoms in the segment
  }	
  return nAtoms;
}

//This version of GetNumBonds applies the segment filter
//Called by the public GetNumBonds
//Assume that Parse has already been called
RbtInt RbtBaseMolecularFileSource::GetNumBondsWithFilter()
{
  RbtInt nBonds(0);
  //More tricky this one, there's no option but to iterate over all the bonds to find if the atoms
  //are in one of the segments in the filter map
  for (RbtBondListIter iter = m_bondList.begin(); iter != m_bondList.end(); iter++) {
    if ( (m_segmentFilterMap.find( (*iter)->GetAtom1Ptr()->GetSegmentName() ) != m_segmentFilterMap.end()) &&
	 (m_segmentFilterMap.find( (*iter)->GetAtom2Ptr()->GetSegmentName() ) != m_segmentFilterMap.end()) )
      nBonds++;
  }
  return nBonds;
}

//This version of GetAtomList applies the segment filter
//Called by the public GetAtomList
//Assume that Parse has already been called
RbtAtomList RbtBaseMolecularFileSource::GetAtomListWithFilter()
{
  RbtAtomList atomList;
  //Loop over all atoms and check whether they belong to one of the requested segments
  //DM 29 Jul 1999 - set consecutive atom ID's starting from 1, irrespective of ID in file
  RbtInt nAtomId(0);
  for (RbtAtomListIter iter = m_atomList.begin(); iter != m_atomList.end(); iter++) {
    if ( m_segmentFilterMap.find((*iter)->GetSegmentName()) != m_segmentFilterMap.end() ) {
      (*iter)->SetAtomId(++nAtomId);
      atomList.push_back(*iter);//Got a match, so store the atom pointer
    }
  }
  return atomList;
}

//This version of GetBondList applies the segment filter
//Called by the public GetBondList
//Assume that Parse has already been called
RbtBondList RbtBaseMolecularFileSource::GetBondListWithFilter()
{
  RbtBondList bondList;
  //Loop over all bonds and check whether their atoms  belong to one of the requested segments
  for (RbtBondListIter iter = m_bondList.begin(); iter != m_bondList.end(); iter++) {
    if ( (m_segmentFilterMap.find( (*iter)->GetAtom1Ptr()->GetSegmentName() ) != m_segmentFilterMap.end()) &&
	 (m_segmentFilterMap.find( (*iter)->GetAtom2Ptr()->GetSegmentName() ) != m_segmentFilterMap.end()) )
      bondList.push_back(*iter);//Got a match, so store the bond pointer;
  }
  return bondList;
}

//Moved from RbtMOL2FileSource - sets partial charges on a per-residue basis
//Automatically determines protonation states of ASP/GLU/HIS using Mandatory and Forbidden
//atoms in spParamSource
void RbtBaseMolecularFileSource::SetupPartialIonicGroups(RbtAtomList& atoms, RbtParameterFileSourcePtr spParamSource)
{
  const RbtString _MANDATORY("MANDATORY");
  const RbtString _FORBIDDEN("FORBIDDEN");
  if (atoms.empty()) {
    cout << "WARNING SetupPartialIonicGroups: Empty atom list" << endl;
    return;
  }
  RbtAtomPtr leadAtom = atoms.front();
  RbtString subunitName = leadAtom->GetSubunitName();
  RbtString match
    = leadAtom->GetSegmentName() + ":"
    + subunitName + "_"
    + leadAtom->GetSubunitId() + ":";
  if (Rbt::GetNumMatchingAtoms(atoms,match) != atoms.size()) {
    cout << "WARNING SetupPartialIonicGroups: Inconsistent subunit names in atom list headed by "
	 << leadAtom->GetFullAtomName() << endl;
    return;
  }
  RbtStringList resList(spParamSource->GetSectionList());
  if (std::find(resList.begin(), resList.end(),subunitName) == resList.end()) {
    //cout << "INFO SetupPartialIonicGroups: No section for residue " << subunitName << endl;
    return;
  }

  spParamSource->SetSection(subunitName);
  RbtStringList atList = spParamSource->GetParameterList();
  if (std::find(atList.begin(), atList.end(),_MANDATORY) != atList.end()) {
    RbtString mandatory = spParamSource->GetParameterValueAsString(_MANDATORY);
    RbtStringList mandAtoms = Rbt::ConvertDelimitedStringToList(mandatory);
    RbtInt nPresent = Rbt::GetNumMatchingAtoms(atoms,mandAtoms);
    if (nPresent != mandAtoms.size()) {
#ifdef _DEBUG
      cout << "INFO SetupPartialIonicGroups: Only " << nPresent << " out of " << mandAtoms.size()
      	   << " mandatory atoms present in atom list headed by "
      	   << leadAtom->GetFullAtomName() << endl;
#endif //_DEBUG
      return;
    }
    std::remove(atList.begin(), atList.end(),_MANDATORY);
  }
  if (std::find(atList.begin(), atList.end(),_FORBIDDEN) != atList.end()) {
    RbtString forbidden = spParamSource->GetParameterValueAsString(_FORBIDDEN);
    RbtStringList forbAtoms = Rbt::ConvertDelimitedStringToList(forbidden);
    RbtInt nPresent = Rbt::GetNumMatchingAtoms(atoms,forbAtoms);
    if (nPresent > 0) {
#ifdef _DEBUG
      cout << "INFO SetupPartialIonicGroups: " << nPresent
      	   << " forbidden atoms present in atom list headed by "
      	   << leadAtom->GetFullAtomName() << endl;
#endif //_DEBUG
      return;
    }
    std::remove(atList.begin(), atList.end(),_FORBIDDEN);
  }
  
  for (RbtStringListConstIter aIter = atList.begin(); aIter != atList.end(); aIter++) {
    RbtDouble partialCharge(spParamSource->GetParameterValue(*aIter));//Get the partial charge value
#ifdef _DEBUG
    cout << endl << "Trying to match " << *aIter << endl;
#endif //_DEBUG
    //Find the atoms which match the specifier string
    RbtAtomList selectedAtoms = Rbt::GetMatchingAtomList(atoms,*aIter);
    //Now we've got the matching atoms, set the group charge on each atom
    for (RbtAtomListIter iter = selectedAtoms.begin(); iter != selectedAtoms.end(); iter++) {
      (*iter)->SetGroupCharge(partialCharge);
#ifdef _DEBUG
      cout << "INFO Setting group charge on " << (*iter)->GetFullAtomName() << " to " << partialCharge << endl;
#endif //_DEBUG
    }
  }
}
