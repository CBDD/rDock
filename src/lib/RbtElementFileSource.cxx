/*This file is part of Rdock.

    Rdock is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Rdock is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Rdock.  If not, see <http://www.gnu.org/licenses/>.*/

#include <sstream>
using std::istringstream;
using std::ostringstream;
using std::ends;

#include "RbtElementFileSource.h"
#include "RbtFileError.h"

//Constructors
//RbtElementFileSource::RbtElementFileSource(const char* fileName) :
//  RbtBaseFileSource(fileName)
//{
//  _RBTOBJECTCOUNTER_CONSTR_("RbtElementFileSource");
//}

RbtElementFileSource::RbtElementFileSource(const RbtString& fileName) :
  RbtBaseFileSource(fileName),m_dHBondRadiusIncr(0.0),m_dImplicitRadiusIncr(0.0)
{
  _RBTOBJECTCOUNTER_CONSTR_("RbtElementFileSource");
}

//Destructor
RbtElementFileSource::~RbtElementFileSource()
{
  ClearElementDataCache();
  _RBTOBJECTCOUNTER_DESTR_("RbtElementFileSource");
}


////////////////////////////////////////
//Public methods
RbtString RbtElementFileSource::GetTitle()
{
  Parse();
  return m_strTitle;
}
  
RbtString RbtElementFileSource::GetVersion()
{
  Parse();
  return m_strVersion;
}

RbtUInt RbtElementFileSource::GetNumElements()
{
  Parse();
  return m_elementNameMap.size();
}


//List of element names
RbtStringList RbtElementFileSource::GetElementNameList()
{
  Parse();
  RbtStringList elementNameList;
  elementNameList.reserve(m_elementNameMap.size());
  for (RbtStringElementDataMapConstIter iter = m_elementNameMap.begin(); iter != m_elementNameMap.end(); iter++)
    elementNameList.push_back((*iter).first);
  return elementNameList;
}

//List of atomic numbers
RbtIntList RbtElementFileSource::GetAtomicNumberList()
{
  Parse();
  RbtIntList atomicNumberList;
  atomicNumberList.reserve(m_atomicNumberMap.size());
  for (RbtIntElementDataMapConstIter iter = m_atomicNumberMap.begin(); iter != m_atomicNumberMap.end(); iter++)
    atomicNumberList.push_back((*iter).first);
  return atomicNumberList;
}


//Get element data for a given element name, throws error if not found
RbtElementData RbtElementFileSource::GetElementData(const RbtString& strElementName) throw (RbtError)
{
  Parse();
  RbtStringElementDataMapIter iter = m_elementNameMap.find(strElementName);
  if (iter != m_elementNameMap.end())
    return (*iter).second;
  else
    throw RbtFileMissingParameter(_WHERE_,"Element " +strElementName+" not found in "+GetFileName());
}

//Get element data for a given atomic number, throws error if not found
RbtElementData RbtElementFileSource::GetElementData(RbtInt nAtomicNumber) throw (RbtError)
{
  Parse();
  RbtIntElementDataMapIter iter = m_atomicNumberMap.find(nAtomicNumber);
  if (iter != m_atomicNumberMap.end())
    return (*iter).second;
  else {
    ostringstream ostr;
    ostr << "Atomic number " << nAtomicNumber << " not found in " << GetFileName() << ends;
    RbtString strError(ostr.str());
    //delete ostr.str();
    throw RbtFileMissingParameter(_WHERE_,strError);
  }
}

//Check if given element name is present
RbtBool RbtElementFileSource::isElementNamePresent(const RbtString& strElementName)
{
  Parse();
  RbtStringElementDataMapIter iter = m_elementNameMap.find(strElementName);
  if (iter != m_elementNameMap.end())
    return true;
  else
    return false;
}

//Check if given atomic number is present
RbtBool RbtElementFileSource::isAtomicNumberPresent(RbtInt nAtomicNumber)
{
  Parse();
  RbtIntElementDataMapIter iter = m_atomicNumberMap.find(nAtomicNumber);
  if (iter != m_atomicNumberMap.end())
    return true;
  else
    return false;
}

//Get vdw radius increment for hydrogen-bonding donors
RbtDouble RbtElementFileSource::GetHBondRadiusIncr()
{
  Parse();
  return m_dHBondRadiusIncr;
}

//Get vdw radius increment for atoms with implicit hydrogens
RbtDouble RbtElementFileSource::GetImplicitRadiusIncr()
{
  Parse();
  return m_dImplicitRadiusIncr;
}

//Private methods
//Pure virtual in RbtBaseFileSource - needs to be defined here
void RbtElementFileSource::Parse() throw (RbtError)
{
  const RbtString strRbtKey = "RBT_ELEMENT_FILE_V1.01";
  const RbtString strTitleKey = "TITLE ";
  const RbtString strVersionKey = "VERSION ";
  const RbtString strElementKey = "ELEMENT";
  const RbtString strHBondRadiusKey = "HBOND_RADIUS_INCREMENT";
  const RbtString strImplicitRadiusKey = "IMPLICIT_RADIUS_INCREMENT";

  //Only parse if we haven't already done so
  if (!m_bParsedOK) {
    ClearElementDataCache();//Clear current cache
    Read();//Read the file

    try {
      RbtFileRecListIter fileIter = m_lineRecs.begin();
      RbtFileRecListIter fileEnd = m_lineRecs.end();

      //////////////////////////////////////////////////////////
      //1. Check for RBT string on line 1
      if ( (*fileIter).find(strRbtKey) != 0)
	throw RbtFileParseError(_WHERE_,"Missing "+strRbtKey+" string in " + GetFileName());

      //////////////////////////////////////////////////////////
      //2. Parse the rest of file
      while (++fileIter != fileEnd) {
	//Ignore blank lines and comment lines
	if ( ((*fileIter).length() == 0) || ((*fileIter).at(0) == '#') ) {
	  continue;
	}
	//Check for Title record
	else if ( (*fileIter).find(strTitleKey) == 0) {
	  m_strTitle = *fileIter;
	  m_strTitle.erase(0,strTitleKey.length());
	}
	//Check for Version record
	else if ( (*fileIter).find(strVersionKey) == 0) {
	  m_strVersion = *fileIter;
	  m_strVersion.erase(0,strVersionKey.length());
	}
	//Check for HBond radius increment
	else if ( (*fileIter).find(strHBondRadiusKey) == 0) {
	  RbtString strDummy;
	  istringstream istr((*fileIter).c_str());
	  istr >> strDummy >> m_dHBondRadiusIncr;
	}
	//Check for Implicit radius increment
	else if ( (*fileIter).find(strImplicitRadiusKey) == 0) {
	  RbtString strDummy;
	  istringstream istr((*fileIter).c_str());
	  istr >> strDummy >> m_dImplicitRadiusIncr;
	}
	//Check for element record
	else if ( (*fileIter).find(strElementKey) == 0) {
	  RbtElementData elemData;
	  RbtString strDummy;
	  istringstream istr((*fileIter).c_str());
	  istr >> strDummy
	       >> elemData.atomicNo
	       >> elemData.element
	       >> elemData.minVal
	       >> elemData.maxVal
	       >> elemData.commonVal
	       >> elemData.mass
	       >> elemData.vdwRadius;
	  //Store the data in two maps, one indexed by element name, the other by atomic number
	  m_elementNameMap[elemData.element] = elemData;
	  m_atomicNumberMap[elemData.atomicNo] = elemData;
	}
      }
      //////////////////////////////////////////////////////////
      //If we get this far everything is OK
      m_bParsedOK = true;
    }
    
    catch (RbtError& error) {
      ClearElementDataCache();
      throw;//Rethrow the RbtError
    }
  }
}

void RbtElementFileSource::ClearElementDataCache()
{
  m_elementNameMap.clear();
  m_atomicNumberMap.clear();
  m_strTitle = "";
  m_strVersion = "";
  m_dHBondRadiusIncr = 0.0;
  m_dImplicitRadiusIncr = 0.0;
}
