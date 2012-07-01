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

#include "RbtCharmmTypesFileSource.h"

//Constructors
RbtCharmmTypesFileSource::RbtCharmmTypesFileSource(const char* fileName) :
  RbtBaseFileSource(fileName)
{
  _RBTOBJECTCOUNTER_CONSTR_("RbtCharmmTypesFileSource");
}

RbtCharmmTypesFileSource::RbtCharmmTypesFileSource(const RbtString fileName) :
  RbtBaseFileSource(fileName)
{
  _RBTOBJECTCOUNTER_CONSTR_("RbtCharmmTypesFileSource");
}

//Destructor
RbtCharmmTypesFileSource::~RbtCharmmTypesFileSource()
{
  ClearTypesCache();
  _RBTOBJECTCOUNTER_DESTR_("RbtCharmmTypesFileSource");
}


////////////////////////////////////////
//Public methods
RbtInt RbtCharmmTypesFileSource::GetNumTypes()
{
  Parse();
  return m_typesList.size();
}

RbtCharmmTypeList RbtCharmmTypesFileSource::GetTypeList()
{
  Parse();
  return m_typesList;
}


//Private methods
//Pure virtual in RbtBaseFileSource - needs to be defined here
void RbtCharmmTypesFileSource::Parse() throw (RbtError)
{
  //Expected string constants in MASSES.RTF file
  const RbtString strMassKey("MASS");

  //Only parse if we haven't already done so
  if (!m_bParsedOK) {
    ClearTypesCache();//Clear current cache
    Read();//Read the file

    try {
      //Very simple parser - loop over all lines looking for string "MASS" in first 4 chars
      for ( RbtFileRecListIter fileIter = m_lineRecs.begin(); fileIter != m_lineRecs.end(); fileIter++) {
	if ( (*fileIter).substr(0,4) == strMassKey) {
	  //We have a match so read in the line and store in the types list
	  CharmmType chrmType;
	  RbtString strDummy;
	  istrstream istr((*fileIter).c_str());
	  istr >> strDummy
	       >> chrmType.nAtomType
	       >> chrmType.strAtomType
	       >> chrmType.mass
	       >> chrmType.element;
	  //				  >> strDummy
	  //				  >> chrmType.comment;
	  //DM 8 Dec 1998 - store the entire line in the comment field as reading via the input stream
	  //was truncating the comment at the first whitespace character
	  chrmType.comment = *fileIter;
	  m_typesList.push_back(chrmType);
#ifdef _DEBUG
	  //	  cout << "Type #" << chrmType.nAtomType << " = " << chrmType.strAtomType << endl;
#endif //_DEBUG
	}
      }
      //////////////////////////////////////////////////////////
      //If we get this far everything is OK
      m_bParsedOK = true;
    }
    
    catch (RbtError& error) {
      ClearTypesCache();
      throw;//Rethrow the RbtError
    }
  }
}

void RbtCharmmTypesFileSource::ClearTypesCache()
{
  m_typesList.clear();
}
