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

#include "RbtCrdFileSource.h"
#include "RbtFileError.h"

//Constructors
RbtCrdFileSource::RbtCrdFileSource(const char* fileName) :
  RbtBaseMolecularFileSource(fileName,"CRD_FILE_SOURCE") //Call base class constructor
{
  _RBTOBJECTCOUNTER_CONSTR_("RbtCrdFileSource");
}

RbtCrdFileSource::RbtCrdFileSource(const RbtString fileName) :
  RbtBaseMolecularFileSource(fileName,"CRD_FILE_SOURCE") //Call base class constructor
{
  _RBTOBJECTCOUNTER_CONSTR_("RbtCrdFileSource");
}


//Default destructor
RbtCrdFileSource::~RbtCrdFileSource()
{
  _RBTOBJECTCOUNTER_DESTR_("RbtCrdFileSource");
}

void RbtCrdFileSource::Parse() throw (RbtError)
{
  //Only parse if we haven't already done so
  if (!m_bParsedOK) {
    ClearMolCache();//Clear current cache
    Read();//Read the file

    try {
      RbtFileRecListIter fileIter = m_lineRecs.begin();
      RbtFileRecListIter fileEnd = m_lineRecs.end();

      //////////////////////////////////////////////////////////
      //2a Read title lines and store them
      while ( ((*fileIter).at(0) == '*') && //first char is a "*"
	      ((*fileIter).size() > 1) && //line length > 1
	      (fileIter != fileEnd) ) { //not at end of file
	m_titleList.push_back(*fileIter++);
      }

      //2b Check we ended on a line consisting of a single *
      if ( (fileIter == fileEnd) ||
	   ( (*fileIter).at(0) != '*'))
	throw RbtFileParseError(_WHERE_,"Missing title record terminator (*) in " + GetFileName());

      //////////////////////////////////////////////////////////
      //3a. Read number of atoms
      RbtInt nAtomRec;
      fileIter++;
      istrstream((*fileIter).c_str()) >> nAtomRec;
      
      //3b ...and store them
      fileIter++;
      m_atomList.reserve(nAtomRec);//Allocate enough memory for the vector

      RbtInt nAtomId; //original atom number in CRD file
      RbtString strSegmentName; //segment name in PSF file
      RbtInt nSubunitId; //subunit(residue) number in CRD file
      //(residues are numbered consecutively across all segments.
      //Note that we don't store the residue number in RbtAtom)
      RbtString strSubunitId; //subunit(residue) ID in CRD file
      RbtString strSubunitName; //subunit(residue) name in CRD file
      RbtString strAtomName; //atom name from CRD file
      RbtCoord coord; //X,Y,Z coords
      //NOTE: we do not read the weighting array value at present

      while ( (m_atomList.size() < nAtomRec) && (fileIter != fileEnd)) {
	istrstream((*fileIter++).c_str()) >> nAtomId
					  >> nSubunitId
					  >> strSubunitName
					  >> strAtomName
					  >> coord.x
					  >> coord.y
					  >> coord.z
					  >> strSegmentName
					  >> strSubunitId;

	//Construct a new atom (constructor only accepts the 2D params)
	RbtAtomPtr spAtom(new RbtAtom(nAtomId,
				      0,//Atomic number undefined (gets set in CheckAtomParams)
				      strAtomName,
				      strSubunitId,
				      strSubunitName,
				      strSegmentName
				      )
			  );


	//Now set the 3-D params we have
	spAtom->SetCoords(coord);

	m_atomList.push_back(spAtom);
	m_segmentMap[strSegmentName]++;//increment atom count in segment map
      }

      //3c ..and check we read them all before reaching the end of the file
      if (m_atomList.size() != nAtomRec)
	throw RbtFileParseError(_WHERE_,"Incomplete atom records in " + GetFileName());

      //////////////////////////////////////////////////////////
      //If we get this far everything is OK
      m_bParsedOK = true;
    }

    catch (RbtError& error) {
      ClearMolCache();//Clear the molecular cache so we don't return incomplete atom and bond lists
      throw;//Rethrow the RbtError
    }
  }
}
