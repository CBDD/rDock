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

#include "RbtNonBondedGrid.h"
#include "RbtFileError.h"

//Static data members
RbtString RbtNonBondedGrid::_CT("RbtNonBondedGrid");

////////////////////////////////////////
//Constructors/destructors
//Construct a NXxNYxNZ grid running from gridMin at gridStep resolution
RbtNonBondedGrid::RbtNonBondedGrid(const RbtCoord& gridMin, const RbtCoord& gridStep,
				   RbtUInt NX, RbtUInt NY, RbtUInt NZ, RbtUInt NPad) :
  RbtBaseGrid(gridMin,gridStep,NX,NY,NZ,NPad)
{
	CreateMap();
  _RBTOBJECTCOUNTER_CONSTR_("RbtNonBondedGrid");
}

//Constructor reading params from binary stream
RbtNonBondedGrid::RbtNonBondedGrid(istream& istr) : RbtBaseGrid(istr)
{
	CreateMap();
  OwnRead(istr);
  _RBTOBJECTCOUNTER_CONSTR_("RbtNonBondedGrid");
}

//Default destructor
RbtNonBondedGrid::~RbtNonBondedGrid()
{
  ClearAtomLists();
  _RBTOBJECTCOUNTER_DESTR_("RbtNonBondedGrid");
}

//Copy constructor
RbtNonBondedGrid::RbtNonBondedGrid(const RbtNonBondedGrid& grid) : RbtBaseGrid(grid)
{
	CreateMap();
  CopyGrid(grid);
  _RBTOBJECTCOUNTER_COPYCONSTR_("RbtNonBondedGrid");
}

//Copy constructor taking base class argument
RbtNonBondedGrid::RbtNonBondedGrid(const RbtBaseGrid& grid) : RbtBaseGrid(grid)
{
	CreateMap();
  _RBTOBJECTCOUNTER_COPYCONSTR_("RbtNonBondedGrid");
}

//Copy assignment
RbtNonBondedGrid& RbtNonBondedGrid::operator=(const RbtNonBondedGrid& grid)
{
  if (this != &grid) {
    //Clear the current atom lists
    ClearAtomLists();
    //In this case we need to explicitly call the base class operator=
    RbtBaseGrid::operator=(grid);
    CopyGrid(grid);
  }
  return *this;
}
//Copy assignment taking base class argument
RbtNonBondedGrid& RbtNonBondedGrid::operator=(const RbtBaseGrid& grid)
{
  if (this != &grid) {
    //Clear the current atom lists
    ClearAtomLists();
    //In this case we need to explicitly call the base class operator=
    RbtBaseGrid::operator=(grid);
  }
  return *this;
}

////////////////////////////////////////
//Virtual functions for reading/writing grid data to streams in
//text and binary format
//Subclasses should provide their own private OwnPrint,OwnWrite, OwnRead
//methods to handle subclass data members, and override the public
//Print,Write and Read methods


//Text output
void RbtNonBondedGrid::Print(ostream& ostr) const
{
  RbtBaseGrid::Print(ostr);
  OwnPrint(ostr);
}
 
//Binary output
void RbtNonBondedGrid::Write(ostream& ostr) const
{
  RbtBaseGrid::Write(ostr);
  OwnWrite(ostr);
}

//Binary input
void RbtNonBondedGrid::Read(istream& istr)
{
  ClearAtomLists();
  RbtBaseGrid::Read(istr);
  OwnRead(istr);
}

////////////////////////////////////////
//Public methods
////////////////

/////////////////////////
//Get attribute functions
/////////////////////////

const RbtAtomRList& RbtNonBondedGrid::GetAtomList(RbtUInt iXYZ) const
{
  //RbtAtomListMapConstIter iter = m_atomMap.find(iXYZ);
  //if (iter != m_atomMap.end())
  //  return iter->second;
  //else
  //  return m_emptyList;
  
  //DM 6 Nov 2000 - map replaced by vector
  if (isValid(iXYZ)) {
  	return m_atomMap[iXYZ];
	}
	else {
		return m_emptyList;
	}
}

const RbtAtomRList& RbtNonBondedGrid::GetAtomList(const RbtCoord& c) const
{
  //if (isValid(c))
  //  return GetAtomList(GetIXYZ(c));
  //else
  //  return m_emptyList;
  
  //DM 6 Nov 2000 - map replaced by vector
  if (isValid(c)) {
  	return m_atomMap[GetIXYZ(c)];
	}
	else {
		return m_emptyList;
	}
}

/////////////////////////
//Set attribute functions
/////////////////////////
void RbtNonBondedGrid::SetAtomLists(RbtAtom* pAtom, RbtDouble radius)
{
  const RbtCoord& c = pAtom->GetCoords();
  RbtUIntList sphereIndices;
  GetSphereIndices(c,radius,sphereIndices);

  for (RbtUIntListConstIter sphereIter = sphereIndices.begin(); sphereIter != sphereIndices.end(); sphereIter++) {
    //RbtAtomListMapIter mapIter = m_atomMap.find(*sphereIter);
    //if (mapIter != m_atomMap.end()) {
    //  (mapIter->second).push_back(pAtom);
    //}
    //else {
    //  m_atomMap[*sphereIter] = RbtAtomRList(1,pAtom);
    //}
    m_atomMap[*sphereIter].push_back(pAtom);
  }
}

void RbtNonBondedGrid::ClearAtomLists()
{
  //m_atomMap.clear();
  //Clear each atom list separately, without clearing the whole map (vector)
  for (RbtAtomListMapIter iter = m_atomMap.begin(); iter != m_atomMap.end(); iter++) {
    (*iter).clear();
  }
}

void RbtNonBondedGrid::UniqueAtomLists() {
  for (RbtAtomListMapIter iter = m_atomMap.begin(); iter != m_atomMap.end(); iter++) {
    //cout << _CT << ": before = " << (*iter).size();
    std::sort((*iter).begin(),(*iter).end(),Rbt::RbtAtomPtrCmp_Ptr());
    RbtAtomRListIter uniqIter = std::unique((*iter).begin(),(*iter).end());
    (*iter).erase(uniqIter,(*iter).end());
    //cout << "; After = " << (*iter).size() << endl;
  }
}

///////////////////////////////////////////////////////////////////////////
//Protected methods

//Protected method for writing data members for this class to text stream
void RbtNonBondedGrid::OwnPrint(ostream& ostr) const {
  ostr << endl << "Class\t" << _CT << endl;
  ostr << "No. of entries in the map: " << m_atomMap.size() << endl;
  //TO BE COMPLETED - no real need for dumping the atom list info
}

//Protected method for writing data members for this class to binary stream
//(Serialisation)
void RbtNonBondedGrid::OwnWrite(ostream& ostr) const {
  //Write all the data members
  //NO MEANS OF WRITING THE ATOM LISTS IN A WAY WHICH CAN BE READ BACK IN
  //i.e. we are holding pointers to atoms which would need to be recreated
  //What we need is some kind of object database I guess.
  //Note: By not writing the title key here, it enables a successful read
  //of ANY grid subclass file. e.g. An RbtNonBondedGrid can be constructed
  //from an RbtRealGrid output stream, so will have the same grid dimensions
  //The RbtRealGrid data array will simply be ignored
}

//Protected method for reading data members for this class from binary stream
//WARNING: Assumes grid data array has already been created
//and is of the correct size
void RbtNonBondedGrid::OwnRead(istream& istr) throw (RbtError) {
  //Read all the data members
  //NOTHING TO READ - see above
}

///////////////////////////////////////////////////////////////////////////
//Private methods
//Helper function called by copy constructor and assignment operator
void RbtNonBondedGrid::CopyGrid(const RbtNonBondedGrid& grid) {
  //This copies the atom lists, but of course the atoms themselves are not copied
  m_atomMap = grid.m_atomMap;
}

//DM 6 Nov 2000 - create AtomListMap of the appropriate size
void RbtNonBondedGrid::CreateMap() {
	m_atomMap = RbtAtomListMap(GetN());
}

