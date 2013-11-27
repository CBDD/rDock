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

#include "RbtInteractionGrid.h"
#include "RbtFileError.h"
#include "RbtPseudoAtom.h"

//Could be a useful general function
//If pAtom is a pseudo atom, then pushes all the constituent atoms onto list
//else, push pAtom itself onto the list
void RbtInteractionCenter::AccumulateAtomList(const RbtAtom* pAtom, RbtAtomRList& atomList) const {
  if (!pAtom) return;
  RbtPseudoAtom* pPseudo = dynamic_cast<RbtPseudoAtom*>(const_cast<RbtAtom*>(pAtom));
  if (pPseudo) {
    RbtAtomList constituents=pPseudo->GetAtomList();
    std::copy(constituents.begin(),constituents.end(),std::back_inserter(atomList));
  }
  else {
    atomList.push_back(const_cast<RbtAtom*>(pAtom));
  }
}

//Returns list of constituent atoms
//(deconvolutes pseudoatoms into their constituent RbtAtom* lists)
RbtAtomRList RbtInteractionCenter::GetAtomList() const {
  RbtAtomRList atomList;
  AccumulateAtomList(m_pAtom1,atomList);
  AccumulateAtomList(m_pAtom2,atomList);
  AccumulateAtomList(m_pAtom3,atomList);
  return atomList;
}

//An interaction is selected if any of its constituent atoms are selected
//If any of the constituent atoms are pseudo-atoms, then check these also
RbtBool RbtInteractionCenter::isSelected() const {
  RbtAtomRList atomList = GetAtomList();
  return (std::count_if(atomList.begin(),atomList.end(),Rbt::isAtomSelected()) > 0);
}

//Select/deselect the interaction center (selects all constituent atoms)
void Rbt::SelectInteractionCenter::operator() (RbtInteractionCenter* pIC) {
  RbtAtomRList atomList = pIC->GetAtomList();
  std::for_each(atomList.begin(),atomList.end(),Rbt::SelectAtom(b));
}

//Static data members
RbtString RbtInteractionGrid::_CT("RbtInteractionGrid");

////////////////////////////////////////
//Constructors/destructors
//Construct a NXxNYxNZ grid running from gridMin at gridStep resolution
RbtInteractionGrid::RbtInteractionGrid(const RbtCoord& gridMin, const RbtCoord& gridStep,
				   RbtUInt NX, RbtUInt NY, RbtUInt NZ, RbtUInt NPad) :
  RbtBaseGrid(gridMin,gridStep,NX,NY,NZ,NPad)
{
	CreateMap();
  _RBTOBJECTCOUNTER_CONSTR_("RbtInteractionGrid");
}

//Constructor reading params from binary stream
RbtInteractionGrid::RbtInteractionGrid(istream& istr) : RbtBaseGrid(istr)
{
	CreateMap();
  OwnRead(istr);
  _RBTOBJECTCOUNTER_CONSTR_("RbtInteractionGrid");
}

//Default destructor
RbtInteractionGrid::~RbtInteractionGrid()
{
	ClearInteractionLists();
  _RBTOBJECTCOUNTER_DESTR_("RbtInteractionGrid");
}

//Copy constructor
RbtInteractionGrid::RbtInteractionGrid(const RbtInteractionGrid& grid) : RbtBaseGrid(grid)
{
	CreateMap();
  CopyGrid(grid);
  _RBTOBJECTCOUNTER_COPYCONSTR_("RbtInteractionGrid");
}

//Copy constructor taking base class argument
RbtInteractionGrid::RbtInteractionGrid(const RbtBaseGrid& grid) : RbtBaseGrid(grid)
{
		CreateMap();
    _RBTOBJECTCOUNTER_COPYCONSTR_("RbtInteractionGrid");
}

//Copy assignment
RbtInteractionGrid& RbtInteractionGrid::operator=(const RbtInteractionGrid& grid)
{
  if (this != &grid) {
    //Clear the current lists
    ClearInteractionLists();
    //In this case we need to explicitly call the base class operator=
    RbtBaseGrid::operator=(grid);
    CopyGrid(grid);
  }
  return *this;
}
//Copy assignment taking base class argument
RbtInteractionGrid& RbtInteractionGrid::operator=(const RbtBaseGrid& grid)
{
  if (this != &grid) {
    //Clear the current lists
    ClearInteractionLists();
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
void RbtInteractionGrid::Print(ostream& ostr) const
{
  RbtBaseGrid::Print(ostr);
  OwnPrint(ostr);
}
 
//Binary output
void RbtInteractionGrid::Write(ostream& ostr) const
{
  RbtBaseGrid::Write(ostr);
  OwnWrite(ostr);
}

//Binary input
void RbtInteractionGrid::Read(istream& istr)
{
  ClearInteractionLists();
  RbtBaseGrid::Read(istr);
  OwnRead(istr);
}

////////////////////////////////////////
//Public methods
////////////////

/////////////////////////
//Get attribute functions
/////////////////////////

const RbtInteractionCenterList& RbtInteractionGrid::GetInteractionList(RbtUInt iXYZ) const
{
  //RbtInteractionListMapConstIter iter = m_intnMap.find(iXYZ);
  //if (iter != m_intnMap.end())
  //  return iter->second;
  //else
  //  return m_emptyList;
  
  //DM 3 Nov 2000 - map replaced by vector
  if (isValid(iXYZ)) {
  	return m_intnMap[iXYZ];
	}
	else {
		return m_emptyList;
	}
}

const RbtInteractionCenterList& RbtInteractionGrid::GetInteractionList(const RbtCoord& c) const
{
  //if (isValid(c))
  //  return GetInteractionList(GetIXYZ(c));
  //else
  //  return m_emptyList;
  
  //DM 3 Nov 2000 - map replaced by vector
  if (isValid(c)) {
    return m_intnMap[GetIXYZ(c)];
  }
  else {
    //cout << _CT << "::GetInteractionList," << c << " is off grid" << endl;
    return m_emptyList;
  }
}

/////////////////////////
//Set attribute functions
/////////////////////////
void RbtInteractionGrid::SetInteractionLists(RbtInteractionCenter* pIntn, RbtDouble radius)
{
  //Index using atom 1 coords - check if atom 1 is present
  RbtAtom* pAtom1 = pIntn->GetAtom1Ptr();
  if (pAtom1 == NULL)
    return;
    
  //DM 27 Oct 2000 - GetCoords now returns by reference
  const RbtCoord& c = pAtom1->GetCoords();
  RbtUIntList sphereIndices;
  GetSphereIndices(c,radius,sphereIndices);

  for (RbtUIntListConstIter sphereIter = sphereIndices.begin(); sphereIter != sphereIndices.end(); sphereIter++) {
    //RbtInteractionListMapIter mapIter = m_intnMap.find(*sphereIter);
    //if (mapIter != m_intnMap.end()) {
    //  (mapIter->second).push_back(pIntn);
    //}
    //else {
    //  m_intnMap[*sphereIter] = RbtInteractionCenterList(1,pIntn);
    //}
    m_intnMap[*sphereIter].push_back(pIntn);
    
  }
}

void RbtInteractionGrid::ClearInteractionLists()
{
  //m_intnMap.clear();
  //Clear each interaction center list separately, without clearing the whole map (vector)
  for (RbtInteractionListMapIter iter = m_intnMap.begin(); iter != m_intnMap.end(); iter++) {
  	(*iter).clear();
	}
}

void RbtInteractionGrid::UniqueInteractionLists() {
  for (RbtInteractionListMapIter iter = m_intnMap.begin(); iter != m_intnMap.end(); iter++) {
    //cout << _CT << ": before = " << (*iter).size();
    std::sort((*iter).begin(),(*iter).end(),Rbt::InteractionCenterCmp());
    RbtInteractionCenterListIter uniqIter = std::unique((*iter).begin(),(*iter).end());
    (*iter).erase(uniqIter,(*iter).end());
    //cout << "; After = " << (*iter).size() << endl;
  }
}

///////////////////////////////////////////////////////////////////////////
//Protected methods

//Protected method for writing data members for this class to text stream
void RbtInteractionGrid::OwnPrint(ostream& ostr) const {
  ostr << endl << "Class\t" << _CT << endl;
  ostr << "No. of entries in the map: " << m_intnMap.size() << endl;
  //TO BE COMPLETED - no real need for dumping the interaction list info
}

//Protected method for writing data members for this class to binary stream
//(Serialisation)
void RbtInteractionGrid::OwnWrite(ostream& ostr) const {
  //Write all the data members
  //NO MEANS OF WRITING THE INTERACTION LISTS IN A WAY WHICH CAN BE READ BACK IN
  //i.e. we are holding pointers to atoms which would need to be recreated
  //What we need is some kind of object database I guess.
  //Note: By not writing the title key here, it enables a successful read
  //of ANY grid subclass file. e.g. An RbtInteractionGrid can be constructed
  //from an RbtRealGrid output stream, so will have the same grid dimensions
  //The RbtRealGrid data array will simply be ignored
}

//Protected method for reading data members for this class from binary stream
//WARNING: Assumes grid data array has already been created
//and is of the correct size
void RbtInteractionGrid::OwnRead(istream& istr) throw (RbtError) {
  //Read all the data members
  //NOTHING TO READ - see above
  CreateMap();
}

///////////////////////////////////////////////////////////////////////////
//Private methods
//Helper function called by copy constructor and assignment operator
void RbtInteractionGrid::CopyGrid(const RbtInteractionGrid& grid) {
  //This copies the interaction lists, but of course the atoms themselves are not copied
  m_intnMap = grid.m_intnMap;
}

//DM 3 Nov 2000 - create InteractionListMap of the appropriate size
void RbtInteractionGrid::CreateMap() {
	m_intnMap = RbtInteractionListMap(GetN());
}
