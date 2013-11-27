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

#include "RbtNonBondedHHSGrid.h"
#include "RbtAtom.h"
#include "RbtFileError.h"

RbtString	RbtNonBondedHHSGrid::_CT("RbtNonBondedHHSGrid");

RbtNonBondedHHSGrid::RbtNonBondedHHSGrid( const RbtCoord& gridMin,
		const RbtCoord& gridStep,
		RbtUInt NX, RbtUInt NY, RbtUInt NZ, RbtUInt NPad/*=0*/) 
		: RbtBaseGrid(gridMin, gridStep, NX,NY,NZ,NPad) 
{
	CreateMap();
	_RBTOBJECTCOUNTER_CONSTR_("RbtNonBondedHHSGrid");
}

RbtNonBondedHHSGrid::RbtNonBondedHHSGrid(istream& istr) : RbtBaseGrid(istr) 
{
	CreateMap();
	OwnRead(istr);
	_RBTOBJECTCOUNTER_CONSTR_("RbtNonBondedHHSGrid");
}

RbtNonBondedHHSGrid::~RbtNonBondedHHSGrid()
{
	ClearHHSLists();
	_RBTOBJECTCOUNTER_DESTR_("RbtNonBondedHHSGrid");
}

RbtNonBondedHHSGrid::RbtNonBondedHHSGrid(const RbtNonBondedHHSGrid& aGrid) : RbtBaseGrid(aGrid)
{
	CreateMap();
	CopyGrid(aGrid);
	_RBTOBJECTCOUNTER_COPYCONSTR_("RbtNonBondedHHSGrid");
}

RbtNonBondedHHSGrid::RbtNonBondedHHSGrid(const RbtBaseGrid& grid) : RbtBaseGrid(grid)
{
	CreateMap();
	_RBTOBJECTCOUNTER_COPYCONSTR_("RbtNonBondedHHSGrid");
}

RbtNonBondedHHSGrid& RbtNonBondedHHSGrid::operator=(const RbtNonBondedHHSGrid& grid)
{
	if (this != &grid) {
		//Clear the current atom lists
		ClearHHSLists();
		//In this case we need to explicitly call the base class operator=
		RbtBaseGrid::operator=(grid);
		CopyGrid(grid);
	}
	return *this;
}
//Copy assignment taking base class argument
RbtNonBondedHHSGrid& RbtNonBondedHHSGrid::operator=(const RbtBaseGrid& grid)
{
	if (this != &grid) {
		//Clear the current atom lists
		ClearHHSLists();
		//In this case we need to explicitly call the base class operator=
		RbtBaseGrid::operator=(grid);
	}
	return *this;
}

//Text output
void RbtNonBondedHHSGrid::Print(ostream& ostr) const
{
  RbtBaseGrid::Print(ostr);
  OwnPrint(ostr);
}
 
//Binary output
void RbtNonBondedHHSGrid::Write(ostream& ostr) const
{
  RbtBaseGrid::Write(ostr);
  OwnWrite(ostr);
}

//Binary input
void RbtNonBondedHHSGrid::Read(istream& istr) 
{
  ClearHHSLists();
  RbtBaseGrid::Read(istr);
  OwnRead(istr);
}

const HHS_SolvationRList& RbtNonBondedHHSGrid::GetHHSList(RbtUInt iXYZ) const
{
  if (isValid(iXYZ)) {
  	return m_hhsMap[iXYZ];
	} else {
		return m_emptyList;
	}
}

const HHS_SolvationRList& RbtNonBondedHHSGrid::GetHHSList(const RbtCoord& c) const
{
  if (isValid(c)) {
  	return m_hhsMap[GetIXYZ(c)];
	} else {
		return m_emptyList;
	}
}

void RbtNonBondedHHSGrid::SetHHSLists(HHS_Solvation* pHHS, RbtDouble radius)
{
  const RbtCoord& c = (pHHS->GetAtom())->GetCoords();
  RbtUIntList sphereIndices;
  GetSphereIndices(c,radius,sphereIndices);

  for (RbtUIntListConstIter sphereIter = sphereIndices.begin(); 
		  sphereIter != sphereIndices.end(); 
		  sphereIter++) {
    m_hhsMap[*sphereIter].push_back(pHHS);
  }
}

void RbtNonBondedHHSGrid::ClearHHSLists()
{
  //m_atomMap.clear();
  //Clear each atom list separately, without clearing the whole map (vector)
  for (HHS_SolvationListMapIter iter = m_hhsMap.begin(); 
		  iter != m_hhsMap.end(); 
		  iter++) {
  	(*iter).clear();
	}
}

void RbtNonBondedHHSGrid::OwnPrint(ostream& ostr) const {
  ostr << endl << "Class\t" << _CT << endl;
  ostr << "No. of entries in the map: " << m_hhsMap.size() << endl;
}


void RbtNonBondedHHSGrid::OwnWrite(ostream& ostr) const {
}

void RbtNonBondedHHSGrid::OwnRead(istream& istr) throw (RbtError) {
}

void RbtNonBondedHHSGrid::CopyGrid(const RbtNonBondedHHSGrid& grid) {
  m_hhsMap = grid.m_hhsMap;
}

void RbtNonBondedHHSGrid::CreateMap() {
	m_hhsMap = HHS_SolvationListMap(GetN());
}

