/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtBaseGrid.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include <algorithm> //for min, max, count
#include <iomanip>

#include "RbtBaseGrid.h"
#include "RbtFileError.h"

//Static data members
RbtString RbtBaseGrid::_CT("RbtBaseGrid");

////////////////////////////////////////
//Constructors/destructors
//Construct a NXxNYxNZ grid running from gridMin at gridStep resolution
RbtBaseGrid::RbtBaseGrid(const RbtCoord& gridMin, const RbtCoord& gridStep,
		       RbtUInt NX, RbtUInt NY, RbtUInt NZ, RbtUInt NPad) :
  m_step(gridStep),m_NX(NX),m_NY(NY),m_NZ(NZ),m_N(NX*NY*NZ),
  m_SX(NY*NZ),m_SY(NZ),m_SZ(1),m_NPad(NPad)
{
  //Set the logical (real-world) dimensions of the grid
  SetGridMin(gridMin);
  _RBTOBJECTCOUNTER_CONSTR_("RbtBaseGrid");
}

//Constructor reading params from binary stream
RbtBaseGrid::RbtBaseGrid(istream& istr) {
	OwnRead(istr);
  _RBTOBJECTCOUNTER_CONSTR_("RbtBaseGrid");
}

//Default destructor
RbtBaseGrid::~RbtBaseGrid()
{
  _RBTOBJECTCOUNTER_DESTR_("RbtBaseGrid");
}

//Copy constructor
RbtBaseGrid::RbtBaseGrid(const RbtBaseGrid& grid)
{
  CopyGrid(grid);
  _RBTOBJECTCOUNTER_COPYCONSTR_("RbtBaseGrid");
}

//Copy assignment
RbtBaseGrid& RbtBaseGrid::operator=(const RbtBaseGrid& grid)
{
  if (this != &grid) {
    CopyGrid(grid);
  }
  return *this;
}


//Insertion operator (primarily for debugging)
ostream& operator<<(ostream& s, const RbtBaseGrid& grid)
{
  grid.Print(s);
  return s;
}

////////////////////////////////////////
//Virtual functions for reading/writing grid data to streams in
//text and binary format
//Subclasses should provide their own private OwnPrint,OwnWrite, OwnRead
//methods to handle subclass data members, and override the public
//Print,Write and Read methods

//These implementations are not strictly necessary for the base class
//as base class objects cannot be instantiated

//Text output
void RbtBaseGrid::Print(ostream& ostr) const
{
  OwnPrint(ostr);
}
 
//Binary output
void RbtBaseGrid::Write(ostream& ostr) const
{
  OwnWrite(ostr);
}

//Binary input
void RbtBaseGrid::Read(istream& istr)
{
  OwnRead(istr);
}

////////////////////////////////////////
//Public methods
////////////////

/////////////////////////
//Set attribute functions
/////////////////////////


//Set grid min without changing the grid step
void RbtBaseGrid::SetGridMin(const RbtCoord& gridMin)
{
  //Determine the integral min and max grid point coords (in multiples of grid step from origin)
  m_nXMin = floor(gridMin.x/m_step.x + 0.5);
  m_nYMin = floor(gridMin.y/m_step.y + 0.5);
  m_nZMin = floor(gridMin.z/m_step.z + 0.5);

  //Max integral coords are just min+N-1
  m_nXMax = m_nXMin+m_NX-1;
  m_nYMax = m_nYMin+m_NY-1;
  m_nZMax = m_nZMin+m_NZ-1;

  //Now determine the min and max real-world coords which are covered by the grid
  //Include a border of +/-0.5*gridStep so that grid points lie in the center of each grid interval
  m_min = (RbtCoord(m_nXMin,m_nYMin,m_nZMin)-0.5)*m_step;
  m_max = (RbtCoord(m_nXMax,m_nYMax,m_nZMax)+0.5)*m_step;

  //Recalc the pad coords as we've moved the grid
  SetPad(m_NPad);
}

//Set grid step
//The integral min and max grid point coords (in multiples of grid step from origin) are left unchanged
//Only the real-world min and max coords are updated
void RbtBaseGrid::SetGridStep(const RbtVector& gridStep)
{
  m_step = gridStep;
  //Now determine the min and max real-world coords which are covered by the grid
  //Include a border of +/-0.5*gridStep so that grid points lie in the center of each grid interval
  m_min = (RbtCoord(m_nXMin,m_nYMin,m_nZMin)-0.5)*m_step;
  m_max = (RbtCoord(m_nXMax,m_nYMax,m_nZMax)+0.5)*m_step;

  //Recalc the pad coords as we've moved the grid
  SetPad(m_NPad);
}


//Set a pad region (cuboid inside the grid which defines the valid coords)
//Coords outside the cuboid are deemed invalid
void RbtBaseGrid::SetPad(RbtUInt NPad)
{
  m_NPad = NPad;
  m_padMin = m_min + m_step*m_NPad;
  m_padMax = m_max - m_step*m_NPad;
}

//Returns list of real-world coordinates for given set of iXYZ indices
RbtCoordList RbtBaseGrid::GetCoordList(const RbtUIntSet& iXYZSet) const
{
  RbtCoordList coordList;
  coordList.reserve(iXYZSet.size());
  for (RbtUIntSetConstIter iter = iXYZSet.begin(); iter != iXYZSet.end(); iter++)
    coordList.push_back(GetCoord(*iter));
  return coordList;
}

//DM 17 May 1999 - returns the set of valid grid points within a sphere of a given center and radius
//DM 17 Jul 2000 - use vector<RbtUInt> and return by reference, for performance boost
void RbtBaseGrid::GetSphereIndices(const RbtCoord& c, RbtDouble radius, RbtUIntList& sIndices) const
{
  sIndices.clear();//Clear the return list
  
  //Iterate over the cube defined by the coord and the radius
  //without exceeding the limits of the pad region
  RbtCoord cubeMin = Rbt::Max(c - radius,m_padMin);
  RbtCoord cubeMax = Rbt::Min(c + radius,m_padMax);

  //Convert to array indices
  //DM 16 Apr 1999 - check again for indices out of range
  //If cubeMax == m_padMax, iMaxX,Y,Z would be one too high
  RbtUInt iMinX = std::max(GetIX(cubeMin),m_NPad+1);
  RbtUInt iMinY = std::max(GetIY(cubeMin),m_NPad+1);
  RbtUInt iMinZ = std::max(GetIZ(cubeMin),m_NPad+1);
  RbtUInt iMaxX = std::min(GetIX(cubeMax),m_NX-m_NPad);
  RbtUInt iMaxY = std::min(GetIY(cubeMax),m_NY-m_NPad);
  RbtUInt iMaxZ = std::min(GetIZ(cubeMax),m_NZ-m_NPad);

  //Determine if points are inside sphere by checking x^2 + y^2 + z^2 <= rad^2
  RbtDouble rad2 = radius*radius;

  RbtDouble rX;//X-distance to sphere center
  RbtUInt iX;
  for (iX = iMinX, rX = GetXCoord(iMinX)-c.x; iX <= iMaxX; iX++, rX += m_step.x) {
    RbtDouble rX2 = rX*rX;
    //Only bother with the Y loop if we haven't already exceeded rad2
    if (rX2 <= rad2) {
      RbtDouble rY; //Y-distance to sphere center
      RbtUInt iY;
      for (iY = iMinY, rY = GetYCoord(iMinY)-c.y; iY <= iMaxY; iY++, rY += m_step.y) {
        RbtDouble rXY2 = rX2 + rY*rY;
        //Only bother with the Z loop if we haven't already exceeded rad2
        if (rXY2 <= rad2) {
          RbtDouble rZ;
          RbtUInt iZ;
          RbtUInt iXYZ;
          for (iZ = iMinZ, iXYZ = GetIXYZ(iX,iY,iMinZ), rZ = GetZCoord(iMinZ)-c.z;
	             iZ <= iMaxZ; iZ++, iXYZ++, rZ += m_step.z) {
            RbtDouble rXYZ2 = rXY2 + rZ*rZ;
            if (rXYZ2 <= rad2) {
              sIndices.push_back(iXYZ);
            }
          }
        }
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////
//Protected methods

//Protected method for writing data members for this class to text stream
void RbtBaseGrid::OwnPrint(ostream& ostr) const {
  ostr << "Class\t" << _CT << endl;  
  ostr << "m_min\t" << m_min << endl;
  ostr << "m_max\t" << m_max << endl;
  ostr << "m_step\t" << m_step << endl;
  ostr << "m_padMin\t" << m_padMin << endl;
  ostr << "m_padMax\t" << m_padMax << endl;
  ostr << "m_NX,m_NY,m_NZ\t" << m_NX << "," << m_NY << "," << m_NZ << endl;
  ostr << "m_N\t" << m_N << endl;
  ostr << "m_SX,m_SY,m_SZ\t" << m_SX << "," << m_SY << "," << m_SZ << endl;
  ostr << "m_NPad\t" << m_NPad << endl;
  ostr << "m_nXMin,m_nXMax\t" << m_nXMin << "," << m_nXMax << endl;
  ostr << "m_nYMin,m_nYMax\t" << m_nYMin << "," << m_nYMax << endl;
  ostr << "m_nZMin,m_nZMax\t" << m_nZMin << "," << m_nZMax << endl;
}

//Protected method for writing data members for this class to binary stream
//(Serialisation)
void RbtBaseGrid::OwnWrite(ostream& ostr) const {
	//Write the class name as a title so we can check the authenticity of streams
	//on read
	const char* const gridTitle = _CT.c_str();
	RbtInt length = strlen(gridTitle);
	Rbt::WriteWithThrow(ostr, (const char*) &length, sizeof(length));
	Rbt::WriteWithThrow(ostr, gridTitle, length);

	//Write all the data members
	//The RbtCoord::Write method doesn't throw an error, but we can live with that for now
	m_min.Write(ostr);
	m_max.Write(ostr);
	m_step.Write(ostr);
	m_padMin.Write(ostr);
	m_padMax.Write(ostr);
	Rbt::WriteWithThrow(ostr, (const char*) &m_NX, sizeof(m_NX));
	Rbt::WriteWithThrow(ostr, (const char*) &m_NY, sizeof(m_NY));
	Rbt::WriteWithThrow(ostr, (const char*) &m_NZ, sizeof(m_NZ));
  Rbt::WriteWithThrow(ostr, (const char*) &m_N, sizeof(m_N));
	Rbt::WriteWithThrow(ostr, (const char*) &m_SX, sizeof(m_SX));
	Rbt::WriteWithThrow(ostr, (const char*) &m_SY, sizeof(m_SY));
	Rbt::WriteWithThrow(ostr, (const char*) &m_SZ, sizeof(m_SZ));
	Rbt::WriteWithThrow(ostr, (const char*) &m_NPad, sizeof(m_NPad));
  Rbt::WriteWithThrow(ostr, (const char*) &m_nXMin, sizeof(m_nXMin));
	Rbt::WriteWithThrow(ostr, (const char*) &m_nXMax, sizeof(m_nXMax));
	Rbt::WriteWithThrow(ostr, (const char*) &m_nYMin, sizeof(m_nYMin));
	Rbt::WriteWithThrow(ostr, (const char*) &m_nYMax, sizeof(m_nYMax));
	Rbt::WriteWithThrow(ostr, (const char*) &m_nZMin, sizeof(m_nZMin));
	Rbt::WriteWithThrow(ostr, (const char*) &m_nZMax, sizeof(m_nZMax));
}

//Protected method for reading data members for this class from binary stream
void RbtBaseGrid::OwnRead(istream& istr) throw (RbtError) {
	//Read title
	RbtInt length;
	Rbt::ReadWithThrow(istr, (char*) &length, sizeof(length));
	char* gridTitle = new char [length+1];
	Rbt::ReadWithThrow(istr, gridTitle, length);
	//Add null character to end of string
	gridTitle[length] = '\0';
	//Compare title with class name
	RbtBool match = (_CT == gridTitle);
	delete [] gridTitle;
	if (!match) {
		throw RbtFileParseError(_WHERE_,"Invalid title string in " + _CT + "::Read()");
	}

	//Read all the data members
	//The RbtCoord::Read method doesn't throw an error, but we can live with that for now
	m_min.Read(istr);
	m_max.Read(istr);
	m_step.Read(istr);
	m_padMin.Read(istr);
	m_padMax.Read(istr);
	Rbt::ReadWithThrow(istr, (char*) &m_NX, sizeof(m_NX));
	Rbt::ReadWithThrow(istr, (char*) &m_NY, sizeof(m_NY));
	Rbt::ReadWithThrow(istr, (char*) &m_NZ, sizeof(m_NZ));
	Rbt::ReadWithThrow(istr, (char*) &m_N, sizeof(m_N));
	Rbt::ReadWithThrow(istr, (char*) &m_SX, sizeof(m_SX));
	Rbt::ReadWithThrow(istr, (char*) &m_SY, sizeof(m_SY));
	Rbt::ReadWithThrow(istr, (char*) &m_SZ, sizeof(m_SZ));
	Rbt::ReadWithThrow(istr, (char*) &m_NPad, sizeof(m_NPad));
	Rbt::ReadWithThrow(istr, (char*) &m_nXMin, sizeof(m_nXMin));
	Rbt::ReadWithThrow(istr, (char*) &m_nXMax, sizeof(m_nXMax));
	Rbt::ReadWithThrow(istr, (char*) &m_nYMin, sizeof(m_nYMin));
	Rbt::ReadWithThrow(istr, (char*) &m_nYMax, sizeof(m_nYMax));
	Rbt::ReadWithThrow(istr, (char*) &m_nZMin, sizeof(m_nZMin));
	Rbt::ReadWithThrow(istr, (char*) &m_nZMax, sizeof(m_nZMax));
}


///////////////////////////////////////////////////////////////////////////
//Private methods

//Helper function called by copy constructor and assignment operator
//Base class doesn't store any pointers so this is a straightforward
//memberwise copy
void RbtBaseGrid::CopyGrid(const RbtBaseGrid& grid) {
  m_min = grid.m_min;
  m_max = grid.m_max;
  m_step = grid.m_step;
  m_padMin = grid.m_padMin;
  m_padMax = grid.m_padMax;
  m_NX = grid.m_NX;
  m_NY = grid.m_NY;
  m_NZ = grid.m_NZ;
  m_N = grid.m_N;
  m_SX = grid.m_SX;
  m_SY = grid.m_SY;
  m_SZ = grid.m_SZ;
  m_NPad = grid.m_NPad;
  m_nXMin = grid.m_nXMin;
  m_nXMax = grid.m_nXMax;
  m_nYMin = grid.m_nYMin;
  m_nYMax = grid.m_nYMax;
  m_nZMin = grid.m_nZMin;
  m_nZMax = grid.m_nZMax;
}


