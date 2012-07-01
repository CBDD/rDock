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

#include <algorithm> //for min, max, count
#include <iomanip>
#include <cstring> 
using std::setw;

#include "RbtRealGrid.h"
#include "RbtFileError.h"

//Static data members
RbtString RbtRealGrid::_CT("RbtRealGrid");

////////////////////////////////////////
//Constructors/destructors
//Construct a NXxNYxNZ grid running from gridMin at gridStep resolution
RbtRealGrid::RbtRealGrid(const RbtCoord& gridMin, const RbtCoord& gridStep,
		                     RbtUInt NX, RbtUInt NY, RbtUInt NZ, RbtUInt NPad) :
  RbtBaseGrid(gridMin,gridStep,NX,NY,NZ,NPad),m_tol(0.001),
  m_grid(NULL),m_data(NULL)
{
  CreateArrays();
  //Initialise the grid to zero
  SetAllValues(0.0);

  _RBTOBJECTCOUNTER_CONSTR_("RbtRealGrid");
}

//Constructor reading params from binary stream
RbtRealGrid::RbtRealGrid(istream& istr) :
                        RbtBaseGrid(istr),m_grid(NULL),m_data(NULL) {
  //Base class constructor has already read the grid dimensions
  //etc, so all we have to do here is created the array
  //and read in the grid values
  CreateArrays();
  OwnRead(istr);
  _RBTOBJECTCOUNTER_CONSTR_("RbtRealGrid");
}

//Default destructor
RbtRealGrid::~RbtRealGrid()
{
  ClearArrays();
  _RBTOBJECTCOUNTER_DESTR_("RbtRealGrid");
}

//Copy constructor
RbtRealGrid::RbtRealGrid(const RbtRealGrid& grid) :
                         RbtBaseGrid(grid),m_grid(NULL),m_data(NULL)
{
  //Base class constructor has already been called
  //so we just need to create the array and copy the array values
  CreateArrays();
  CopyGrid(grid);
  _RBTOBJECTCOUNTER_COPYCONSTR_("RbtRealGrid");
}

//Copy constructor taking a base class argument
//Sets up the grid dimensions, then creates an empty data array
RbtRealGrid::RbtRealGrid(const RbtBaseGrid& grid) :
                         RbtBaseGrid(grid),m_grid(NULL),m_data(NULL)
{
  CreateArrays();
  SetAllValues(0.0);
  _RBTOBJECTCOUNTER_COPYCONSTR_("RbtRealGrid");
}

//Copy assignment
RbtRealGrid& RbtRealGrid::operator=(const RbtRealGrid& grid)
{
  if (this != &grid) {
    //Clear the current grid
    ClearArrays();
    //In this case we need to explicitly call the base class operator=
    RbtBaseGrid::operator=(grid);
    CreateArrays();  
    CopyGrid(grid);
  }
  return *this;
}
//Copy assignment taking a base class argument
//Sets up the grid dimensions, then creates an empty data array
RbtRealGrid& RbtRealGrid::operator=(const RbtBaseGrid& grid)
{
  if (this != &grid) {
    //Clear the current grid
    ClearArrays();
    //In this case we need to explicitly call the base class operator=
    RbtBaseGrid::operator=(grid);
    CreateArrays();  
    SetAllValues(0.0);
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
void RbtRealGrid::Print(ostream& ostr) const
{
  RbtBaseGrid::Print(ostr);
  OwnPrint(ostr);
}
 
//Binary output
void RbtRealGrid::Write(ostream& ostr) const
{
  RbtBaseGrid::Write(ostr);
  OwnWrite(ostr);
}

//Binary input
void RbtRealGrid::Read(istream& istr)
{
  //Clear the current grid before reading the new grid dimensions
  ClearArrays();
  RbtBaseGrid::Read(istr);//Base class read
  //Now create the new grid and read in the grid values
  CreateArrays();
  OwnRead(istr);
}



////////////////////////////////////////
//Public methods
////////////////

/////////////////////////
//Get/Set value functions
/////////////////////////

//DM 20 Jul 2000 - get values smoothed by trilinear interpolation
//D. Oberlin and H.A. Scheraga, J. Comp. Chem. (1998) 19, 71.
RbtDouble RbtRealGrid::GetSmoothedValue(const RbtCoord& c) const
{
  const RbtCoord& gridMin = GetGridMin();
  const RbtVector& gridStep = GetGridStep();
  RbtDouble rx = 1.0 / gridStep.x;//reciprocal of grid step (x)
  RbtDouble ry = 1.0 / gridStep.y;//reciprocal of grid step (y)
  RbtDouble rz = 1.0 / gridStep.z;//reciprocal of grid step (z)
  //Get lower left corner grid point
  //(not necessarily the nearest grid point as returned by GetIX() etc)
  //Need to shift the int(..) argument by half a grid step
  RbtUInt iX = int(rx * (c.x - gridMin.x) - 0.5) + 1;
  RbtUInt iY = int(ry * (c.y - gridMin.y) - 0.5) + 1;
  RbtUInt iZ = int(rz * (c.z - gridMin.z) - 0.5) + 1;
#ifdef _DEBUG
  cout << "GetSmoothedValue" << c << "\tiX,iY,iZ=" << iX << "\t" << iY << "\t" << iZ << endl;
#endif //_DEBUG
  //Check this point (iX,iY,iZ) and (iX+1,iY+1,iZ+1) are all in bounds
  //else return the unsmoothed GetValue(c)
  if (!isValid(iX,iY,iZ) || !isValid(iX+1,iY+1,iZ+1)) {
#ifdef _DEBUG
    cout << "Out of bounds" << endl;
#endif //_DEBUG
    return GetValue(c);
  }
  //p is the vector relative to the lower left corner
  RbtVector p = c - GetCoord(iX,iY,iZ);
#ifdef _DEBUG
  cout << "p=" << p << endl;
#endif //_DEBUG
  //Set up B0 and B1 for each of x,y,z axes
  //RbtDoubleList bx(2);
  //RbtDoubleList by(2);
  //RbtDoubleList bz(2);
  //bx[1] = rx * p.x;
  //bx[0] = 1.0 - bx[1];
  //by[1] = ry * p.y;
  //by[0] = 1.0 - by[1];
  //bz[1] = rz * p.z;
  //bz[0] = 1.0 - bz[1];
  
  RbtDouble val(0.0);
  //DM 3/5/2005 - fully unroll this loop
  RbtDouble bx1 = rx * p.x;
  RbtDouble bx0 = 1.0 - bx1;
  RbtDouble by1 = ry * p.y;
  RbtDouble by0 = 1.0 - by1;
  RbtDouble bz1 = rz * p.z;
  RbtDouble bz0 = 1.0 - bz1;
  RbtDouble bx0by0 = bx0 * by0;
  RbtDouble bx0by1 = bx0 * by1;
  RbtDouble bx1by0 = bx1 * by0;
  RbtDouble bx1by1 = bx1 * by1;
  val += m_grid[iX][iY][iZ]       * bx0by0 * bz0;
  val += m_grid[iX][iY][iZ+1]     * bx0by0 * bz1;
  val += m_grid[iX][iY+1][iZ]     * bx0by1 * bz0;
  val += m_grid[iX][iY+1][iZ+1]   * bx0by1 * bz1;
  val += m_grid[iX+1][iY][iZ]     * bx1by0 * bz0;
  val += m_grid[iX+1][iY][iZ+1]   * bx1by0 * bz1;
  val += m_grid[iX+1][iY+1][iZ]   * bx1by1 * bz0;
  val += m_grid[iX+1][iY+1][iZ+1] * bx1by1 * bz1;
  //for (RbtUInt i = 0; i < 2; i++) {
  //  for (RbtUInt j = 0; j < 2; j++) {
  //    for (RbtUInt k = 0; k < 2; k++) {
  //      val += m_grid[iX+i][iY+j][iZ+k]*bx[i]*by[j]*bz[k];
  //    }
  //  }
  //}
  return val;   
}

//Set all grid points to the given value
void RbtRealGrid::SetAllValues(RbtDouble val)
{
  for (RbtUInt i = 0; i < GetN(); i++) {
    m_data[i]=val;
  }
}

//Replaces all grid points between oldValMin and oldValMax with newVal
void RbtRealGrid::ReplaceValueRange(RbtDouble oldValMin, RbtDouble oldValMax, RbtDouble newVal)
{
  for (RbtUInt i = 0; i < GetN(); i++) {
    float d = m_data[i];
    if ( (d >= oldValMin) && (d < oldValMax))
      m_data[i] = newVal;
  }
}

//Set all grid points within radius of coord to the given value
//If bOverwrite is false, does not replace non-zero values
//If bOverwrite is true, all grid points are set to the new value
void RbtRealGrid::SetSphere(const RbtCoord& c, RbtDouble radius, RbtDouble val, RbtBool bOverwrite)
{
  RbtUIntList sphereIndices;
  GetSphereIndices(c,radius,sphereIndices);
  SetValues(sphereIndices,val,bOverwrite);
}

//Set all grid points with radii between rad1 and rad2 from coord to the given value
//If bOverwrite is false, does not replace non-zero values
//If bOverwrite is true, all grid points are set the new value
void RbtRealGrid::SetSurface(const RbtCoord& c, RbtDouble innerRad, RbtDouble outerRad, RbtDouble val, RbtBool bOverwrite)
{
  RbtUIntList oIndices;//List of sphere points for outer sphere
  RbtUIntList iIndices;//List of sphere points for outer sphere
  RbtUIntList sIndices;//List of points for surface region (points in outer sphere, but not inner sphere)
  GetSphereIndices(c,outerRad,oIndices);
  GetSphereIndices(c,innerRad,iIndices);
  //std::sort(oIndices.begin(),oIndices.end());
  //std::sort(iIndices.begin(),iIndices.end());
  std::set_difference(oIndices.begin(),oIndices.end(),iIndices.begin(),iIndices.end(),std::back_inserter(sIndices));
  SetValues(sIndices,val,bOverwrite);
}

//Sets all grid points with value=oldValue, which are adjacent to those with value=adjacentValue, to value=newValue
//+/- tolerance is applied to oldValue and adjacentValue
void RbtRealGrid::CreateSurface(RbtDouble oldVal, RbtDouble adjVal, RbtDouble newVal)
{
  //Iterate over the cuboid defined by the pad coords
  RbtUInt iMinX = GetPad()+1;
  RbtUInt iMinY = GetPad()+1;
  RbtUInt iMinZ = GetPad()+1;
  RbtUInt iMaxX = GetNX()-GetPad();
  RbtUInt iMaxY = GetNY()-GetPad();
  RbtUInt iMaxZ = GetNZ()-GetPad();

  for (RbtUInt iX = iMinX; iX <= iMaxX; iX++) {
    for (RbtUInt iY = iMinY; iY <= iMaxY; iY++) {
      for (RbtUInt iZ = iMinZ; iZ <= iMaxZ; iZ++) {
        //We have a match with oldVal
	      if (fabs(m_grid[iX][iY][iZ] - oldVal) < m_tol) {
	        //Check the six adjacent points for a match with adjVal
          if ( ( (iX > iMinX) && (fabs(m_grid[iX-1][iY][iZ] - adjVal) < m_tol)) ||
               ( (iX < iMaxX) && (fabs(m_grid[iX+1][iY][iZ] - adjVal) < m_tol)) ||
	             ( (iY > iMinY) && (fabs(m_grid[iX][iY-1][iZ] - adjVal) < m_tol)) ||
	             ( (iY < iMaxY) && (fabs(m_grid[iX][iY+1][iZ] - adjVal) < m_tol)) ||
	             ( (iZ > iMinZ) && (fabs(m_grid[iX][iY][iZ-1] - adjVal) < m_tol)) ||
	             ( (iZ < iMaxZ) && (fabs(m_grid[iX][iY][iZ+1] - adjVal) < m_tol)) )
            m_grid[iX][iY][iZ] = newVal;
	      }
      }
    }
  }
}

//DM 16 Apr 1999 - helper function for determining solvent accessible regions
//Returns true if any of the grid points within a sphere around the central coord have the specified value
RbtBool RbtRealGrid::isValueWithinSphere(const RbtCoord& c, RbtDouble radius, RbtDouble val)
{
  RbtUIntList sphereIndices;
  GetSphereIndices(c,radius,sphereIndices);
  return isValueWithinList(sphereIndices, val);
}

//Sets all grid points with value=oldValue, which have no grid points with value=adjacentValue within a sphere
//of given radius, to value=newValue
//+/- tolerance is applied to oldValue and adjacentValue
void RbtRealGrid::SetAccessible(RbtDouble radius, RbtDouble oldVal, RbtDouble adjVal, RbtDouble newVal, RbtBool bCenterOnly)
{
  //Iterate over the cuboid defined by the pad coords
  RbtUInt iMinX = GetPad()+1;
  RbtUInt iMinY = GetPad()+1;
  RbtUInt iMinZ = GetPad()+1;
  RbtUInt iMaxX = GetNX()-GetPad();
  RbtUInt iMaxY = GetNY()-GetPad();
  RbtUInt iMaxZ = GetNZ()-GetPad();
  
  //Work out the maximum no. of grid points in the sphere and reserve enough space in the indices vector.
  //Actually, this is a considerable overestimate (no. of points in the enclosing cube)
  RbtUIntList sphereIndices;
  RbtUInt nMax = (int(radius/GetGridStep().x)+1)*(int(radius/GetGridStep().y)+1)*(int(radius/GetGridStep().z)+1);
  sphereIndices.reserve(nMax);
  for (RbtUInt iX = iMinX; iX <= iMaxX; iX++) {
    for (RbtUInt iY = iMinY; iY <= iMaxY; iY++) {
      for (RbtUInt iZ = iMinZ; iZ <= iMaxZ; iZ++) {
	//We have a match with oldVal
	if (fabs(m_grid[iX][iY][iZ] - oldVal) < m_tol) {
	  RbtCoord c = GetCoord(iX,iY,iZ);          
	  //Check the sphere around this grid point
          GetSphereIndices(c,radius,sphereIndices);
	  if (!isValueWithinList(sphereIndices,adjVal)) {
	    if (bCenterOnly)
	      m_grid[iX][iY][iZ] = newVal;//Set just the center grid point
	    else
	      //SetValues(sphereIndices,newVal,false);//Set all grid points in the sphere
	      SetValues(sphereIndices,newVal,true);//DM 3 April 2002 - now overwrite ALL gridpoints with newVal
	  }	      
	}
      }
    }
  }
}


/////////////////////////
//Statistical functions
/////////////////////////

//Returns number of occurrences of a given value range
RbtUInt RbtRealGrid::CountRange(RbtDouble valMin, RbtDouble valMax) const
{
  RbtUInt n(0);
  for (RbtUInt i = 0; i < GetN(); i++) {
    float d = m_data[i];
    if ( (d >= valMin) && (d < valMax))
      n++;
  }
  return n;
}

//Min/max values
RbtDouble RbtRealGrid::MinValue() const
{
  float fMin = m_data[0];
  for (RbtUInt i = 0; i < GetN(); i++) {
    if (m_data[i] < fMin)
      fMin = m_data[i];
  }
  return fMin;
}

RbtDouble RbtRealGrid::MaxValue() const
{
  float fMax = m_data[0];
  for (RbtUInt i = 0; i < GetN(); i++) {
    if (m_data[i] > fMax)
      fMax = m_data[i];
  }
  return fMax;
}

//iXYZ index of grid point with minimum value
RbtUInt RbtRealGrid::FindMinValue() const
{
  RbtUInt iMin = 0;
  for (RbtUInt i = 0; i < GetN(); i++) {
    if (m_data[i] < m_data[iMin])
      iMin = i;
  }
  return iMin; 
}

//iXYZ index of grid point with maximum value
RbtUInt RbtRealGrid::FindMaxValue() const
{
  RbtUInt iMax = 0;
  for (RbtUInt i = 0; i < GetN(); i++) {
    if (m_data[i] > m_data[iMax])
      iMax = i;
  }
  return iMax; 
}


/////////////////////////
//Output functions
/////////////////////////

//Dump grid in a format readable by Insight
void RbtRealGrid::PrintInsightGrid(ostream& s) const
{
  s << "RBT FFT GRID" << endl;
  s << "(1F15.10)" << endl;
  s.precision(3);
  s.setf(ios_base::fixed,ios_base::floatfield);
  s.setf(ios_base::right,ios_base::adjustfield);
  s << setw(8) << GetGridStep().x*(GetNX()-1)
    << setw(8) << GetGridStep().y*(GetNY()-1)
    << setw(8) << GetGridStep().z*(GetNZ()-1)
    << setw(8) << 90.0
    << setw(8) << 90.0
    << setw(8) << 90.0
    << endl;
  //s.setf(0,ios_base::floatfield);
  s << setw(5) << GetNX()-1
    << setw(5) << GetNY()-1
    << setw(5) << GetNZ()-1
    << endl;
  s << setw(5) << 1
    << setw(5) << GetnXMin() << setw(5) << GetnXMax()
    << setw(5) << GetnYMin() << setw(5) << GetnYMax()
    << setw(5) << GetnZMin() << setw(5) << GetnZMax()
    << endl;
  s.precision(10);
  //s.setf(ios_base::fixed,ios_base::floatfield);
  //s.setf(ios_base::right,ios_base::adjustfield);
  //Insight expects data in [1][1][1],[2][1][1]..[NX][1][1] order
  for (RbtUInt iZ = 1; iZ <= GetNZ(); iZ++) {
    for (RbtUInt iY = 1; iY <= GetNY(); iY++) {
      for (RbtUInt iX = 1; iX <= GetNX(); iX++) {
        s << setw(15) << m_grid[iX][iY][iZ] << endl;
      }
    }
  }
}


///////////////////////////////////////////////////////////////////////////
//Protected methods

//Protected method for writing data members for this class to text stream
void RbtRealGrid::OwnPrint(ostream& ostr) const {
  ostr << "Class\t" << _CT << endl;  
  //Iterate over all grid points
  float tol = GetTolerance();
  for (RbtUInt iX = 1; iX <= GetNX(); iX++) {
    ostr << endl << endl << "Plane iX=" << iX << endl;
    for (RbtUInt iY = 1; iY <= GetNY(); iY++) {
      for (RbtUInt iZ = 1; iZ <= GetNZ(); iZ++) {
	      float f = m_grid[iX][iY][iZ];
	      ostr << ( (f < -tol) ? '-' : (f > tol) ? '+': '.');
      }
      ostr << endl;
    }
  }
}

//Protected method for writing data members for this class to binary stream
//(Serialisation)
void RbtRealGrid::OwnWrite(ostream& ostr) const {
  //Write the class name as a title so we can check the authenticity of streams
  //on read
  const char* const gridTitle = _CT.c_str();
  RbtInt length = strlen(gridTitle);
  Rbt::WriteWithThrow(ostr, (const char*) &length, sizeof(length));
  Rbt::WriteWithThrow(ostr, gridTitle, length);
  
  //Write all the data members
  Rbt::WriteWithThrow(ostr, (const char*) &m_tol, sizeof(m_tol));
  for (RbtUInt i = 0; i < GetN(); i++) {
    Rbt::WriteWithThrow(ostr, (const char*) &m_data[i], sizeof(m_data[i]));
  }
}

//Protected method for reading data members for this class from binary stream
//WARNING: Assumes grid data array has already been created
//and is of the correct size
void RbtRealGrid::OwnRead(istream& istr) throw (RbtError) {
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
  Rbt::ReadWithThrow(istr, (char*) &m_tol, sizeof(m_tol));
  for (RbtUInt i = 0; i < GetN(); i++) {
    Rbt::ReadWithThrow(istr, (char*) &m_data[i], sizeof(m_data[i]));
  }
}


///////////////////////////////////////////////////////////////////////////
//Private methods

//DM 17 Jul 2000 - analogous to isValueWithinSphere but iterates over arbitrary set
//of IXYZ indices. Private method as there is no error checking on iXYZ values out of bounds
RbtBool RbtRealGrid::isValueWithinList(const RbtUIntList& iXYZList, RbtDouble val)
{
  for (RbtUIntListConstIter iter = iXYZList.begin(); iter != iXYZList.end(); iter++) {
	  if (fabs(m_data[*iter] - val) < m_tol) {
		  return true;
		}
  }
  return false; //Got all the way without finding a match
}

//Analogous to SetSphere but iterates over arbitrary set
//of IXYZ indices. Private method as there is no error checking on iXYZ values out of bounds
//If bOverwrite is false, does not replace non-zero values
//If bOverwrite is true, all grid points are set the new value
void RbtRealGrid::SetValues(const RbtUIntList& iXYZList, RbtDouble val, RbtBool bOverwrite)
{
  for (RbtUIntListConstIter iter = iXYZList.begin(); iter != iXYZList.end(); iter++) {
	  if (bOverwrite || (fabs(m_data[*iter]) < m_tol)) {
		  m_data[*iter] = val;
		}
  }
}

void RbtRealGrid::CreateArrays() {
  if (m_grid != NULL) {//Clear existing grid
    ClearArrays();
  }
  int nX = GetNX();
  int nY = GetNY();
  int nZ = GetNZ();
  m_grid = new float**[nX+1];
  m_grid[1] = new float*[nX*nY+1];
  m_grid[1][1] = new float[nX*nY*nZ+1];

  for (int iY = 2; iY <= nY; iY++) {
    m_grid[1][iY] = m_grid[1][iY-1] + nZ;
  }
  for (int iX = 2; iX <= nX; iX++) {
    m_grid[iX] = m_grid[iX-1] + nY;
    m_grid[iX][1] = m_grid[iX-1][1] + nY*nZ;
    for (int iY = 2; iY <= nY; iY++) {
      m_grid[iX][iY] = m_grid[iX][iY-1] + nZ;
    }
  }
  m_data = &m_grid[1][1][1];
}

void RbtRealGrid::ClearArrays() {
  delete m_grid[1][1];
  delete m_grid[1];
  delete m_grid;
  m_grid = NULL;
  m_data = NULL;
}

//Helper function called by copy constructor and assignment operator
//Need to create a new deep copy of the data array
//Gets called after array has been created, and base class copy has been done
void RbtRealGrid::CopyGrid(const RbtRealGrid& grid) {
  m_tol = grid.m_tol;
  float* gridData = grid.GetGridData();
  for (RbtUInt i = 0; i < GetN(); i++) {
    m_data[i] = gridData[i];
  }
}


