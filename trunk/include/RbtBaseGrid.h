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

//Interface to all grid classes. Provides conversion methods between different
//addressing methods: (iXYZ),(iX,iY,iZ), and (x,y,z) coords. Underlying data
//values are managed by subclasses.

#ifndef _RBTBASEGRID_H_
#define _RBTBASEGRID_H_

#include "RbtConfig.h"
#include "RbtCoord.h"

class RbtBaseGrid
{
 public:
	//Class type string
	static RbtString _CT;

  ////////////////////////////////////////
  //Constructors/destructors
  //Construct a NXxNYxNZ grid running from gridMin at gridStep resolution
  RbtBaseGrid(const RbtCoord& gridMin, const RbtVector& gridStep,
	     RbtUInt NX, RbtUInt NY, RbtUInt NZ, RbtUInt NPad=0);
  //Constructor reading all params from binary stream
  RbtBaseGrid(istream& istr);
  RbtBaseGrid(const RbtBaseGrid&);//Copy constructor
  RbtBaseGrid& operator=(const RbtBaseGrid&);//Copy assignment
  virtual ~RbtBaseGrid(); //Default destructor

  //Friend functions
  //Insertion operator (primarily for debugging)
  friend ostream& operator<<(ostream& s, const RbtBaseGrid& grid);
  

  ////////////////////////////////////////
  //Virtual functions for reading/writing grid data to streams in
  //text and binary format
  //Subclasses should provide their own private OwnPrint,OwnWrite, OwnRead
  //methods to handle subclass data members, and override the public
  //Print,Write and Read methods
	virtual void Print(ostream& ostr) const;//Text output
	virtual void Write(ostream& ostr) const;//Binary output (serialisation)
	virtual void Read(istream& istr);//Binary input, replaces existing grid

  ////////////////////////////////////////
  //Public methods
  ////////////////

  /////////////////////////
  //Get attribute functions
  /////////////////////////

  RbtUInt GetNX() const {return m_NX;}
  RbtUInt GetNY() const {return m_NY;}
  RbtUInt GetNZ() const {return m_NZ;}
  RbtUInt GetN() const {return m_N;}
  RbtUInt GetStrideX() const {return m_SX;}
  RbtUInt GetStrideY() const {return m_SY;}
  RbtUInt GetStrideZ() const {return m_SZ;}
  RbtInt GetnXMin() const {return m_nXMin;}
  RbtInt GetnYMin() const {return m_nYMin;}
  RbtInt GetnZMin() const {return m_nZMin;}
  RbtInt GetnXMax() const {return m_nXMax;}
  RbtInt GetnYMax() const {return m_nYMax;}
  RbtInt GetnZMax() const {return m_nZMax;}

  const RbtCoord& GetGridMin() const {return m_min;}
  const RbtCoord& GetGridMax() const {return m_max;}
  const RbtVector& GetGridStep() const {return m_step;}
  RbtCoord GetGridCenter() const {return (GetGridMax()+GetGridMin())/2.0;}
  RbtCoord GetGridSize() const {return GetGridMax()-GetGridMin();}
  RbtUInt GetPad() const {return m_NPad;}
  const RbtCoord& GetPadMin() const {return m_padMin;}
  const RbtCoord& GetPadMax() const {return m_padMax;}


  /////////////////////////
  //Set attribute functions
  /////////////////////////

  //Set grid min without changing the grid step
  void SetGridMin(const RbtCoord& gridMin);
  //Translate the grid by the given vector
  void TranslateGrid(const RbtVector& vec) {SetGridMin(GetGridMin()+vec);}
  //Set the center of the grid to the given coord
  void SetGridCenter(const RbtCoord& gridCenter) {TranslateGrid(gridCenter-GetGridCenter());}
  //Set grid step
  void SetGridStep(const RbtVector& gridStep);
  //Set a pad region (cuboid inside the grid which defines the valid coords)
  void SetPad(RbtUInt NPad=0);


  /////////////////////////
  //Bounds checking and index conversion routines
  /////////////////////////

  //When checking coords, take account of pad region
  RbtBool isValid(const RbtCoord& c) const {return (c >= m_padMin) && (c < m_padMax);}
  //When checking 3-D indices, take account of pad region
  RbtBool isValid(RbtUInt iX, RbtUInt iY, RbtUInt iZ) const {return (iX >= m_NPad+1) && (iX <= m_NX-m_NPad) &&
							       (iY >= m_NPad+1) && (iY <= m_NY-m_NPad) &&
							       (iZ >= m_NPad+1) && (iZ <= m_NZ-m_NPad);}
  //When checking 1-D index, don't take account of pad region
  RbtBool isValid(RbtUInt iXYZ) const {return (iXYZ < m_N);}

  //Return iX,iY,iZ indices for given coord
  RbtUInt GetIX(const RbtCoord& c) const {return int((c.x-m_min.x)/m_step.x)+1;}
  RbtUInt GetIY(const RbtCoord& c) const {return int((c.y-m_min.y)/m_step.y)+1;}
  RbtUInt GetIZ(const RbtCoord& c) const {return int((c.z-m_min.z)/m_step.z)+1;}

  //Return iX,iY,iZ indices for given coord
  RbtUInt GetIX(RbtDouble x) const {return int((x-m_min.x)/m_step.x)+1;}
  RbtUInt GetIY(RbtDouble y) const {return int((y-m_min.y)/m_step.y)+1;}
  RbtUInt GetIZ(RbtDouble z) const {return int((z-m_min.z)/m_step.z)+1;}
  
  //Return iX,iY,iZ indices for given iXYZ index
  RbtUInt GetIX(RbtUInt iXYZ) const {return iXYZ/m_SX+1;}
  RbtUInt GetIY(RbtUInt iXYZ) const {return (iXYZ%m_SX)/m_SY+1;}
  RbtUInt GetIZ(RbtUInt iXYZ) const {return (iXYZ%m_SY)/m_SZ+1;}

  //Return iXYZ index for given iX,iY,iZ indices
  RbtUInt GetIXYZ(RbtUInt iX, RbtUInt iY, RbtUInt iZ) const {return (iX-1)*m_SX+(iY-1)*m_SY+(iZ-1)*m_SZ;}
  //Return iXYZ index for given coord
  RbtUInt GetIXYZ(const RbtCoord& c) const {return GetIXYZ(GetIX(c),GetIY(c),GetIZ(c));}

  //Returns real-world coordinate for given iX,iY,iZ indices
  RbtCoord GetCoord(RbtUInt iX, RbtUInt iY, RbtUInt iZ) const {
    return (RbtCoord(m_nXMin,m_nYMin,m_nZMin)+RbtCoord(iX-1,iY-1,iZ-1))*m_step;
  }
  //Returns real-world X coordinate for given iX index
  RbtDouble GetXCoord(RbtUInt iX) const {return (m_nXMin+(RbtInt)iX-1)*m_step.x;}
  //Returns real-world Y coordinate for given iY index
  RbtDouble GetYCoord(RbtUInt iY) const {return (m_nYMin+(RbtInt)iY-1)*m_step.y;}
  //Returns real-world Z coordinate for given iZ index
  RbtDouble GetZCoord(RbtUInt iZ) const {return (m_nZMin+(RbtInt)iZ-1)*m_step.z;}

  //Returns real-world coordinate for given iXYZ index
  RbtCoord GetCoord(RbtUInt iXYZ) const {return GetCoord(GetIX(iXYZ),GetIY(iXYZ),GetIZ(iXYZ));}
  //Returns list of real-world coordinates for given set of iXYZ indices
  RbtCoordList GetCoordList(const RbtUIntSet& iXYZSet) const;

  //DM 17 May 1999 - returns the set of valid grid points within a sphere of a given center and radius
  //DM 17 Jul 2000 - use vector<RbtUInt> and return by reference, for performance boost
  void GetSphereIndices(const RbtCoord& c, RbtDouble radius, RbtUIntList& sIndices) const;



 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////

  //Protected method for writing data members for this class to text stream
	void OwnPrint(ostream& ostr) const;
  //Protected method for writing data members for this class to binary stream
	void OwnWrite(ostream& ostr) const;
	//Protected method for reading data members for this class from binary stream
	void OwnRead(istream& istr) throw (RbtError);
	

 private:
  ////////////////////////////////////////
  //Private methods
  /////////////////
  RbtBaseGrid(); //Disable default constructor
  
  //Helper function called by copy constructor and assignment operator
  void CopyGrid(const RbtBaseGrid&);

 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////


 private:
  ////////////////////////////////////////
  //Private data
  //////////////
  RbtUInt m_NX;//Number of grid points in X direction (fixed in constructor)
  RbtUInt m_NY;//Number of grid points in Y direction (fixed in constructor)
  RbtUInt m_NZ;//Number of grid points in Z direction (fixed in constructor)
  RbtUInt m_N;//Total number of grid points (fixed in constructor)
  RbtUInt m_SX;//Stride of X
  RbtUInt m_SY;//Stride of Y
  RbtUInt m_SZ;//Stride of Z
  RbtCoord m_min;//Minimum grid coords (real-world units)
  RbtCoord m_max;//Maximum grid coords (real-world units)
  RbtVector m_step;//Grid step in X,Y,Z (real-world units)
  RbtUInt m_NPad;//Defines a zero-padded border around the grid
  RbtCoord m_padMin;//Minimum pad coords (used to define an active region of the grid)
  RbtCoord m_padMax;//Maximum pad coords (used to define an active region of the grid)
  RbtInt m_nXMin;//Minimum X grid point (in units of grid step)
  RbtInt m_nXMax;//Maximum X grid point (in units of grid step)
  RbtInt m_nYMin;//Minimum Y grid point (in units of grid step)
  RbtInt m_nYMax;//Maximum Y grid point (in units of grid step)
  RbtInt m_nZMin;//Minimum Z grid point (in units of grid step)
  RbtInt m_nZMax;//Maximum Z grid point (in units of grid step)
};

//Useful typedefs
typedef SmartPtr<RbtBaseGrid> RbtBaseGridPtr;//Smart pointer

#endif //_RBTBASEGRID_H_
