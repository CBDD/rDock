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

//Manages an array of floats representing a 3-D lattice of grid points

#ifndef _RBTREALGRID_H_
#define _RBTREALGRID_H_

#include "RbtBaseGrid.h"

class RbtRealGrid : public RbtBaseGrid
{
 public:
	//Class type string
	static RbtString _CT;

  ////////////////////////////////////////
  //Constructors/destructors
  //Construct a NXxNYxNZ grid running from gridMin at gridStep resolution
  RbtRealGrid(const RbtCoord& gridMin, const RbtCoord& gridStep,
	     RbtUInt NX, RbtUInt NY, RbtUInt NZ, RbtUInt NPad=0);

  //Constructor reading all params from binary stream
  RbtRealGrid(istream& istr);

  ~RbtRealGrid(); //Default destructor

  //Copy constructor
  RbtRealGrid(const RbtRealGrid&);
  //Copy constructor taking a base class argument
  //Sets up the grid dimensions, then creates an empty data array
  RbtRealGrid(const RbtBaseGrid&);
  //Copy assignment
  RbtRealGrid& operator=(const RbtRealGrid&);
  //Copy assignment taking a base class argument
  //Sets up the grid dimensions, then creates an empty data array
  RbtRealGrid& operator=(const RbtBaseGrid&);

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
  float* GetGridData() const {return m_data;}


  /////////////////////////
  //Get/Set value functions
  /////////////////////////

  //Get/Set tolerance used for assessing value equality
  RbtDouble GetTolerance() const {return m_tol;}
  void SetTolerance(RbtDouble tol) {m_tol=tol;}

  //Get/Set single grid point value with bounds checking
  RbtDouble GetValue(const RbtCoord& c) const {return isValid(c) ? m_grid[GetIX(c)][GetIY(c)][GetIZ(c)] : 0.0;}
  RbtDouble GetValue(RbtUInt iX, RbtUInt iY, RbtUInt iZ) const {return isValid(iX,iY,iZ) ? m_grid[iX][iY][iZ] : 0.0;}
  RbtDouble GetValue(RbtUInt iXYZ) const {return isValid(iXYZ) ? m_data[iXYZ] : 0.0;}
  
  //DM 20 Jul 2000 - get values smoothed by trilinear interpolation
  //D. Oberlin and H.A. Scheraga, J. Comp. Chem. (1998) 19, 71.
  RbtDouble GetSmoothedValue(const RbtCoord& c) const;

  void SetValue(const RbtCoord& c, RbtDouble val) {if (isValid(c)) m_grid[GetIX(c)][GetIY(c)][GetIZ(c)] = val;}
  void SetValue(RbtUInt iX, RbtUInt iY, RbtUInt iZ, RbtDouble val) {if (isValid(iX,iY,iZ)) m_grid[iX][iY][iZ] = val;}
  void SetValue(RbtUInt iXYZ, RbtDouble val) {if (isValid(iXYZ)) m_data[iXYZ] = val;}

  //Set all grid points to the given value
  void SetAllValues(RbtDouble val);
  //Replaces all grid points with a given value (+/- tolerance) with newVal
  void ReplaceValue(RbtDouble oldVal, RbtDouble newVal) {ReplaceValueRange(oldVal-m_tol,oldVal+m_tol,newVal);}
  //Replaces all grid points between oldValMin and oldValMax with newVal
  void ReplaceValueRange(RbtDouble oldValMin, RbtDouble oldValMax, RbtDouble newVal);

  //Set all grid points within radius of coord to the given value
  //If bOverwrite is false, does not replace non-zero values
  //If bOverwrite is true, all grid points are set the new value
  void SetSphere(const RbtCoord& c, RbtDouble radius, RbtDouble val, RbtBool bOverwrite=true);

  //Set all grid points with radii between rad1 and rad2 from coord to the given value
  //If bOverwrite is false, does not replace non-zero values
  //If bOverwrite is true, all grid points are set the new value
  void SetSurface(const RbtCoord& c, RbtDouble innerRad, RbtDouble outerRad, RbtDouble val, RbtBool bOverwrite=true);

  //Sets all grid points with value=oldValue, which are adjacent to those with value=adjacentValue, to value=newValue
  //+/- tolerance is applied to oldValue and adjacentValue
  void CreateSurface(RbtDouble oldVal, RbtDouble adjVal, RbtDouble newVal);

  //DM 16 Apr 1999 - helper function for determining solvent accessible regions
  //Returns true if any of the grid points within a sphere around the central coord have the specified value
  RbtBool isValueWithinSphere(const RbtCoord& c, RbtDouble radius, RbtDouble val);
  //Sets all grid points with value=oldValue, which have no grid points with value=adjacentValue within a sphere
  //of given radius, to value=newValue
  //+/- tolerance is applied to oldValue and adjacentValue
  //If bCenterOnly is true, just the center of the sphere is set the newValue
  //If bCenterOnly is false, all grid points in the sphere are set to the newValue
  void SetAccessible(RbtDouble radius, RbtDouble oldVal, RbtDouble adjVal, RbtDouble newVal, RbtBool bCenterOnly=true);


  /////////////////////////
  //Statistical functions
  /////////////////////////

  //Returns number of occurrences of a given value (+/- tolerance)
  RbtUInt Count(RbtDouble val) const {return CountRange(val-m_tol,val+m_tol);}
  //Returns number of occurrences of a given value range
  RbtUInt CountRange(RbtDouble valMin, RbtDouble valMax) const;

  //Min/max values
  RbtDouble MinValue() const;
  RbtDouble MaxValue() const;

  RbtUInt FindMinValue() const;//iXYZ index of grid point with minimum value
  RbtUInt FindMaxValue() const;//iXYZ index of grid point with maximum value


  /////////////////////////
  //I/O functions
  /////////////////////////

  //Dump grid in a format readable by Insight
  void PrintInsightGrid(ostream& s) const;


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
  RbtRealGrid(); //Disable default constructor

  //DM 17 Jul 2000 - analogous to isValueWithinSphere but iterates over arbitrary set
  //of IXYZ indices. Private method as there is no error checking on iXYZ values out of bounds
  RbtBool isValueWithinList(const RbtUIntList& iXYZList, RbtDouble val);
  //Analogous to SetSphere but iterates over arbitrary set
  //of IXYZ indices. Private method as there is no error checking on iXYZ values out of bounds
  //If bOverwrite is false, does not replace non-zero values
  //If bOverwrite is true, all grid points are set the new value
  void SetValues(const RbtUIntList& iXYZList, RbtDouble val, RbtBool bOverwrite=true);

  void CreateArrays();
  void ClearArrays();

  //Helper function called by copy constructor and assignment operator
  void CopyGrid(const RbtRealGrid&);

 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////


 private:
  ////////////////////////////////////////
  //Private data
  //////////////
  float*** m_grid;//3-D array, accessed as m_grid[i][j][k], indicies from 1
  float* m_data;//1-D view of same 3-D array, accessed as m_data[i], index from 0
  RbtDouble m_tol;//Tolerance for comparing grid values;
};

//Useful typedefs
typedef SmartPtr<RbtRealGrid> RbtRealGridPtr;//Smart pointer
typedef vector<RbtRealGridPtr> RbtRealGridList;
typedef RbtRealGridList::iterator RbtRealGridListIter;
typedef RbtRealGridList::const_iterator RbtRealGridListConstIter;

#endif //_RBTREALGRID_H_
