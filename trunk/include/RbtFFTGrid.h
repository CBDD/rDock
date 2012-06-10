/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtFFTGrid.h#4 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Subclass of RbtRealGrid with extra methods for finding peaks
//The original FFT functionality has long since been abandoned

#ifndef _RBTFFTGRID_H_
#define _RBTFFTGRID_H_

#include "RbtRealGrid.h"

//struct for holding FFT Peak info
class RbtFFTPeak
{
 public:
  //Constructor
  RbtFFTPeak() : index(0),coord(0.0,0.0,0.0),volume(0),height(0.0) {}

  //Data attributes - public to avoid need for accessor functions
  RbtUInt index;//Location of peak maximum - iXYZ index into grid
  RbtCoord coord;//Location of peak maximum - real world vector/coord corresponding to displacement
  RbtDouble height;//Maximum value of peak
  RbtUInt volume;//Number of grid points in peak (=points.size())
  RbtUIntSet points;//Set of iXYZ peak indices
};

typedef SmartPtr<RbtFFTPeak> RbtFFTPeakPtr;//Smart pointer
//Multi-map of key=peak height,value=FFTPeakPtr
typedef multimap<RbtDouble,RbtFFTPeakPtr> RbtFFTPeakMap;
typedef RbtFFTPeakMap::iterator RbtFFTPeakMapIter;
typedef RbtFFTPeakMap::reverse_iterator RbtFFTPeakMapRIter;
typedef RbtFFTPeakMap::const_iterator RbtFFTPeakMapConstIter;
typedef RbtFFTPeakMap::const_reverse_iterator RbtFFTPeakMapConstRIter;


class RbtFFTGrid : public RbtRealGrid
{
 public:
	//Class type string
	static RbtString _CT;

  ////////////////////////////////////////
  //Constructors/destructors
  //Construct a NXxNYxNZ grid running from gridMin at gridStep resolution
  RbtFFTGrid(const RbtCoord& gridMin, const RbtCoord& gridStep,
	     RbtUInt NX, RbtUInt NY, RbtUInt NZ, RbtUInt NPad=0);

  //Constructor reading all params from binary stream
  RbtFFTGrid(istream& istr);

  ~RbtFFTGrid(); //Default destructor

  RbtFFTGrid(const RbtFFTGrid&);//Copy constructor
  //Copy constructors taking base class arguments
  RbtFFTGrid(const RbtRealGrid&);
  RbtFFTGrid(const RbtBaseGrid&);
  
  RbtFFTGrid& operator=(const RbtFFTGrid&);//Copy assignment
  //Copy assignment taking base class arguments
  RbtFFTGrid& operator=(const RbtRealGrid&);
  RbtFFTGrid& operator=(const RbtBaseGrid&);

  ////////////////////////////////////////
  //Public methods
  ////////////////

  ////////////////////////////////////////
  //Virtual functions for reading/writing grid data to streams in
  //text and binary format
  //Subclasses should provide their own private OwnPrint,OwnWrite, OwnRead
  //methods to handle subclass data members, and override the public
  //Print,Write and Read methods
	virtual void Print(ostream& ostr) const;//Text output
	virtual void Write(ostream& ostr) const;//Binary output (serialisation)
	virtual void Read(istream& istr);//Binary input, replaces existing grid

  //Find the coords of all (separate) peaks above the threshold value
  //whose volumes are not less than minVol
  //Returns a map of RbtFFTPeaks
  RbtFFTPeakMap FindPeaks(RbtDouble threshold, RbtUInt minVol=1) const;
  //Returns the grid point with the maximum value in RbtFFTPeak format
  //Just a wrapper around FindMaxValue() (see below)
  RbtFFTPeak FindMaxPeak() const;


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
  RbtFFTGrid(); //Disable default constructor

  //Helper function called by copy constructor and assignment operator
  void CopyGrid(const RbtFFTGrid&);


 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////


 private:
  ////////////////////////////////////////
  //Private data
  //////////////
};

//Useful typedefs
typedef SmartPtr<RbtFFTGrid> RbtFFTGridPtr;//Smart pointer

#endif //_RBTFFTGRID_H_
