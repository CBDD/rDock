/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtFFTGrid.cxx#4 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include <algorithm> //for min, max, count
#include <iomanip>
#include <queue>

#include "RbtFFTGrid.h"
#include "RbtFileError.h"

//Static data members
RbtString RbtFFTGrid::_CT("RbtFFTGrid");

////////////////////////////////////////
//Constructors/destructors
//Construct a NXxNYxNZ grid running from gridMin at gridStep resolution
RbtFFTGrid::RbtFFTGrid(const RbtCoord& gridMin, const RbtCoord& gridStep,
		       RbtUInt NX, RbtUInt NY, RbtUInt NZ, RbtUInt NPad) :
  RbtRealGrid(gridMin,gridStep,NX,NY,NZ,NPad)
{
  _RBTOBJECTCOUNTER_CONSTR_("RbtFFTGrid");
}

//Constructor reading params from binary stream
RbtFFTGrid::RbtFFTGrid(istream& istr) :
                        RbtRealGrid(istr) {
  OwnRead(istr);
  _RBTOBJECTCOUNTER_CONSTR_("RbtFFTGrid");
}

//Default destructor
RbtFFTGrid::~RbtFFTGrid()
{
  _RBTOBJECTCOUNTER_DESTR_("RbtFFTGrid");
}

//Copy constructor
RbtFFTGrid::RbtFFTGrid(const RbtFFTGrid& grid) : RbtRealGrid(grid)
{
  CopyGrid(grid);
  _RBTOBJECTCOUNTER_COPYCONSTR_("RbtFFTGrid");
}

//Copy constructors taking base class arguments
RbtFFTGrid::RbtFFTGrid(const RbtRealGrid& grid) : RbtRealGrid(grid)
{
  _RBTOBJECTCOUNTER_COPYCONSTR_("RbtFFTGrid");
}

RbtFFTGrid::RbtFFTGrid(const RbtBaseGrid& grid) : RbtRealGrid(grid)
{
  _RBTOBJECTCOUNTER_COPYCONSTR_("RbtFFTGrid");
}

//Copy assignment
RbtFFTGrid& RbtFFTGrid::operator=(const RbtFFTGrid& grid)
{
  if (this != &grid) {
    RbtRealGrid::operator=(grid);
    CopyGrid(grid);
  }
  return *this;
}
//Copy assignment taking base class arguments
RbtFFTGrid& RbtFFTGrid::operator=(const RbtRealGrid& grid)
{
  if (this != &grid) {
    RbtRealGrid::operator=(grid);
  }
  return *this;
}
//Copy assignment taking base class arguments
RbtFFTGrid& RbtFFTGrid::operator=(const RbtBaseGrid& grid)
{
  if (this != &grid) {
    RbtRealGrid::operator=(grid);
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
void RbtFFTGrid::Print(ostream& ostr) const
{
  RbtRealGrid::Print(ostr);
  OwnPrint(ostr);
}
 
//Binary output
void RbtFFTGrid::Write(ostream& ostr) const
{
  RbtRealGrid::Write(ostr);
  OwnWrite(ostr);
}

//Binary input
void RbtFFTGrid::Read(istream& istr)
{
  RbtRealGrid::Read(istr);
  OwnRead(istr);
}

//Find the coords of all (separate) peaks above the threshold value
//whose volumes are not less than minVol
//Returns a map of RbtFFTPeaks
RbtFFTPeakMap RbtFFTGrid::FindPeaks(RbtDouble threshold, RbtUInt minVol) const
{
  RbtFFTPeakMap peakMap;//Initialise the return set

  //First compile a list of all points higher than the threshold
  RbtUIntSet stillToProcess;
  float* nrData = GetGridData();

  //DM 21 Jan 2000 - take account of tolerance when assessing threshold
  threshold -= GetTolerance();
  
  for (RbtUInt i=0; i < GetN(); i++) {
    if (nrData[i] > threshold)
      stillToProcess.insert(i);
  }

#ifdef _DEBUG
  cout << stillToProcess.size() << " data points found higher than  " << threshold << endl;
#endif //_DEBUG

  //Repeat while we still have data points to process
  while (!stillToProcess.empty()) {

    //Start a new peak set
    RbtUIntSet currentPeak;
    //Queue of points to be added to current peak
    std::queue<RbtUInt> toAddToPeak;

    //Seed the queue with the first unprocessed data point
    RbtUIntSetIter iter = stillToProcess.begin();
    RbtUInt iXYZ0 = *iter;
    RbtUInt peakPos = iXYZ0;
    float peakHeight = nrData[peakPos];
    toAddToPeak.push(iXYZ0);
    stillToProcess.erase(iter);

#ifdef _DEBUG
    //cout << "Seeding new peak at point " << iXYZ0 << endl;
#endif //_DEBUG
    
    //Repeat while the peak keeps growing
    while (!toAddToPeak.empty()) {
      //Take a new point from the front of the pending queue
      RbtUInt iXYZ0 = toAddToPeak.front();
      RbtUInt iX0 = GetIX(iXYZ0);
      RbtUInt iY0 = GetIY(iXYZ0);
      RbtUInt iZ0 = GetIZ(iXYZ0);
      toAddToPeak.pop();
      //Check if new point is higher than the current maximum
      float f = nrData[iXYZ0];
      if (f > peakHeight) {
        peakHeight = f;
	      peakPos = iXYZ0;
      }
      //Add the point to the current peak
      currentPeak.insert(iXYZ0);
#ifdef _DEBUG
      //cout << "Popping point " << iXYZ0 << endl;
#endif //_DEBUG

      //Now check if any neighbours of the point can also be added to the queue

      //Is X+1 within range ?
      if ( iX0 < GetNX()) {
        RbtUInt iXYZ1 = iXYZ0+GetStrideX();//Add stride(X) to get to X+1
        //Is X+1 in the "still to process" set
        RbtUIntSetIter iter = stillToProcess.find(iXYZ1);
        if (iter != stillToProcess.end()) {
          toAddToPeak.push(iXYZ1);
	        stillToProcess.erase(iter);
#ifdef _DEBUG
	        //cout << "Pushing X+1 point " << iXYZ1 << endl;
#endif //_DEBUG
        }
      }
      //Is X-1 within range ?
      if ( iX0 > 1 ) {
        RbtUInt iXYZ1 = iXYZ0-GetStrideX();
        //Is X-1 in the "still to process" set
        RbtUIntSetIter iter = stillToProcess.find(iXYZ1);
        if (iter != stillToProcess.end()) {
          toAddToPeak.push(iXYZ1);
          stillToProcess.erase(iter);
#ifdef _DEBUG
          //cout << "Pushing X-1 point " << iXYZ1 << endl;
#endif //_DEBUG
        }
      }
      //Is Y+1 within range ?
      if ( iY0 < GetNY()) {
        RbtUInt iXYZ1 = iXYZ0+GetStrideY();
        //Is Y+1 in the "still to process" set
        RbtUIntSetIter iter = stillToProcess.find(iXYZ1);
        if (iter != stillToProcess.end()) {
          toAddToPeak.push(iXYZ1);
          stillToProcess.erase(iter);
#ifdef _DEBUG
	        //cout << "Pushing Y+1 point " << iXYZ1 << endl;
#endif //_DEBUG
        }
      }
      //Is Y-1 within range ?
      if ( iY0 > 1 ) {
        RbtUInt iXYZ1 = iXYZ0-GetStrideY();
        //Is Y-1 in the "still to process" set
        RbtUIntSetIter iter = stillToProcess.find(iXYZ1);
        if (iter != stillToProcess.end()) {
          toAddToPeak.push(iXYZ1);
          stillToProcess.erase(iter);
#ifdef _DEBUG
          //cout << "Pushing Y-1 point " << iXYZ1 << endl;
#endif //_DEBUG
        }
      }
      //Is Z+1 within range ?
      if ( iZ0 < GetNZ()) {
        RbtUInt iXYZ1 = iXYZ0+GetStrideZ();
        //Is Z+1 in the "still to process" set
        RbtUIntSetIter iter = stillToProcess.find(iXYZ1);
        if (iter != stillToProcess.end()) {
          toAddToPeak.push(iXYZ1);
          stillToProcess.erase(iter);
#ifdef _DEBUG
          //cout << "Pushing Z+1 point " << iXYZ1 << endl;
#endif //_DEBUG
        }
      }
      //Is Z-1 within range ?
      if ( iZ0 > 1 ) {
        RbtUInt iXYZ1 = iXYZ0-GetStrideZ();
        //Is Z-1 in the "still to process" set
        RbtUIntSetIter iter = stillToProcess.find(iXYZ1);
        if (iter != stillToProcess.end()) {
          toAddToPeak.push(iXYZ1);
          stillToProcess.erase(iter);
#ifdef _DEBUG
          //cout << "Pushing Z-1 point " << iXYZ1 << endl;
#endif //_DEBUG
        }
      }
    }

    //Store current peak in peak map if volume is not less than minVol
    if (currentPeak.size() >= minVol) {
      RbtFFTPeakPtr spPeak(new RbtFFTPeak);
      spPeak->index = peakPos;
      spPeak->height = peakHeight;
      spPeak->coord = GetCoord(peakPos);
      spPeak->volume = currentPeak.size();
      spPeak->points = currentPeak;
      peakMap.insert(std::pair<RbtDouble, RbtFFTPeakPtr>(peakHeight,spPeak));
#ifdef _DEBUG
      cout.precision(3);
      cout.setf(ios_base::fixed,ios_base::floatfield);
      cout << "Finished growing peak #" << peakMap.size() << ", pos=" << spPeak->coord
	   << ", height=" << spPeak->height << ", volume=" << spPeak->volume << endl << endl;
#endif //_DEBUG
    }
  }
  return peakMap;
}


//Returns the grid point with the maximum value in RbtFFTPeak format
//Just a wrapper around FindMaxValue() (see below)
RbtFFTPeak RbtFFTGrid::FindMaxPeak() const
{
  RbtFFTPeak maxPeak;
  RbtUInt peakPos = FindMaxValue();
  maxPeak.index = peakPos;
  maxPeak.height = GetValue(peakPos);
  maxPeak.coord = GetCoord(peakPos);
  maxPeak.volume = 1;
  maxPeak.points.insert(peakPos);
  return maxPeak;
}


///////////////////////////////////////////////////////////////////////////
//Protected methods

//Protected method for writing data members for this class to text stream
void RbtFFTGrid::OwnPrint(ostream& ostr) const {
  ostr << "Class\t" << _CT << endl;  
}

//Protected method for writing data members for this class to binary stream
//(Serialisation)
void RbtFFTGrid::OwnWrite(ostream& ostr) const {
  //Write the class name as a title so we can check the authenticity of streams
  //on read
  const char* const gridTitle = _CT.c_str();
  RbtInt length = strlen(gridTitle);
  Rbt::WriteWithThrow(ostr, (const char*) &length, sizeof(length));
  Rbt::WriteWithThrow(ostr, gridTitle, length);
}

//Protected method for reading data members for this class from binary stream
void RbtFFTGrid::OwnRead(istream& istr) throw (RbtError) {
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
}


///////////////////////////////////////////////////////////////////////////
//Private methods

//Helper function called by copy constructor and assignment operator
//No longer required
void RbtFFTGrid::CopyGrid(const RbtFFTGrid& grid) {
}

