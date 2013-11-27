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

//Class for writing RbtModels to PSF files

#ifndef _RBTPSFFILESINK_H_
#define _RBTPSFFILESINK_H_

#include "RbtBaseMolecularFileSink.h"

class RbtPsfFileSink : public RbtBaseMolecularFileSink
{
 public:
  ////////////////////////////////////////
  //Constructors/destructors
  //RbtPsfFileSink(const char* fileName, RbtModelPtr spModel);
  RbtPsfFileSink(const RbtString& fileName, RbtModelPtr spModel);

  virtual ~RbtPsfFileSink(); //Default destructor


  ////////////////////////////////////////
  //Public methods
  ////////////////
  //Override public methods from RbtBaseFileSink
  virtual void Render() throw (RbtError);

 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////

 private:
  ////////////////////////////////////////
  //Private methods
  /////////////////

  RbtPsfFileSink(); //Disable default constructor
  RbtPsfFileSink(const RbtPsfFileSink&);//Copy constructor disabled by default
  RbtPsfFileSink& operator=(const RbtPsfFileSink&);//Copy assignment disabled by default


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
typedef SmartPtr<RbtPsfFileSink> RbtPsfFileSinkPtr;//Smart pointer

#endif //_RBTPSFFILESINK_H_
