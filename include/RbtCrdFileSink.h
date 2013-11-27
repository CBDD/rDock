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

//Class for writing RbtModels to CRD files

#ifndef _RBTCRDFILESINK_H_
#define _RBTCRDFILESINK_H_

#include "RbtBaseMolecularFileSink.h"

class RbtCrdFileSink : public RbtBaseMolecularFileSink
{
 public:
  ////////////////////////////////////////
  //Constructors/destructors
  //RbtCrdFileSink(const char* fileName, RbtModelPtr spModel);
  RbtCrdFileSink(const RbtString& fileName, RbtModelPtr spModel);

  virtual ~RbtCrdFileSink(); //Default destructor


  ////////////////////////////////////////
  //Public methods
  ////////////////
  //
  //Derived classes should override
  //Controls whether file type can support the writing of multiple conformations/models to a single file
  virtual RbtBool isMultiConfSupported() {return true;};
  //
  ////////////////////////////////////////
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

  RbtCrdFileSink(); //Disable default constructor
  RbtCrdFileSink(const RbtCrdFileSink&);//Copy constructor disabled by default
  RbtCrdFileSink& operator=(const RbtCrdFileSink&);//Copy assignment disabled by default


 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////

 private:
  ////////////////////////////////////////
  //Private data
  //////////////
  RbtUInt m_numAtomsLineRec;
};

//Useful typedefs
typedef SmartPtr<RbtCrdFileSink> RbtCrdFileSinkPtr;//Smart pointer

#endif //_RBTCRDFILESINK_H_
