/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtCrdFileSink.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
