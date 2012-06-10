/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtPsfFileSink.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
