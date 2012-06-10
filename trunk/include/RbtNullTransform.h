/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtNullTransform.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Null transform which does nothing except fire off any requests to the
//scoring function

#ifndef _RBTNULLTRANSFORM_H_
#define _RBTNULLTRANSFORM_H_

#include "RbtBaseTransform.h"

class RbtNullTransform : public RbtBaseTransform
{
 public:
	//Static data member for class type
	static RbtString _CT;	

  ////////////////////////////////////////
  //Constructors/destructors
 	RbtNullTransform(const RbtString& strName = "NULL");
	virtual ~RbtNullTransform();

  ////////////////////////////////////////
  //Public methods
  ////////////////	
  
	virtual void Update(RbtSubject* theChangedSubject);//Does nothing

 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////
  virtual void Execute();//Does nothing

 private:
  ////////////////////////////////////////
  //Private methods
  /////////////////
  RbtNullTransform(const RbtNullTransform&);//Copy constructor disabled by default
  RbtNullTransform& operator=(const RbtNullTransform&);//Copy assignment disabled by default

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
typedef SmartPtr<RbtNullTransform> RbtNullTransformPtr;//Smart pointer

#endif //_RBTNULLTRANSFORM_H_
