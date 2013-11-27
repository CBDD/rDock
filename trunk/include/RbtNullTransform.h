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
