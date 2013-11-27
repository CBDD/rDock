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

//Randomises all rotable bonds of the ligand, within specified +/- range
//No bump check on final conformation

#ifndef _RBTRANDLIGTRANSFORM_H_
#define _RBTRANDLIGTRANSFORM_H_

#include "RbtBaseUniMolTransform.h"
#include "RbtRand.h"

class RbtRandLigTransform : public RbtBaseUniMolTransform
{
 public:
  //Static data member for class type
  static RbtString _CT;
  //Parameter names
  static RbtString _TORS_STEP;
  
  ////////////////////////////////////////
  //Constructors/destructors
  RbtRandLigTransform(const RbtString& strName = "RANDLIG");
  virtual ~RbtRandLigTransform();

  ////////////////////////////////////////
  //Public methods
  ////////////////	
  
 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////
  virtual void SetupTransform();//Called by Update when model has changed
  virtual void Execute();

 private:
  ////////////////////////////////////////
  //Private methods
  /////////////////
  RbtRandLigTransform(const RbtRandLigTransform&);//Copy constructor disabled by default
  RbtRandLigTransform& operator=(const RbtRandLigTransform&);//Copy assignment disabled by default

 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////


 private:
  ////////////////////////////////////////
  //Private data
  //////////////
  RbtRand& m_rand;//keep a reference to the singleton random number generator
  RbtBondList m_rotableBonds;
};

//Useful typedefs
typedef SmartPtr<RbtRandLigTransform> RbtRandLigTransformPtr;//Smart pointer

#endif //_RBTRANDLIGTRANSFORM_H_
