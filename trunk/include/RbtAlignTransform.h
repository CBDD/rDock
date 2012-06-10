/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtAlignTransform.h#2 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Aligns ligand with the principal axes of one of the active site cavities

#ifndef _RBTALIGNTRANSFORM_H_
#define _RBTALIGNTRANSFORM_H_

#include "RbtBaseBiMolTransform.h"
#include "RbtRand.h"
#include "RbtCavity.h"

class RbtAlignTransform : public RbtBaseBiMolTransform
{
 public:
	//Static data member for class type
	static RbtString _CT;
	//Parameter names
	static RbtString _COM;
	static RbtString _AXES;	

  ////////////////////////////////////////
  //Constructors/destructors
 	RbtAlignTransform(const RbtString& strName = "ALIGN");
	virtual ~RbtAlignTransform();

  ////////////////////////////////////////
  //Public methods
  ////////////////	
  
 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////
	virtual void SetupReceptor();//Called by Update when receptor is changed
	virtual void SetupLigand();//Called by Update when ligand is changed
	virtual void SetupTransform();//Called by Update when either model has changed
  virtual void Execute();

 private:
  ////////////////////////////////////////
  //Private methods
  /////////////////
  RbtAlignTransform(const RbtAlignTransform&);//Copy constructor disabled by default
  RbtAlignTransform& operator=(const RbtAlignTransform&);//Copy assignment disabled by default

 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////


 private:
  ////////////////////////////////////////
  //Private data
  //////////////
  RbtRand& m_rand;//keep a reference to the singleton random number generator
  RbtCavityList m_cavities;//List of active site cavities to choose from
  RbtIntList m_cumulSize;//Cumulative sizes, for weighted probabilities
  RbtInt m_totalSize;//Total size of all cavities
};

//Useful typedefs
typedef SmartPtr<RbtAlignTransform> RbtAlignTransformPtr;//Smart pointer

#endif //_RBTALIGNTRANSFORM_H_
