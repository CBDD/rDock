/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtTransformFactory.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Factory class for creating transform objects. Transforms are specified by
//string names (as defined in each class's static _CT string)
//Directly analogous to RbtSFFactory for creating scoring function objects

#ifndef _RBTTRANSFORMFACTORY_H_
#define _RBTTRANSFORMFACTORY_H_

#include "RbtConfig.h"
#include "RbtTransformAgg.h"
#include "RbtParameterFileSource.h"

class RbtTransformFactory
{
	//Parameter name which identifies a scoring function definition 
	static RbtString _TRANSFORM;
	public:
	////////////////////////////////////////
	//Constructors/destructors
      
	RbtTransformFactory(); //Default constructor

	virtual ~RbtTransformFactory(); //Default destructor


	////////////////////////////////////////
	//Public methods
	////////////////
	
	//Creates a single transform object of type strTransformClass, and name strName
	//e.g. strTransformClass = RbtSimAnnTransform
	virtual RbtBaseTransform* Create(const RbtString& strTransformClass, const RbtString& strName) throw (RbtError);
	
	//Creates an aggregate transform from a parameter file source
	//Each component transform is in a named section, which should minimally contain a TRANSFORM parameter
	//whose value is the class name to instantiate
	//strTransformClasses contains a comma-delimited list of transform class names to instantiate
	//If strTransformClasses is empty, all named sections in spPrmSource are scanned for valid transform definitions
	//Transform parameters and scoring function requests are set from the list of parameters in each named section
	virtual RbtTransformAgg* CreateAggFromFile(RbtParameterFileSourcePtr spPrmSource,
																			const RbtString& strName,
																			const RbtString& strTransformClasses = RbtString()) throw (RbtError);
      

	protected:
	////////////////////////////////////////
	//Protected methods
	///////////////////
      
      
	private:
	////////////////////////////////////////
	//Private methods
	/////////////////
      
	RbtTransformFactory(const RbtTransformFactory&);//Copy constructor disabled by default
      
	RbtTransformFactory& operator=(const RbtTransformFactory&);//Copy assignment disabled by default
      
      
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
typedef SmartPtr<RbtTransformFactory> RbtTransformFactoryPtr;//Smart pointer

#endif //_RBTTRANSFORMFACTORY_H_
