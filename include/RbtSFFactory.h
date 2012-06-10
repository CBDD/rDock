/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtSFFactory.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Factory class for creating scoring function objects. Scoring functions are
//specified by string names (as defined in each class's static _CT string)

#ifndef _RBTSFFACTORY_H_
#define _RBTSFFACTORY_H_

#include "RbtConfig.h"
#include "RbtSFAgg.h"
#include "RbtParameterFileSource.h"

class RbtSFFactory
{
	//Parameter name which identifies a scoring function definition 
	static RbtString _SF;
	public:
	////////////////////////////////////////
	//Constructors/destructors
      
	RbtSFFactory(); //Default constructor

	virtual ~RbtSFFactory(); //Default destructor


	////////////////////////////////////////
	//Public methods
	////////////////
	
	//Creates a single scoring function object of type strSFClass, and name strName
	//e.g. strSFClass = RbtHBondIntnSF
	virtual RbtBaseSF* Create(const RbtString& strSFClass, const RbtString& strName) throw (RbtError);
	
	//Creates an aggregate scoring function from a parameter file source
	//Each component SF is in a named section, which should minimally contain a SCORING_FUNCTION parameter
	//whose value is the class name to instantiate
	//strSFClasses contains a comma-delimited list of SF class names to instantiate
	//If strSFClasses is empty, all named sections in spPrmSource are scanned for valid scoring function definitions
	//SF parameters are set from the list of parameters in each named section
	virtual RbtSFAgg* CreateAggFromFile(RbtParameterFileSourcePtr spPrmSource,
																			const RbtString& strName,
																			const RbtString& strSFClasses = RbtString()) throw (RbtError);
      

	protected:
	////////////////////////////////////////
	//Protected methods
	///////////////////
      
      
	private:
	////////////////////////////////////////
	//Private methods
	/////////////////
      
	RbtSFFactory(const RbtSFFactory&);//Copy constructor disabled by default
      
	RbtSFFactory& operator=(const RbtSFFactory&);//Copy assignment disabled by default
      
      
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
typedef SmartPtr<RbtSFFactory> RbtSFFactoryPtr;//Smart pointer

#endif //_RBTSFFACTORY_H_
