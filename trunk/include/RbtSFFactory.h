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
