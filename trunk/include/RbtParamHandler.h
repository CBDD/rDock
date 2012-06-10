/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtParamHandler.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Base implementation class for managing a collection of named parameters
//Any class requiring parameter handling can derive from RbtParamHandler
//Parameters are stored as RbtVariants (double, string or stringlist)
//Only derived classes can add and delete parameters from the collection

#ifndef _RBTPARAMHANDLER_H_
#define _RBTPARAMHANDLER_H_

#include "RbtConfig.h"
#include "RbtVariant.h"

class RbtParamHandler
{
	public:
		////////////////////////////////////////
		//Constructors/destructors
		virtual ~RbtParamHandler(); //Default destructor


		////////////////////////////////////////
		//Public methods
		////////////////
		//Get number of stored parameters
		RbtUInt GetNumParameters() const;
		//Get a named parameter, throws error if name not found
		RbtVariant GetParameter(const RbtString& strName) const throw (RbtBadArgument);
		//Check if named parameter is present
		RbtBool isParameterValid(const RbtString& strName) const;
		//Get list of all parameter names
		RbtStringList GetParameterNames() const;
		//Get list of all parameter
		RbtStringVariantMap GetParameters() const;
		
		//Set named parameter to new value, throws error if name not found
		void SetParameter(const RbtString& strName, const RbtVariant& vValue) throw (RbtBadArgument);

    //Virtual function for dumping parameters to an output stream
  	//Called by operator <<
  	virtual void Print(ostream& s) const;
		
	protected:
		////////////////////////////////////////
		//Protected methods
		///////////////////
		RbtParamHandler(); //Default constructor
		
		//Only derived classes can mess with the parameter list
		void AddParameter(const RbtString& strName, const RbtVariant& vValue);
		void DeleteParameter(const RbtString& strName);
		void ClearParameters();
		//DM 25 Oct 2000 - ParameterUpdated is invoked whenever SetParameter is called
		//to allow derived classes to manage a data member which tracks the param value
		//Useful for performance purposes as there is quite an overhead in finding a string
		//in a map, then converting from a Variant to the native datatype
		//Base class version does nothing
		virtual void ParameterUpdated(const RbtString& strName) {};
	
	private:
		////////////////////////////////////////
		//Private methods
		/////////////////
	
		RbtParamHandler(const RbtParamHandler&);//Copy constructor disabled by default
		RbtParamHandler& operator=(const RbtParamHandler&);//Copy assignment disabled by default
	
	
	protected:
		////////////////////////////////////////
		//Protected data
		////////////////

	
	private:
		////////////////////////////////////////
		//Private data
		//////////////
		RbtStringVariantMap m_parameters;	
};

////////////////////////////////////////
//Non-member functions

//Insertion operator (primarily for debugging)
ostream& operator<<(ostream& s, const RbtParamHandler& ph);

#endif //_RBTPARAMHANDLER_H_
