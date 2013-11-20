/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtParamHandler.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtParamHandler.h"

RbtParamHandler::RbtParamHandler() {
  _RBTOBJECTCOUNTER_CONSTR_("RbtParamHandler");
}

RbtParamHandler::~RbtParamHandler() {
  _RBTOBJECTCOUNTER_DESTR_("RbtParamHandler");
}

////////////////////////////////////////
//Public methods
////////////////
//Get number of stored parameters
RbtUInt RbtParamHandler::GetNumParameters() const {return m_parameters.size();}

//Get a named parameter, throws error if name not found
RbtVariant RbtParamHandler::GetParameter(const RbtString& strName) const throw (RbtBadArgument)
{
	RbtStringVariantMapConstIter iter = m_parameters.find(strName);
	if (iter == m_parameters.end()) {
		throw RbtBadArgument(_WHERE_,"Undefined parameter "+strName);
	}
	else
		return (*iter).second;
}

//Check if named parameter is present
RbtBool RbtParamHandler::isParameterValid(const RbtString& strName) const
{
	return m_parameters.find(strName) != m_parameters.end();
}

//Get list of all parameter names
RbtStringList RbtParamHandler::GetParameterNames() const
{
	RbtStringList nameList;
	for (RbtStringVariantMapConstIter iter = m_parameters.begin(); iter != m_parameters.end(); iter++) {
		nameList.push_back((*iter).first);
	}
	return nameList;
}

//Get map of all parameters
RbtStringVariantMap RbtParamHandler::GetParameters() const
{
	return m_parameters;
}

//Set named parameter to new value, throws error if name not found
void RbtParamHandler::SetParameter(const RbtString& strName, const RbtVariant& vValue) throw (RbtBadArgument)
{
	if (!isParameterValid(strName)) {
		throw RbtBadArgument(_WHERE_,"Undefined parameter "+strName);
	}
	else {
		m_parameters[strName] = vValue;
		//DM 25 Oct 2000 - notify derived class that parameter has changed
		ParameterUpdated(strName);
	}
}

////////////////////////////////////////
//Protected methods
///////////////////
//Only derived classes can mess with the parameter list
void RbtParamHandler::AddParameter(const RbtString& strName, const RbtVariant& vValue)
{
	m_parameters[strName] = vValue;
	//DM 25 Oct 2000 - notify derived class that parameter has changed
	//DM 12 Apr 2002 - no need to call ParameterUpdated here
	//Causes problems during construction of derived classes with virtual base classes
	//Better solution is for derived classes to synchronise the initial values of private
	//data members themselves during construction
	//ParameterUpdated(strName);
}

void RbtParamHandler::DeleteParameter(const RbtString& strName)
{
	m_parameters.erase(strName);
}

void RbtParamHandler::ClearParameters()
{
	m_parameters.clear();
}

//Virtual function for dumping parameters to an output stream
//Called by operator <<
void RbtParamHandler::Print(ostream& s) const {
	for (RbtStringVariantMapConstIter iter = m_parameters.begin(); iter != m_parameters.end(); iter++) {
		cout << (*iter).first << "\t" << (*iter).second << endl;
	}
}

////////////////////////////////////////
//Non-member functions

//Insertion operator (primarily for debugging)
ostream& operator<<(ostream& s, const RbtParamHandler& ph) {
	ph.Print(s);
  return s;
}

