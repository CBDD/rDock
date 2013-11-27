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

//Simple struct to represent a request message
//Stores a request ID (int) and an arbitrary length list of variant params
//Base class cannot be instantiated.
//Derived class constructors should set the ID and the param list accordingly

#ifndef _RBTREQUEST_H_
#define _RBTREQUEST_H_

#include "RbtVariant.h"

typedef RbtInt RbtRequestID;

//Only check request assertions in debug build
#ifdef _NDEBUG
const RbtBool REQ_CHECK = false;
#else
const RbtBool REQ_CHECK = true;
#endif //_NDEBUG

class RbtRequest
{
	public:
	////////////////////////////////////////
	//Constructors/destructors
	virtual ~RbtRequest() {}
	
	////////////////////////////////////////
	//Public methods
	////////////////
	RbtRequestID GetID() const {return m_id;}
	RbtUInt GetNumParameters() const {return m_parameters.size();}
	RbtVariantList GetParameters() const {return m_parameters;}

	////////////////////////////////////////
	//Protected methods
	////////////////
	//Protected constructor to be used by derived classes to set request ID
	RbtRequest(RbtRequestID id) : m_id(id) {}
	//AddParameter can be used by derived classes to add parameters to list
	void AddParameter(const RbtVariant& v) {
		m_parameters.push_back(v);
	}
	
	private:
	////////////////////////////////////////
	//Private data
	//////////////
	RbtRequestID m_id;
	RbtVariantList m_parameters;
};

//Useful typedefs
typedef SmartPtr<RbtRequest> RbtRequestPtr;
typedef vector<RbtRequestPtr> RbtRequestList;
typedef RbtRequestList::iterator RbtRequestListIter;
typedef RbtRequestList::const_iterator RbtRequestListConstIter;

#endif //_RBTREQUEST_H_
