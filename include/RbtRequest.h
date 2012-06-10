/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtRequest.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
