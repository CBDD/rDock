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

//Scoring function specific requests

#ifndef _RBTSFREQUEST_H_
#define _RBTSFREQUEST_H_

#include "RbtRequest.h"

//Request IDs
const RbtRequestID ID_REQ_SF_ENABLE = 1;
const RbtRequestID ID_REQ_SF_DISABLE = 2;
const RbtRequestID ID_REQ_SF_PARTITION = 3;
const RbtRequestID ID_REQ_SF_SETPARAM = 4;

class RbtSFEnableRequest : public RbtRequest
{
	public:
	//Request to enable a particular named scoring function
	RbtSFEnableRequest(const RbtString& sfName) : RbtRequest(ID_REQ_SF_ENABLE) {
		AddParameter(sfName);
	}
};

class RbtSFDisableRequest : public RbtRequest
{
	public:
	//Request to disable a particular named scoring function
	RbtSFDisableRequest(const RbtString& sfName) : RbtRequest(ID_REQ_SF_DISABLE) {
		AddParameter(sfName);
	}
};

class RbtSFPartitionRequest : public RbtRequest
{
	public:
	//Request to partition all scoring functions to a particular distance
	RbtSFPartitionRequest(RbtDouble dist) : RbtRequest(ID_REQ_SF_PARTITION) {
		AddParameter(dist);
	}
	//Request to partition a particular named scoring function to a particular distance
	RbtSFPartitionRequest(const RbtString& sfName, RbtDouble dist) : RbtRequest(ID_REQ_SF_PARTITION) {
		AddParameter(sfName);
		AddParameter(dist);
	}
};

class RbtSFSetParamRequest : public RbtRequest
{
	public:
	//Request to set a named parameter of all scoring functions to a new value (e.g. TRACE)
	RbtSFSetParamRequest(const RbtString& paramName, const RbtVariant& paramValue) :
		RbtRequest(ID_REQ_SF_SETPARAM) {
		AddParameter(paramName);
		AddParameter(paramValue);
	}
	//Request to set a named parameter of a named scoring function to a new value
	RbtSFSetParamRequest(const RbtString& sfName, const RbtString& paramName, const RbtVariant& paramValue) :
		RbtRequest(ID_REQ_SF_SETPARAM) {
		AddParameter(sfName);
		AddParameter(paramName);
		AddParameter(paramValue);
	}
};


#endif //_RBTSFREQUEST_H_
