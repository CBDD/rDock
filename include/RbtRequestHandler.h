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

//Abstract class for handling requests. Classes who want to handle requests
//should derive from RbtRequestHandler and override HandleRequest.
//Useful for implementing a Chain of Responsibility in aggregate containers

#ifndef _RBTREQUESTHANDLER_H_
#define _RBTREQUESTHANDLER_H_

#include "RbtRequest.h"

class RbtRequestHandler
{
	public:
	////////////////////////////////////////
	//Constructors/destructors
	RbtRequestHandler() {};
	virtual ~RbtRequestHandler() {};
	
	////////////////////////////////////////
	//Public methods
	////////////////
	virtual void HandleRequest(RbtRequestPtr spRequest) {};
};

#endif //_RBTREQUESTHANDLER_H_
