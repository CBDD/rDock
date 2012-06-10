/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtRequestHandler.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
