/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtBaseObject.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Base class for all workspace observers which have named parameters and
//handle arbitrary requests (i.e. scoring functions and transforms)

#ifndef _RBTBASEOBJECT_H_
#define _RBTBASEOBJECT_H_

#include "RbtConfig.h"
#include "RbtParamHandler.h"
#include "RbtObserver.h"
#include "RbtRequestHandler.h"

class RbtWorkSpace;//forward definition

class RbtBaseObject : public RbtParamHandler, public RbtObserver, public RbtRequestHandler
{
	public:
	
	//Class type string
	static RbtString _CT;
	//Parameter names
	static RbtString _CLASS;
	static RbtString _NAME;
	static RbtString _ENABLED;
	static RbtString _TRACE;//DM 1 Mar 2002 - move from RbtBaseTransform

	////////////////////////////////////////
	//Constructors/destructors
	virtual ~RbtBaseObject();
	
	////////////////////////////////////////
	//Public methods
	////////////////
	//Class name (e.g. RbtConstSF)
	RbtString GetClass() const;
	//Short name (e.g. HBOND)
	RbtString GetName() const;
	void SetName(const RbtString&);
	//Fully qualified name (should be overridden by subclasses which can be aggregated
	//to prefix the name with the parent's name)
	virtual RbtString GetFullName() const;
	void Enable();
	void Disable();
	RbtBool isEnabled() const;

	//DM 1 Mar 2002 - move from RbtBaseTransform
	//get/set the trace level for debugging
	RbtInt GetTrace() const;
	void SetTrace(RbtInt);

	//WorkSpace handling methods
	//Register scoring function with a workspace
	//Base class version just registers itself
	virtual void Register(RbtWorkSpace*);
	//Unregister with a workspace
	//Base class version just unregisters itself
	virtual void Unregister();
	//Get workspace pointer
	RbtWorkSpace* GetWorkSpace() const;
	
	//Override Observer method
	//Notify observer that subject is about to be deleted
	virtual void Deleted(RbtSubject* theDeletedSubject);

	//Request Handling method
	virtual void HandleRequest(RbtRequestPtr spRequest);
	
	//Virtual function for dumping parameters to an output stream
	//Called by operator <<
	virtual void Print(ostream& s) const;
		
	protected:
	////////////////////////////////////////
	//Protected methods
	///////////////////
 	RbtBaseObject(const RbtString& strClass, const RbtString& strName);
 	RbtBaseObject();
	
	//DM 25 Oct 2000 - track changes to parameter values in local data members
	//ParameterUpdated is invoked by RbtParamHandler::SetParameter
	void ParameterUpdated(const RbtString& strName);
	
	private:
	////////////////////////////////////////
	//Private methods
	/////////////////
 	RbtBaseObject(const RbtBaseObject&);//Copy constructor disabled by default      
	RbtBaseObject& operator=(const RbtBaseObject&);//Copy assignment disabled by default

	protected:
	////////////////////////////////////////
	//Protected data
	////////////////

      
	private:
	////////////////////////////////////////
	//Private data
	//////////////
	RbtWorkSpace* m_workspace;
	RbtBool m_enabled;
	RbtInt m_trace;
};

//Useful typedefs
typedef vector<RbtBaseObject*> RbtBaseObjectList;//Vector of smart pointers
typedef RbtBaseObjectList::iterator RbtBaseObjectListIter;
typedef RbtBaseObjectList::const_iterator RbtBaseObjectListConstIter;

#endif //_RBTBASEOBJECT_H_
