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

#include "RbtBaseObject.h"
#include "RbtSFRequest.h"
#include "RbtWorkSpace.h"

//Static data members
RbtString RbtBaseObject::_CT("RbtBaseObject");
RbtString RbtBaseObject::_CLASS("CLASS");
RbtString RbtBaseObject::_NAME("NAME");
RbtString RbtBaseObject::_ENABLED("ENABLED");
RbtString RbtBaseObject::_TRACE("TRACE");

////////////////////////////////////////
//Constructors/destructors
RbtBaseObject::RbtBaseObject(const RbtString& strClass, const RbtString& strName)
  : m_workspace(NULL),m_enabled(true),m_trace(0) {
  //Add parameters
  AddParameter(_CLASS,strClass);
  AddParameter(_NAME,strName);
  AddParameter(_ENABLED,m_enabled);
  AddParameter(_TRACE,m_trace);//DM 1 Mar 2002 - move from RbtBaseTransform
#ifdef _DEBUG
  cout << _CT << " parameterised constructor for " << strClass << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

//Dummy default constructor for virtual base subclasses
//Should never get called
RbtBaseObject::RbtBaseObject() {
#ifdef _DEBUG
  cout << "WARNING: " << _CT << " default constructor" << endl;
#endif //_DEBUG
  //_RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtBaseObject::~RbtBaseObject() {
	Unregister();//Unregister object from workspace
#ifdef _DEBUG
	cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}


////////////////////////////////////////
//Public methods
////////////////

//Class name (e.g. RbtConstSF)
RbtString RbtBaseObject::GetClass() const {return GetParameter(_CLASS);}

RbtString RbtBaseObject::GetName() const {return GetParameter(_NAME);}
void RbtBaseObject::SetName(const RbtString& strName) {SetParameter(_NAME,strName);}
//Fully qualified name (should be overridden by subclasses which can be aggregated
//to prefix the name with the parent's name)
RbtString RbtBaseObject::GetFullName() const{return GetName();}

void RbtBaseObject::Enable() {SetParameter(_ENABLED,true);}
void RbtBaseObject::Disable() {SetParameter(_ENABLED,false);}
RbtBool RbtBaseObject::isEnabled() const {return m_enabled;}

RbtInt RbtBaseObject::GetTrace() const {return m_trace;}
void RbtBaseObject::SetTrace(RbtInt trace) {SetParameter(_TRACE,trace);}

//WorkSpace handling methods
//Register scoring function with a workspace
//Base class version just registers itself
void RbtBaseObject::Register(RbtWorkSpace* pWorkSpace)
{
	//First unregister with current workspace
	Unregister();
	m_workspace = pWorkSpace;
#ifdef _DEBUG
	cout << _CT << "::Register(): Registering " << GetName() << " with workspace" << endl;
#endif //_DEBUG
	m_workspace->Attach(this);
}

//Unregister with a workspace
//Base class version just unregisters itself
void RbtBaseObject::Unregister()
{
	if (m_workspace) {
#ifdef _DEBUG
		cout << _CT << "::Unregister(): Unregistering " << GetName() << " from workspace" << endl;
#endif //_DEBUG
		m_workspace->Detach(this);
	}
	m_workspace = NULL;
}

//Get workspace pointer
RbtWorkSpace* RbtBaseObject::GetWorkSpace() const {return m_workspace;}
	
//Override Observer method
//Notify observer that subject is about to be deleted
void RbtBaseObject::Deleted(RbtSubject* theDeletedSubject) {
	if (theDeletedSubject == m_workspace) {
#ifdef _DEBUG
		cout << _CT << "::Deleted(): " << GetName() << " received notice that subject is to be deleted " << endl;
#endif //_DEBUG
		Unregister();
	}
}

//Request Handling method
void RbtBaseObject::HandleRequest(RbtRequestPtr spRequest) {
  //Base class can handle ENABLE, DISABLE and SETPARAM
  RbtVariantList params = spRequest->GetParameters();
  switch (spRequest->GetID()) {
    
    //Enable object
    //params[0] is fully-qualified object name
  case ID_REQ_SF_ENABLE:
    Assert<RbtAssert>(!REQ_CHECK||params.size()==1);
    if (!params.empty() && (params[0].String() == GetFullName())) {
      if (m_trace > 2) {
	cout << _CT << "::HandleRequest: Enabling " << GetFullName() << endl;
      }
      Enable();
    }
    break;
			
    //Disable object
    //params[0] is fully-qualified object name
  case ID_REQ_SF_DISABLE:
    Assert<RbtAssert>(!REQ_CHECK||params.size()==1);
    if (!params.empty() && (params[0].String() == GetFullName())) {
      if (m_trace > 2) {
	cout << _CT << "::HandleRequest: Disabling " << GetFullName() << endl;
      }
      Disable();
    }
    break;
    
    //Set parameter has two flavours
  case ID_REQ_SF_SETPARAM:
    //Either:
    //params[0] is fully-qualified object name
    //params[1] is parameter name
    //params[2] is parameter value
    //Check if object has this named parameter before calling SetParameter
    if ( (params.size()==3) && (params[0].String() == GetFullName()) && isParameterValid(params[1])) {
      if (m_trace > 2) {
	cout << _CT << "::HandleRequest: Setting parameter "
	     << params[1] << " to " << params[2] << " for " << GetFullName() << endl;
      }
      SetParameter(params[1],params[2]);
    }
    //Or:
    //params[0] is parameter name
    //params[1] is parameter value
    //Check if object has this named parameter before calling SetParameter
    else if ( (params.size()==2) && isParameterValid(params[0])) {
      if (m_trace > 2) {
	cout << _CT << "::HandleRequest: Setting parameter "
	     << params[0] << " to " << params[1] << " for " << GetFullName() << endl;
      }
      SetParameter(params[0],params[1]);
    }
    break;
    
    //Ignore all other requests
  default:
    if (m_trace > 2) {
      cout << _CT << "::HandleRequest: " << GetFullName() << " is ignoring request" << endl;
    }
    break;
  }		
}
//Virtual function for dumping parameters to an output stream
//Called by operator <<
void RbtBaseObject::Print(ostream& s) const {
	cout << endl << GetFullName() << endl;
	RbtParamHandler::Print(s);
}

//DM 25 Oct 2000 - track changes to parameter values in local data members
//ParameterUpdated is invoked by RbtParamHandler::SetParameter
void RbtBaseObject::ParameterUpdated(const RbtString& strName) {
  if (strName == _ENABLED) {
    m_enabled = GetParameter(_ENABLED);
  }
  else if (strName == _TRACE) {
    m_trace = GetParameter(_TRACE);
  }
  else {
    RbtParamHandler::ParameterUpdated(strName);
  }
}
