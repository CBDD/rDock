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

//Abstract Subject class for Subject-Observer design pattern. Subject will
//notify registered observers of any change of state of the Subject, by
//invoking Update method on each observer. It is up to the observer to
//interrogate subject as to what has changed.
//Main use is for RbtWorkspace which manages the model list and notifies
//scoring functions and transforms of any changes.
//
//Design considerations:
//Constructor is protected to prevent instantiation of base class.
//No smart pointers are used for storing observers.
//Subject is NOT responsible for deleting observers in destructor, but will
//notify registered observers of impending destruction by invoking Deleted
//method on each observer.

#ifndef _RBTSUBJECT_H_
#define _RBTSUBJECT_H_

#include "RbtConfig.h"
#include "RbtObserver.h"

class RbtSubject
{
public:
	////////////////////////////////////////
	//Constructors/destructors
	virtual ~RbtSubject();

	////////////////////////////////////////
	//Public methods
	////////////////
	virtual void Attach(RbtObserver*) throw (RbtError);
	virtual void Detach(RbtObserver*) throw (RbtError);
	virtual void Notify();
	
protected:
	////////////////////////////////////////
	//Protected methods
	///////////////////
	RbtSubject();	//Disable constructor
	
private:
	////////////////////////////////////////
	//Private methods
	/////////////////

protected:
	////////////////////////////////////////
	//Protected data
	////////////////

	
private:
	////////////////////////////////////////
	//Private data
	//////////////
	RbtObserverList m_observers;      
};

#endif //_RBTSUBJECT_H_
