/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtSubject.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
