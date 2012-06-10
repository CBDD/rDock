/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtObserver.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Abstract Observer class for Subject-Observer design pattern, see RbtSubject.h
//Design considerations:
//Subclasses should override the Update and Deleted methods, which are invoked
//by the subject when the subject changes state or is about to be deleted,
//respectively.
//Observer subclasses also need to store a pointer to the concrete subject
//subclass (e.g. RbtWorkSpace), either via Register and Unregister methods or
//via the constructor, so that subject subclass methods can be called to get
//details of the change in state. No pointer to the base subject object is
//stored in the base observer object.

#ifndef _RBTOBSERVER_H_
#define _RBTOBSERVER_H_

#include <vector>

class RbtSubject;

class RbtObserver
{
	public:
	////////////////////////////////////////
	//Constructors/destructors
	virtual ~RbtObserver(); //Default destructor


	////////////////////////////////////////
	//Public methods
	////////////////
	//PURE VIRTUAL - DERIVED CLASSES MUST OVERRIDE
	//Notify observer that subject has changed
	virtual void Update(RbtSubject* theChangedSubject) = 0;
	//Notify observer that subject is about to be deleted
	virtual void Deleted(RbtSubject* theDeletedSubject) = 0;
	
	
	protected:
	////////////////////////////////////////
	//Protected methods
	///////////////////
	RbtObserver();
	
	
	private:
	////////////////////////////////////////
	//Private methods
	/////////////////
	RbtObserver(const RbtObserver&);//Copy constructor disabled by default
	RbtObserver& operator=(const RbtObserver&);//Copy assignment disabled by default
	
	
	protected:
	////////////////////////////////////////
	//Protected data
	////////////////

	
	private:
	////////////////////////////////////////
	//Private data
	//////////////	
};

//Useful typedefs
typedef std::vector<RbtObserver*> RbtObserverList;
typedef RbtObserverList::iterator RbtObserverListIter;
typedef RbtObserverList::const_iterator RbtObserverListConstIter;

#endif //_RBTOBSERVER_H_
