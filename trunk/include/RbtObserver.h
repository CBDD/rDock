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
