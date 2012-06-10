/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtTransformAgg.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Transform aggregate class. Manages collection of child transforms

#ifndef _RBTTRANSFORMAGG_H_
#define _RBTTRANSFORMAGG_H_

#include "RbtBaseTransform.h"

//Only check transform aggregate assertions in debug build
#ifdef _NDEBUG
const RbtBool TRANSFORMAGG_CHECK = false;
#else
const RbtBool TRANSFORMAGG_CHECK = true;
#endif //_NDEBUG

//Useful typedefs

class RbtTransformAgg : public RbtBaseTransform
{
	public:
	//Static data member for class type (i.e. "RbtTransformAgg")
	static RbtString _CT;
	
	////////////////////////////////////////
	//Constructors/destructors
 	RbtTransformAgg(const RbtString& strName = "DOCK");
	virtual ~RbtTransformAgg();


	////////////////////////////////////////
	//Public methods
	////////////////

	//Aggregate handling methods
	virtual void Add(RbtBaseTransform*) throw (RbtError);
	virtual void Remove(RbtBaseTransform*) throw (RbtError);
	virtual RbtBool isAgg() const;
	virtual RbtUInt GetNumTransforms() const;
	virtual RbtBaseTransform* GetTransform(RbtUInt iTransform) const throw (RbtError);
	
	//WorkSpace handling methods
	//Register scoring function with a workspace
	//Aggregate version registers all children, but NOT itself
	//(Aggregates are just containers, and have no need for model information
	virtual void Register(RbtWorkSpace*);
	//Unregister with a workspace
	//Aggregate version unregisters all children, but NOT itself
	virtual void Unregister();
	
	//Override RbtObserver pure virtual
	//Notify observer that subject has changed
	//Does nothing in RbtTransformAgg as aggregates do not require updating
	virtual void Update(RbtSubject* theChangedSubject);
	
	//Request Handling method
	virtual void HandleRequest(RbtRequestPtr spRequest);
	
  //Virtual function for dumping transform details to an output stream
  //Called by operator <<
  virtual void Print(ostream& s) const;

	protected:
	////////////////////////////////////////
	//Protected methods
	///////////////////
	//Actually apply the transform
	//Aggregate version loops over all child transforms
  virtual void Execute();

	private:
	////////////////////////////////////////
	//Private methods
	/////////////////
 	RbtTransformAgg(const RbtTransformAgg&);//Copy constructor disabled by default      
	RbtTransformAgg& operator=(const RbtTransformAgg&);//Copy assignment disabled by default
                  
	protected:
	////////////////////////////////////////
	//Protected data
	////////////////

      
	private:
	////////////////////////////////////////
	//Private data
	//////////////
	RbtBaseTransformList m_transforms;
};

//Useful typedefs
typedef SmartPtr<RbtTransformAgg> RbtTransformAggPtr;//Smart pointer
typedef vector<RbtTransformAggPtr> RbtTransformAggList;//Vector of smart pointers
typedef RbtTransformAggList::iterator RbtTransformAggListIter;
typedef RbtTransformAggList::const_iterator RbtTransformAggListConstIter;

#endif //_RBTTRANSFORMAGG_H_
