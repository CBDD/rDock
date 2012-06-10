/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtSFAgg.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Scoring Function aggregate class. Manages collection of child scoring
//functions.

#ifndef _RBTSFAGG_H_
#define _RBTSFAGG_H_

#include "RbtBaseSF.h"

//Only check SF aggregate assertions in debug build
#ifdef _NDEBUG
const RbtBool SFAGG_CHECK = false;
#else
const RbtBool SFAGG_CHECK = true;
#endif //_NDEBUG

//Useful typedefs

class RbtSFAgg : public RbtBaseSF
{
	public:
	//Static data member for class type (i.e. "RbtSFAgg")
	static RbtString _CT;
	
	////////////////////////////////////////
	//Constructors/destructors
 	RbtSFAgg(const RbtString& strName = "SCORE");
	virtual ~RbtSFAgg();


	////////////////////////////////////////
	//Public methods
	////////////////

  //Returns all child component scores as a string-variant map
  //Key = fully qualified component name, value = weighted score
  //(for saving in a Model's data fields)
  virtual void ScoreMap(RbtStringVariantMap& scoreMap) const;

	//Aggregate handling methods
	virtual void Add(RbtBaseSF*) throw (RbtError);
	virtual void Remove(RbtBaseSF*) throw (RbtError);
	virtual RbtBool isAgg() const;
	virtual RbtUInt GetNumSF() const;
	virtual RbtBaseSF* GetSF(RbtUInt iSF) const throw (RbtError);
	
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
	//Does nothing in RbtSFAgg as aggregates do not require updating
	virtual void Update(RbtSubject* theChangedSubject);
	
	//Request Handling method
	virtual void HandleRequest(RbtRequestPtr spRequest);
	
  //Virtual function for dumping scoring function details to an output stream
  //Called by operator <<
  virtual void Print(ostream& s) const;

	protected:
	////////////////////////////////////////
	//Protected methods
	///////////////////
  virtual RbtDouble RawScore() const;

	private:
	////////////////////////////////////////
	//Private methods
	/////////////////
 	RbtSFAgg(const RbtSFAgg&);//Copy constructor disabled by default      
	RbtSFAgg& operator=(const RbtSFAgg&);//Copy assignment disabled by default
                  
	protected:
	////////////////////////////////////////
	//Protected data
	////////////////

      
	private:
	////////////////////////////////////////
	//Private data
	//////////////
	RbtBaseSFList m_sf;
	RbtInt m_nNonHLigandAtoms;//for normalised scores (score / non-H ligand atoms)
};

//Useful typedefs
typedef SmartPtr<RbtSFAgg> RbtSFAggPtr;//Smart pointer
typedef vector<RbtSFAggPtr> RbtSFAggList;//Vector of smart pointers
typedef RbtSFAggList::iterator RbtSFAggListIter;
typedef RbtSFAggList::const_iterator RbtSFAggListConstIter;

#endif //_RBTSFAGG_H_
