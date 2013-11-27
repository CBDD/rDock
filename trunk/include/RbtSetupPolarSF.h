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

//Pseudo-scoring function who's only role is to setup the receptor local
//neighbour density property and store in each atom's User1Value
//This SF should be registered with the workspace BEFORE any SF which needs
//this property (e.g. HBOND, IONIC) to ensure that User1Value is calculated
//before being used. Scoring function is disabled by default.

#ifndef _RBTSETUPPOLARSF_H_
#define _RBTSETUPPOLARSF_H_

#include "RbtBaseInterSF.h"

class RbtSetupPolarSF : public RbtBaseInterSF
{
	public:
	//Class type string
	static RbtString _CT;
	//Parameter names
	static RbtString _RADIUS;
	static RbtString _NORM;
	static RbtString _POWER;
	static RbtString _CHGFACTOR;
	//DM 14 Nov 2001 - relative strength of guanidinium intns
	static RbtString _GUANFACTOR;
  
  RbtSetupPolarSF(const RbtString& strName = "SETUP_POLAR");
	virtual ~RbtSetupPolarSF();
	
	protected:
	virtual void SetupReceptor();
	virtual void SetupLigand();
	virtual void SetupSolvent();
	virtual void SetupScore();
  	virtual RbtDouble RawScore() const;

	private:
	void SetupAtomList(RbtAtomList& atomList,
						const RbtAtomList& neighbourList,
						RbtInt traceTriggerLevel);
};

#endif //_RBTSETUPPOLARSF_H_
