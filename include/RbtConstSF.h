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

//Scoring function to penalise ligand and solvent binding
//Ligand binding RawScore is always 1
//Solvent binding RawScore is:
//SOLVENT_PENALTY * the number of enabled solvent models
#ifndef _RBTCONSTSF_H_
#define _RBTCONSTSF_H_

#include "RbtBaseInterSF.h"

class RbtConstSF : public RbtBaseInterSF
{
	public:
	//Static data member for class type (i.e. "RbtConstSF")
	static RbtString _CT;
	static RbtString _SOLVENT_PENALTY;

 	RbtConstSF(const RbtString& strName = "CONST");
	virtual ~RbtConstSF();

	virtual void ScoreMap(RbtStringVariantMap& scoreMap) const;

	protected:
  	virtual void SetupReceptor() {};
  	virtual void SetupLigand() {};
  	virtual void SetupScore() {};
  	virtual RbtDouble RawScore() const;
  	void ParameterUpdated(const RbtString& strName);

	private:
	//The original constant score for ligand binding
	RbtDouble InterScore() const;
	//The solvent binding penalty
	RbtDouble SystemScore() const;
	RbtDouble m_solventPenalty;
};

#endif //_RBTCONSTSF_H_
