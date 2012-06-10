/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtConstSF.h#4 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
