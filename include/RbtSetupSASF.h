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

#ifndef _RBTSETUPSASF_H_
#define _RBTSETUPSASF_H_

#include "RbtBaseInterSF.h"

// after Hasel, Hendrickson and Still
struct RbtSimpleSolvationParameters {
	RbtDouble	p;
	RbtDouble	r;
};

class RbtSetupSASF : public RbtBaseInterSF
{
		RbtAtomList			theLigandList;	// ligand typing
		RbtAtomList			theReceptorList;// receptor typing
		
	public:
		static RbtString _CT;

		RbtSetupSASF(const RbtString& strName="SETUP_SA");
		~RbtSetupSASF();

	protected:
		virtual void 		SetupReceptor();
		virtual void 		SetupLigand();
		virtual void 		SetupScore();
		virtual RbtDouble	RawScore() const;

		void SetupReceptorSATypes(void)	throw (RbtError);
		void SetupLigandSATypes(void)		throw (RbtError);
		
};
#endif //_RBTSETUPSASF_H_

