/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtSetupSASF.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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

