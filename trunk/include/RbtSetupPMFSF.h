/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtSetupPMFSF.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

#ifndef _RBTSETUPPMFSF_H_
#define _RBTSETUPPMFSF_H_

#include "RbtBaseInterSF.h"

class RbtSetupPMFSF : public RbtBaseInterSF
{
		RbtAtomList			theLigandList;	// ligand typing
		RbtAtomList			theReceptorList;// receptor typing
		
	public:
		static RbtString _CT;

		RbtSetupPMFSF(const RbtString& strName="SETUP_PMF");
		~RbtSetupPMFSF();

	protected:
		virtual void 		SetupReceptor();
		virtual void 		SetupLigand();
		virtual void 		SetupScore();
		virtual RbtDouble	RawScore() const;

		void SetupReceptorPMFTypes(void)	throw (RbtError);
		void SetupLigandPMFTypes(void)		throw (RbtError);
		/**
		 * Routines to figure out complex PMF types
		 */
		RbtPMFType	GetPMFfor_lC(RbtAtomPtr); // for ligands
		RbtPMFType	GetPMFfor_lN(RbtAtomPtr);
		RbtPMFType	GetPMFfor_lO(RbtAtomPtr);
		RbtPMFType	GetPMFfor_lS(RbtAtomPtr);

		RbtPMFType	GetPMFfor_rC(RbtAtomPtr); // for receptor
		RbtPMFType	GetPMFfor_rN(RbtAtomPtr);
		RbtPMFType	GetPMFfor_rO(RbtAtomPtr);
		RbtPMFType	GetPMFfor_rS(RbtAtomPtr);

		RbtBool		IsChargedNitrogen(RbtAtomPtr);	// true for guanidino or other charged
		
};
#endif //_RBTSETUPPMFSF_H_

