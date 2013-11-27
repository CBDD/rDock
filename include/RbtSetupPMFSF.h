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

