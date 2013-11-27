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

//Indexed-grid-based PMF intermolecular scoring function

#ifndef _RBTPMFIDXSF_H_
#define _RBTPMFIDXSF_H_

#include "RbtBaseIdxSF.h"
#include "RbtBaseInterSF.h"
#include "RbtPMF.h"
#include "RbtAnnotationHandler.h"
// for PMF pseudogrids where 
// X is distance, X and Y are receptor and ligand
// distances respectively 
#include "RbtRealGrid.h"	

/** RbtPMFIdxSF class for PMF scoring.
 */
class RbtPMFIdxSF : public RbtBaseInterSF, public RbtBaseIdxSF, public virtual RbtAnnotationHandler 
{
	// by default "private" members
		/**
		 * vGrid: vector of grids. every PMF type needs a grid
		 */
		RbtAtomList						theLigandList;		// ligand typing
		RbtAtomList						theReceptorList;	// recepor typing
		// regular lists
		RbtAtomRList					theLigandRList;		// ligand
		RbtAtomRList					theReceptorRList;	// recepor
		vector<RbtRealGridPtr>			theTypeGrid;		// grids for PMF values for different atom types in receptor
		vector<RbtPMFType>				theLigandTypes;		// type values in theTypeGrid
		RbtNonBondedGridPtr				theSurround;		// atoms arond a gridpoint
		RbtRealGridPtr					thePMFGrid;			// grid for X-distance Y 
															// this is the representation of the PMFs
		RbtRealGridPtr					theSlopeGrid;		// grid to store values where the plateaus starts
		
	public:
		RbtPMFIdxSF(const RbtString& strName="PMF"); /**< The only one constructor */
		virtual ~RbtPMFIdxSF(); /**< The virtual destructor */
		
		/**
		 * This string serves as a name for the weight factor
		 * of the PMF score. Should be read from the .prm
		 * file.
		 */
		static RbtString _CT;
		static RbtString _PMFDIR; /**< String pointing to the dir with PMF tables*/
		static RbtString _CC_CUTOFF;
		static RbtString _SLOPE;
		/**
		 *  Observer functionalities pure virtual, so it must be
		 *  overloaded.
		 */
		virtual void 		Update(RbtSubject* theChangedSubject);


	protected:
		/**
		 * Setting up receptor means atom typing for all the
		 * receptor atoms and create RbtNonBondedGrid for NN
		 */
		virtual void SetupReceptor();
		/**
		 * determine PMF types for all ligand atoms
		 */
		virtual void SetupLigand();
		/**
		 * Setting up scores 
		 */
		virtual void SetupScore();
		/** 
		 * RawScore(void) returning with the actual PMF score
		 */
		virtual RbtDouble RawScore(void) const;
		/**
		 * Estimate value for short distances instead of using plateau in PMFs
		 */
		RbtDouble GetLinearCloseRangeValue(RbtDouble aDist, RbtPMFType aRecType, RbtPMFType aLigType) const;
};

#endif // _RBTPMFIDXSF_H_
