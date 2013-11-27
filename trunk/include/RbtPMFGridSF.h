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

#ifndef _RBTPMFGRIDSF_H_
#define _RBTPMFGRIDSF_H_

#include "RbtBaseInterSF.h"
#include "RbtRealGrid.h"

class RbtPMFGridSF : public RbtBaseInterSF 
{
		RbtBool					m_bSmoothed;
		RbtAtomList				theLigandList;			// vector to store the ligand
		vector<RbtPMFType>		theTypeList;			// store PMF used types here
		vector<RbtRealGridPtr>	theGrids;				// grids with PMF data
	
		void				ReadGrids(istream& istr) throw (RbtError);	

	public:
		static				RbtString	_CT;		// class name
		static				RbtString	_GRID;		// filename extension (.grd)
		static				RbtString	_SMOOTHED;	// controls wether to smooth the grid values
		
							RbtPMFGridSF(const RbtString& strName = "PMFGRID");
		virtual				~RbtPMFGridSF();
	protected:
		virtual void		SetupReceptor();
		virtual void		SetupLigand();
		virtual void		SetupScore() {};
		virtual RbtDouble	RawScore() const;
		RbtUInt				GetCorrectedType(RbtPMFType aType) const;
};

#endif	// _RBTPMFGRIDSF_H_

