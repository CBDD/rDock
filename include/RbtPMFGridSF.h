/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtPMFGridSF.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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

