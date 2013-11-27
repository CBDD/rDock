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

//Base class for indexed-grid-based scoring functions
//Provides protected methods for subclasses to use to create indexing grids

#ifndef _RBTBASEIDXSF_H_
#define _RBTBASEIDXSF_H_

#include "RbtBaseSF.h"
#include "RbtInteractionGrid.h"
#include "RbtNonBondedGrid.h"
#include "RbtNonBondedHHSGrid.h"

class RbtBaseIdxSF : public virtual RbtBaseSF
{
	public:
		//Class type string
		static RbtString _CT;
		//Parameter names
		static RbtString _GRIDSTEP;
		static RbtString _BORDER;

		////////////////////////////////////////
		//Constructors/destructors
		virtual ~RbtBaseIdxSF();

		////////////////////////////////////////
		//Public methods
		////////////////
		RbtDouble GetGridStep() const;
		void SetGridStep(RbtDouble step);
		RbtDouble GetBorder() const;
		void SetBorder(RbtDouble border);

	protected:
		////////////////////////////////////////
		//Protected methods
		///////////////////
		RbtBaseIdxSF();

		RbtInteractionGridPtr	CreateInteractionGrid() const;
		RbtNonBondedGridPtr		CreateNonBondedGrid() const;
		RbtNonBondedHHSGridPtr	CreateNonBondedHHSGrid() const;
		RbtDouble 				GetMaxError() const;
		//DM 12 Apr 2002
		//Returns the maximum range of the scoring function,
		//corrected for max grid error, and grid border around docking site
		//This should be used by subclasses for selecting the receptor atoms to index
		//GetCorrectedRange() = GetRange() + GetMaxError() + GetBorder()
		RbtDouble GetCorrectedRange() const;

		//As this has a virtual base class we need a separate OwnParameterUpdated
		//which can be called by concrete subclass ParameterUpdated methods
		//See Stroustrup C++ 3rd edition, p395, on programming virtual base classes
		void OwnParameterUpdated(const RbtString& strName);

	private:
		////////////////////////////////////////
		//Private methods
		/////////////////

	protected:
		////////////////////////////////////////
		//Protected data
		////////////////

	private:
		////////////////////////////////////////
		//Private data
		//////////////
		RbtDouble m_gridStep;
		RbtDouble m_border;
};

#endif //_RBTBASEIDXSF_H_
