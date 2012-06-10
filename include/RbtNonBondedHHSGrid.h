/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtNonBondedHHSGrid.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

#ifndef _RBTNONBONDEDHHSGRID_H_
#define _RBTNONBONDEDHHSGRID_H_

#include "RbtBaseGrid.h"
#include "RbtSATypes.h"

// typedefs are in RbtSATypes.h 
// rest of the class architecture from RbtNonBondedGrid.h
// but instead RbtAtom HHS_Solvation is used as operand
 
class RbtNonBondedHHSGrid : public RbtBaseGrid
{
	public:
		static RbtString _CT;
		RbtNonBondedHHSGrid(const RbtCoord& gridMin, const RbtCoord& gridStep,
				         RbtUInt NX, RbtUInt NY, RbtUInt NZ, RbtUInt NPad=0);
		RbtNonBondedHHSGrid(istream& istr);
		~RbtNonBondedHHSGrid();

		RbtNonBondedHHSGrid(const RbtNonBondedHHSGrid&);
		RbtNonBondedHHSGrid(const RbtBaseGrid&);
		RbtNonBondedHHSGrid& operator=(const RbtNonBondedHHSGrid&);
		RbtNonBondedHHSGrid& operator=(const RbtBaseGrid&);

		virtual	void Print(ostream& ostr) const;
		virtual	void Write(ostream& ostr) const;
		virtual	void Read(istream& istr);

		const	HHS_SolvationRList&	GetHHSList(RbtUInt iXYZ) const; 
		const	HHS_SolvationRList&	GetHHSList(const RbtCoord& c) const; 

				void 				SetHHSLists(HHS_Solvation* pHHS, RbtDouble radius);
				void				ClearHHSLists(void);
				
	protected:
		void	OwnPrint(ostream& ostr) const;
		void	OwnWrite(ostream& ostr) const;
		void	OwnRead(istream& istr) throw (RbtError);
		
	private:
		RbtNonBondedHHSGrid();

		void	CopyGrid(const RbtNonBondedHHSGrid&);
		void	CreateMap();

		HHS_SolvationListMap		m_hhsMap;
		const HHS_SolvationRList	m_emptyList;
};

typedef	SmartPtr<RbtNonBondedHHSGrid> RbtNonBondedHHSGridPtr;
#endif // _RBTNONBONDEDHHSGRID_H_
