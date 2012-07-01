/*This file is part of Rdock.

    Rdock is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Rdock is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Rdock.  If not, see <http://www.gnu.org/licenses/>.*/

#include "RbtBaseIdxSF.h"
#include "RbtWorkSpace.h"

//Static data members
RbtString RbtBaseIdxSF::_CT("RbtBaseIdxSF");
RbtString RbtBaseIdxSF::_GRIDSTEP("GRIDSTEP");
RbtString RbtBaseIdxSF::_BORDER("BORDER");

RbtBaseIdxSF::RbtBaseIdxSF() : m_gridStep(0.5), m_border(1.0)
{
#ifdef _DEBUG
	cout << _CT << " default constructor" << endl;
#endif //_DEBUG
	//Add parameters
	AddParameter(_GRIDSTEP,m_gridStep);
	AddParameter(_BORDER,m_border);
	_RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtBaseIdxSF::~RbtBaseIdxSF() {
	_RBTOBJECTCOUNTER_DESTR_(_CT);
}

RbtDouble RbtBaseIdxSF::GetGridStep() const {
	return m_gridStep;
}

void RbtBaseIdxSF::SetGridStep(RbtDouble step) {
	SetParameter(_GRIDSTEP,step);
}

RbtDouble RbtBaseIdxSF::GetBorder() const {
	return m_border;
}

void RbtBaseIdxSF::SetBorder(RbtDouble border) {
	SetParameter(_BORDER,border);
}

//DM 10 Apr 2002
//I know, I know, grids should be templated to avoid the need for two different CreateGrid methods...
RbtInteractionGridPtr RbtBaseIdxSF::CreateInteractionGrid() const {
	//Create a grid covering the docking site
	RbtDockingSitePtr spDS = GetWorkSpace()->GetDockingSite();
	if (spDS.Null())
		return RbtInteractionGridPtr();

	//Extend grid by _BORDER, mainly to allow for the possibility of polar hydrogens
	//falling outside of the docking site
	RbtCoord minCoord = spDS->GetMinCoord()-m_border;
	RbtCoord maxCoord = spDS->GetMaxCoord()+m_border;
	RbtVector recepExtent = maxCoord-minCoord;
	RbtVector gridStep(m_gridStep,m_gridStep,m_gridStep);
	RbtUInt nX = int(recepExtent.x/gridStep.x)+1;
	RbtUInt nY = int(recepExtent.y/gridStep.y)+1;
	RbtUInt nZ = int(recepExtent.z/gridStep.z)+1;
	return RbtInteractionGridPtr(new RbtInteractionGrid(minCoord,gridStep,nX,nY,nZ));
}

RbtNonBondedGridPtr RbtBaseIdxSF::CreateNonBondedGrid() const {
	//Create a grid covering the docking site
	RbtDockingSitePtr spDS = GetWorkSpace()->GetDockingSite();
	if (spDS.Null())
		return RbtInteractionGridPtr();

	//Extend grid by _BORDER, mainly to allow for the possibility of polar hydrogens
	//falling outside of the docking site
	RbtCoord minCoord = spDS->GetMinCoord()-m_border;
	RbtCoord maxCoord = spDS->GetMaxCoord()+m_border;
	RbtVector recepExtent = maxCoord-minCoord;
	RbtVector gridStep(m_gridStep,m_gridStep,m_gridStep);
	RbtUInt nX = int(recepExtent.x/gridStep.x)+1;
	RbtUInt nY = int(recepExtent.y/gridStep.y)+1;
	RbtUInt nZ = int(recepExtent.z/gridStep.z)+1;
	return RbtNonBondedGridPtr(new RbtNonBondedGrid(minCoord,gridStep,nX,nY,nZ));
}

RbtNonBondedHHSGridPtr RbtBaseIdxSF::CreateNonBondedHHSGrid() const {
	//Create a grid covering the docking site
	RbtDockingSitePtr spDS = GetWorkSpace()->GetDockingSite();
	if (spDS.Null())
		return RbtNonBondedHHSGridPtr();

	//Extend grid by _BORDER, mainly to allow for the possibility of polar hydrogens
	//falling outside of the docking site
	RbtCoord minCoord = spDS->GetMinCoord()-m_border;
	RbtCoord maxCoord = spDS->GetMaxCoord()+m_border;
	RbtVector recepExtent = maxCoord-minCoord;
	RbtVector gridStep(m_gridStep,m_gridStep,m_gridStep);
	RbtUInt nX = int(recepExtent.x/gridStep.x)+1;
	RbtUInt nY = int(recepExtent.y/gridStep.y)+1;
	RbtUInt nZ = int(recepExtent.z/gridStep.z)+1;
	return RbtNonBondedHHSGridPtr(new RbtNonBondedHHSGrid(minCoord,gridStep,nX,nY,nZ));
}

RbtDouble RbtBaseIdxSF::GetMaxError() const {
	//maxError is half a grid diagonal. This is the tolerance we have to allow when indexing the receptor atoms on the grid
	//When retrieving the nearest neighbours to a ligand atom, we do the lookup on the nearest grid point to the ligand atom.
	//However the actual ligand atom - receptor atom distance may be maxError Angstroms closer than the
	//grid point - receptor atom distance used in the indexing.
	return 0.5*sqrt(3.0)*m_gridStep;
}

//DM 12 Apr 2002
//Returns the maximum range of the scoring function,
//corrected for max grid error, and grid border around docking site
//This should be used by subclasses for selecting the receptor atoms to index
//GetCorrectedRange() = GetRange() + GetMaxError() + GetBorder()
RbtDouble RbtBaseIdxSF::GetCorrectedRange() const {
	return GetRange() + GetMaxError() + m_border;
}

//As this has a virtual base class we need a separate OwnParameterUpdated
//which can be called by concrete subclass ParameterUpdated methods
//See Stroustrup C++ 3rd edition, p395, on programming virtual base classes
void RbtBaseIdxSF::OwnParameterUpdated(const RbtString& strName) {
	if (strName == _GRIDSTEP) {
		m_gridStep = GetParameter(_GRIDSTEP);
	} else if (strName == _BORDER) {
		m_border = GetParameter(_BORDER);
	}
}

