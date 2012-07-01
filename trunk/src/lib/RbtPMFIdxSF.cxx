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

#include "RbtAtom.h"
#include "RbtPMFIdxSF.h"
#include "RbtPMFDirSource.h"
#include "RbtWorkSpace.h"

// the PMFs are defined only in this range in the Muegge database
const RbtDouble cPMFStart	= 0.2;	// closest PMF distance
const RbtDouble cPMFRes		= 0.2;	// PMF resolution
const RbtDouble cPMFEnd		= 12.0;	// farest point in PMFs
const RbtUInt	cPlStart	= 1;	// index value to get where the PMF plateau starts
const RbtUInt	cPlVal		= 2;	// index value to get the PMF plateau value 

RbtString RbtPMFIdxSF::_CT("RbtPMFIdxSF");
RbtString RbtPMFIdxSF::_PMFDIR("PMFDIR");
RbtString RbtPMFIdxSF::_CC_CUTOFF("CC_CUTOFF");
RbtString RbtPMFIdxSF::_SLOPE("SLOPE");

RbtDouble	delta;	// used for linear interpolation

RbtPMFIdxSF::RbtPMFIdxSF(const RbtString& aName)
: RbtBaseSF(_CT,aName)
{
	// see PMF-related .prm files for explanation
	AddParameter(_PMFDIR, "data/pmf");
	AddParameter(_CC_CUTOFF, 6.0);
	AddParameter(_SLOPE, -3.0);
	delta	 				= cPMFRes/2.0;	// half of the PMF grid resoluton: delta for linear interpolation
	// create the PMF pseudogrid
	RbtInt		nTypes		= 37;	// must be changed when we are defining new types :(
	RbtCoord	thePMFGridMin (cPMFStart, 0.0, 0.0);
	RbtCoord	thePMFGridStep(cPMFRes, 1.0, 1.0);
	thePMFGrid = RbtRealGridPtr(new RbtRealGrid(thePMFGridMin,thePMFGridStep,(RbtInt)(cPMFEnd/cPMFRes),nTypes,nTypes));
	// we are building a grid for slopes, so we can get the plateau
	// info similarly like the PMF values
	RbtCoord	theSlopeGridCoords(1.0, 1.0, 1.0);
	theSlopeGrid = RbtRealGridPtr(new RbtRealGrid(theSlopeGridCoords,theSlopeGridCoords,2,// cPlStart & cPlVal	
				nTypes,nTypes));

	// setup the grids for every types
	// first read .pmf files from directory
	RbtPMFDirSource	theSrcDir( Rbt::GetRbtRoot()+"/"+(RbtString)GetParameter(_PMFDIR));
	vector< vector<RbtPMFValue> >		thePMF;			// vector helping to read PMF files
	vector<string>						theFileNames;	// vector storing filenames (needed to find out types in thePMF)
	vector<RbtPMFValue>					theSlopeIndex;	// contains the first value where plateau starts
	theSrcDir.ReadFiles(&thePMF,&theFileNames,&theSlopeIndex);
	// convert vector to grid
	for(int i=0; i<thePMF.size(); i++) {
#ifdef _DEBUG1
		cout << _CT << " " << theFileNames[i]; 
		cout << " " << theFileNames[i].substr(0,2) << " " << Rbt::PMFStr2Type( theFileNames[i].substr(0,2) );
		cout << " " << theFileNames[i].substr(2,2) << " " << Rbt::PMFStr2Type( theFileNames[i].substr(2,2) );
		cout << " " << endl;
#endif	// _DEBUG1

		// iterate through each individual PMF
		for(int j=0; j<thePMF[i].size(); j++) {
			RbtPMFType rType = Rbt::PMFStr2Type( theFileNames[i].substr(0,2));		// receptor type
			RbtPMFType lType = Rbt::PMFStr2Type( theFileNames[i].substr(2,2));		// ligand type
			thePMFGrid->SetValue( thePMFGrid->GetIX((RbtDouble)thePMF[i][j].distance),	// distance
					rType,lType,
					thePMF[i][j].density						// the density value
					);
		}
	}
	// fill slope grid
	for(int i=0; i<theSlopeIndex.size(); i++) {
		RbtPMFType rType = Rbt::PMFStr2Type( theFileNames[i].substr(0,2));		// receptor type
		RbtPMFType lType = Rbt::PMFStr2Type( theFileNames[i].substr(2,2));		// ligand type
		theSlopeGrid->SetValue( cPlStart,		// index 1 for the distance where the		
				rType,lType,	// plateau starts
				theSlopeIndex[i].distance
				);
		theSlopeGrid->SetValue( cPlVal,			// that is the actual value  
				rType,lType,	// can be more than 3.0
				theSlopeIndex[i].density
				);
	}
#ifdef _DEBUG1
	// checking values for type cP only. This test assumes the original
	// Muegge file list where the very first is the cPcP.pmf 
	for(int j=0; j<thePMF[0].size(); j++) {	
		cout << _CT << " " << thePMFGrid->GetValue( thePMFGrid->GetIX(thePMF[0][j].distance),cP,cP) - thePMF[0][j].density << endl; 
	}
	cout << _CT << " " << theFileNames[0] << endl;

	//	 this is a printout of values NOT the grid
	cout << _CT << " number of PMF functions: " << thePMF.size() << endl;
	for(long k=0; k<thePMF.size(); k++) {
		for(int l=0, m=0; l<thePMF[k].size(); l++, m++) {
			cout << "| "<<thePMF[k][l].distance << "  " << thePMF[k][l].density << "| ";
			if(5==m) {
				m=0;
				cout << endl;
			}
		}
		cout << endl;
		cout << endl;
	}
#endif //_DEBUG1
#ifdef _DEBUG1
	// other test for grid. RbtPMFType enum starts with 1.
	for (int i=1; i<PMF_UNDEFINED; i++) {
		for(int j=1; j<PMF_UNDEFINED; j++) {
			cout << "PMF Pair: "<<Rbt::PMFType2Str((RbtPMFType)i)<<Rbt::PMFType2Str((RbtPMFType)j)<<" *********"<<endl; 
			for(RbtDouble grid = cPMFStart; grid<cPMFEnd; grid+=cPMFRes) {
				RbtUInt     idx     = thePMFGrid->GetIX(grid);
				cout <<Rbt::PMFType2Str((RbtPMFType)i)<<Rbt::PMFType2Str((RbtPMFType)j);
				cout <<" "<<grid <<" "<<idx<<" "<<thePMFGrid->GetValue(idx,(RbtPMFType)i,(RbtPMFType)j)<<endl;

			}
		}
	}
#endif//_DEBUG1

	_RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtPMFIdxSF::~RbtPMFIdxSF() 
{
#ifdef _DEBUG
	cout << _CT << " PMF destructor" << endl;
#endif //_DEBUG

	_RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtPMFIdxSF::Update(RbtSubject* theChangedSubject)
{
#ifdef _DEBUG
	cout << _CT << " PMF Update" << endl;
#endif //_DEBUG
	RbtBaseInterSF::Update(theChangedSubject);
}

/**
 * SetupReceptor: 
 * determine PMF atom types for all receptor atoms
 * create RbtNonBondedGrid of nearest neighbours
 */
void RbtPMFIdxSF::SetupReceptor()
{
#ifdef _DEBUG
	cout << _CT << " PMF SetupReceptor" << endl;
#endif //_DEBUG
	// clear to be on the safe side
	theReceptorList.clear();
	theReceptorRList.clear();
	theSurround	= RbtNonBondedGridPtr();
	if (GetReceptor().Null()) {
		cout << _CT << "WARNING: no receptor defined. " << endl;
		return;
	}
	else {	// load PMFs from table files
		theSurround				= CreateNonBondedGrid();
		RbtDockingSitePtr spDS	= GetWorkSpace()->GetDockingSite();
		RbtDouble range			= GetRange()+GetMaxError();
		RbtAtomList atomList	= spDS->GetAtomList(GetReceptor()->GetAtomList(),0.0,GetCorrectedRange());
		theReceptorList			= Rbt::GetAtomList(atomList,std::not1(Rbt::isAtomicNo_eq(1)));
		//cout << "Receptor list size: "<<theReceptorList.size()<<endl;
		// create non-bonded grid to get the atoms around each gridpoint
		for(RbtAtomListIter sIter	= theReceptorList.begin(); sIter != theReceptorList.end(); ++sIter) {
			theSurround->SetAtomLists(*sIter,range); // get surround at cutoff+gridstep*sqrt(3)/2
			// initialise cumulative PMF values for annotation
			(*sIter)->SetUser2Value(0.0);
		}
	}
	// transform smartpointers into regular ones
	std::copy( theReceptorList.begin(), theReceptorList.end(), std::back_inserter(theReceptorRList));
}

void RbtPMFIdxSF::SetupLigand()
{
#ifdef _DEBUG
	cout << _CT << " PMF SetupLigand" << endl;
#endif //_DEBUG

	theLigandList	= Rbt::GetAtomList(GetLigand()->GetAtomList(),std::not1(Rbt::isAtomicNo_eq(1)));

	// clean old list
	theLigandRList.clear();
	// transform smartpointers into regular ones
	std::copy( theLigandList.begin(), theLigandList.end(), std::back_inserter(theLigandRList));
}

void RbtPMFIdxSF::SetupScore()
{
#ifdef _DEBUG
	cout << _CT << " PMF SetupScore" << endl;
#endif //_DEBUG
}


RbtDouble RbtPMFIdxSF::RawScore() const 
{
	//	return 0.0;
#ifdef _DEBUG
	cout << _CT << " PMF RawScore " << endl;
#endif //_DEBUG
	RbtDouble	theScore = 0.0;

	// check for existence of  atom list grid
	if (theSurround.Null()) {
		cout << _CT << "No index grid" << endl;
		return theScore;
	}
	// enable/disable annotations
	RbtBool bAnnotate = isAnnotationEnabled();

	// for all ligand atoms:
	for(RbtAtomRListConstIter lIter=theLigandRList.begin(); lIter!=theLigandRList.end(); ++lIter ) {
		const RbtCoord& ligCoord		= (*lIter)->GetCoords();
		// get receptor atoms that are within the PMF radius - if there are any
		const RbtAtomRList& rAtomList	= theSurround->GetAtomList(ligCoord);
		if(rAtomList.empty())
			continue;
		//cout << _CT << " loop over the receptor atoms around the ligand atom" << endl;
		for (RbtAtomRListConstIter rIter=rAtomList.begin(); rIter!=rAtomList.end(); rIter++) {
			// get distance of the atom
			RbtDouble theDist	= Rbt::Length((*rIter)->GetCoords(),ligCoord);
			RbtDouble i_score;	// interpolated score

			if(theDist > GetRange() ) 	// skip distances out of a given distance
				continue;
			//cout << _CT << " theDist " << theDist << endl;
			// to make code readable introduce some local variables
			// hoping cc can optimize them out
			const RbtPMFType	rType	= (*rIter)->GetPMFType();
			//cout << _CT << " receptor type " << Rbt::PMFType2Str(rType) << endl;	
			const RbtPMFType	lType	= (*lIter)->GetPMFType();
			// optimal distance for C-C interactions is
			// under 6A. Note NC is the next item in RbtPMFType
			// after the carbon types
			if(theDist>(RbtDouble)GetParameter(_CC_CUTOFF) && rType<NC && lType<NC)
				continue;
			// if we are in the plateau region
			else if(theDist < theSlopeGrid->GetValue(1,rType, lType)) {
				i_score = GetLinearCloseRangeValue(theDist,rType,lType);
#ifdef _DEBUG1
				cout << "Plateau reached: "<<theDist<<" "<<theSlopeGrid->GetValue(1,rType, lType);
				cout <<" "<<Rbt::PMFType2Str(rType)<<" "<<Rbt::PMFType2Str(lType)<< endl; 
#endif // _DEBUG1
			} else {	
				//if(HH==rType || HL==lType || HL==rType || HH==lType)
				//	cout << Rbt::PMFType2Str(rType)<<Rbt::PMFType2Str(lType)<<" should not happen" << endl;
				//cout << _CT << " ligand type " << Rbt::PMFType2Str(lType) << endl;	
				//RbtUInt		idx		= thePMFGrid->GetIX(theDist);
				//RbtDouble score		= thePMFGrid->GetValue( idx, rType, lType );
				//theScore += score;
				// make a linear interpolation
				RbtUInt   inf_idx	= thePMFGrid->GetIX(theDist-delta);	
				RbtDouble inf_score = thePMFGrid->GetValue( inf_idx, rType, lType );
				RbtUInt   sup_idx	= thePMFGrid->GetIX(theDist+delta);	
				RbtDouble sup_score = thePMFGrid->GetValue( sup_idx, rType, lType );
				// now calculate the distances from the gridpoints
				RbtDouble inf_d     = (theDist-thePMFGrid->GetXCoord(inf_idx))/cPMFRes;
				RbtDouble sup_d     = 1.0-inf_d;
				// weight the score with the distances from gridpoints
				i_score	= inf_score*sup_d+sup_score*inf_d;
			}
			// store (increment) contribution of receptor atom 
			(*rIter)->SetUser2Value( (*rIter)->GetUser2Value() + i_score );
			theScore += i_score;
		}
	}	
#if 0
	// check PMF scores: the two values should be the same
	RbtDouble receptorPMFscore = 0.0;
	for(RbtAtomRListConstIter sIter	= theReceptorRList.begin(); sIter != theReceptorRList.end(); ++sIter) {
		receptorPMFscore += (*sIter)->GetUser2Value();
	}
	cout << _CT << " Receptor score is: " << receptorPMFscore << endl;
	cout << _CT << " PMF score is     : " << theScore << endl;
#endif
	
	// save annotation (if needed)
	// since there is no ligand side here we are assigning the the very first ligand atom
	// also no distance 
	if(bAnnotate) {
		for(RbtAtomRListConstIter sIter	= theReceptorRList.begin(); sIter != theReceptorRList.end(); ++sIter) {
			RbtAtom* lAtom	= (*theLigandRList.begin());
			RbtAnnotationPtr spAnnotation( 
					new RbtAnnotation (
						lAtom,						// ligand
						(*sIter),					// receptor
						0.0,						// distance
						(*sIter)->GetUser2Value() )	// cumulative receptor PMF contribution
						);
			AddAnnotation(spAnnotation);
			(*sIter)->SetUser2Value(0.0);			// clean old values
		}
	}
	
	return theScore;
}

RbtDouble RbtPMFIdxSF::GetLinearCloseRangeValue(RbtDouble aDist, RbtPMFType aRecType,RbtPMFType aLigType)
	const
{
	RbtDouble thePlateauStart 	= theSlopeGrid->GetValue(cPlStart,aRecType, aLigType);
	RbtDouble thePlateauVal		= theSlopeGrid->GetValue(cPlVal,aRecType, aLigType);
	RbtDouble theSlope			= (RbtDouble)GetParameter(_SLOPE);

	return theSlope*aDist-theSlope*thePlateauStart+thePlateauVal;
}
