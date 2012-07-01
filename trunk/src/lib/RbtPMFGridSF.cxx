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

#include "RbtPMFGridSF.h"
#include "RbtFileError.h"
#include "RbtWorkSpace.h"

RbtString RbtPMFGridSF::_CT("RbtPMFGridSF");
RbtString RbtPMFGridSF::_GRID("GRID");
RbtString RbtPMFGridSF::_SMOOTHED("SMOOTHED");

RbtPMFGridSF::RbtPMFGridSF(const RbtString& strName)
  : RbtBaseSF(_CT,strName),m_bSmoothed(true) 
{
  AddParameter(_GRID,".grd");
  AddParameter(_SMOOTHED,m_bSmoothed);
  
  cout << _CT << " parameterised constructor" << endl;
  
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtPMFGridSF::~RbtPMFGridSF()
{
  cout << _CT << " destructor" << endl;
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtPMFGridSF::SetupReceptor()
{
  cout <<_CT<<"::SetupReceptor()"<<endl;
  theGrids.clear();
  
  if (GetReceptor().Null())
    return;
  
  RbtString strWSName = GetWorkSpace()->GetName();
  
  RbtString strSuffix = GetParameter(_GRID);
  RbtString strFile = Rbt::GetRbtFileName("data/grids",strWSName+strSuffix);
  cout << _CT << " Reading PMF grid from " <<strFile<<endl;
  ifstream istr(strFile.c_str(),ios_base::in|ios_base::binary);
  ReadGrids(istr);
  istr.close();
}
//Determine PMF grid type for each atom
void RbtPMFGridSF::SetupLigand()
{
	theLigandList.clear();
	if (GetLigand().Null())
		return;

	// get the  non-H atoms only
	theLigandList = Rbt::GetAtomList(GetLigand()->GetAtomList(),std::not1(Rbt::isAtomicNo_eq(1)));
#ifdef _DEBUG
	cout << _CT << "::SetupLigand(): #ATOMS = " << theLigandList.size() << endl;
#endif //_DEBUG
}

RbtDouble RbtPMFGridSF::RawScore() const
{
	RbtDouble theScore = 0.0;
	if(theGrids.empty())	// if no grids defined
		return theScore;

	//Loop over all ligand atoms
	RbtAtomListConstIter iter = theLigandList.begin();
	if (m_bSmoothed) { 
		for (RbtUInt i = 0; iter != theLigandList.end(); iter++,i++) {
			RbtUInt		theType	= GetCorrectedType((*iter)->GetPMFType());
			RbtDouble	score 	= theGrids[theType-1]->GetSmoothedValue((*iter)->GetCoords());
			theScore			+= score;
		}
	} else { 
		for (RbtUInt i = 0; iter != theLigandList.end(); iter++,i++) {
			RbtUInt		theType	= GetCorrectedType((*iter)->GetPMFType());
			RbtDouble	score	= theGrids[theType-1]->GetValue((*iter)->GetCoords());
			theScore			+= score;
		}
	}
	return theScore;

}

void RbtPMFGridSF::ReadGrids(istream& istr) throw (RbtError) 
{
	cout << "**************************************************************"<<endl;
	cout << "                    Reading grid..."<<endl;
	cout << "**************************************************************"<<endl;
	cout << _CT <<"::ReadGrids(istream&)"<<endl;
	theGrids.clear();
  
	//Read header string
	RbtInt	length;
	Rbt::ReadWithThrow(istr, (char*) &length, sizeof(length));
	char*	header	= new char [length+1];
	Rbt::ReadWithThrow(istr, header, length);
	//Add null character to end of string
	header[length]	= '\0';
	//Compare title with 
	RbtBool	match	= (_CT == header);
	delete [] header;
	if (!match) {
		throw RbtFileParseError(_WHERE_,"Invalid title string in " + _CT + "::ReadGrids()");
	}
	
	//Now read number of grids
	RbtInt	nGrids;
	Rbt::ReadWithThrow(istr, (char*) &nGrids, sizeof(nGrids));
	cout << "Reading " << nGrids << " grids..." << endl;
	theGrids.reserve(nGrids);
	for (RbtInt i = CF; i <= nGrids; i++) {
		cout << "Grid# " << i << " ";
		// Read type
		RbtPMFType theType;
		Rbt::ReadWithThrow(istr, (char*) &theType, sizeof(theType));
		cout << "type "<< Rbt::PMFType2Str(theType);
		//Now we can read the grid
		RbtRealGridPtr spGrid(new RbtRealGrid(istr));
		cout << " done"<<endl;
		theGrids.push_back(spGrid);
	}
}

// since there is no  HH,HL,Fe,V,Mn grid, we have to correct
RbtUInt RbtPMFGridSF::GetCorrectedType(RbtPMFType aType) const
{
	if(aType<HL)
		return (RbtUInt) aType;
	else if (aType<Mn)
		return (RbtUInt) aType - 1;	// HL and HH but list starts at 0 still PMFType starts at 1
	else if (aType<Fe)
		return (RbtUInt) aType - 2; // HL,HH,Mn (3-1)
	else 
		return (RbtUInt) aType - 3; // HL,HH,Mn,Fe - V is the very last in the list
}
