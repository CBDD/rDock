/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtVdwGridSF.cxx#4 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtVdwGridSF.h"
#include "RbtFileError.h"
#include "RbtWorkSpace.h"

//Static data members
RbtString RbtVdwGridSF::_CT("RbtVdwGridSF");
RbtString RbtVdwGridSF::_GRID("GRID");
RbtString RbtVdwGridSF::_SMOOTHED("SMOOTHED");
	
//NB - Virtual base class constructor (RbtBaseSF) gets called first,
//implicit constructor for RbtBaseInterSF is called second
RbtVdwGridSF::RbtVdwGridSF(const RbtString& strName)
  : RbtBaseSF(_CT,strName),m_bSmoothed(true) {
  //Add parameters
  AddParameter(_GRID,".grd");
  AddParameter(_SMOOTHED,m_bSmoothed);
#ifdef _DEBUG
  cout << _CT << " parameterised constructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtVdwGridSF::~RbtVdwGridSF() {
#ifdef _DEBUG
  cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}


void RbtVdwGridSF::SetupReceptor() {
  m_grids.clear();
  if (GetReceptor().Null())
    return;

  //Trap multiple receptor conformations and flexible OH/NH3 here: this SF does not support them yet
  RbtBool bEnsemble = (GetReceptor()->GetNumSavedCoords() > 1);
  RbtBool bFlexRec = GetReceptor()->isFlexible();
  if (bEnsemble || bFlexRec) {
    RbtString message("Vdw grid scoring function does not support multiple receptor conformations\n");
    message += "or flexible OH/NH3 groups yet";
    throw RbtInvalidRequest(_WHERE_,message);
  }
    
  //Read grids
  //File names are composed of workspace name + grid suffix
  //Reasonably safe to assume that GetWorkSpace is not NULL,
  //as otherwise SetupReceptor would not have been called
  RbtString strWSName = GetWorkSpace()->GetName();
  
  RbtString strSuffix = GetParameter(_GRID);
  RbtString strFile = Rbt::GetRbtFileName("data/grids",strWSName+strSuffix);
  //DM 26 Sep 2000 - ios_base::binary qualifier doesn't appear to be valid
  //with IRIX CC
#ifdef __sgi
  ifstream istr(strFile.c_str(),ios_base::in);
#else
  ifstream istr(strFile.c_str(),ios_base::in|ios_base::binary);
#endif
  ReadGrids(istr);
  istr.close();
}

void RbtVdwGridSF::SetupLigand() {
  m_ligAtomList.clear();
  m_ligAtomTypes.clear();
  if (GetLigand().Null())
    return;

  RbtAtomList tmpList = GetLigand()->GetAtomList();
  //Strip off the smart pointers
  std::copy(tmpList.begin(),tmpList.end(),std::back_inserter(m_ligAtomList));
}

void RbtVdwGridSF::SetupSolvent() {
  RbtModelList solvent = GetSolvent();
  if (!solvent.empty()) {
    RbtString message("Vdw grid scoring function does not support explicit solvent yet\n");
    throw RbtInvalidRequest(_WHERE_,message);
  }
}

void RbtVdwGridSF::SetupScore() {
  //Determine probe grid type for each atom, based on comparing Tripos atom type with probe atom types
  //This needs to be in SetupScore as it is dependent on both the ligand and receptor grid data
  m_ligAtomTypes.clear();
  if (m_ligAtomList.empty())
    return;

  RbtInt iTrace = GetTrace();
  m_ligAtomTypes.reserve(m_ligAtomList.size());
  //Check if we have a grid for the UNDEFINED type:
  //If so, we can use it if a particular atom type grid is missing
  //If not, then we have to throw an error if a particular atom type grid is missing
  RbtBool bHasUndefined = !m_grids[RbtTriposAtomType::UNDEFINED].Null();
  RbtTriposAtomType triposType;

  if (iTrace > 1) {
    if (bHasUndefined) {
      cout << "The grid for the UNDEFINED atom type will be used if a grid is missing for any atom type" << endl;
    }
    else {
      cout << "There is no grid for the UNDEFINED atom type. An error will be thrown if a grid is missing for any atom type" << endl;
    }
  }

  for (RbtAtomRListConstIter iter = m_ligAtomList.begin(); iter != m_ligAtomList.end(); iter++) {
    RbtTriposAtomType::eType aType = (*iter)->GetTriposType();

    //If there is no grid for this atom type, revert to using the UNDEFINED grid if available
    //else throw an error
    if (m_grids[aType].Null()) {
      RbtString strError = "No vdw grid available for " + (*iter)->GetFullAtomName() + " (type " + triposType.Type2Str(aType) + ")";
      if (iTrace > 1) {
	cout << strError << endl;
      }
      if (bHasUndefined) {
	aType = RbtTriposAtomType::UNDEFINED;
      }
      else {
	throw (RbtFileError(_WHERE_,strError));
      }
    }

    m_ligAtomTypes.push_back(aType);
    if (iTrace > 1) {
      cout << "Using grid #" << aType << " for " << (*iter)->GetFullAtomName() << endl;
    }
  }
}

RbtDouble RbtVdwGridSF::RawScore() const {
  RbtDouble score = 0.0;
  
  //Check grids are defined
  if (m_grids.empty())
    return score;
  
  //Loop over all ligand atoms
  RbtAtomRListConstIter aIter = m_ligAtomList.begin();
  RbtTriposAtomTypeListConstIter tIter = m_ligAtomTypes.begin();
  if (m_bSmoothed) {
    for (; aIter != m_ligAtomList.end(); aIter++,tIter++) {
      score += m_grids[*tIter]->GetSmoothedValue((*aIter)->GetCoords());
    }
  }
  else {
    for (; aIter != m_ligAtomList.end(); aIter++,tIter++) {
      score += m_grids[*tIter]->GetValue((*aIter)->GetCoords());
    }
  }
  return score;
}

//Read grids from input stream, checking that header string matches RbtVdwGridSF
void RbtVdwGridSF::ReadGrids(istream& istr) throw (RbtError)
{
  m_grids.clear();
  RbtInt iTrace = GetTrace();
  
  //Read header string
  RbtInt length;
  Rbt::ReadWithThrow(istr, (char*) &length, sizeof(length));
  char* header = new char [length+1];
  Rbt::ReadWithThrow(istr, header, length);
  //Add null character to end of string
  header[length] = '\0';
  //Compare title with 
  RbtBool match = (_CT == header);
  delete [] header;
  if (!match) {
    throw RbtFileParseError(_WHERE_,"Invalid title string in " + _CT + "::ReadGrids()");
  }
  
  //Now read number of grids
  RbtInt nGrids;
  Rbt::ReadWithThrow(istr, (char*) &nGrids, sizeof(nGrids));
  if (iTrace > 0) {
    cout << _CT << ": reading " << nGrids << " grids..." << endl;
  }

  //We actually create a vector of size RbtTriposAtomType::MAXTYPES
  //Each grid is in the file is prefixed by the atom type string (e.g. C.2)
  //This string is converted to the corresponding RbtTriposAtomType::eType
  //and the grid stored at m_grids[eType]
  //This is slightly more transferable as it means grids can be read correctly even if the atom
  //type enums change (as long as the atom type string stay the same)
  //It also means we do not have to have a grid for each and every atom type if we don't want to
  m_grids = RbtRealGridList(RbtTriposAtomType::MAXTYPES);
  for (RbtInt i = 0; i < nGrids; i++) {
    //Read the atom type string
    Rbt::ReadWithThrow(istr, (char*) &length, sizeof(length));
    char* szType = new char [length+1];
    Rbt::ReadWithThrow(istr, szType, length);
    //Add null character to end of string
    szType[length] = '\0';
    RbtString strType(szType);
    delete [] szType;
    RbtTriposAtomType triposType;
    RbtTriposAtomType::eType aType = triposType.Str2Type(strType);
    if (iTrace > 0) {
      cout << "Grid# " << i << "\t" << "atom type=" << strType << " (type #" << aType << ")" << endl;
    }
    //Now we can read the grid
    RbtRealGridPtr spGrid(new RbtRealGrid(istr));
    m_grids[aType] = spGrid;
  }
}

//DM 25 Oct 2000 - track changes to parameter values in local data members
//ParameterUpdated is invoked by RbtParamHandler::SetParameter
void RbtVdwGridSF::ParameterUpdated(const RbtString& strName) {
  //DM 25 Oct 2000 - heavily used params
  if (strName == _SMOOTHED) {
    m_bSmoothed = GetParameter(_SMOOTHED);
  }
  else {
    RbtBaseSF::ParameterUpdated(strName);
  }
}


