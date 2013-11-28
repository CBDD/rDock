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

//XB cctype for check atom in MOL2 is number
//#include <cctype>
#include <sstream>
using std::stringstream;

#include "RbtMOL2FileSource.h"
#include "RbtFileError.h"
#include "RbtModelError.h"
#include "RbtAtomFuncs.h"

RbtString RbtMOL2FileSource::_CT("RbtMOL2FileSource"); 
// record delimiter strings
RbtString RbtMOL2FileSource::_TRIPOS_DELIM("@<TRIPOS>"); 
RbtInt RbtMOL2FileSource::_TRIPOS_DELIM_SIZE(RbtMOL2FileSource::_TRIPOS_DELIM.size());
RbtString RbtMOL2FileSource::_IDS_MOL2_RECDELIM("@<TRIPOS>MOLECULE");

// straightforward ctor 
RbtMOL2FileSource::RbtMOL2FileSource(const RbtString& fileName, RbtBool bImplHydrogens) : 
  RbtBaseMolecularFileSource(fileName,_IDS_MOL2_RECDELIM,"MOL2_FILE_SOURCE"),
  m_bImplHydrogens(bImplHydrogens)
{
  m_spElementData = RbtElementFileSourcePtr(new RbtElementFileSource(Rbt::GetRbtFileName("data","RbtElements.dat")));
  m_spParamSource = RbtParameterFileSourcePtr(new RbtParameterFileSource(Rbt::GetRbtFileName("data/sf","RbtIonicAtoms.prm")));
  m_NL = 0;	// make line counter zero
  // initialize MOLECULE descriptors
  nAtoms = nBonds = nSubstructures = nFeatures = nSets = 0;
  _RBTOBJECTCOUNTER_CONSTR_("RbtMOL2FileSource");
}

RbtMOL2FileSource::~RbtMOL2FileSource()
{
  _RBTOBJECTCOUNTER_DESTR_("RbtMOL2FileSource");
}

// parsing file and setting up atom descriptors
// at the beginning nothing is set, at the end everything should be set correctly
void RbtMOL2FileSource::Parse(void) throw (RbtError) 
{
  // start parsing if we haven't already done so
  if (!m_bParsedOK) { 
    // clear current cache
    ClearMolCache();
    // read the record, setting the delimiter to the beginning
    Read(false);
    // create a switchboard map - each member is a pointer to a function 
    // (a sort of functor, but I prefer to call using a meaningful string instead a number)
    map<RbtString,fcPtr> theSwitchBoard;
    theSwitchBoard["MOLECULE"] = &RbtMOL2FileSource::ParseRecordMOLECULE;
    theSwitchBoard["ATOM"] = &RbtMOL2FileSource::ParseRecordATOM;
    theSwitchBoard["BOND"] = &RbtMOL2FileSource::ParseRecordBOND;
    theSwitchBoard["SUBSTRUCTURE"] = &RbtMOL2FileSource::ParseRecordSUBSTRUCTURE;
    theSwitchBoard["UNSUPPORTED"] = &RbtMOL2FileSource::ParseRecordUNSUPPORTED;
    // the real parser
    try {
      // it would be nice to have a regexp stuff here
      // but the patterns are rather simple so it is not that bad
      // --------------------------------------------------------
      // declare the callback switched parser
      void (RbtMOL2FileSource::*theCurrentParser)(const RbtString&);
      // first assume the parser function is that that handles the MOLECULE record
      theCurrentParser	= &RbtMOL2FileSource::ParseRecordMOLECULE;
      // than get lines one-by-one
      for(RbtFileRecListIter lineIter = m_lineRecs.begin();
	  lineIter < m_lineRecs.end();
	  ++lineIter ) {
	//Ignore blank lines and comment lines
	if ( (*lineIter).empty() || ((*lineIter).at(0) == '#')) {
	  //cout << "Skipping blank/comment record" << endl;
	  continue;
	}
	// check wether is it a record field delimiter (tag like: @<TRIPOS>BLAHBLAH )line
	RbtString theNewTag = GetMOL2Tag( (*lineIter) );
	if ( !theNewTag.empty() ) {
	  // check key existence
	  map<RbtString,fcPtr>::const_iterator i = theSwitchBoard.find(theNewTag);
	  if(i == theSwitchBoard.end()){
	    //cout << "INFO Skipping unsupported record: " << theNewTag << endl;
	    theCurrentParser = theSwitchBoard["UNSUPPORTED"];
	  }
	  else {
	    // cout << "Reading " << theNewTag << endl;
	    theCurrentParser	= theSwitchBoard[theNewTag];
	  }
	} else { // call the parser
	  (*this.*theCurrentParser)( (*lineIter) );
	}
      }
      
      // we are ready with parsing
      // let's assign chain ids to atoms
      //and sort atomlist by substructure ID
      m_atomList.clear();
      for(map<RbtInt,RbtAtomList>::const_iterator iter = m_ssAtoms.begin(); iter != m_ssAtoms.end(); iter++) {
	RbtInt subst_id = iter->first;
	RbtAtomList ssAtomList = iter->second;
	//cout << "Found " << ssAtomList.size() << " atoms in SUBSTRUCTURE#" << subst_id << endl;
	std::copy(ssAtomList.begin(),ssAtomList.end(),std::back_inserter(m_atomList));
	RbtMOL2SubstructureMapConstIter ssIter = m_ssInfo.find(subst_id);
	if (ssIter != m_ssInfo.end()) {
	  RbtMOL2Substructure ss = ssIter->second;
	  RbtString subst_name = ss.GetName();
	  RbtString chain = ss.GetChain();
	  RbtString sub_type = ss.GetType();
	  RbtString sID,sName;
	  GetSSIDandName(subst_name,subst_id,sID,sName);
	  for (RbtAtomListIter aIter = ssAtomList.begin(); aIter != ssAtomList.end(); aIter++) {
	    RbtString atom_sID = (*aIter)->GetSubunitId();
	    RbtString atom_sName = (*aIter)->GetSubunitName();
	    if ( (sID != atom_sID) || (sName != atom_sName)) {
	      cout << _CT << ": WARNING Mismatch between SUBSTRUCTURE (" << sName << sID << ") and ATOM ("
		   << atom_sName << atom_sID << ") info" << endl;
	    }
	    if ( (atom_sName == "UNK") && (sub_type != "UNK") ) {
	      cout << _CT << ": INFO Could use sub_type (" << sub_type << ") to update SubunitName for "
		   << (*aIter)->GetFullAtomName() << endl;
	    }
	    (*aIter)->SetSegmentName(chain);
	    //cout << "Updating segment name to " << chain << " for " << (*aIter)->GetFullAtomName() << endl;
	  }
	}
	else {
	  //cout << _CT << ": INFO No SUBSTRUCTURE record found for subst_id=" << subst_id << endl;
	}
	//Set the group charges on a residue-by-residue basis
	//so we can check for presence or absence of all required atoms
	//e.g. in HIS need both polar Hs (HD1/HE2) to be present
	SetupPartialIonicGroups(ssAtomList,m_spParamSource);
      }

      if (m_atomList.size() != nAtoms) {
	cout << _CT << ": WARNING Mismatch in #ATOMS read: Expected=" << nAtoms << "; Actual = " << m_atomList.size() << endl;
      }
      if (m_bondList.size() != nBonds) {
	cout << _CT << ": WARNING Mismatch in #BONDS read: Expected=" << nBonds << "; Actual = " << m_bondList.size() << endl;
      }
      if (m_ssInfo.size() != nSubstructures) {
	cout << _CT << ": WARNING Mismatch in #SUBSTRUCTURES read: Expected=" << nSubstructures << "; Actual = " << m_ssInfo.size() << endl;
      }
      
      SetupAtomParams();

      for (RbtAtomListConstIter iter = m_atomList.begin(); iter != m_atomList.end(); ++iter) {
	m_segmentMap[(*iter)->GetSegmentName()]++;//increment atom count in segment map
	RbtTriposAtomType::eType tt_file = (*iter)->GetTriposType();
	RbtTriposAtomType::eType tt_auto = m_typer(*iter,true);
	if (tt_file == RbtTriposAtomType::UNDEFINED) {
	  cout << _CT << ": WARNING Undefined Tripos type for " << (*iter)->GetFullAtomName() << endl;
	}
	if (tt_file != tt_auto) {
#ifdef _DEBUG
	  cout << _CT << ": INFO Disagreement in Tripos atom types for " << (*iter)->GetFullAtomName();
	  cout << ": File=" << m_typer.Type2Str(tt_file) << "; Auto=" << m_typer.Type2Str(tt_auto) << endl;
#endif //_DEBUG
	}
      }

      m_bParsedOK = true;
    } 
    catch (RbtError& error) {
      ClearMolCache();
      throw; 
    }
  }
}

// @<TRIPOS>MOLECULE
// six data lines (after tag)
// ------------------------------
// 1 name
// ------------------------------
// 2 numer of 
// atoms [ bonds substructures features sets ]
// ------------------------------
// 3 type 
// [SMALL|PROTEIN|NUCLEIC_ACID|SACCHARIDE]
// ------------------------------
// 4 type of charges
// [NO_CHARGES|DEL_RE|GASTEIGER|GAST_HUCK|HUCKEL|PULLMAN|GAUSS80_CHARGES
// |AMPAC_CHARGES|MULLIKEN_CHARGES|DICT_CHARGES|MMFF94_CHARGES|USER_CHARGES]
// ------------------------------
// 5 status_bits internal SYBYL stuff, should never be set by the user 
// ------------------------------
// 6 one comment line 
// ------------------------------
//
// it is the m_titleList that is set up here only
// 
void RbtMOL2FileSource::ParseRecordMOLECULE(const RbtString& aLine)
{
#ifdef _DEBUG_
  cout << "MOLECULE " << aLine << "   " << m_NL << endl; 
#endif // _DEBUG_
  
  switch(m_NL){
  case 0:		// first line
    m_titleList.push_back(aLine);	// put it into title
    m_NL++;							// increment record line counter
    break;
  case 1:		// second line
    ParseCountFields(aLine);
    m_NL++;
    break;
  default:	// do nothing
    break;
  }
}

void RbtMOL2FileSource::ParseRecordATOM(const RbtString& aLine) {
  RbtStringList tokens;
  Tokenize(aLine,tokens);
  // there must be at least 6 fields, rest is optional
  if(tokens.size() < 6)
    throw RbtFileParseError(_WHERE_, "Corrupted MOL2 file: not enough fields in ATOM ");

  //Compulsory fields
  RbtInt atom_id = atoi(tokens[0].c_str());
  RbtString atom_name = tokens[1];
  RbtCoord atom_coord;
  atom_coord.x = atof(tokens[2].c_str());
  atom_coord.y = atof(tokens[3].c_str());
  atom_coord.z = atof(tokens[4].c_str());
  RbtString atom_type = tokens[5];

  //Optional fields.
  RbtInt subst_id = (tokens.size() > 6) ? atoi(tokens[6].c_str()) : 1;
  RbtString subst_name = (tokens.size() > 7 ) ? tokens[7] : "****";
  RbtString sID,sName;
  GetSSIDandName(subst_name,subst_id,sID,sName);
  RbtDouble charge = (tokens.size() > 8) ? atof(tokens[8].c_str()) : 0.0;
// XB reweighting parameters
//	RbtDouble wxb = (tokens.size() > 9) ? atof(tokens[9].c_str()) : 1.0;
//XB mod for only doing if number, not characters	
//	RbtDouble wxb;
//  if (tokens.size() > 9 && isdigit(tokens[9][0])){
//	  wxb = atof(tokens[9].c_str());
//	} else {
//	  wxb = 1.0;
//	}
//end XB reweighting parameters

  //Derived atom params (some may be updated later)
  RbtTriposAtomType::eType tt = m_typer.Str2Type(atom_type);
  RbtInt atomic_number = m_typer.Type2AtomicNo(tt);
  //Workaround for metals such as Zn, Mn, Mg
  //We should add Tripos types for these metals, but the most important thing is
  //to pick up the atomic number (can do this from the atom type string)
  if (tt == RbtTriposAtomType::UNDEFINED) {
    RbtString el = atom_type.substr(0,2);
    cout << _CT << ": INFO Attempting to identify element by 1st 2 chars of FFType string: " << el << endl;
    if (m_spElementData->isElementNamePresent(el)) {
      atomic_number = m_spElementData->GetElementData(el).atomicNo;
      cout << _CT << ": INFO Atomic number = " << atomic_number << endl;
    }
    else {
      cout << _CT << ": INFO Element not found" << endl;
    }
  }
  RbtElementData elementData = m_spElementData->GetElementData(atomic_number);
  RbtAtom::eHybridState hybrid_state =  (RbtAtom::eHybridState) m_typer.Type2Hybrid(tt);
  RbtAtomPtr newAtom(new RbtAtom(atom_id,
				 atomic_number,
				 atom_name,
				 sID,// subunit id undef 
				 sName,// subunit name undef
				 "UNK"));// segment name undefined at this point (need SUBSTRUCTURE records)
  newAtom->SetCoords(atom_coord);
  newAtom->SetFFType(atom_type);//Original Tripos type in string format
  newAtom->SetTriposType(tt);//Same type, but stored as an enum (may be updated later)
  newAtom->SetPartialCharge(charge);
  newAtom->SetHybridState(hybrid_state);
  newAtom->SetAtomicMass(elementData.mass);
//  newAtom->SetReweight(wxb); //XB

  m_atomList.push_back(newAtom);
  m_ssAtoms[subst_id].push_back(newAtom);
  //cout << "ParseRecordATOM: " << *newAtom << endl;
}

void RbtMOL2FileSource::ParseRecordBOND(const RbtString& aLine) {
  RbtStringList tokens;
  Tokenize(aLine,tokens);
  // there must be at least 4 fields, rest is optional
  if(tokens.size() < 4)
    throw RbtFileParseError(_WHERE_, "Corrupted MOL2 file: not enough fields in BOND");

  //Compulsory fields
  RbtInt bond_id = atoi(tokens[0].c_str()); 
  RbtInt origin_bond_id	= atoi(tokens[1].c_str());	
  RbtInt target_bond_id = atoi(tokens[2].c_str());
  RbtString bond_type = tokens[3];

  //In rDock we only have integer bond order
  //Aromatics should be defined as alternate double/single bonds as in MDL mol files
  //Can't do much here so define everything as single, double or triple
  //We have the hybridisation state of each atoms already (from the Tripos type)
  //so are not so dependent on bond orders to deduce hybridisation
  RbtInt bond_order = (bond_type == "3") ? 3 : (bond_type == "2") ? 2 : 1;
  //Check that atom IDs are in range
  //This assumes that the Bond records are read after the Atom records (should be always the case)
  if ( (origin_bond_id < 1) || (origin_bond_id > m_atomList.size()) ||
       (target_bond_id < 1) || (target_bond_id > m_atomList.size()) )
    throw RbtFileParseError(_WHERE_, "atom_id out of range in BOND record\n"+aLine);

  //Decrement the atom indices as the atoms are numbered from zero in our atom vector
  RbtAtomPtr spAtom1(m_atomList[origin_bond_id-1]);
  RbtAtomPtr spAtom2(m_atomList[target_bond_id-1]);
  RbtBondPtr spBond( new RbtBond(bond_id,spAtom1,spAtom2,bond_order));
  m_bondList.push_back(spBond);
  //cout << "ParseRecordBOND: " << *spBond << endl;
}

void RbtMOL2FileSource::ParseRecordSUBSTRUCTURE(const RbtString& aLine) {
  RbtStringList	tokens;
  Tokenize(aLine,tokens);
  // there must be at least 3 fields, rest is optional
  if(tokens.size() < 3)
    throw RbtFileParseError(_WHERE_, "Corrupted MOL2 file: not enough fields in SUBSTRUCTURE");

  //Compulsory fields
  RbtInt subst_id = atoi(tokens[0].c_str());
  RbtString subst_name = tokens[1];
  RbtInt root_atom = atoi(tokens[2].c_str());

  //Optional fields	
  RbtString chain = ( (tokens.size() > 5) && (tokens[5] != "****") ) ? tokens[5] : "UNK";
  RbtString sub_type = ( (tokens.size() > 6) && (tokens[6] != "****") ) ? tokens[6] : "UNK"; 

  RbtMOL2Substructure ss(subst_name,root_atom,chain,sub_type);
  m_ssInfo[subst_id] = ss;
  //cout << "ParseRecordSUBSTRUCTURE: " << subst_id << "," << subst_name << "," << root_atom << ","
  //     << chain << "," << sub_type << endl;
}

void RbtMOL2FileSource::ParseRecordUNSUPPORTED(const RbtString& aLine)
{
  //	do nothing...
}

RbtString RbtMOL2FileSource::GetMOL2Tag(const RbtString& aLine) {
  // template based tokenizer is at boost.org, but we have to wait until 
  // they are ready to provide Makefile instead of the crap jam they have
  // so, look for a line starting with "@<TRIPOS>"
  if( aLine.find(_TRIPOS_DELIM) != string::npos) {
    // tokenize remainder line in case of garbage at the rest of line
    RbtStringList tokens;
    Tokenize(aLine,tokens);
    //get rid of the "@<TRIPOS>" part 
    if(tokens[0].size() <= _TRIPOS_DELIM_SIZE) {
      cout << "Corrupted MOL2: no tag identifier after @<TRIPOS>." << endl;
      // why the heck is a segfault here?
      throw RbtFileParseError(_WHERE_,"Corrupted MOL2: no tag identifier after @<TRIPOS>.");
    }
    else {
      tokens[0].erase(0,_TRIPOS_DELIM_SIZE);
      return tokens[0];	// return with the very first string
    }
  } else
    // return empty string if it is not a MOL2 tag
    return "";	
}

void RbtMOL2FileSource::ParseCountFields(const RbtString& aLine)
{
  RbtStringList	tokens;
  Tokenize(aLine,tokens);
  if(tokens.size()<1) 
    throw RbtFileParseError(_WHERE_,"Corrupted MOL2: no atom numbers in MOLECULE record.");
  // number of atoms:
  nAtoms = atoi(tokens[0].c_str());
  if(tokens.size()>1)	// we have bonds number as well
    nBonds = atoi(tokens[1].c_str());
  if(tokens.size()>2)	// substructures
    nSubstructures = atoi(tokens[2].c_str());
  if(tokens.size()>3)	
    nFeatures = atoi(tokens[3].c_str());
  if(tokens.size()>4)	// substructures
    nSets = atoi(tokens[4].c_str());
  //cout << "nAtoms " << nAtoms << endl;
  //cout << "nBonds " << nBonds << endl;
  //cout << "nSubstructures " << nSubstructures << endl;
  //cout << "nFeatures " << nFeatures << endl;
  //cout << "nSets "<< nSets << endl;
}

//Correct all the atom attributes that could not be set explicitly when
//parsing the atom records
void RbtMOL2FileSource::SetupAtomParams() throw (RbtError)
{
  FixImplicitHydrogenCount();
  if (m_bImplHydrogens) {
    RemoveNonPolarHydrogens();
    //Quick fix: remove all "Lone Pair" atoms as well
    RbtAtomList lpList = Rbt::GetAtomList(m_atomList,Rbt::isFFType_eq("LP"));
    for (RbtAtomListIter iter = lpList.begin(); iter != lpList.end(); iter++) {
      //cout << "INFO Removing Lone Pair " << (*iter)->GetFullAtomName() << endl;
      RemoveAtom(*iter);
    }
    RenumberAtomsAndBonds();
  }
  FixHybridState();
  FixTriposTypes();
  SetupVdWRadii();
  RenumberAtomsAndBonds();
}

//We don't have the benefit of accurate bond orders (as we haven't bothered to separate
//the aromatic bonds into alternate double/single formal bonds)
//but we do have accurate hybridisation states from the Tripos type lookup.
//Can detect implicit hydrogens in the input file by a simple check of the
//number of connections to each carbon
void RbtMOL2FileSource::FixImplicitHydrogenCount() {
  //Get list of all carbons
  RbtAtomList cList = Rbt::GetAtomList(m_atomList,Rbt::isAtomicNo_eq(6));
  for (RbtAtomListIter iter = cList.begin(); iter != cList.end(); iter++) {
    RbtAtom::eHybridState hyb = (*iter)->GetHybridState();
    //Compare actual and expected bond count based on hybridisation state
    //Asssume SP3 (4 bonds) if hyb state is undefined
    RbtInt actual = (*iter)->GetNumBonds();
    RbtInt expected = (hyb==RbtAtom::SP3) ? 4 : (hyb==RbtAtom::SP2) ? 3 : (hyb==RbtAtom::AROM) ? 3 : (hyb==RbtAtom::SP) ? 2 : 4;
    RbtInt nImplH = expected-actual;
    //Pick up fragmented aromatic rings here
    //Cannot have more than one explicit or implicit hydrogen to an aromatic carbon
    if ( (hyb == RbtAtom::AROM) && (nImplH > 0) ) {
      RbtInt nExplH = (*iter)->GetCoordinationNumber(1);
      if ( (nImplH + nExplH) > 1) {
	cout << _CT << ": INFO Too few bonds detected to aromatic carbon " << (*iter)->GetFullAtomName() << "; #hydrogens capped at 1" << endl;
	nImplH = 1 - nExplH;
      }
    }
    if (nImplH > 0) {
      (*iter)->SetNumImplicitHydrogens(nImplH);
    }
    else if (nImplH < 0) {
      cout << _CT << ": INFO Too many bonds detected to " << (*iter)->GetFullAtomName() << endl;
    }
  }   
}

//Convert {O,S}_SP3 to TRI, if and only if the atom is bonded to pi-atom
void RbtMOL2FileSource::FixHybridState() {
  RbtAtomList sp3AtomList = Rbt::GetAtomList(m_atomList,Rbt::isHybridState_eq(RbtAtom::SP3));
  sp3AtomList = Rbt::GetAtomList(sp3AtomList,std::not1(Rbt::isAtomPosCharged()));
  
  for (RbtAtomListIter iter = sp3AtomList.begin(); iter != sp3AtomList.end(); iter++) {
    RbtAtomList bondedAtomList;
    switch ((*iter)->GetAtomicNo()) {
    case 8://O
    case 16://S
      bondedAtomList = Rbt::GetBondedAtomList(*iter);
      if (Rbt::FindAtom(bondedAtomList,Rbt::isPiAtom()) != bondedAtomList.end()) {
	(*iter)->SetHybridState(RbtAtom::TRI);
	//cout << "Changing " << (*iter)->GetFullAtomName() << " from SP3 to TRI" << endl;
      }
      break;
    default:
      break;
    }
  }
}

//Correct the Tripos types read from file, for extended carbons and polar hydrogens
void RbtMOL2FileSource::FixTriposTypes() {
  RbtAtomList cList = Rbt::GetAtomList(m_atomList,Rbt::isAtomicNo_eq(6));
  for (RbtAtomListIter iter = cList.begin(); iter != cList.end(); iter++) {
    RbtTriposAtomType::eType t = (*iter)->GetTriposType();
    RbtTriposAtomType::eType xt = t;
    RbtInt nImplH = (*iter)->GetNumImplicitHydrogens();
    if (nImplH > 0) {
      switch (t) {
      case RbtTriposAtomType::C_3:
	xt = (nImplH == 3) ? RbtTriposAtomType::C_3_H3 :
	  (nImplH == 2) ? RbtTriposAtomType::C_3_H2 :
	  (nImplH == 1) ? RbtTriposAtomType::C_3_H1 : RbtTriposAtomType::C_3;
	break;
      case RbtTriposAtomType::C_2:
	xt = (nImplH == 2) ? RbtTriposAtomType::C_2_H2 :
	  (nImplH == 1) ? RbtTriposAtomType::C_2_H1 : RbtTriposAtomType::C_2;
	break;
      case RbtTriposAtomType::C_ar:
	xt = (nImplH == 1) ? RbtTriposAtomType::C_ar_H1 : RbtTriposAtomType::C_ar;
	break;
      case RbtTriposAtomType::C_1:
	xt = (nImplH == 1) ? RbtTriposAtomType::C_1_H1 : RbtTriposAtomType::C_1;
	break;
      default:
	break;
      }
      if (xt != t) {
	//	cout << "Correcting Tripos type for " << (*iter)->GetFullAtomName() << " from "
      	//     << m_typer.Type2Str(t) << " to " << m_typer.Type2Str(xt) << endl;
      	(*iter)->SetTriposType(xt);
      }
    }
  }
  RbtAtomList hpList = Rbt::GetAtomList(m_atomList,Rbt::isAtomHBondDonor());
  for (RbtAtomListIter iter = hpList.begin(); iter != hpList.end(); iter++) {
    RbtTriposAtomType::eType t = (*iter)->GetTriposType();
    RbtTriposAtomType::eType xt = RbtTriposAtomType::H_P;
    if (xt != t) {
      //cout << "Correcting Tripos type for " << (*iter)->GetFullAtomName() << " from "
      //   << m_typer.Type2Str(t) << " to " << m_typer.Type2Str(xt) << endl;
      (*iter)->SetTriposType(xt);
    }
  }
}


void RbtMOL2FileSource::RemoveNonPolarHydrogens() {
  RbtAtomList cList = Rbt::GetAtomList(m_atomList,Rbt::isAtomicNo_eq(6));
  for (RbtAtomListIter cIter = cList.begin(); cIter != cList.end(); cIter++) {
    //Get list of all bonded hydrogens
    RbtAtomList hList = Rbt::GetAtomList(Rbt::GetBondedAtomList(*cIter),Rbt::isAtomicNo_eq(1));
    RbtInt nH = hList.size();
    if (nH > 0) {
      for (RbtAtomListIter hIter = hList.begin(); hIter != hList.end(); hIter++) {
	RemoveAtom(*hIter);
      }
      //Adjust number of implicit hydrogens
      (*cIter)->SetNumImplicitHydrogens((*cIter)->GetNumImplicitHydrogens()+nH);
      //cout << "Removing " << nH << " hydrogens from " << (*cIter)->GetFullAtomName() << endl;
    }
  }
  //Remove orphan hydrogens (from fragmented residues for e.g.)
  RbtAtomList hList = Rbt::GetAtomList(m_atomList,Rbt::isAtomicNo_eq(1));
  hList = Rbt::GetAtomList(hList,Rbt::isCoordinationNumber_eq(0));
  for (RbtAtomListIter hIter = hList.begin(); hIter != hList.end(); hIter++) {
    cout << _CT << ": INFO Removing orphan hydrogen " << (*hIter)->GetFullAtomName() << endl;
    RemoveAtom(*hIter);
  }
}

//Defines nominal vdW radius (used for cavity mapping)
//correcting for extended atoms and H-bond donor hydrogens
void RbtMOL2FileSource::SetupVdWRadii() {
  //Radius increment for atoms with implicit hydrogens
  //DM 22 Jul 1999 - only increase the radius for sp3 atoms with implicit hydrogens
  //For sp2 and aromatic, leave as is
  Rbt::isHybridState_eq bIsSP3(RbtAtom::SP3);
  RbtDouble dImplRadIncr = m_spElementData->GetImplicitRadiusIncr();
  //Element data for hydrogen
  RbtElementData elHData = m_spElementData->GetElementData(1);
  //Radius increment and predicate for H-bonding hydrogens
  Rbt::isAtomHBondDonor bIsHBondDonor;
  RbtDouble dHBondRadius = elHData.vdwRadius + m_spElementData->GetHBondRadiusIncr();
  
  for (RbtAtomListIter iter = m_atomList.begin(); iter != m_atomList.end(); iter++) {
    //Get the element data for this atom
    RbtInt nAtomicNo = (*iter)->GetAtomicNo();
    RbtElementData elData = m_spElementData->GetElementData(nAtomicNo);
    RbtDouble vdwRadius = elData.vdwRadius;
    RbtInt nImplH = (*iter)->GetNumImplicitHydrogens();
    //Adjust atomic mass and vdw radius for atoms with implicit hydrogens
    if (nImplH > 0) {
      (*iter)->SetAtomicMass(elData.mass+(nImplH*elHData.mass));//Adjust atomic mass
      if (bIsSP3(*iter)) {
	vdwRadius += dImplRadIncr;//adjust vdw radius (for sp3 implicit atoms only)
      }
    }
    //Adjust vdw radius for H-bonding hydrogens
    else if (bIsHBondDonor(*iter)) {
      vdwRadius = dHBondRadius;
    }
    //Finally we can set the radius
    (*iter)->SetVdwRadius(vdwRadius);
#ifdef _DEBUG
    //cout << (*iter)->GetFullAtomName() << ": #H=" << nImplH
    //<< "; vdwR=" << (*iter)->GetVdwRadius()
    //<< "; mass=" << (*iter)->GetAtomicMass() << endl;
#endif //_DEBUG
  }
}

//Splits the Tripos MOL2 subst_name field into appropriate residue name (sName) and residue ID (sID) components
//If subst_name is **** or in the style <148>, then substitutes default name (UNK) and default ID (subst_id)
void RbtMOL2FileSource::GetSSIDandName(const RbtString& subst_name, RbtInt subst_id, RbtString& sID, RbtString& sName) {
  //subst_names beginning with * or < do not follow the normal pattern of ALA99 etc
  const RbtString nonstandard("*<");
  const RbtString numeric("0123456789");
  const RbtString defaultName("UNK");
  const RbtString defaultID = RbtVariant(subst_id);//Quick way to convert from int to string

  if ( subst_name.empty() || (subst_name.find_first_of(nonstandard) == 0)) {
    sID = defaultID;
    sName = defaultName;
    return;
  }
  //Find first numeric char
  RbtString::size_type nFirst = subst_name.find_first_of(numeric);
  //Split into name and ID:
  //Name is up to the beginning of 1st numeric char
  //ID is from the first numeric char onwards (may have non-numeric suffix as well e.g 100B)
  sID = (nFirst != RbtString::npos) ? subst_name.substr(nFirst) : defaultID;
  sName = (nFirst > 0) ? subst_name.substr(0,nFirst) : defaultName;
  return;
}

// simple tokenizer that splits at whitspace
void RbtMOL2FileSource::Tokenize(const RbtString& aString, RbtStringList& aTokensBuf) {
  RbtString buf;
  stringstream ss(aString);
  while(ss >> buf)
    aTokensBuf.push_back(buf);
}
