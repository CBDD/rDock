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

#include <functional> //For STL function objects
using std::ptr_fun;
using std::bind2nd;

#include "RbtAtom.h"
#include "RbtPseudoAtom.h"
#include "RbtBond.h"
#include "RbtModel.h"

///////////////////////////////////////////////
//Constructors / destructors
///////////////////////////////////////////////

//Default constructor
//Initialise 2-D attributes to suitable defaults
//Initialise 3-D attributes to zero/null
RbtAtom::RbtAtom() :
  m_nAtomId(0),
  m_nAtomicNo(6),
  m_strAtomName("C"),
  m_strSubunitId("1"),
  m_strSubunitName("RES"),
  m_strSegmentName("SEG1"),
  m_eState(UNDEFINED),//DM 8 Dec 1998 Changed from SP3 to UNDEFINED
  m_nHydrogens(0),
  m_nFormalCharge(0),
  m_coord(0.0,0.0,0.0),
  m_dPartialCharge(0.0),
  m_dGroupCharge(0.0),
  m_dAtomicMass(0.0),
  m_dVdwRadius(0.0),
  m_strFFType(""),
  m_pModel(NULL),
  m_bCyclic(false),
  m_bSelected(false),
  m_bUser1(false),
  m_dUser1(1.0),//DM 27 Jul 2000 - initialise user values to 1 as they are commonly used as weightings
  m_dUser2(1.0),
  m_nPMFType(PMF_UNDEFINED),
  m_triposType(RbtTriposAtomType::UNDEFINED)
{
  _RBTOBJECTCOUNTER_CONSTR_("RbtAtom");
}

//Constructor supplying all 2-D parameters
//Initialise 3-D attributes to zero/null
RbtAtom::RbtAtom(RbtInt nAtomId,
		 RbtInt nAtomicNo /*= 6*/,
		 RbtString strAtomName /*= "C"*/,
		 RbtString strSubunitId /*= "1"*/,
		 RbtString strSubunitName /*= "RES"*/,
		 RbtString strSegmentName /*= "SEG1"*/,
		 eHybridState eState /*= UNDEFINED*/,
		 RbtUInt nHydrogens /*= 0*/,
		 RbtInt nFormalCharge /*= 0.0*/
		 ) :
  m_nAtomId(nAtomId),
  m_nAtomicNo(nAtomicNo),
  m_strAtomName(strAtomName),
  m_strSubunitId(strSubunitId),
  m_strSubunitName(strSubunitName),
  m_strSegmentName(strSegmentName),
  m_eState(eState),
  m_nHydrogens(nHydrogens),
  m_nFormalCharge(nFormalCharge),
  m_coord(0.0,0.0,0.0),
  m_dPartialCharge(0.0),
  m_dGroupCharge(0.0),
  m_dAtomicMass(0.0),
  m_dVdwRadius(0.0),
  m_strFFType(""),
  m_pModel(NULL),
  m_bCyclic(false),
  m_bSelected(false),
  m_bUser1(false),
  m_dUser1(1.0),//DM 27 Jul 2000 - initialise user values to 1 as they are commonly used as weightings
  m_dUser2(1.0),
  m_nPMFType(PMF_UNDEFINED),
  m_triposType(RbtTriposAtomType::UNDEFINED)
{
  _RBTOBJECTCOUNTER_CONSTR_("RbtAtom");
}

//Default destructor
RbtAtom::~RbtAtom()
{
  _RBTOBJECTCOUNTER_DESTR_("RbtAtom");
}

//Copy constructor
RbtAtom::RbtAtom(const RbtAtom& atom)
{
  m_nAtomId = atom.m_nAtomId;
  m_nAtomicNo = atom.m_nAtomicNo;
  m_strAtomName = atom.m_strAtomName;
  m_strSubunitId = atom.m_strSubunitId;
  m_strSubunitName = atom.m_strSubunitName;
  m_strSegmentName = atom.m_strSegmentName;
  m_eState = atom.m_eState;
  m_nHydrogens = atom.m_nHydrogens;
  m_nFormalCharge = atom.m_nFormalCharge;
  m_coord = atom.m_coord;
  m_dPartialCharge = atom.m_dPartialCharge;
  m_dGroupCharge = atom.m_dGroupCharge;
  m_dAtomicMass = atom.m_dAtomicMass;
  m_dVdwRadius = atom.m_dVdwRadius;
  m_strFFType = atom.m_strFFType;
  m_savedCoords = atom.m_savedCoords;
  m_bUser1 = atom.m_bUser1;
  m_dUser1 = atom.m_dUser1;
  m_dUser2 = atom.m_dUser2;
  m_nPMFType = atom.m_nPMFType;
  m_triposType = atom.m_triposType;
  //Copied atoms no longer belong to the model so set to NULL here
  SetModelPtr(NULL);
  //Copied atoms no longer belong to the bonds so erase the bond map
  ClearBondMap();
  //Set the cyclic flag to false as the atom isn't bonded to anything
  SetCyclicFlag(false);
  //Reset the selection flag also
  SetSelectionFlag(false);
  _RBTOBJECTCOUNTER_COPYCONSTR_("RbtAtom");
}

//Copy assignment
RbtAtom& RbtAtom::operator=(const RbtAtom& atom)
{
  if (this != &atom) { //beware of self-assignment
    m_nAtomId = atom.m_nAtomId;
    m_nAtomicNo = atom.m_nAtomicNo;
    m_strAtomName = atom.m_strAtomName;
    m_strSubunitId = atom.m_strSubunitId;
    m_strSubunitName = atom.m_strSubunitName;
    m_strSegmentName = atom.m_strSegmentName;
    m_eState = atom.m_eState;
    m_nHydrogens = atom.m_nHydrogens;
    m_nFormalCharge = atom.m_nFormalCharge;
    m_coord = atom.m_coord;
    m_dPartialCharge = atom.m_dPartialCharge;
    m_dGroupCharge = atom.m_dGroupCharge;
    m_dAtomicMass = atom.m_dAtomicMass;
    m_dVdwRadius = atom.m_dVdwRadius;
    m_strFFType = atom.m_strFFType;
    m_savedCoords = atom.m_savedCoords;
    m_bUser1 = atom.m_bUser1;
    m_dUser1 = atom.m_dUser1;
    m_dUser2 = atom.m_dUser2;
    m_nPMFType = atom.m_nPMFType;
    m_triposType = atom.m_triposType;
    //Copied atoms no longer belong to the model so set to NULL here
    SetModelPtr(NULL);
    //Copied atoms no longer belong to the bonds so erase the bond map
    ClearBondMap();
    //Set the cyclic flag to false as the atom isn't bonded to anything
    SetCyclicFlag(false);
    //Reset the selection flag also
    SetSelectionFlag(false);
  }
  return *this;
}
 

///////////////////////////////////////////////
//Stream functions
///////////////////////////////////////////////

//Insertion operator (primarily for debugging)
ostream& operator<<(ostream& s, const RbtAtom& atom)
{
  return atom.Print(s);
}

//Virtual function for dumping atom details to an output stream
//Derived classes (e.g. pseudoatom) can override if required
ostream& RbtAtom::Print(ostream& s) const
{
  //Get owning model name (if any)
  RbtString strModelName = "Orphan";
  if (m_pModel != NULL)
    strModelName = m_pModel->GetName();

  RbtTriposAtomType triposType;

  return s << "(" << strModelName << ")"
	   << ", ID=" << m_nAtomId
	   << ", AtNo=" << m_nAtomicNo
	   << ", " << m_strSegmentName
	   << "," << m_strSubunitName
	   << "," << m_strSubunitId
	   << "," << m_strAtomName
	   << ", Hybrid=" << Rbt::ConvertHybridStateToString(m_eState)
	   << ", #Hyd=" << m_nHydrogens
	   << ", GrpChg=" << m_dGroupCharge
	   << ", Tripos=" << triposType.Type2Str(GetTriposType())
//	   << ", FFType=" << m_strFFType
	   << ", Cyclic=" << m_bCyclic;
}



///////////////////////////////////////////////
//Public accessor functions
//2-D attributes
///////////////////////////////////////////////

RbtBool RbtAtom::GetEnabled() const {
  return (m_pModel) ? m_pModel->GetEnabled() : true;
}
  
//Returns composite of segment name, subunit id and name, and atom name
RbtString RbtAtom::GetFullAtomName() const
{
  return m_strSegmentName + ":" + m_strSubunitName + "_" + m_strSubunitId + ":" + m_strAtomName;
}

//DM 04 Dec 1998  Add functions to handle bond map
//Add a bond to the bond map - returns true if OK, false if this atom is not a member of the bond
//or if bond already added
RbtBool RbtAtom::AddBond(RbtBond* pBond)
{
  //Case 1 - this atom is the first atom in the bond
  //(key=pBond, value = true) as this atom is atom 1
  if (pBond->GetAtom1Ptr().Ptr() == this) {
    pair<RbtBondMapIter,RbtBool> p = m_bondMap.insert(make_pair(pBond,true));
    return p.second;//p.second is true if bond inserted, false if bond was already present
  }
  //Case 2 - this atom is the second atom in the bond
  //(key=pBond, value = false) as this atom is atom 2
  else if (pBond->GetAtom2Ptr().Ptr() == this) {
    pair<RbtBondMapIter,RbtBool> p = m_bondMap.insert(make_pair(pBond,false));
    return p.second;//p.second is true if bond inserted, false if bond was already present
  }
  else
    return false;
}

//Remove a bond from the bond map - returns true if OK, false if bond not present in map
RbtBool RbtAtom::RemoveBond(RbtBond* pBond)
{
  if (m_bondMap.erase(pBond) > 0) //erase returns the number of element deleted
    return true;
  else
    return false;
}

//Returns number of cyclic bonds in map
RbtUInt RbtAtom::GetNumCyclicBonds() const
{
  RbtUInt nCount(0);
  for (RbtBondMapConstIter iter = m_bondMap.begin(); iter != m_bondMap.end(); iter++)
    if ((*iter).first->GetCyclicFlag())
      nCount++;
  return nCount;
}

//Returns bond map with only cyclic bonds included
RbtBondMap RbtAtom::GetCyclicBondMap() const
{
  RbtBondMap cyclicBondMap;
  for (RbtBondMapConstIter iter = m_bondMap.begin(); iter != m_bondMap.end(); iter++)
    if ((*iter).first->GetCyclicFlag())
      cyclicBondMap.insert(*iter);
  return cyclicBondMap;
}

//Returns total formal bond order around atom
RbtInt RbtAtom::GetTotalFormalBondOrder() const
{
  RbtInt nTot(0);
  for (RbtBondMapConstIter iter = m_bondMap.begin(); iter != m_bondMap.end(); iter++)
    nTot += (*iter).first->GetFormalBondOrder();
  return nTot;
}

//Returns max formal bond order for all bonds to atom
RbtInt RbtAtom::GetMaxFormalBondOrder() const
{
  RbtInt nMax(0);
  for (RbtBondMapConstIter iter = m_bondMap.begin(); iter != m_bondMap.end(); iter++)
    nMax = std::max(nMax,(*iter).first->GetFormalBondOrder());
  return nMax;
}

//Clears the bond map - should only need to be called by the copy constructors, hence private
void RbtAtom::ClearBondMap()
{
  m_bondMap.clear();
}


///////////////////////////////////////////////
//Other public methods
///////////////////////////////////////////////

//DM 08 Feb 1999 - all saved coords are now saved in a map<RbtUInt,RbtCoord>
//map key=0 is reserved for the default SaveCoords and RevertCoords
void RbtAtom::SaveCoords(RbtUInt coordNum)
{
  m_savedCoords[coordNum] = m_coord;
}

void RbtAtom::RevertCoords(RbtUInt coordNum) throw (RbtError)
{
  RbtUIntCoordMapConstIter iter = m_savedCoords.find(coordNum);
  if (iter != m_savedCoords.end()) {
    m_coord = (*iter).second;
  }
  else
    throw RbtInvalidRequest(_WHERE_,"RevertCoords failed on atom " + GetFullAtomName());
}



//DM 04 Dec 1998  Now we have the bond map, we can easily provide coordination numbers
//This version returns the total number of coordinated atoms (includes implicit hydrogens)
//Equivalent to GetNumBonds()+GetNumImplicitHydrogens()
RbtUInt RbtAtom::GetCoordinationNumber() const
{
  return m_bondMap.size() + m_nHydrogens;
}

//This version returns the number of coordinated atoms of a given element
//Note: for hydrogens, the implicit atoms are excluded
//so GetCoordinationNumber(1) will return the number of explicit hydrogens
RbtUInt RbtAtom::GetCoordinationNumber(RbtInt nAtomicNo) const
{
  return Rbt::GetNumAtoms(Rbt::GetBondedAtomList(m_bondMap),Rbt::isAtomicNo_eq(nAtomicNo));
}

//This version returns the number of coordinated atoms of a given force field type
RbtUInt RbtAtom::GetCoordinationNumber(const RbtString& strFFType) const
{
  return Rbt::GetNumAtoms(Rbt::GetBondedAtomList(m_bondMap),Rbt::isFFType_eq(strFFType));
}

//This version returns the number of coordinated atoms of a given hybridisation state
RbtUInt RbtAtom::GetCoordinationNumber(eHybridState e) const
{
  return Rbt::GetNumAtoms(Rbt::GetBondedAtomList(m_bondMap),Rbt::isHybridState_eq(e));
}

//////////////////////////////////////////
// Non-member functions (in Rbt namespace)
//////////////////////////////////////////

//DM 23 Apr 1999 - provide custom comparison function
//to bond map, so that map is sorted by bond ID, not by pointer
RbtBool Rbt::RbtBondPCmp_BondId::operator()(RbtBond* pBond1, RbtBond* pBond2) const {
  return pBond1->GetBondId() < pBond2->GetBondId();
}

////////////////////////////////////////////
//Conversion functions
///////////////////////

//Converts hybridisation state enum to a string
RbtString Rbt::ConvertHybridStateToString(RbtAtom::eHybridState eState)
{
  switch (eState) {
  case RbtAtom::SP:
    return "SP";
  case RbtAtom::SP2:
    return "SP2";
  case RbtAtom::SP3:
    return "SP3";
  case RbtAtom::AROM:
    return "AROM";
  case RbtAtom::TRI:
    return "TRI";
  default:
    return "";
  }
}

//Convert formal charge to a string (e.g. +, -, ++, --, +3, -3 etc)
RbtString Rbt::ConvertFormalChargeToString(RbtInt nCharge)
{
  switch (nCharge) {
  case 0:
    return "";
  case 1:
    return "+";
  case -1:
    return "-";
  case 2:
    return "++";
  case -2:
    return "--";
  default:
    //For higher charges, return as +3,-3 etc
    ostrstream ostr;
    ostr.setf(ios_base::showpos);
    ostr << nCharge << ends;
    RbtString strCharge(ostr.str());
    delete ostr.str();
    return strCharge;
  }
}


////////////////////////////////////////////
//Predicates (for use by STL algorithms)
////////////////////////////////////////////

//Is atom a H-Bond Acceptor ?
//Checks atomic number for O(8), N(7) or S(16)
//DM 15 Feb 1999 - screens out positively charged atoms
//Checks nitrogen hybridisation state and coordination number for available lone pairs
RbtBool Rbt::isAtomHBondAcceptor::operator() (const RbtAtom* pAtom) const
{
  //Screen out positively charged atoms
  //if (spAtom->GetFormalCharge() > 0.0)
  //  return false;

  const RbtInt nAtomicNo = pAtom->GetAtomicNo();
  const RbtAtom::eHybridState hybrid = pAtom->GetHybridState();
  const RbtUInt nCoord = pAtom->GetCoordinationNumber();
  switch (nAtomicNo) {
  case 7: //N
    //Only include nitrogens which have available lone pairs
    switch (hybrid) {
    case RbtAtom::SP3:
      return (nCoord < 4);
    case RbtAtom::SP2:
    case RbtAtom::TRI:
    case RbtAtom::AROM:
      return (nCoord < 3);
    case RbtAtom::SP:
      return (nCoord < 2);
    default:
      return false;
    }
  case 8:
    switch (hybrid) {
    case RbtAtom::SP3:
      return (nCoord < 3);
    case RbtAtom::SP2:
    case RbtAtom::TRI:
    case RbtAtom::AROM:
      //Remove terminal oxygens bonded to nitrogen (e.g. Nitro oxygens)
      return ((nCoord > 1) || (pAtom->GetCoordinationNumber(7) == 0));
    default:
      return false;
    }
  case 16:
    //16 May 2003 (DM) Sulphur - only include terminal S_SP2, e.g. in C=S
    //Previously, all S were included
    switch (hybrid) {
    case RbtAtom::SP2:
      return (nCoord < 2);
    default:
      return false;
    }
  default:
    return false;
  }
}

//Is atom a H-Bond Donor ?
//Checks 1) is it a hydrogen, 2) does it make exactly one bond, 3) is the bond to O, N or S ?
RbtBool Rbt::isAtomHBondDonor::operator() (const RbtAtom* pAtom) const
{
  return ( (pAtom->GetAtomicNo()==1) && (pAtom->GetNumBonds()==1) &&
	   ( (pAtom->GetCoordinationNumber(8)==1) ||
	     (pAtom->GetCoordinationNumber(7)==1) ||
	     (pAtom->GetCoordinationNumber(16)==1) ) );
}

//Is atom planar ?
//Checks if 1) atom makes 2 bonds (in which case must be planar) or
//2) hybridisation state is SP2, AROM or TRI
RbtBool Rbt::isAtomPlanar::operator() (const RbtAtom* pAtom) const
{
  const RbtAtom::eHybridState hybrid = pAtom->GetHybridState();
  //17 Dec 1998 (DM) GetCoordinationNumber works for extended atoms
  // as it includes any implicit hydrogens
  //  return ( (spAtom->GetNumBonds() == 2) ||
  return ( (pAtom->GetCoordinationNumber() == 2) ||
	   (hybrid == RbtAtom::SP2) ||
	   (hybrid == RbtAtom::AROM) ||
	   (hybrid == RbtAtom::TRI) );
}

//Is atom a pi-atom ?
//SP2,TRI or AROM
RbtBool Rbt::isPiAtom::operator() (const RbtAtom* pAtom) const
{
  const RbtAtom::eHybridState hybrid = pAtom->GetHybridState();
  return ( (hybrid == RbtAtom::SP2) ||
	   (hybrid == RbtAtom::AROM) ||
	   (hybrid == RbtAtom::TRI) );
}

//DM 29 Jul 1999
//Is atom within a given distance of any coord in the coord list
RbtBool Rbt::isAtomNearCoordList::operator() (const RbtAtom* pAtom) const
{
  //DM 27 Oct 2000 - GetCoords now returns by reference
  const RbtCoord& atomCoord = pAtom->GetCoords();
  for (RbtCoordListConstIter iter = cl.begin(); iter != cl.end(); iter++) {
    if (Rbt::Length2(atomCoord-(*iter)) <= r2)
      return true;
  }
  return false;
}


//Is atom2 1-2 connected (i.e. bonded) to atom1 ?
Rbt::isAtom_12Connected::isAtom_12Connected(RbtAtom* pAtom) : pAtom1(pAtom)
{
  bondedAtomList1 = GetBondedAtomList(pAtom1);
}

RbtBool Rbt::isAtom_12Connected::operator() (RbtAtom* pAtom2) const
{
  //Check if atom1 and atom2 are the same
  if (pAtom1 == pAtom2)
    return false;

  //Checks if atom2 appears in the bonded atom list for atom1.
  //STL find_if will do the trick
  //Need to use bind2nd to convert binary predicate to unary
  return std::find_if(bondedAtomList1.begin(),bondedAtomList1.end(),
		      std::bind2nd(Rbt::isAtomPtr_eq(),pAtom2)) != bondedAtomList1.end(); 
}

//Is atom2 1-3 connected (i.e. via a bond angle) to atom1 ?
Rbt::isAtom_13Connected::isAtom_13Connected(RbtAtom* pAtom) : pAtom1(pAtom)
{
  bondedAtomList1 = GetBondedAtomList(pAtom1);
}

RbtBool Rbt::isAtom_13Connected::operator() (RbtAtom* pAtom2) const
{
  //Check if atom1 and atom2 are the same
  if (pAtom1 == pAtom2)
    return false;

  RbtAtomList bondedAtomList2 = GetBondedAtomList(pAtom2);
  //Check if any atoms appear in both bonded atom lists
  //STL find_first_of will do the trick.
  return std::find_first_of(bondedAtomList1.begin(),bondedAtomList1.end(),
			    bondedAtomList2.begin(),bondedAtomList2.end(),
			    Rbt::isAtomPtr_eq()) !=  bondedAtomList1.end();
}


//Is atom's subunit an RNA-type (i.e. A,G,C or U)
//DM 4 Jan 1999 - check for 3 letter residue names also
//Note: strictly speaking these are used for DNA, not RNA, but they often crop up in RNA
RbtBool Rbt::isAtomRNA::operator() (const RbtAtom* pAtom) const
{
  RbtString strSubunitName = pAtom->GetSubunitName();
  return ( (strSubunitName == "A") || (strSubunitName == "ADE") ||
	   (strSubunitName == "G") || (strSubunitName == "GUA") ||
	   (strSubunitName == "C") || (strSubunitName == "CYT") ||
	   (strSubunitName == "U") || (strSubunitName == "URI") );
}


//DM 21 Jul 1999 Is atom lipophilic ?
//DM 16 May 2003 Total rewrite to be much more comprehensive
RbtBool Rbt::isAtomLipophilic::operator() (RbtAtom* pAtom) const
{
  //Remove the obvious polar stuff
  if ( isIonic(pAtom) || isHBD(pAtom) || isHBA(pAtom) || isMetal(pAtom) || isO(pAtom) || isN(pAtom)) return false;

  const RbtInt nAtomicNo = pAtom->GetAtomicNo();
  RbtAtomList bondedAtomList;
  RbtAtomListConstIter iter;

  switch (nAtomicNo) {
  case 1://H
    return true;
  case 6://C
  case 16://S
    if (isSP3(pAtom))
      return true;
    bondedAtomList = Rbt::GetBondedAtomList(pAtom);
    //Remove C,S bonded to O_SP2 (gets rid of carbonyl, sulphone, sulphoxide, sulphonamide etc)
    for (iter = bondedAtomList.begin(); iter != bondedAtomList.end(); iter++) {
      //Check this way so that predicate works for PSF/CRD files (formal bond order not set)
      if (isO(*iter) && isSP2(*iter))
	return false;
    }
    return true;
  case 9://F
    return false;
  case 17://Cl
  case 35://Br
  case 53://I
    return true;
  default:
    return false;
  }
}

//DM 24 Jan 2001
//Checks for common metal ions by atomic number (Na,Mg,K->Zn)
Rbt::isAtomMetal::isAtomMetal()  {}
RbtBool Rbt::isAtomMetal::operator() (const RbtAtom* pAtom) const {
  RbtInt atNo = pAtom->GetAtomicNo();
  return (atNo == 11) || (atNo == 12) || ((atNo >= 19) && (atNo <= 30));
}

//Select the intramolecular flexible interactions to this atom (sets SelectionFlag = true)
void Rbt::SelectFlexAtoms::operator() (RbtAtom* pAtom) {
  RbtModel* pModel = pAtom->GetModelPtr();
  if (pModel) {
    pModel->SelectFlexAtoms(pAtom);
  }
}

////////////////////////////////////////////
//Atom list functions (implemented as STL algorithms)
////////////////////////////////////////////

void Rbt::SetAtomSelectionFlags(RbtAtomList& atomList,RbtBool bSelected) {
  Rbt::SelectAtom select(bSelected);
  std::for_each(atomList.begin(),atomList.end(),select);
}

void Rbt::InvertAtomSelectionFlags(RbtAtomList& atomList) {
  Rbt::InvertSelectAtom invert;
  std::for_each(atomList.begin(),atomList.end(),invert);
}

void Rbt::SetAtomCyclicFlags(RbtAtomList& atomList,RbtBool bCyclic) {
  Rbt::CyclicAtom cyclic(bCyclic);
  std::for_each(atomList.begin(),atomList.end(),cyclic);
}

//BondedAtoms
//Helper function to return the atom pointer for the "other" atom in the bond
//i.e. return (bondPair.second) ? bondPair.first->GetAtom2Ptr() : bondPair.first->GetAtom1Ptr();
RbtAtomPtr Rbt::GetBondedAtomPtr(pair<RbtBond*,RbtBool> bondBoolPair)
{
  //If this atom is first in the bond (bondPair.second == true), get atom 2 ptr, else get atom 1 ptr 
  return (bondBoolPair.second) ? (bondBoolPair.first)->GetAtom2Ptr() : (bondBoolPair.first)->GetAtom1Ptr();
}

RbtUInt Rbt::GetNumBondedAtoms(const RbtBondMap& bondMap)
{
  return bondMap.size();//This function is not really needed, but provide for consistency
}

RbtAtomList Rbt::GetBondedAtomList(const RbtBondMap& bondMap)
{
  RbtAtomList atomList;
  std::transform(bondMap.begin(),bondMap.end(),std::back_inserter(atomList),std::ptr_fun(&Rbt::GetBondedAtomPtr));
  //DM 23 Apr 1999 - sort the bonded atom list by atom ID so we return the atoms in the same order every time
  //The problem with the bond map is that the key is RbtBond* which is not guaranteed to increase in order
  //of construction.
  std::sort(atomList.begin(),atomList.end(),RbtAtomPtrCmp_AtomId());
  return atomList;
}

//If spAtom is a regular RbtAtom, these two functions just behave like the two above
//If spAtom can be dynamically_casted to an RbtPseudoAtom, these return the constituent atom list
//for the pseudoatom
RbtUInt Rbt::GetNumBondedAtoms(const RbtAtom* pAtom)
{
  //Attempt the dynamic cast (on the regular pointer)
  RbtPseudoAtom* pPseudoAtom = dynamic_cast<RbtPseudoAtom*>(const_cast<RbtAtom*>(pAtom));
  if (pPseudoAtom != NULL)
    return pPseudoAtom->GetNumAtoms();
  else
    return Rbt::GetNumBondedAtoms(pAtom->GetBondMap());
}

RbtAtomList Rbt::GetBondedAtomList(const RbtAtom* pAtom)
{
  //Attempt the dynamic cast (on the regular pointer)
  RbtPseudoAtom* pPseudoAtom = dynamic_cast<RbtPseudoAtom*>(const_cast<RbtAtom*>(pAtom));
  if (pPseudoAtom != NULL)
    return pPseudoAtom->GetAtomList();
  else
    return Rbt::GetBondedAtomList(pAtom->GetBondMap());
}

//Binary



//DM 1 Feb 1999
//Attempts to match atoms in atomList2 with those in atomList1, where the match is performed
//using isAtom_eq (tests subunit ID, subunit name and atom name for equality) and NOT isAtomPtr_eq
//(which tests the RbtAtom* pointers themselves for equality, i.e. same objects)
//Returns list of RbtAtom smart pointers to atoms in atomList1 for which a match is found.
RbtUInt Rbt::GetNumMatchingAtoms(const RbtAtomList& atomList1, const RbtAtomList& atomList2)
{
  RbtUInt nCount(0);
  for (RbtAtomListConstIter iter = atomList2.begin(); iter != atomList2.end(); iter++) {
    nCount += GetNumAtoms(atomList1,std::bind2nd(Rbt::isAtom_eq(),*iter));
  }
  return nCount;
}

RbtAtomList Rbt::GetMatchingAtomList(const RbtAtomList& atomList1, const RbtAtomList& atomList2)
{
  RbtAtomList matchingAtomList;
  for (RbtAtomListConstIter iter = atomList2.begin(); iter != atomList2.end(); iter++) {
    RbtAtomList tempList = GetAtomList(atomList1,std::bind2nd(Rbt::isAtom_eq(),*iter));
    std::copy(tempList.begin(),tempList.end(),std::back_inserter(matchingAtomList));
  }
  return matchingAtomList;
}  

//DM 30 Mar 1999
//These versions match atoms against a full atom name specifier of the form given by GetFullAtomName():
//<segment name>:<subunit name>_<subunit id>:<atom name>  e.g. A:U_23:O4
//All fields are optional, and if missing, will match all atoms
//So, for example, these are all valid full names:
//N7     (equiv. to ::N7) - matches N7's in all subunits, all segments
//A::O4                   - matches O4's in all subunits in segment A
//U:O4 (equiv. to :U:O4)  - matches O4's in subunits named U in all segments
//U_23:O4                 - matches O4's in subunit U23 in all segments
//_23:                    - matches all atoms in subunit ID=23 in all segments
RbtUInt Rbt::GetNumMatchingAtoms(const RbtAtomList& atomList, const RbtString& strFullName)
{
  return Rbt::GetMatchingAtomList(atomList,strFullName).size();
}

RbtAtomList Rbt::GetMatchingAtomList(const RbtAtomList& atomList, const RbtString& strFullName)
{
  RbtAtomList matchingAtoms = atomList;

  //Split the name into it's constituent subunit and atom names
  RbtStringList componentList = Rbt::ConvertDelimitedStringToList(strFullName,":");
  RbtUInt idx(0);//Index into constituent list
  //Switch on how many constituent names are specified
  switch (componentList.size()) {      

  case 3://we have segment, subunit and atom names
    if (!componentList[idx].empty()) {
      matchingAtoms = Rbt::GetAtomList(matchingAtoms,Rbt::isSegmentName_eq(componentList[idx]));
#ifdef _DEBUG
      //cout << "Matching segment name=" << componentList[idx] << ", #atoms=" << matchingAtoms.size() << endl;
#endif //_DEBUG
    }
    idx++;

  case 2://Fall-through!! We have subunit and atom names
    //Check if the subunit name also specifies a subunit ID
    if (!componentList[idx].empty()) {
      RbtStringList subunitList =  Rbt::ConvertDelimitedStringToList(componentList[idx],"_");
      switch (subunitList.size()) {
      case 2://Subunit ID and name are specified
	if (!subunitList[1].empty()) {
	  matchingAtoms = Rbt::GetAtomList(matchingAtoms,Rbt::isSubunitId_eq(subunitList[1]));
#ifdef _DEBUG
	  //cout << "Matching subunit id=" << subunitList[1] << ", #atoms=" << matchingAtoms.size() << endl;
#endif //_DEBUG
	}
      case 1://Fall-through!! Only the subunit name is specified
	if (!subunitList[0].empty()) {
	  matchingAtoms = Rbt::GetAtomList(matchingAtoms,Rbt::isSubunitName_eq(subunitList[0]));
#ifdef _DEBUG
	  //cout << "Matching subunit name=" << subunitList[0] << ", #atoms=" << matchingAtoms.size() << endl;
#endif //_DEBUG
	}
	break;
      default:
#ifdef _DEBUG
	//cout << "Invalid subunit string in " << strFullName << endl;
#endif //_DEBUG
	break;
      }
    }
    idx++;
    
  case 1://Fall-through!! We only have an atom name specifier
    if (!componentList[idx].empty()) {
      matchingAtoms = Rbt::GetAtomList(matchingAtoms,Rbt::isAtomName_eq(componentList[idx]));
#ifdef _DEBUG
      //cout << "Matching atom name=" << componentList[idx] << ", #atoms=" << matchingAtoms.size() << endl;
#endif //_DEBUG
    }
    break;

  default:
#ifdef _DEBUG
    //cout << "Too many colons (:) in " << strFullName << endl;
#endif //_DEBUG
    break;
  }

  return matchingAtoms;
}

//DM 15 Apr 1999 - as above, but match against a list of full atom name specifiers
//Returns total list (i.e. all matches OR'd). Does not remove duplicates.
RbtUInt Rbt::GetNumMatchingAtoms(const RbtAtomList& atomList, const RbtStringList& fullNameList)
{
  return Rbt::GetMatchingAtomList(atomList,fullNameList).size();
}

RbtAtomList Rbt::GetMatchingAtomList(const RbtAtomList& atomList, const RbtStringList& fullNameList)
{
  RbtAtomList matchingAtoms;
  for (RbtStringListConstIter iter = fullNameList.begin(); iter != fullNameList.end(); iter++) {
    //Find the atoms which match the specifier string
    RbtAtomList selectedAtoms = Rbt::GetMatchingAtomList(atomList,*iter);
    //Add to the total list of matching atoms
    std::copy(selectedAtoms.begin(),selectedAtoms.end(),std::back_inserter(matchingAtoms));
  }
  return matchingAtoms;
}

/////////////////////////
//DM 07 Jan 1999
//Actions on atom lists
/////////////////////////

//Save coords by number for all atoms in the list
void Rbt::SaveAtomCoords(const RbtAtomList& atomList, RbtUInt coordNum)
{
  for (RbtAtomListConstIter iter = atomList.begin(); iter != atomList.end(); iter++) {
    RbtAtomPtr spAtom(*iter);
    spAtom->SaveCoords(coordNum);
  }
}

//Revert to numbered coords for all atoms in the list
void Rbt::RevertAtomCoords(const RbtAtomList& atomList, RbtUInt coordNum)
{
  for (RbtAtomListConstIter iter = atomList.begin(); iter != atomList.end(); iter++) {
    RbtAtomPtr spAtom(*iter);
    spAtom->RevertCoords(coordNum);
  }
}

//Returns total atomic mass (molecular weight) for all atoms in the list
RbtDouble Rbt::GetTotalAtomicMass(const RbtAtomList& atomList)
{
  return std::accumulate(atomList.begin(),atomList.end(),0.0,AccumAtomicMass);
}

//Returns center of mass of atoms in the list
RbtCoord Rbt::GetCenterOfMass(const RbtAtomList& atomList)
{
  //Accumulate sum of mass*coord
  RbtCoord com = std::accumulate(atomList.begin(),atomList.end(),RbtCoord(),AccumMassWeightedCoords);
  //Divide by total mass
  com /= Rbt::GetTotalAtomicMass(atomList);
  return com;
}

//DM 20 May 1999 - returns the coords of all atoms in the list
RbtCoordList Rbt::GetCoordList(const RbtAtomList& atomList)
{
  RbtCoordList coordList;
  std::transform(atomList.begin(),atomList.end(),std::back_inserter(coordList),ExtractAtomCoord);
  return coordList;
}

//DM 09 Aug 2001 - returns coordList via argument
void Rbt::GetCoordList(const RbtAtomList& atomList, RbtCoordList& coordList)
{
  coordList.clear();
  std::transform(atomList.begin(),atomList.end(),std::back_inserter(coordList),ExtractAtomCoord);
}

//Streams an atom list in Quanta CSD file format (for easy color coding of selected atoms in Quanta)
//nFormat  = 0 => Receptor atom format: "zone A:1 # pick O5T = col 2"
//nFormat != 0 => Ligand atom format "pick N1 .and. segm H = col 2"
//ostream should have been opened before calling QuantaCSDFileDump
//DM 30 Jul 1999 - added segment name to receptor atom format
void Rbt::PrintQuantaCSDFormat(const RbtAtomList& atomList, ostream& s, RbtInt nColor, RbtInt nFormat)
{
  for (RbtAtomListConstIter iter = atomList.begin(); iter != atomList.end(); iter++) {
    RbtAtomPtr spAtom(*iter);
    switch (nFormat) {
    case 0:
      s << "zone " << spAtom->GetSegmentName() << ":" << spAtom->GetSubunitId() << " # pick "
	<< spAtom->GetAtomName() << " = col " << nColor << endl;
      break;
    default:
      s << "pick " << spAtom->GetAtomName() << " .and. segm " << spAtom->GetSegmentName()
	<< " = col " << nColor << endl;
    }
  }
}

//DM 25 Feb 1999
//Modified DM 6 Apr 1999 to operate using the GroupCharge attribute
//Scans an atom list and sets the group charge to zero on any
//cationic-ionic pairs which are 1-2 or 1-3 connected
//(e.g. on -OOC.CR.NH3+ would zero the C and N)
void Rbt::RemoveZwitterions(RbtAtomList& atomList)
{
  RbtAtomList cationicAtomList = Rbt::GetAtomList(atomList, Rbt::isAtomCationic());
  RbtAtomList anionicAtomList = Rbt::GetAtomList(atomList, Rbt::isAtomAnionic());

  //Progressively partition the cationic atom list into those that are 1-3 connected or less
  //to one of the anionic atoms, and those that aren't
  //c14begin points to the first cationic atom which is 1-4 connected or above to all of the anionic atoms
  //All cationic atoms up to c14begin can be removed
  RbtAtomListIter c14begin = cationicAtomList.begin();
  for (RbtAtomListConstIter aIter = anionicAtomList.begin(); aIter != anionicAtomList.end(); aIter++) {
    //Need to partition in two stages: firstly 1-2 connected, then 1-3 connected
    c14begin = std::stable_partition(c14begin,
				     cationicAtomList.end(),
				     Rbt::isAtom_12Connected(*aIter));
    c14begin = std::stable_partition(c14begin,
				     cationicAtomList.end(),
				     Rbt::isAtom_13Connected(*aIter));
  }

  //Now repeat the other way round, checking the anionic atom list against ALL the cationic atoms
  //(including those we know are due to be removed)
  //We need to do it in this way to cope with situations like (O-)-(N+)-(O-)
  //We couldn't simply delete the first N,O pair we find as this would leave the second O
  //with nothing to match.
  //
  RbtAtomListIter a14begin = anionicAtomList.begin();
  for (RbtAtomListConstIter cIter = cationicAtomList.begin(); cIter != cationicAtomList.end(); cIter++) {
    //Need to partition in two stages: firstly 1-2 connected, then 1-3 connected
    a14begin = std::stable_partition(a14begin,
				     anionicAtomList.end(),
				     Rbt::isAtom_12Connected(*cIter));
    a14begin = std::stable_partition(a14begin,
				     anionicAtomList.end(),
				     Rbt::isAtom_13Connected(*cIter));
  }

  //Now we can safely neutralise the 1-2 and 1-3 atoms
  for (RbtAtomListIter iter = cationicAtomList.begin(); iter != c14begin; iter++) {
    (*iter)->SetGroupCharge(0.0);
#ifdef _DEBUG
    cout << "RemoveZwitterions: Neutralising cation " << (*iter)->GetFullAtomName() << endl;
#endif //_DEBUG
  }

  for (RbtAtomListIter iter = anionicAtomList.begin(); iter != a14begin; iter++) {
    (*iter)->SetGroupCharge(0.0);
#ifdef _DEBUG
    cout << "RemoveZwitterions: Neutralising anion " << (*iter)->GetFullAtomName() << endl;
#endif //_DEBUG
  }
}
