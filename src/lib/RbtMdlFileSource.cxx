/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtMdlFileSource.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include <iomanip>

#include "RbtMdlFileSource.h"
#include "RbtFileError.h"
#include "RbtModelError.h"
#include "RbtPlane.h"
#include "RbtAtomFuncs.h"

RbtMdlFileSource::RbtMdlFileSource(const RbtString& fileName, RbtBool bPosIonisable, RbtBool bNegIonisable, RbtBool bImplHydrogens) :
  RbtBaseMolecularFileSource(fileName,IDS_MDL_RECDELIM,"MDL_FILE_SOURCE"), //Call base class constructor
  m_bPosIonisable(bPosIonisable),m_bNegIonisable(bNegIonisable),m_bImplHydrogens(bImplHydrogens)
{
  //Open an Element data source
  m_spElementData = RbtElementFileSourcePtr(new RbtElementFileSource(Rbt::GetRbtFileName("data","RbtElements.dat")));
  _RBTOBJECTCOUNTER_CONSTR_("RbtMdlFileSource");
}


//Default destructor
RbtMdlFileSource::~RbtMdlFileSource()
{
  _RBTOBJECTCOUNTER_DESTR_("RbtMdlFileSource");
}

void RbtMdlFileSource::Parse() throw (RbtError)
{
  //Only parse if we haven't already done so
  if (!m_bParsedOK) {
    ClearMolCache();//Clear current cache
    Read();//Read the current record

    try {
      RbtFileRecListIter fileIter = m_lineRecs.begin();
      RbtFileRecListIter fileEnd = m_lineRecs.end();


      //////////////////////////////////////////////////////////
      //1a. Store title lines (first 3)..
      RbtInt nTitleRec = 3;
      m_titleList.reserve(nTitleRec);//Allocate enough memory for the vector
      while ( (m_titleList.size() < nTitleRec) && (fileIter != fileEnd)) {
	m_titleList.push_back(*fileIter++);
      }

      //1b ..and check we read them all before reaching the end of the file
      if (m_titleList.size() != nTitleRec)
	throw RbtFileParseError(_WHERE_,"Incomplete title records in " + GetFileName());

      //1c DM 20 Dec 1999 - throw model error if record is 2D instead of 3D
      //DM 08 Aug 2000 - check for "3" at position 21 as well as "3D" at position 20
      if ((m_titleList[1].find("3",21) != 21) && (m_titleList[1].find("3D",20) != 20)) {
	throw RbtModelError(_WHERE_,"2D record in " + GetFileName());
      }

      //////////////////////////////////////////////////////////
      //2. Read number of atoms and bonds
      RbtInt nAtomRec;
      RbtInt nBondRec;
      if (fileIter != fileEnd) {
	//The SD file format only uses a field width of 3 to store nAtoms, nBonds
	//so for values over 99 the two fields coalesce.
	//Workaround is to insert a space between the two fields (or use sscanf)
	if ( (*fileIter).size() > 3)
	  (*fileIter).insert(3," ");
	istrstream istr((*fileIter++).c_str());
	istr >> nAtomRec >> nBondRec;
#ifdef _DEBUG
	//cout << nAtomRec << " atoms, " << nBondRec << " bonds" << endl;
#endif //_DEBUG
      }
      else
	throw RbtFileParseError(_WHERE_,"Missing atom and bond information in " + GetFileName());

      //////////////////////////////////////////////////////////
      //3a. Create and store atoms
      m_atomList.reserve(nAtomRec);//Allocate enough memory for the vector

      RbtInt nAtomId(0);
      RbtCoord coord;//X,Y,Z coords
      RbtString strElementName; //element name
      RbtInt nMassDiff; //mass difference
      RbtInt nFormalCharge; //formal charge (MDL stores in a funny way, see below)
      RbtString strSegmentName("H"); //constant (i.e. one segment, one residue)
      RbtString strSubunitId("1"); //constant
      RbtString strSubunitName("MOL"); //constant

      //DM 24 Mar 1999 - allocate a C-string for streaming the atom name into
      RbtInt lenAtomName(10);
      char* szAtomName(new char[lenAtomName+1]);

      while ( (m_atomList.size() < nAtomRec) && (fileIter != fileEnd)) {
	istrstream istr((*fileIter++).c_str());
	istr >> coord.x
	     >> coord.y
	     >> coord.z
	     >> strElementName
	     >> nMassDiff
	     >> nFormalCharge;

	//Look up the element data
	RbtElementData elData = m_spElementData->GetElementData(strElementName);
	RbtInt nAtomicNo = elData.atomicNo;

	//Correct the formal charge (for non-zero charges, MDL stores as 4-charge)
	//i.e. -3 is stored as 7, -2 as 6, -1 as 5, 0 as 0 or 4
	//+1 as 3, +2 as 2, +3 as 1
	if (nFormalCharge > 0)
	  nFormalCharge = 4-nFormalCharge;

	//Compose the atom name from element+atomID (i.e. C1, N2, C3 etc)
	nAtomId++;
	ostrstream ostr(szAtomName,lenAtomName);
	ostr << strElementName << nAtomId << ends;
	RbtString strAtomName(szAtomName);

	//Construct a new atom (constructor only accepts the 2D params)
	RbtAtomPtr spAtom(new RbtAtom(nAtomId,
				      nAtomicNo,
				      strAtomName,
				      strSubunitId,
				      strSubunitName,
				      strSegmentName,
				      RbtAtom::UNDEFINED,//Hybridisation state
				      0,//Num implicit hydrogens
				      nFormalCharge) //Formal charge
			  );

	//Now set the remaining 2-D and 3-D params
	spAtom->SetCoords(coord);
	spAtom->SetPartialCharge(0.0);
	spAtom->SetGroupCharge(nFormalCharge);//Default is group charge=formal charge (gets corrected later)
	spAtom->SetAtomicMass(elData.mass+nMassDiff);

	m_atomList.push_back(spAtom);
	m_segmentMap[strSegmentName]++;//increment atom count in segment map
      }

      delete [] szAtomName;//Tidy up the C-string

      //3b ..and check we read them all before reaching the end of the file
      if (m_atomList.size() != nAtomRec)
	throw RbtFileParseError(_WHERE_,"Incomplete atom records in " + GetFileName());


      //////////////////////////////////////////////////////////
      //4a. Create and store bonds
      m_bondList.reserve(nBondRec);//Allocate enough memory for the vector

      RbtInt nBondId(0);
      RbtUInt idxAtom1;
      RbtUInt idxAtom2;
      RbtInt nBondOrder;

      while ( (m_bondList.size() < nBondRec) && (fileIter != fileEnd)) {
	//The SD file format only uses a field width of 3 to store atom1,atom2
	//so for values over 99 the two fields coalesce.
	//Workaround is to insert a space between the two fields (or use sscanf)
	if ( (*fileIter).size() > 3)
	  (*fileIter).insert(3," ");
	istrstream istr((*fileIter++).c_str());
	istr >> idxAtom1
	     >> idxAtom2
	     >> nBondOrder;
	if ( (idxAtom1 > nAtomRec) || (idxAtom2 > nAtomRec) ) {//Check for indices in range
	  throw RbtFileParseError(_WHERE_,"Atom index out of range in bond records in " + GetFileName());
	}
	RbtAtomPtr spAtom1(m_atomList[idxAtom1-1]);//Decrement the atom index as the atoms are numbered from zero in our atom vector
	RbtAtomPtr spAtom2(m_atomList[idxAtom2-1]);//Decrement the atom index as the atoms are numbered from zero in our atom vector
	RbtBondPtr spBond(new RbtBond(++nBondId,//Store a nominal bond ID starting from 1
				      spAtom1,
				      spAtom2,
				      nBondOrder
				      )
			  );
	m_bondList.push_back(spBond);
      }

      //4b ..and check we read them all before reaching the end of the file
      if (m_bondList.size() != nBondRec)
	throw RbtFileParseError(_WHERE_,"Incomplete bond records in " + GetFileName());

      //DM 12 May 1999 - read data records (if any)
      for ( ; fileIter != fileEnd; fileIter++) {
	if ((*fileIter).find(">") == 0) { //Found a data record
	  RbtString::size_type ob = (*fileIter).find("<");//First open bracket
	  RbtString::size_type cb = (*fileIter).rfind(">");//Last closed bracket
	  if ( (ob != RbtString::npos) && (cb != RbtString::npos) ) {
	    RbtString fieldName = (*fileIter).substr(ob+1,cb-ob-1);//Data field name
	    RbtStringList sl;//String list for storing data value
	    while ( (++fileIter != fileEnd) && !(*fileIter).empty()) {
	      sl.push_back(*fileIter);
	    }
	    m_dataMap[fieldName] = RbtVariant(sl);
	  }
	}
      }

      //DM 18 May 1999 - check for <Name> record
      //If not present, define <Name> as the first title record
      if ( (m_dataMap.find("Name") == m_dataMap.end()) &&
	   !m_titleList.empty() ) {
	m_dataMap["Name"] = RbtVariant(m_titleList.front());
      }

      //Setup the atomic params not stored in the file (e.g. hybridisation state etc)
      SetupAtomParams();

      //DM 4 June 1999 - set up distinct segment names for each molecular fragment in the record
      SetupSegmentNames();

      //////////////////////////////////////////////////////////
      //If we get this far everything is OK
      m_bParsedOK = true;
    }

    catch (RbtError& error) {
      ClearMolCache();//Clear the molecular cache so we don't return incomplete atom and bond lists
      throw;//Rethrow the RbtError
    }

  }
}

//Sets up all the atomic attributes that are not explicitly stored in the PSF file
void RbtMdlFileSource::SetupAtomParams() throw (RbtError)
{
  //In the SD/MOL file we have:
  //  Element names
  //  X,Y,Z coords
  //  Formal bond orders
  //  Formal charges
  //
  //We need to determine:
  //  Hybridisation state - sp3,sp2,sp,trigonal planar
  //  Valency of each atom - is it in range? - use to assign implicit hydrogens if necessary
  //  Corrected vdW radii on H-bonding hydrogens and extended atoms
  //  Simple "forcefield" type - e.g. Osp3-
  //  Interaction group charges

  //DM 2 Aug 1999 - remove all non-polar hydrogens if required
  if (m_bImplHydrogens) {
    RemoveNonPolarHydrogens();
    RenumberAtomsAndBonds();//Ensure atom and bond IDs are consecutive prior to autoprotonation
  }

  SetupHybridState();
  SetupTheRest();
  SetCyclicFlags();//DM 8 Feb 2000 - needed for SetupIonicGroups()
  if (m_bPosIonisable)
    SetupPosIonisableGroups();
  if (m_bNegIonisable)
    SetupNegIonisableGroups();
  SetupIonicGroups();
  Rbt::RemoveZwitterions(m_atomList);//DM 6 Apr 1999 - might as well neutralise any zwitterions here
  RenumberAtomsAndBonds();//DM 30 Oct 2000 - clean up atom and bond list numbering
}

//Helper functions for SetupAtomParams()

//Sets hybridisation state and checks for Nsp3,Osp3 or Ssp3 next to sp2 special cases
void RbtMdlFileSource::SetupHybridState() throw (RbtError)
{
  //////////////////
  //FIRST PASS
  //////////////////
  for (RbtAtomListIter iter = m_atomList.begin(); iter != m_atomList.end(); iter++) {
    //Determine hybridisation state from the maximum bond order to this atom
    //i.e. if it makes a triple bond it is SP, double bond=SP2, single bond=SP3,
    //no bonds=UNDEFINED
    switch((*iter)->GetMaxFormalBondOrder()) {
    case 3:
      (*iter)->SetHybridState(RbtAtom::SP);
      break;
    case 2:
      (*iter)->SetHybridState(RbtAtom::SP2);
      break;
    case 1:
      (*iter)->SetHybridState(RbtAtom::SP3);
      break;
    case 0:
      (*iter)->SetHybridState(RbtAtom::UNDEFINED);
      break;
    default:
      //throw RbtModelError(_WHERE_, (*iter)->GetAtomName()+" makes a bond with an unknown bond order");
      cout << RbtModelError(_WHERE_, (*iter)->GetAtomName()+" makes a bond with an unknown bond order") << endl;
      (*iter)->SetHybridState(RbtAtom::UNDEFINED);
      break;
    }
  }
  //////////////////
  //SECOND PASS - check for Nsp3, Osp3 and Ssp3 next to a pi atom. Convert to trigonal planar
  //For Nsp3, also check the actual planarity of the atom, based on the distance from the nitrogen
  //to the plane of the substituent. If non-planar, then leave as Nsp3.
  //Needed for amides, pyrroles, furans and thiophenes.
  //Trigonal planar atoms are sp2 hybridised with two electrons in a pi-orbital, and only
  //make formal single bonds.
  //(as opposed to SP2 atoms which make a formal double bond with only one electron in the pi-orbital)
  //////////////////
  //
  //Get a list of all the sp3 atoms which are not positively charged
  RbtAtomList sp3AtomList = Rbt::GetAtomList(m_atomList,Rbt::isHybridState_eq(RbtAtom::SP3));
  sp3AtomList = Rbt::GetAtomList(sp3AtomList,std::not1(Rbt::isAtomPosCharged()));

  for (RbtAtomListIter iter = sp3AtomList.begin(); iter != sp3AtomList.end(); iter++) {
    RbtAtomList bondedAtomList;
    switch ((*iter)->GetAtomicNo()) {
    case 7://N
    case 8://O
    case 16://S
      bondedAtomList = Rbt::GetBondedAtomList(*iter);
      //Is is bonded to a pi-atom ?
      if (Rbt::FindAtom(bondedAtomList,Rbt::isPiAtom()) != bondedAtomList.end()) {
	//Determine distance from atom to plane of bonded atoms (only nitrogen will have 3 substituents)
	if (bondedAtomList.size() == 3) {
	  const RbtCoord& c0 = (*iter)->GetCoords();
	  RbtVector v1 = bondedAtomList[0]->GetCoords()-c0;//Vector to subst 1
	  RbtVector v2 = bondedAtomList[1]->GetCoords()-c0;//Vector to subst 2
	  RbtVector v3 = bondedAtomList[2]->GetCoords()-c0;//Vector to subst 3
	  //Take unit vectors along each bond vector to eliminate the effect of different bond lengths
	  RbtPlane p = RbtPlane(c0+v1.Unit(),c0+v2.Unit(),c0+v3.Unit());//Plane of three substs
	  //This is the signed distance from the atom to the plane
	  RbtDouble dist = fabs(Rbt::DistanceFromPointToPlane(c0,p));
#ifdef _DEBUG
	  cout << "Distance from " << (*iter)->GetAtomName() << " to plane of substituents=" << dist << endl;
#endif //_DEBUG
	  if (dist > 0.1) {//Not planar
#ifdef _DEBUG
	    cout << "Not planar" << endl;
#endif //_DEBUG
	    break;
	  }
	}
	(*iter)->SetHybridState(RbtAtom::TRI);
#ifdef _DEBUG
	cout << "Changing " << (*iter)->GetAtomName() << " from SP3 to TRI" << endl;
#endif //_DEBUG
      }
      break;
    default:
      break;
    }
  }
}


//Checks for valency within range
//Adds implicit hydrogens to C,N,S
//Defines vdW radius, correcting for extended atoms and H-bond donor hydrogens
//Defines formal "force-field" type string
void RbtMdlFileSource::SetupTheRest() throw (RbtError)
{
  //Radius increment for atoms with implicit hydrogens
  //DM 22 Jul 1999 - only increase the radius for sp3 atoms with implicit hydrogens
  //For sp2 and aromatic, leave as is
  Rbt::isHybridState_eq bIsSP3(RbtAtom::SP3);
  RbtDouble dImplRadIncr = m_spElementData->GetImplicitRadiusIncr();
  //Radius increment and predicate for H-bonding hydrogens
  RbtDouble dHBondRadIncr = m_spElementData->GetHBondRadiusIncr();
  Rbt::isAtomHBondDonor bIsHBondDonor;
  //Element data for hydrogen
  RbtElementData elHData = m_spElementData->GetElementData(1);

  for (RbtAtomListIter iter = m_atomList.begin(); iter != m_atomList.end(); iter++) {
    //Get the element data for this atom
    RbtInt nAtomicNo = (*iter)->GetAtomicNo();
    RbtElementData elData = m_spElementData->GetElementData(nAtomicNo);
    RbtDouble vdwRadius = elData.vdwRadius;
    //Determine valency
    RbtInt nFormalCharge = (*iter)->GetFormalCharge();
    RbtInt nTotalFormalBondOrder = (*iter)->GetTotalFormalBondOrder();
    RbtInt nValency = nTotalFormalBondOrder-nFormalCharge;
    //Is it within range ?
    RbtInt nMinVal = elData.minVal;
    RbtInt nMaxVal = elData.maxVal;
    //If it is too low, then add implicit hydrogens but only for C,N,or S
    if (nValency < nMinVal) {
      switch (nAtomicNo) {
      case 6:
	//case 7:
      case 16:
	(*iter)->SetNumImplicitHydrogens(nMinVal-nValency);
	(*iter)->SetAtomicMass(elData.mass+(nMinVal-nValency)*elHData.mass);//Adjust atomic mass
	if (bIsSP3(*iter)) {
	  vdwRadius += dImplRadIncr;//adjust vdw radius (for sp3 implicit atoms only)
#ifdef _DEBUG
	  cout << "Adding " << nMinVal-nValency << " implicit hydrogens to " << (*iter)->GetAtomName() << endl;
#endif //_DEBUG
	}
	break;
      default:
	//DM 26 Jul 1999 - don't complain about valency errors
	//throw RbtModelError(_WHERE_, (*iter)->GetAtomName()+" makes too few bonds");
	cout << "**WARNING** " << RbtModelError(_WHERE_, (*iter)->GetAtomName()+" makes too few bonds") << endl;
	break;
      }
    }
    //DM 26 Jul 1999 - don't complain about valency errors
    else if (nValency > nMaxVal) {
      //throw RbtModelError(_WHERE_, (*iter)->GetAtomName()+" makes too many bonds");
      cout << "**WARNING** " << RbtModelError(_WHERE_, (*iter)->GetAtomName()+" makes too many bonds") << endl;
    }

    //Set the vdw radius, first correcting for H-bonding hydrogens
    if (bIsHBondDonor(*iter))
      vdwRadius += dHBondRadIncr;
    (*iter)->SetVdwRadius(vdwRadius);

    //Now set a nominal "force-field" type of the style OSP3- etc
    RbtString strFFType = elData.element;
    strFFType += "_";
    strFFType += Rbt::ConvertHybridStateToString((*iter)->GetHybridState());
    strFFType += Rbt::ConvertFormalChargeToString((*iter)->GetFormalCharge());
    (*iter)->SetFFType(strFFType);
  }
}

//Protonates all amines, guanidines and imidazoles
void RbtMdlFileSource::SetupPosIonisableGroups()
{
  //Useful predicates
  Rbt::isFFType_eq bIsN_SP3("N_SP3");
  Rbt::isFFType_eq bIsN_SP2("N_SP2");
  Rbt::isFFType_eq bIsN_TRI("N_TRI");
  Rbt::isFFType_eq bIsC_SP2("C_SP2");

  //1. Amines
  //Protonate all Nsp3 atoms
  //DM 14 May 1999 - unless adjacent to a pi-atom
  //Most will be N_TRI so won't be protonated anyway, but some will be non-planar
  //so are still N_SP3. However we don't want to protonate these guys as the pKa
  //for e.g. Ph-NH3+ is much lower than for R(aliphatic)-NH3+.
  //This will also prevent e.g. SO2NHR getting protonated.
  RbtAtomList nsp3AtomList = Rbt::GetAtomList(m_atomList,bIsN_SP3);
  //DM 04 Aug 1999 - only consider Nsp3 atoms which are bonded to at least one Csp3
  //to avoid protonating hydrazines for e.g. This is the definition of an amine in any case
  nsp3AtomList = Rbt::GetAtomList(nsp3AtomList,std::not1(Rbt::isCoordinationNumber_eq(0,"C_SP3")));

  for (RbtAtomListIter iter = nsp3AtomList.begin(); iter != nsp3AtomList.end(); iter++) {
    RbtAtomList bondedAtomList = Rbt::GetBondedAtomList(*iter);
    //Only protonate if not bonded to pi-atom ?
    if (Rbt::FindAtom(bondedAtomList,Rbt::isPiAtom()) == bondedAtomList.end()) {
      //#ifdef _DEBUG
      cout << "Protonating neutral amine " << (*iter)->GetAtomName() << endl;
      //#endif //_DEBUG
      AddHydrogen(*iter);
    }
#ifdef _DEBUG
    else {
      cout << "Amine adjacent to pi-atom, not protonated: " << (*iter)->GetAtomName() << endl;
    }
#endif //_DEBUG
  }
  
  //2. Guanidine, imidazole, amidine
  //Guanidine
  // N_TRI    N_TRI
  //     `    '
  //      C_SP2
  //       ||
  //       N_SP2
  // - protonate the N_SP2
  //
  //Imidazole
  //    C_SP2-N_SP2
  //    //     `
  //   C_SP2   C_SP2
  //     `     '
  //      N_TRI
  //        |
  //        H
  // - protonate the N_SP2
  //

  //Loop over all C_SP2
  RbtAtomList csp2AtomList = Rbt::GetAtomList(m_atomList,bIsC_SP2);
  for (RbtAtomListIter iter = csp2AtomList.begin(); iter != csp2AtomList.end(); iter++) {
    RbtAtomList bondedAtomList = Rbt::GetBondedAtomList(*iter);

    RbtAtomListIter nsp2Iter;
    RbtAtomListIter ntriIter;
    RbtAtomList nsp2BondedAtomList;
    RbtAtomList ntriBondedAtomList;

    //////////////////
    //switch on #N_SP2
    //////////////////
    switch (Rbt::GetNumAtoms(bondedAtomList,bIsN_SP2)) {

    case 1://1 x N_SP2
      nsp2Iter = Rbt::FindAtom(bondedAtomList,bIsN_SP2);

      //////////////////
      //switch on #N_TRI
      //////////////////
      switch (Rbt::GetNumAtoms(bondedAtomList,bIsN_TRI)) {

      //////////////////////////////////
      //1 x N_SP2, 2 x N_TRI = guanidine
      //////////////////////////////////
      case 2:
	if (!(*iter)->GetCyclicFlag()) {
	  //#ifdef _DEBUG
	  cout << "Protonating neutral guanidine " << (*nsp2Iter)->GetAtomName() << endl;
	  //#endif //_DEBUG
	  AddHydrogen(*nsp2Iter);
	}
#ifdef _DEBUG
	else {
	  cout << "Cyclic guanidine-like fragment - ignored " << (*nsp2Iter)->GetAtomName() << endl;
	}
#endif //_DEBUG
	break;

      //////////////////////////////////
      //1 x N_SP2, 1 x N_TRI = possible imidazole, also need to check for 5 membered ring
      //////////////////////////////////
      case 1:
	ntriIter = Rbt::FindAtom(bondedAtomList,bIsN_TRI);
	nsp2BondedAtomList = Rbt::GetBondedAtomList(*nsp2Iter);
	ntriBondedAtomList = Rbt::GetBondedAtomList(*ntriIter);
#ifdef _DEBUG
	cout << "Possible imidazole/amidine, found 1xN_SP2, 1xN_TRI bonded to " << (*iter)->GetAtomName() << endl;
#endif //_DEBUG
	//IMIDAZOLE - Check that the N_SP2 is bonded to 2 x C_SP2 and that
	//the N_TRI is bonded to 2 x C_SP2 and 1 hydrogen
	if ( (Rbt::GetNumAtoms(nsp2BondedAtomList,bIsC_SP2) == 2) &&
	     (Rbt::GetNumAtoms(ntriBondedAtomList,bIsC_SP2) == 2) &&
	     (Rbt::GetNumAtoms(ntriBondedAtomList,Rbt::isAtomicNo_eq(1)) == 1)) {
#ifdef _DEBUG
	  cout << "Possible imidazole, bonding requirements for N_SP2 and N_TRI met" << endl;
#endif //_DEBUG
	  //Now check if any of the atoms bonded to the N_TRI are 1-2 connected to any
	  //of the atoms bonded to the N_SP2. If so, it is a 5-membered ring
	  //DM 25 Jul 2002 - also check whether the atoms are bridgeheads or not
	  //We don't want to protonate imidazoles which are part of larger fused ring systems
	  //(can have very different pKa's)
	  RbtAtomList atoms12Conn;
	  for (RbtAtomListConstIter iter2 = ntriBondedAtomList.begin(); iter2 != ntriBondedAtomList.end(); iter2++) {
	    RbtAtomList tmpList = Rbt::GetAtomList(nsp2BondedAtomList,Rbt::isAtom_12Connected(*iter2));
	    std::copy(tmpList.begin(), tmpList.end(), std::back_inserter(atoms12Conn));
	  }
	  if (atoms12Conn.size() == 1) {
#ifdef _DEBUG
	    cout << "5-membered imidazole ring found" << endl;
#endif //_DEBUG
	    //Remove bridgehead atoms (fused rings)
	    atoms12Conn = Rbt::GetAtomList(atoms12Conn,std::not1(Rbt::isAtomBridgehead()));
	    if (atoms12Conn.size() == 1) {//Imidazole!!
	      cout << "Protonating neutral imidazole " << (*nsp2Iter)->GetAtomName() << endl;
	      AddHydrogen(*nsp2Iter);
	    }
#ifdef _DEBUG
	    else {
	      cout << "Fused imidazole - ignored " << (*nsp2Iter)->GetAtomName() << endl;
	    }
#endif //_DEBUG
	  }
#ifdef _DEBUG
	  else {
	    cout << "Not a 5-membered ring" << endl;
	  }
#endif //_DEBUG
	}
	//AMIDINE Check that the N_SP2 is bonded to 1 hydrogen and the N_TRI is bonded to 2 hydrogens
	else if ((Rbt::GetNumAtoms(nsp2BondedAtomList,Rbt::isAtomicNo_eq(1)) == 1) && 
		 (Rbt::GetNumAtoms(ntriBondedAtomList,Rbt::isAtomicNo_eq(1)) == 2)) {
	  //#ifdef _DEBUG
	  cout << "Protonating neutral amidine " << (*nsp2Iter)->GetAtomName() << endl;
	  //#endif //_DEBUG
	  AddHydrogen(*nsp2Iter);
	}

#ifdef _DEBUG
	else {
	  cout << "Bonding requirements for N_SP2 and N_TRI not met" << endl;
	}
#endif //_DEBUG
	break;
	//////////////////////
	//End of case #N_TRI=1
	//////////////////////
	
      default:
	break;
      }
      //////////////////////
      //End of switch #N_TRI
      //////////////////////
      break;
      //////////////////////
      //End of case #N_SP2=1
      //////////////////////

    default:
      break;
    }
    //////////////////////
    //End of switch #N_SP2
    //////////////////////
  }
  //End of loop over all C_SP2
}

//Deprotonate all carboxylic acids, phosphates, phosphonates, sulphates, sulphonates
void RbtMdlFileSource::SetupNegIonisableGroups()
{
  //Useful predicates
  Rbt::isFFType_eq bIsO_TRI("O_TRI");
  Rbt::isFFType_eq bIsO_TRIM("O_TRI-");
  Rbt::isFFType_eq bIsO_SP2("O_SP2");
  Rbt::isFFType_eq bIsC_SP2("C_SP2");
  Rbt::isFFType_eq bIsS_SP2("S_SP2");
  Rbt::isFFType_eq bIsP_SP2("P_SP2");

  //Compile list of all C_SP2, S_SP2 and P_SP2
  RbtAtomList csp2AtomList = Rbt::GetAtomList(m_atomList,bIsC_SP2);
  RbtAtomList ssp2AtomList = Rbt::GetAtomList(m_atomList,bIsS_SP2);
  RbtAtomList psp2AtomList = Rbt::GetAtomList(m_atomList,bIsP_SP2);
  RbtAtomList acidAtomList = csp2AtomList;
  acidAtomList.insert(acidAtomList.end(),ssp2AtomList.begin(),ssp2AtomList.end());
  acidAtomList.insert(acidAtomList.end(),psp2AtomList.begin(),psp2AtomList.end());
  
  //Loop over all potential acid centers
  for (RbtAtomListIter iter = acidAtomList.begin(); iter != acidAtomList.end(); iter++) {
    RbtAtomList bondedAtomList = Rbt::GetBondedAtomList(*iter);
    //Check if there is at least one O_SP2 bonded (i.e.C=O, S=O, P=O)
    //DM 27 Apr 1999 - check there is already an O_TRI- present. We only want to deprotonate one oxygen.
    if ( (Rbt::GetNumAtoms(bondedAtomList,bIsO_SP2) > 0) &&
	 (Rbt::GetNumAtoms(bondedAtomList,bIsO_TRIM) == 0 ) ) {
      //Now get the list of all O_TRI which are bonded (i.e. C-O-, S-O-, P-O-)
      RbtAtomList otriAtomList = Rbt::GetAtomList(bondedAtomList,bIsO_TRI);
      for (RbtAtomListIter otriIter = otriAtomList.begin(); otriIter != otriAtomList.end(); otriIter++) {
	//Is there a hydrogen bonded to the O_TRI ?
	RbtAtomList otriBondedAtomList = Rbt::GetBondedAtomList(*otriIter);
	RbtAtomListIter hIter = Rbt::FindAtom(otriBondedAtomList,Rbt::isAtomicNo_eq(1));
	if (hIter != otriBondedAtomList.end()) {
	  //#ifdef _DEBUG
	  cout << "Removing " << (*hIter)->GetAtomName() << " from " << (*otriIter)->GetAtomName()
	       << ", acid center=" << (*iter)->GetAtomName() << endl;
	  //#endif //_DEBUG
	  RemoveAtom(*hIter);
	  //Adjust the attributes on the O_TRI
	  (*otriIter)->SetFormalCharge(-1);
	  (*otriIter)->SetGroupCharge(-1.0);
	  (*otriIter)->SetFFType("O_TRI-");
	  break;//DM 19 Apr 1999 - only deprotonate the first O-H on any acid group
	}
      }
    }
  }
}

//Adds a proton to N_SP3 and N_SP2 atoms
void RbtMdlFileSource::AddHydrogen(RbtAtomPtr spAtom) throw (RbtError)
{
  Rbt::isFFType_eq bIsNSP3("N_SP3");
  Rbt::isFFType_eq bIsNSP2("N_SP2");
  
  //Element data for hydrogen
  RbtElementData elHData = m_spElementData->GetElementData(1);
  //Radius increment for H-bonding hydrogens
  RbtDouble dHBondRadIncr = m_spElementData->GetHBondRadiusIncr();
  
  //Attributes for new atom (to be defined)
  RbtCoord hCoord;//New coord
  RbtDouble hRadius;//Vdw radius

  //Attributes for heavy atom
  const RbtCoord& c0 = spAtom->GetCoords();
  RbtAtomList bondedAtomList = Rbt::GetBondedAtomList(spAtom);

  ////////////////
  //Nsp3
  //We determine the plane formed by the three substituents and add the
  //hydrogen along a vector normal to the plane (and in the other direction)
  //    H
  //    |
  //    N
  //   '|`
  ////////////////
  if (bIsNSP3(spAtom)) {
    if (bondedAtomList.size() != 3)
      throw RbtModelError(_WHERE_, spAtom->GetAtomName()+" has inconsistent atom type");
    RbtDouble dNHLength = 1.00;//N-H bond length
    RbtVector v1 = bondedAtomList[0]->GetCoords()-c0;//Vector to subst 1
    RbtVector v2 = bondedAtomList[1]->GetCoords()-c0;//Vector to subst 2
    RbtVector v3 = bondedAtomList[2]->GetCoords()-c0;//Vector to subst 3
    //Take unit vectors along each bond vector to eliminate the effect of different bond lengths
    RbtPlane p = RbtPlane(c0+v1.Unit(),c0+v2.Unit(),c0+v3.Unit());//Plane of three substs
    //This is the signed distance from the nitrogen to the plane
    RbtDouble dist = Rbt::DistanceFromPointToPlane(c0,p);
    //Take sign so we go in the opposite direction to the other substituents
    RbtInt iSign = (dist>0) ? 1 : -1;
    hCoord = c0+(iSign*dNHLength*p.VNorm()); //New H coordinate
    hRadius = elHData.vdwRadius+dHBondRadIncr; //New H vdw radius
    //Adjust the attributes on the Nitrogen
    spAtom->SetFormalCharge(1);
    spAtom->SetGroupCharge(1.0);
    spAtom->SetFFType("N_SP3+");
  }

  ////////////////
  //Nsp2
  //We determine the average of the bond vectors to the two substituents and add the
  //hydrogen along this vector in the other direction
  //   H
  //   |
  //   N
  //  ' `
  ////////////////
  else if (bIsNSP2(spAtom)) {
    if (bondedAtomList.size() != 2)
      throw RbtModelError(_WHERE_, spAtom->GetAtomName()+" has inconsistent atom type");
    RbtDouble dNHLength = 1.00;//N-H bond length
    RbtVector v1 = bondedAtomList[0]->GetCoords()-c0;//Vector to subst 1
    RbtVector v2 = bondedAtomList[1]->GetCoords()-c0;//Vector to subst 2
    //Take unit vectors along each bond vector to eliminate the effect of different bond lengths
    RbtVector vMean = (v1.Unit()+v2.Unit())/2.0;
    hCoord = c0-(dNHLength*vMean.Unit()); //New H coordinate
    hRadius = elHData.vdwRadius+dHBondRadIncr; //New H vdw radius
    //Adjust the attributes on the Nitrogen
    spAtom->SetFormalCharge(1);
    spAtom->SetGroupCharge(1.0);
    spAtom->SetFFType("N_SP2+");
  }

  ////////////////
  //Don't consider other element types for now
  ////////////////
  else
    return;

  //Construct the new hydrogen atom (constructor only accepts the 2D params)
  RbtInt nAtomId = m_atomList.size()+1;
  ostrstream ostr;
  ostr << "H" << nAtomId << ends;
  RbtString strAtomName(ostr.str());
  delete ostr.str();
  RbtAtomPtr spHAtom(new RbtAtom(nAtomId,
				 1,//nAtomicNo,
				 strAtomName,
				 "1",//strSubunitId,
				 "MOL",//strSubunitName,
				 "H",//strSegmentName,
				 RbtAtom::SP3,//Hybridisation state
				 0,//Num implicit hydrogens
				 0) //Formal charge
		     );

  //Now set the remaining 2-D and 3-D params
  spHAtom->SetCoords(hCoord);
  spHAtom->SetPartialCharge(0.0);
  spHAtom->SetGroupCharge(0.0);
  spHAtom->SetAtomicMass(elHData.mass);
  spHAtom->SetVdwRadius(hRadius);
  m_atomList.push_back(spHAtom);

  //Construct a new bond to the hydrogen
  RbtInt nBondId = m_bondList.size()+1;
  RbtBondPtr spBond(new RbtBond(nBondId,spAtom,spHAtom,1));//Bond order = 1
  m_bondList.push_back(spBond);
}



//SetupIonicGroups sets the "group charge" attribute of the atoms. The idea is that the "group charge" is
//file-format independent, and eliminates the need for the scoring function setup to have any knowledge of
//the input file format. For example, PSF files represent COO- as OC-C-OC (both oxygens charged) whereas
//SD files represent it as the formal resonance structure O=C-O-. For RiboDock we actually want the negative
//charge on the central carbon!
//
//It is envisaged each molecular file source will have a version of SetupIonicGroups to convert from the native
//representation to the RiboDock representation. The FormalCharge and PartialCharge attributes should be left
//untouched to allow the model to be rewritten back in the same format.
void RbtMdlFileSource::SetupIonicGroups()
{
  //Group charges have been initialised to the formal charge by Parse() so we only need
  //to check for special cases here.
  SetupNSP3Plus();
  SetupNSP2Plus();
  SetupOTRIMinus();

}

//Helper function to set up protonated amines
void RbtMdlFileSource::SetupNSP3Plus()
{
  Rbt::isFFType_eq bIsN_SP3p("N_SP3+");
  Rbt::isAtomicNo_eq bIsH(1);

  //Divide up the charge equally amongst the hydrogens bonded to the nitrogen
  RbtAtomList nsp3pAtomList = Rbt::GetAtomList(m_atomList,bIsN_SP3p);
  for (RbtAtomListIter iter = nsp3pAtomList.begin(); iter != nsp3pAtomList.end(); iter++) {
    //Get the list of hydrogens
    RbtAtomList hbdAtomList = Rbt::GetAtomList(Rbt::GetBondedAtomList(*iter),bIsH);
    RbtDouble nCharge = (*iter)->GetGroupCharge();//Total charge
    RbtDouble pCharge = nCharge/(hbdAtomList.size());//Partial charge
    (*iter)->SetGroupCharge(0.0);//Neutralise the nitrogen
    for (RbtAtomListIter hIter = hbdAtomList.begin(); hIter != hbdAtomList.end(); hIter++) {
      (*hIter)->SetGroupCharge(pCharge);//Charge up the hydrogens
#ifdef _DEBUG
      cout << "Transferring charge of " << pCharge << " from " << (*iter)->GetAtomName()
	   << " to " << (*hIter)->GetAtomName() << endl;
#endif //_DEBUG
    }
  }
}

//Helper function to set up guanidines, imidazoles, amidines
//DM 23 Nov 2000 - also handles protonated pyridines correctly (transfers +1 charge from N to H)
void RbtMdlFileSource::SetupNSP2Plus()
{
  //Useful predicates
  Rbt::isFFType_eq bIsN_SP2p("N_SP2+");
  Rbt::isFFType_eq bIsN_TRI("N_TRI");
  Rbt::isFFType_eq bIsC_SP2("C_SP2");
  Rbt::isAtomicNo_eq bIsH(1);

  //Loop over all N_SP2+
  RbtAtomList nsp2PlusAtomList = Rbt::GetAtomList(m_atomList,bIsN_SP2p);
  for (RbtAtomListIter iter = nsp2PlusAtomList.begin(); iter != nsp2PlusAtomList.end(); iter++) {
    //Now get the list of all bonded C_SP2
    RbtAtomList csp2AtomList = Rbt::GetAtomList(Rbt::GetBondedAtomList(*iter),bIsC_SP2);
#ifdef _DEBUG
    cout << csp2AtomList.size() << " CSP2 atoms found bonded to NSP2+" << endl;
#endif //_DEBUG
    //Eliminate bridgehead carbons, so we can guarantee we handle imidazoles correctly
    csp2AtomList = Rbt::GetAtomList(csp2AtomList,std::not1(Rbt::isAtomBridgehead()));
#ifdef _DEBUG
    cout << csp2AtomList.size() << " are non-bridgehead" << endl;
#endif //_DEBUG
    //There may be cases where the guanidinium and imidazole bonding patterns are both present
    //(on different C_SP2 atoms) so we want to check for both. Let the guanidinium moiety take precedence
    RbtAtomListIter guanIter = Rbt::FindAtom(csp2AtomList,Rbt::isCoordinationNumber_eq(2,"N_TRI"));
    RbtAtomListIter imidIter = Rbt::FindAtom(csp2AtomList,Rbt::isCoordinationNumber_eq(1,"N_TRI"));
    RbtAtomListIter csp2Iter = csp2AtomList.end();//This will point to the C_SP2 atom to transfer the charge to
    if (guanIter != csp2AtomList.end())
      csp2Iter = guanIter;
    else if (imidIter != csp2AtomList.end())
      csp2Iter = imidIter;

#ifdef _DEBUG
    cout << "Checking N_SP2+ atom " << (*iter)->GetAtomName() << endl;
    if (guanIter != csp2AtomList.end())
      cout << (*guanIter)->GetAtomName() << " is a guanidinium-like carbon" << endl;
    if (imidIter != csp2AtomList.end())
      cout << (*imidIter)->GetAtomName() << " is a imidazole-like carbon" << endl;    
#endif //_DEBUG
    
    //DM 23 Nov 2000 - move this code out of the following if {} block, to handle protonated pyridines
    //Get the list of hydrogens on the attached nitrogens
    //Start with the hydrogen on the N_SP2+
    RbtAtomList hbdAtomList = Rbt::GetAtomList(Rbt::GetBondedAtomList(*iter),bIsH);
#ifdef _DEBUG
    cout << "Found " << hbdAtomList.size() << " hydrogens bonded to " << (*iter)->GetAtomName() << endl;
#endif //_DEBUG

    //DM 23 Nov 2000 - this bit handles guanidines, imidazoles, amidines
    if (csp2Iter != csp2AtomList.end()) {
      //Now add the hydrogens on the (1 or 2) N_TRI atoms
      RbtAtomList ntriAtomList = Rbt::GetAtomList(Rbt::GetBondedAtomList(*csp2Iter),bIsN_TRI);
#ifdef _DEBUG
      cout << "Found " << ntriAtomList.size() << " nitrogens bonded to " << (*csp2Iter)->GetAtomName() << endl;
#endif //_DEBUG
      for (RbtAtomListConstIter nIter = ntriAtomList.begin(); nIter != ntriAtomList.end(); nIter++) {
	RbtAtomList hAtomList = Rbt::GetAtomList(Rbt::GetBondedAtomList(*nIter),bIsH);
#ifdef _DEBUG
	cout << "Found " << hAtomList.size() << " hydrogens bonded to " << (*nIter)->GetAtomName() << endl;
#endif //_DEBUG
	hbdAtomList.insert(hbdAtomList.end(),hAtomList.begin(),hAtomList.end());
      }
      hbdAtomList.push_back(*csp2Iter);//Add the central carbon to the list of charged atoms
    }

    //Divide up the charge
    RbtInt nH = hbdAtomList.size();
    //DM 23 Nov 2000 - only transfer the charge if there is at least one attached H
    if (nH > 0) {
      RbtDouble nCharge = (*iter)->GetGroupCharge();
      RbtDouble pCharge = nCharge/nH;
      (*iter)->SetGroupCharge(0.0);//N_SP2+
      for (RbtAtomListIter hIter = hbdAtomList.begin(); hIter != hbdAtomList.end(); hIter++) {
	(*hIter)->SetGroupCharge(pCharge);//H
#ifdef _DEBUG
	cout << "Transferring charge of " << pCharge << " from " << (*iter)->GetAtomName()
	     << " to " << (*hIter)->GetAtomName() << endl;
#endif //_DEBUG
      }
    }
  }
}

//Helper function to set up deprotonated acids
void RbtMdlFileSource::SetupOTRIMinus()
{
  Rbt::isFFType_eq bIsO_TRIm("O_TRI-");
  Rbt::isFFType_eq bIsO_SP2("O_SP2");
  Rbt::isFFType_eq bIsC_SP2("C_SP2");
  Rbt::isFFType_eq bIsS_SP2("S_SP2");
  Rbt::isFFType_eq bIsP_SP2("P_SP2");

  //Compile list of all C_SP2, S_SP2 and P_SP2
  RbtAtomList csp2AtomList = Rbt::GetAtomList(m_atomList,bIsC_SP2);
  RbtAtomList ssp2AtomList = Rbt::GetAtomList(m_atomList,bIsS_SP2);
  RbtAtomList psp2AtomList = Rbt::GetAtomList(m_atomList,bIsP_SP2);
  RbtAtomList acidAtomList = csp2AtomList;
  acidAtomList.insert(acidAtomList.end(),ssp2AtomList.begin(),ssp2AtomList.end());
  acidAtomList.insert(acidAtomList.end(),psp2AtomList.begin(),psp2AtomList.end());
  
  //Look for acid centers which have at least one bonded O_SP2 and at least one
  //bonded O_TRI-
  //e.g. C=O
  //     |
  //     O-
  //Divide total charge equally over the oxygens

  //Loop over all potential acid centers
  for (RbtAtomListConstIter iter = acidAtomList.begin(); iter != acidAtomList.end(); iter++) {
    //Get the list of bonded O_SP2 and O_TRI-
    RbtAtomList osp2AtomList = Rbt::GetAtomList(Rbt::GetBondedAtomList(*iter),bIsO_SP2);
    RbtAtomList otrimAtomList = Rbt::GetAtomList(Rbt::GetBondedAtomList(*iter),bIsO_TRIm);
    if ((osp2AtomList.size() > 0) && (otrimAtomList.size() > 0)) {
      //Combine O_SP2 and O_TRI- lists into one
      RbtAtomList oAtomList = osp2AtomList;
      oAtomList.insert(oAtomList.end(),otrimAtomList.begin(),otrimAtomList.end());
      //Total charge is the number of O_TRI-
      RbtDouble nCharge = otrimAtomList.size();//Total charge
      //Partial charge is total / divided by number of oxygens (O_SP2 + O_TRI-)
      RbtDouble pCharge = -nCharge/(oAtomList.size());//Partial charge
      for (RbtAtomListIter oIter = oAtomList.begin(); oIter != oAtomList.end(); oIter++) {
	(*oIter)->SetGroupCharge(pCharge);//Charge up the oxygens
#ifdef _DEBUG
	cout << "Transferring charge of " << pCharge << " from " << (*iter)->GetAtomName()
	   << " to " << (*oIter)->GetAtomName() << endl;
#endif //_DEBUG
      }
    }
  }
}

//DM 4 June 1999 - search for the separate fragments in the record
//and set distinct segment names for each one
//Segment names are initially H1,H2,H3...
//except that the largest segment is renamed back to H
void RbtMdlFileSource::SetupSegmentNames()
{
  RbtString strLargestSegName;//Name of the largest segment
  RbtInt nMaxSize(0);//Size of the largest segment
  RbtInt nSeg;
  RbtAtomListIter seed;
  for (nSeg=1,seed = m_atomList.begin();
       seed != m_atomList.end();
       nSeg++,seed = Rbt::FindAtom(m_atomList,Rbt::isSegmentName_eq("H"))) {
    //New segment name (H1, H2 etc)
    ostrstream ostr;
    ostr << "H" << nSeg << ends;
    RbtString strSegName(ostr.str());
    delete ostr.str();
    //Temporary atom list containing atoms to be processed
    //Note: this is a true list (not a vector) as we will be making numerous insertions and deletions
    (*seed)->SetSegmentName(strSegName);
    RbtAtomTrueList pendingAtomList;
    pendingAtomList.push_back(*seed);
    //While we still have atoms to process
    //DM 12 May 2000 - correct nSize for ANSI scoping
    RbtInt nSize(0);
    for (nSize = 0; !pendingAtomList.empty(); nSize++) {
      //Take the last atom from the list and remove it
      RbtAtomPtr spAtom = pendingAtomList.back();
      pendingAtomList.pop_back();
      //Get the list of bonded atoms which haven't yet been processed (i.e. those whose segName is still H)
      RbtAtomList atomsToAdd = Rbt::GetAtomList(Rbt::GetBondedAtomList(spAtom),Rbt::isSegmentName_eq("H"));
      //Set the segment name for each atom, and add to the queue
      for (RbtAtomListIter iter = atomsToAdd.begin(); iter != atomsToAdd.end(); iter++) {
	(*iter)->SetSegmentName(strSegName);
	pendingAtomList.push_back(*iter);
      }
    }
#ifdef _DEBUG
    cout << "Segment " << strSegName << " contains " << nSize << " atoms" << endl;
#endif //_DEBUG

    //Keep track of the largest segment
    if (nSize > nMaxSize) {
      nMaxSize = nSize;
      strLargestSegName = strSegName;
    }
  }

  //Rename the largest segment back to H
#ifdef _DEBUG
  cout << "Largest segment is " << strLargestSegName << endl;
#endif //_DEBUG
  RbtAtomList largestAtomList = Rbt::GetAtomList(m_atomList,Rbt::isSegmentName_eq(strLargestSegName));
  for (RbtAtomListIter iter = largestAtomList.begin(); iter != largestAtomList.end(); iter++) {
    (*iter)->SetSegmentName("H");
  }

  //Now update the segment map
  m_segmentMap.clear();
  for (RbtAtomListIter iter = m_atomList.begin(); iter != m_atomList.end(); iter++) {
    m_segmentMap[(*iter)->GetSegmentName()]++;
  }
}

//DM 2 Aug 1999 - remove all non-polar hydrogens
void RbtMdlFileSource::RemoveNonPolarHydrogens()
{
  //Get list of all hydrogens
  RbtAtomList hList = Rbt::GetAtomList(m_atomList,Rbt::isAtomicNo_eq(1));
  //Get list of all hydrogens bonded to carbon
  RbtAtomList implHList = Rbt::GetAtomList(hList,Rbt::isCoordinationNumber_eq(1,6));

  //Remove them one by one (not very efficient)
  for (RbtAtomListIter iter = implHList.begin(); iter != implHList.end(); iter++) {
    RemoveAtom(*iter);
  }
#ifdef _DEBUG
  for (RbtBondListConstIter bIter = m_bondList.begin(); bIter != m_bondList.end(); bIter++) {
    cout << "Bond ID=" << (*bIter)->GetBondId() << " ("
	 << (*bIter)->GetAtom1Ptr()->GetAtomName() << "-"
	 << (*bIter)->GetAtom2Ptr()->GetAtomName() << ")"<< endl;
  }
#endif //_DEBUG
}

//DM 8 Feb 2000 - setup atom and bond cyclic flags (previously in RbtModel)
//Set the atom and bond cyclic flags for all atoms and bonds in the source
void RbtMdlFileSource::SetCyclicFlags()
{
  Rbt::SetAtomAndBondCyclicFlags(m_atomList,m_bondList);
}
