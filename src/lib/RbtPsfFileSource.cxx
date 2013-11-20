/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtPsfFileSource.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtElementFileSource.h"
#include "RbtPsfFileSource.h"
#include "RbtFileError.h"

//Constructors
//RbtPsfFileSource::RbtPsfFileSource(const char* fileName) :
//  RbtBaseMolecularFileSource(fileName,"PSF_FILE_SOURCE") //Call base class constructor
//{
//  _RBTOBJECTCOUNTER_CONSTR_("RbtPsfFileSource");
//}

RbtPsfFileSource::RbtPsfFileSource(const RbtString& fileName,
				   const RbtString& strMassesFile,
				   RbtBool bImplHydrogens) :
  RbtBaseMolecularFileSource(fileName,"PSF_FILE_SOURCE"), //Call base class constructor
  m_bImplHydrogens(bImplHydrogens)
{
  //Open a Charmm data source for converting atom type number to the more friendly string variety
  m_spCharmmData = RbtCharmmDataSourcePtr(new RbtCharmmDataSource(strMassesFile));
  //DM 13 Jun 2000 - new separate prm file for receptor ionic atom definitions
  m_spParamSource = RbtParameterFileSourcePtr(new RbtParameterFileSource(Rbt::GetRbtFileName("data/sf","RbtIonicAtoms.prm")));
  //Open an Element data source
  m_spElementData = RbtElementFileSourcePtr(new RbtElementFileSource(Rbt::GetRbtFileName("data","RbtElements.dat")));
  _RBTOBJECTCOUNTER_CONSTR_("RbtPsfFileSource");
}


//Default destructor
RbtPsfFileSource::~RbtPsfFileSource()
{
  _RBTOBJECTCOUNTER_DESTR_("RbtPsfFileSource");
}

void RbtPsfFileSource::Parse() throw (RbtError)
{
  //Expected string constants in PSF files
  const RbtString strPsfKey("PSF");
  const RbtString strTitleKey("!NTITLE");
  const RbtString strAtomKey("!NATOM");
  //DM 19 Feb 2002 - InsightII writes PSF files with !NBONDS: string rather than !NBOND:
  //so need to cope with either
  const RbtString strBondKey("!NBOND:");
  //const RbtString strBondKey("!NBOND");
  RbtString strKey;

  //Only parse if we haven't already done so
  if (!m_bParsedOK) {
    ClearMolCache();//Clear current cache
    Read();//Read the file

    try {
      RbtFileRecListIter fileIter = m_lineRecs.begin();
      RbtFileRecListIter fileEnd = m_lineRecs.end();


      //////////////////////////////////////////////////////////
      //1. Check for PSF string on line 1
      if ( (fileIter == fileEnd) || ((*fileIter).substr(0,3) != strPsfKey))
	throw RbtFileParseError(_WHERE_,"Missing "+strPsfKey+" string in " + GetFileName());
    
      //////////////////////////////////////////////////////////
      //2a Read number of title lines and check for correct title key...
      RbtInt nTitleRec;
      fileIter+=2;
      istrstream((*fileIter).c_str()) >> nTitleRec >> strKey;
      if (strKey != strTitleKey)
	throw RbtFileParseError(_WHERE_,"Missing "+strTitleKey+" string in " + GetFileName());

      //2b ...and store them
      fileIter++;
      m_titleList.reserve(nTitleRec);//Allocate enough memory for the vector
      while ( (m_titleList.size() < nTitleRec) && (fileIter != fileEnd)) {
	m_titleList.push_back(*fileIter++);
      }

      //2c ..and check we read them all before reaching the end of the file
      if (m_titleList.size() != nTitleRec)
	throw RbtFileParseError(_WHERE_,"Incomplete title records in " + GetFileName());

      //////////////////////////////////////////////////////////
      //3a. Read number of atoms and check for correct atom key...
      RbtInt nAtomRec;
      fileIter++;
      istrstream((*fileIter).c_str()) >> nAtomRec >> strKey;
      if (strKey != strAtomKey)
	throw RbtFileParseError(_WHERE_,"Missing "+strAtomKey+" string in " + GetFileName());
      
      //3b ...and store them
      fileIter++;
      m_atomList.reserve(nAtomRec);//Allocate enough memory for the vector

      RbtInt nAtomId; //original atom number in PSF file
      RbtString strSegmentName; //segment name in PSF file
      RbtString strSubunitId; //subunit(residue) ID in PSF file
      RbtString strSubunitName; //subunit(residue) name in PSF file
      RbtString strAtomName; //atom name from PSF file
      RbtString strFFType; //force field type from PSF file
      //RbtInt nFFType; //force field type from PSF file (integer)
      RbtDouble dPartialCharge; //partial charge
      RbtDouble dAtomicMass; //atomic mass from PSF file

      while ( (m_atomList.size() < nAtomRec) && (fileIter != fileEnd)) {
	istrstream istr((*fileIter++).c_str());
	istr >> nAtomId
	     >> strSegmentName
	     >> strSubunitId
	     >> strSubunitName
	     >> strAtomName
	     >> strFFType
	     >> dPartialCharge
	     >> dAtomicMass;

	//Construct a new atom (constructor only accepts the 2D params)
	RbtAtomPtr spAtom(new RbtAtom(nAtomId,
				      0,//Atomic number undefined (gets set below)
				      strAtomName,
				      strSubunitId,
				      strSubunitName,
				      strSegmentName
				      )
			  );

	//Now set the remaining 2-D and 3-D params
	//Use the Charmm data source to provide some params not in the PSF file
	//
	//DM 4 Jan 1999 - check if force field type is numeric (parm22) or string (parm19)
	RbtInt nFFType = atoi(strFFType.c_str());
	if (nFFType > 0)
	  strFFType = m_spCharmmData->AtomTypeString(nFFType);//Convert atom type from int to string
	spAtom->SetFFType(strFFType);
	spAtom->SetAtomicNo(m_spCharmmData->AtomicNumber(strFFType));//Now we can set the atomic number correctly
	spAtom->SetNumImplicitHydrogens(m_spCharmmData->ImplicitHydrogens(strFFType));//Set # of implicit hydrogens
	spAtom->SetHybridState(m_spCharmmData->HybridState(strFFType)); //DM 8 Dec 1998 Set hybridisation state
	//spAtom->SetVdwRadius(m_spCharmmData->VdwRadius(strFFType));
	//spAtom->SetFormalCharge(m_spCharmmData->FormalCharge(strFFType));
	spAtom->SetGroupCharge(0.0);
	spAtom->SetPartialCharge(dPartialCharge);
	spAtom->SetAtomicMass(dAtomicMass);

	m_atomList.push_back(spAtom);
	m_segmentMap[strSegmentName]++;//increment atom count in segment map
      }

      //3c ..and check we read them all before reaching the end of the file
      if (m_atomList.size() != nAtomRec)
	throw RbtFileParseError(_WHERE_,"Incomplete atom records in " + GetFileName());

      //////////////////////////////////////////////////////////
      //4a. Read number of bonds and check for correct bond key...
      RbtInt nBondRec;
      fileIter++;
      //istrstream((*fileIter).c_str()) >> nBondRec >> strKey;
      istrstream((*fileIter).c_str()) >> nBondRec >> strKey;
      //cout << "strKey "<<strKey << " strBondKey " << strBondKey << endl;
      //if(strBondKey.compare(strKey,0,6)) {// 6 for "!NBOND" old style API
      RbtInt iCmp = strBondKey.compare(0,5,strKey);
      //cout << "iCmp = " << iCmp << endl;
      if(iCmp>0) {// 6 for "!NBOND" new API
	throw RbtFileParseError(_WHERE_,"Missing "+strBondKey+" string in " + GetFileName());
      }
      
      //4b ...and store them
      fileIter++;
      m_bondList.reserve(nBondRec);//Allocate enough memory for the vector

      RbtInt nBondId(0);
      RbtUInt idxAtom1;
      RbtUInt idxAtom2;
      while ( (m_bondList.size() < nBondRec) && (fileIter != fileEnd)) {
	istrstream istr((*fileIter++).c_str());

	//Read bonds slightly differently to atoms as we have 4 bonds per line
	//Current method assumes atoms are numbered consectively from 1
	//May be better to store a map rather than a vector ??
	while (istr >> idxAtom1 >> idxAtom2) {
	  if ( (idxAtom1 > nAtomRec) || (idxAtom2 > nAtomRec) ) {//Check for indices in range
	    throw RbtFileParseError(_WHERE_,"Atom index out of range in bond records in " + GetFileName());
	  }
	  RbtAtomPtr spAtom1(m_atomList[idxAtom1-1]);//Decrement the atom index as the atoms are numbered from zero in our atom vector
	  RbtAtomPtr spAtom2(m_atomList[idxAtom2-1]);//Decrement the atom index as the atoms are numbered from zero in our atom vector
	  RbtBondPtr spBond(new RbtBond(++nBondId,//Store a nominal bond ID starting from 1
					spAtom1,
					spAtom2
					)
			    );
	  m_bondList.push_back(spBond);
	}
      }
      //4c ..and check we read them all before reaching the end of the file
      if (m_bondList.size() != nBondRec)
	throw RbtFileParseError(_WHERE_,"Incomplete bond records in " + GetFileName());

      //Setup the atomic params not stored in the file
      SetupAtomParams();

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
void RbtPsfFileSource::SetupAtomParams() throw (RbtError)
{
  //In the PSF file we have:
  //  Force field type
  //  X,Y,Z coords
  //  Partial charges
  //  Atomic mass
  //
  //CharmmDataSource has already provided (based on force field type):
  //  Element type
  //  Hybridisation state - sp3,sp2,sp,trigonal planar
  //  Number of implicit hydrogens

  //We need to define:
  //  Corrected vdW radii on H-bonding hydrogens and extended atoms
  //  Interaction group charges

  //DM 16 Feb 2000 - remove all non-polar hydrogens if required
  if (m_bImplHydrogens)
    RemoveNonPolarHydrogens();
  SetupVdWRadii();
  SetupPartialIonicGroups();//Partial ionic groups (N7 etc)
  RenumberAtomsAndBonds();//DM 30 Oct 2000 - tidy up atom and bond numbering
}

//Defines vdW radius, correcting for extended atoms and H-bond donor hydrogens
void RbtPsfFileSource::SetupVdWRadii() throw (RbtError)
{
  //Radius increment for atoms with implicit hydrogens
  //DM 22 Jul 1999 - only increase the radius for sp3 atoms with implicit hydrogens
  //For sp2 and aromatic, leave as is
  Rbt::isHybridState_eq bIsSP3(RbtAtom::SP3);
  Rbt::isHybridState_eq bIsTri(RbtAtom::TRI);
  Rbt::isCoordinationNumber_eq bTwoBonds(2);

  RbtDouble dImplRadIncr = m_spElementData->GetImplicitRadiusIncr();
  //Element data for hydrogen
  RbtElementData elHData = m_spElementData->GetElementData(1);
  //Radius increment and predicate for H-bonding hydrogens
  Rbt::isAtomHBondDonor bIsHBondDonor;
  RbtDouble dHBondRadius = elHData.vdwRadius + m_spElementData->GetHBondRadiusIncr();

  for (RbtAtomListIter iter = m_atomList.begin(); iter != m_atomList.end(); iter++) {
    //Get the element data for this atom
    RbtInt nAtomicNo = (*iter)->GetAtomicNo();
    //DM 04 Dec 2003 - correction to hybrid state for Ntri
    //Problem is Charmm type 34 is used ambiguously for "Nitrogen in 5-membered ring"
    //Some should be Ntri (3 bonds), others should Nsp2 (2 bonds)
    //masses.rtf defines hybrid state as TRI so switch here if necessary
    if ( (nAtomicNo == 7) && bIsTri(*iter) && bTwoBonds(*iter)) {
      (*iter)->SetHybridState(RbtAtom::SP2);
      //cout << "Switch from N_tri to N_sp2: " << (*iter)->GetFullAtomName() << endl;
    }
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


void RbtPsfFileSource::SetupPartialIonicGroups() {
  //Break the atom list up into substructures and call the base class SetupPartialIonicGroups
  //to assign group charges by residue
  //Completely general code, does not assume that all atoms in each substructure are continguous
  //
  //Deselect all atoms
  std::for_each(m_atomList.begin(),m_atomList.end(),Rbt::SelectAtom(false));
  //The loop increment searches for the next unselected atom
  for (RbtAtomListIter iter = m_atomList.begin();
       iter != m_atomList.end();
       iter = std::find_if(iter+1,m_atomList.end(),std::not1(Rbt::isAtomSelected()))) {
    //Copy all atoms that belong to same substructure as head atom
    RbtAtomList ss;
    std::copy_if(iter,m_atomList.end(),std::back_inserter(ss),isSS_eq(*iter));
    //cout << "Psf SetupPartialIonicGroups: SS from " << ss.front()->GetFullAtomName()
    //	 << " to " << ss.back()->GetFullAtomName() << " (" << ss.size() << ") atoms" << endl;
    //Assign group charges for this residue
    RbtBaseMolecularFileSource::SetupPartialIonicGroups(ss,m_spParamSource);
    //Mark each atom in substructure as having been processed
    std::for_each(ss.begin(),ss.end(),Rbt::SelectAtom(true));
  }
}

//DM 16 Feb 2000 - remove all non-polar hydrogens, adjusts number of implicit hydrogens
void RbtPsfFileSource::RemoveNonPolarHydrogens()
{
  //Get list of all carbons
  RbtAtomList cList = Rbt::GetAtomList(m_atomList,Rbt::isAtomicNo_eq(6));

  for (RbtAtomListIter cIter = cList.begin(); cIter != cList.end(); cIter++) {
    //Get list of all bonded hydrogens
    RbtAtomList hList = Rbt::GetAtomList(Rbt::GetBondedAtomList(*cIter),Rbt::isAtomicNo_eq(1));
    RbtInt nImplH = hList.size();
    if (nImplH > 0) {
      for (RbtAtomListIter hIter = hList.begin(); hIter != hList.end(); hIter++) {
				RemoveAtom(*hIter);
      }
      //Adjust number of implicit hydrogens
      (*cIter)->SetNumImplicitHydrogens((*cIter)->GetNumImplicitHydrogens()+nImplH);
#ifdef _DEBUG
      cout << "Removing " << nImplH << " hydrogens from " << (*cIter)->GetFullAtomName() << endl;
#endif //_DEBUG
    }
  }
}

