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

//SD file dump, listing number of atoms, bonds, interaction centers etc

#include <algorithm> //for min,max
#include <iomanip>
using namespace std;
#include "RbtMdlFileSource.h"
#include "RbtMdlFileSink.h"

const RbtString EXEVERSION = " ($Id: //depot/dev/client3/rdock/2013.1/src/exe/rblist.cxx#4 $)";

namespace Rbt {
  class isAtomNitroN : public std::unary_function<RbtAtomPtr,RbtBool> {
    //Useful predicates
    Rbt::isFFType_eq bIsN_SP2P;
    Rbt::isFFType_eq bIsO_SP2;
    Rbt::isFFType_eq bIsO_SP3M;
    
    //DM 11 Jul 2001 - checks for nitro nitrogen (NSP2+ bonded to OSP3- and OSP2)
  public:
    explicit isAtomNitroN() : bIsN_SP2P("N_SP2+"),bIsO_SP2("O_SP2"),bIsO_SP3M("O_SP3-") {}
    RbtBool operator() (RbtAtomPtr spAtom) const {
      if (!bIsN_SP2P(spAtom))
	return false;
      RbtAtomList bondedAtomList = Rbt::GetBondedAtomList(spAtom);
      return (Rbt::GetNumAtoms(bondedAtomList,bIsO_SP2) == 1) && (Rbt::GetNumAtoms(bondedAtomList,bIsO_SP3M) == 1);
    }
  };
}

void ListAtoms(const RbtAtomList& atomList, RbtString strTitle, RbtInt nMaxToPrint = 20)
{
  RbtInt nSize = atomList.size();
  RbtInt nToPrint = std::min(nMaxToPrint,nSize);
  cout << endl << "1st " << nToPrint << " " << strTitle << " (out of " << nSize << ")" << endl;
  for (RbtAtomListConstIter iter = atomList.begin(); (iter != (atomList.begin()+nMaxToPrint)) && (iter != atomList.end()); iter++)
    cout << (**iter) << endl;
}
  
void ListBonds(const RbtBondList& bondList, RbtString strTitle, RbtInt nMaxToPrint = 20)
{
  RbtInt nSize = bondList.size();
  RbtInt nToPrint = std::min(nMaxToPrint,nSize);
  cout << endl << "1st " << nToPrint << " " << strTitle << " (out of " << nSize << ")" << endl;
  for (RbtBondListConstIter iter = bondList.begin(); (iter != (bondList.begin()+nMaxToPrint)) && (iter != bondList.end()); iter++)
    cout << (**iter) << endl;
}

//DM 7 Jun 1999 - set primary and secondary amide bonds to 180.0 deg
void CheckAmideBonds(RbtModelPtr spModel)
{
  //Useful predicates
  Rbt::isFFType_eq bIsN_TRI("N_TRI");
  Rbt::isFFType_eq bIsC_SP2("C_SP2");
  Rbt::isFFType_eq bIsO_SP2("O_SP2");
  Rbt::isAtomicNo_eq bIsH(1);

  RbtBondList amideBondList = Rbt::GetBondList(spModel->GetBondList(),Rbt::isBondAmide());
  //DM 9 Oct 2000 - remove cyclic bonds from list
  amideBondList = Rbt::GetBondList(amideBondList,std::not1(Rbt::isBondCyclic()));
  RbtInt nAMIDE = amideBondList.size();
  cout << nAMIDE << " amide bonds:" << endl;
  for (RbtBondListConstIter bIter = amideBondList.begin(); bIter != amideBondList.end(); bIter++) {
    RbtAtomPtr spAtom2 = (*bIter)->GetAtom1Ptr();
    RbtAtomPtr spAtom3 = (*bIter)->GetAtom2Ptr();

    //Find N and C
    RbtAtomPtr spN = (bIsN_TRI(spAtom2)) ? spAtom2 : spAtom3;
    RbtAtomPtr spC = (bIsC_SP2(spAtom2)) ? spAtom2 : spAtom3;

    //Find bonded O_SP2 and H (if any)
    RbtAtomList bondedAtomsN = Rbt::GetBondedAtomList(spN);
    RbtAtomList bondedAtomsC = Rbt::GetBondedAtomList(spC);
    RbtAtomListIter oIter = Rbt::FindAtom(bondedAtomsC,bIsO_SP2);
    RbtAtomListIter hIter = Rbt::FindAtom(bondedAtomsN,bIsH);

    if ( (oIter != bondedAtomsC.end()) && (hIter != bondedAtomsN.end())) {
      RbtAtomPtr spH = (*hIter);
      RbtAtomPtr spO = (*oIter);
      RbtDouble phi = Rbt::BondDihedral(spH,spN,spC,spO);
      cout << spH->GetAtomName() << "(" << spH->GetFFType() << ") - "
	   << spN->GetAtomName() << "(" << spN->GetFFType() << ") - "
	   << spC->GetAtomName() << "(" << spC->GetFFType() << ") - "
	   << spO->GetAtomName() << "(" << spO->GetFFType() << "): "
	   << phi << " deg" << endl;
      RbtDouble deltaPhi = 180.0 - phi;
      cout << "Rotating bond by " << deltaPhi << " deg" << endl;
      spModel->RotateBond(*bIter,deltaPhi);
      phi = Rbt::BondDihedral(spH,spN,spC,spO);
      cout << "Dihedral is now " << phi << " deg" << endl;
    }
    else {
      cout << spN->GetAtomName() << "(" << spN->GetFFType() << ") - "
	   << spC->GetAtomName() << "(" << spC->GetFFType() << "): TERTIARY" << endl;
    }
  }
}

int main(int argc,char* argv[])
{
  cout.setf(ios_base::left,ios_base::adjustfield);

  //Strip off the path to the executable, leaving just the file name
  RbtString strExeName(argv[0]);
  RbtString::size_type i = strExeName.rfind("/");
  if (i != RbtString::npos)
    strExeName.erase(0,i+1);

  //Print a standard header
  Rbt::PrintStdHeader(cout,strExeName+EXEVERSION);

  //Default values for optional arguments
  RbtString strInputSDFile;
  RbtString strOutputSDFile;
  RbtBool bPosIonise(false);
  RbtBool bNegIonise(false);
  RbtBool bImplH(true);//if true, read only polar hydrogens from SD file, else read all H's present
  RbtBool bCheckAmides(false);//if true, sets all 2ndry amide bonds to trans
  RbtBool bList(false);

  //Display brief help message if no args
  if (argc == 1) {
    cout << endl << strExeName << " - output interaction center info for ligands in SD file (with optional autoionisation)" << endl;
    cout << endl << "Usage:\t" << strExeName << " -i<InputSDFile> [-o<OutputSDFile>] [-ap] [-an] [-allH]" << endl;
    cout << endl << "Options:\t-i<InputSDFile> - input ligand SD file" << endl;
    cout << "\t\t-o<OutputSDFile> - output SD file with descriptors (default=no output)" << endl;
    cout << "\t\t-ap - protonate all neutral amines, guanidines, imidazoles (default=disabled)" << endl;
    cout << "\t\t-an - deprotonate all carboxylic, sulphur and phosphorous acid groups (default=disabled)" << endl;
    cout << "\t\t-allH - read all hydrogens present (default=polar hydrogens only)" << endl;
    cout << "\t\t-tr - rotate all 2ndry amides to trans (default=leave alone)" << endl;
    cout << "\t\t-l - verbose listing of ligand atoms and rotable bonds (default = compact table format)" << endl;
    return 1;
  }

  //Check command line arguments
  cout << endl << "Command line args:" << endl;
  for (RbtInt iarg = 1; iarg < argc; iarg++) {
    cout << argv[iarg];
    RbtString strArg(argv[iarg]);
    if (strArg.find("-i")==0)
      strInputSDFile = strArg.substr(2);
    else if (strArg.find("-o")==0)
      strOutputSDFile = strArg.substr(2);
    else if (strArg.find("-ap")==0)
      bPosIonise = true;
    else if (strArg.find("-an")==0)
      bNegIonise = true;
    else if (strArg.find("-allH")==0)
      bImplH = false;
    else if (strArg.find("-tr")==0)
      bCheckAmides = true;
    else if (strArg.find("-l")==0)
      bList = true;
    else {
      cout << " ** INVALID ARGUMENT" << endl;
      return 1;
    }
    cout << endl;
  }
  cout << endl;

  RbtBool bWriteLigand = !strOutputSDFile.empty();

  if (!bList) {
    cout << setw(8) << "RECORD"
	 << setw(30)<< "NAME"
	 << setw(8) << "#ATOMS"
	 << setw(8) << "#BONDS"
	 << setw(8) << "MW"
	 << setw(8) << "#ROT"
	 << setw(8) << "#HBD"
	 << setw(8) << "#HBA"
	 << setw(8) << "#AROM"
	 << setw(8) << "#GUAN"
	 << setw(8) << "#LIPOC"
	 << endl;
  }

  try {
    //Prepare the SD file sink for saving the autoionised structures for each ligand
    RbtMolecularFileSinkPtr spMdlFileSink;
    if (bWriteLigand)
      spMdlFileSink = RbtMolecularFileSinkPtr(new RbtMdlFileSink(strOutputSDFile,RbtModelPtr()));

    ///////////////////////////////////
    //MAIN LOOP OVER LIGAND RECORDS
    ///////////////////////////////////
    RbtMolecularFileSourcePtr spMdlFileSource(new RbtMdlFileSource(strInputSDFile,bPosIonise,bNegIonise,bImplH));
    //DM 4 June 1999 - only read the largest segment
    //spMdlFileSource->SetSegmentFilterMap(Rbt::ConvertStringToSegmentMap("H"));

    for (RbtInt nRec=1; spMdlFileSource->FileStatusOK(); spMdlFileSource->NextRecord(), nRec++) {
      RbtError molStatus = spMdlFileSource->Status();
      if (!molStatus.isOK()) {
	cout << molStatus << endl;
	continue;
      }

      //DM 26 Jul 1999 - only read the largest segment (guaranteed to be called H)
      spMdlFileSource->SetSegmentFilterMap(Rbt::ConvertStringToSegmentMap("H"));

      //Create the ligand model
      RbtModelPtr spModel(new RbtModel(spMdlFileSource));
      RbtString strModelName = spModel->GetName();

      if (bList) {
	RbtAtomList atomList = spModel->GetAtomList();
	cout << endl << "Model = " << spModel->GetName() << endl;
	for (RbtAtomListConstIter iter = atomList.begin(); iter != atomList.end(); iter++) {
	  cout << (**iter) << endl;
	}
	//continue;
      }
	
      //DM 11 Nov 1999
      RbtDouble dMolWt = spModel->GetTotalAtomicMass();

      //Determine the interaction center counts
      RbtInt nAtoms= spModel->GetNumAtoms();

      RbtInt nBonds = spModel->GetNumBonds();
      RbtInt nSegs = spModel->GetNumSegments();
      RbtAtomList atomList = spModel->GetAtomList();
      RbtBondList bondList = spModel->GetBondList();


      Rbt::isHybridState_eq bIsArom(RbtAtom::AROM);
      RbtInt nLipoC = Rbt::GetNumAtoms(atomList,Rbt::isAtomLipophilic());
      RbtInt nAromAtoms = Rbt::GetNumAtoms(atomList,bIsArom);//# aromatic atoms
      RbtInt nNHBD = Rbt::GetNumAtoms(atomList,Rbt::isAtomHBondDonor());
      RbtInt nMetal = Rbt::GetNumAtoms(atomList,Rbt::isAtomMetal());
      RbtInt nGuan = Rbt::GetNumAtoms(atomList,Rbt::isAtomGuanidiniumCarbon());
      RbtInt nNHBA = Rbt::GetNumAtoms(atomList,Rbt::isAtomHBondAcceptor());

      //Total +ve and -ve charges
      RbtDouble posChg(0.0);
      RbtDouble negChg(0.0);

      for (RbtAtomListConstIter iter = atomList.begin(); iter != atomList.end(); iter++) {
	RbtDouble chg = (*iter)->GetGroupCharge();
	if (chg > 0.0) {
	  posChg += chg;
	}
	else if (chg < 0.0) {
	  negChg += chg;
	}
      }

      //Count aromatic rings
      RbtAtomListList ringLists = spModel->GetRingAtomLists();
      RbtInt nAromRings(0);//# aromatic rings
      for (RbtAtomListListIter rIter = ringLists.begin(); rIter != ringLists.end(); rIter++) {
	if (Rbt::GetNumAtoms(*rIter,Rbt::isPiAtom()) == (*rIter).size())
	  nAromRings++;
      }

      RbtBondList rotatableBondList = Rbt::GetBondList(bondList,Rbt::isBondRotatable());
      //DM 24 Sep 2001 - remove bonds to terminal NH3+
      rotatableBondList = Rbt::GetBondList(rotatableBondList,std::not1(Rbt::isBondToNH3()));
      //DM 19 Sep 2002 - remove bonds to terminal OH
      rotatableBondList = Rbt::GetBondList(rotatableBondList,std::not1(Rbt::isBondToOH()));
      RbtInt nROT = rotatableBondList.size();

      if (bList) {
	cout << nROT << " rotatable bonds:" << endl;
	for (RbtBondListConstIter bIter = rotatableBondList.begin(); bIter != rotatableBondList.end(); bIter++) {
	  RbtAtomPtr spAtom1 = (*bIter)->GetAtom1Ptr();
	  RbtAtomPtr spAtom2 = (*bIter)->GetAtom2Ptr();
	  cout << spAtom1->GetAtomName() << "(" << spAtom1->GetFFType() << ") - "
	       << spAtom2->GetAtomName() << "(" << spAtom2->GetFFType() << ")" << endl;
	}
      }
      else {
	//Output the ligand info
	cout << setw(8) << nRec
	     << setw(30)<< strModelName.c_str()
	     << setw(8) << nAtoms
	     << setw(8) << nBonds
	     << setw(8) << dMolWt
	     << setw(8) << nROT
	     << setw(8) << nNHBD
	     << setw(8) << nNHBA
	     << setw(8) << nAromRings
	     << setw(8) << nGuan
	     << setw(8) << nLipoC
	     << endl;
      }
      
      //Dump model to MdlFileSink
      if (bWriteLigand) {
	if (bCheckAmides) {
	  CheckAmideBonds(spModel);//DM 7 June 1999 - check and fix amide bonds
	}
	//DM 25 Jul 2002 - do not change ligand coords
	//spModel->AlignPrincipalAxes();//Align principal axes with Cartesian axes

	//DM 11 Nov 1999 Set Data fields for mol.wt, #interaction centers etc
	//NB Molecular weight is calculated for the fragment actually read
	//i.e. largest fragment present if multifragment record.
	spModel->SetDataValue("LIG_MW",RbtVariant(dMolWt));
	spModel->SetDataValue("LIG_NATOMS",RbtVariant(nAtoms));
	spModel->SetDataValue("LIG_NLIPOC",RbtVariant(nLipoC));
	spModel->SetDataValue("LIG_NAROMATOMS",RbtVariant(nAromAtoms));
	spModel->SetDataValue("LIG_NAROMRINGS",RbtVariant(nAromRings));
	spModel->SetDataValue("LIG_NHBD",RbtVariant(nNHBD));
	spModel->SetDataValue("LIG_NMETAL",RbtVariant(nMetal));
	spModel->SetDataValue("LIG_NGUAN",RbtVariant(nGuan));
	spModel->SetDataValue("LIG_NHBA",RbtVariant(nNHBA));
	spModel->SetDataValue("LIG_NROT",RbtVariant(nROT));

	//Atom type percentages
	spModel->SetDataValue("LIG_PERC_LIPOC",RbtVariant(100.0*nLipoC/nAtoms));
	spModel->SetDataValue("LIG_PERC_AROMATOMS",RbtVariant(100.0*nAromAtoms/nAtoms));
	spModel->SetDataValue("LIG_PERC_HBD",RbtVariant(100.0*nNHBD/nAtoms));
	spModel->SetDataValue("LIG_PERC_METAL",RbtVariant(100.0*nMetal/nAtoms));
	spModel->SetDataValue("LIG_PERC_GUAN",RbtVariant(100.0*nGuan/nAtoms));
	spModel->SetDataValue("LIG_PERC_HBA",RbtVariant(100.0*nNHBA/nAtoms));

	//Charges
	spModel->SetDataValue("LIG_POS_CHG",RbtVariant(posChg));
	spModel->SetDataValue("LIG_NEG_CHG",RbtVariant(negChg));
	spModel->SetDataValue("LIG_TOT_CHG",RbtVariant(posChg+negChg));

	spMdlFileSink->SetModel(spModel);
	spMdlFileSink->Render();
      }

    }
    //END OF MAIN LOOP OVER LIGAND RECORDS
    ////////////////////////////////////////////////////
  }
  catch (RbtError& e) {
    cout << e << endl;
  }
  catch (...) {
    cout << "Unknown exception" << endl;
  }

  _RBTOBJECTCOUNTER_DUMP_(cout)

  return 0;
}
