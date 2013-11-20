//Ancillary Software from Vernalis
//$Id: //depot/dev/client3/rdock/2006.1/src/exe/rbtether.cxx#5 $

//Filters a ligand library with a SMARTS query and writes down the
//corresponding tethered atoms in the sd files

#include <strstream>
#include <iomanip>

#include "RbtBiMolWorkSpace.h"
#include "RbtMdlFileSource.h"
#include "RbtMdlFileSink.h"
#include "RbtParameterFileSource.h"
#include "RbtPRMFactory.h"
#include "RbtSmarts.h"
const RbtString EXEVERSION = 
      " ($Id: //depot/dev/client3/rdock/2006.1/src/exe/rbtether.cxx#5 $)";

void print_atoms(RbtAtomList &atoms, ostrstream &ost);
/////////////////////////////////////////////////////////////////////
// MAIN PROGRAM STARTS HERE
/////////////////////////////////////////////////////////////////////

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

  //Command line arguments and default values
  RbtString strLigandMdlFile;
  RbtBool bOutput(false);
  RbtString strRunName;
  RbtString strReferenceSDFile;//Receptor param file
  RbtString strQuery;//Smarts query
  RbtBool bPosIonise(false);
  RbtBool bNegIonise(false);
  RbtBool bImplH(true);

  //Brief help message
  if (argc == 1) 
  {
    cout << endl << "rbtether - marks the tethered atoms for use in docking" 
         << endl;
    cout << endl << "Usage:\trbtether -i<InputSDFile> -o<OutputRoot>" 
         << "-r<ReferenceSDFile> -s<Smarts Query>"
         << " [-ap] [-an] [-allH]" << endl;
    cout << endl << "Options:\t-i<InputSDFile> - input ligand SD file" << endl;
    cout << "\t\t-o<OutputRoot> - root name for output SD file(s)" 
         << "(.sd suffix not required)" << endl;
    cout << "\t\t-r<ReferenceSDFile> - reference SD file "
         << "(with the correct coordinates for the tethered atoms)" << endl;
    cout << "\t\t-s<Smarts Query> - SMARTS query " 
         << "(specifying the atoms to be tethered)" << endl;
    cout << "\t\t-ap - protonate all neutral amines, guanidines, "
         << "imidazoles (default=disabled)" << endl;
    cout << "\t\t-an - deprotonate all carboxylic, sulphur and phosphorous "
         << "acid groups (default=disabled)" << endl;
    cout << "\t\t-allH - read all hydrogens present (default=polar "
         << "hydrogens only)" << endl;
    return 1;
  }

  //Check command line arguments
  cout << endl << "Command line args:" << endl;
  for (RbtInt iarg = 1; iarg < argc; iarg++) 
  {
    cout << argv[iarg];
    RbtString strArg(argv[iarg]);
    if (strArg.find("-i")==0)
      strLigandMdlFile = strArg.substr(2);
    else if (strArg.find("-o")==0) 
    {
      strRunName = strArg.substr(2);
      bOutput = true;
    }
    else if (strArg.find("-r")==0)
      strReferenceSDFile = strArg.substr(2);
    else if (strArg.find("-s")==0)
      strQuery = strArg.substr(2);
    else if (strArg.find("-ap")==0)
      bPosIonise = true;
    else if (strArg.find("-an")==0)
      bNegIonise = true;
    else if (strArg.find("-allH")==0)
      bImplH = false;
    else 
    {
      cout << " ** INVALID ARGUMENT" << endl;
      return 1;
    }
    cout << endl;
  }

  cout << endl;

  //DM 20 Apr 1999 - set the auto-ionise flags
  if (bPosIonise)
    cout << "Automatically protonating positive ionisable groups " 
         << "(amines, imidazoles, guanidines)" << endl;
  if (bNegIonise)
    cout << "Automatically deprotonating negative ionisable groups "
         << "(carboxylic acids, phosphates, sulphates, sulphonates)" << endl;
  if (bImplH)
    cout << "Reading polar hydrogens only from ligand SD file" << endl;
  else
    cout << "Reading all hydrogens from ligand SD file" << endl;

  try {
    //Read the reference SD file
    RbtMolecularFileSourcePtr  spReferenceSD(new RbtMdlFileSource
					     (strReferenceSDFile, bPosIonise, bNegIonise, bImplH));
    RbtModelPtr spReferenceMdl(new RbtModel(spReferenceSD));
    RbtString strSmiles;
    RbtAtomListList tetherAtoms = DT::QueryModel(spReferenceMdl, strQuery, strSmiles);
    cout << endl << "REFERENCE SD FILE: " << strReferenceSDFile << endl;
    if (strQuery.empty()) {
      cout << "No SMARTS query defined - tethering to whole of reference ligand" << endl;
      strQuery = strSmiles;
    }
    cout << "SMILES: " << strSmiles << endl;
    cout << "SMARTS: " << strQuery << endl;
    for (RbtAtomListListIter alli = tetherAtoms.begin(); alli != tetherAtoms.end(); alli++) {
      cout << "Path: ";
      for (RbtAtomListIter ali = alli->begin(); ali != alli->end(); ali++) {
	cout << (*ali)->GetAtomName() << "\t";
      }
      cout << endl;
    }
    if ((argc == 2) && !strReferenceSDFile.empty()) {
      // smile expression found for the reference file. Nothing
      // else to be done
      exit(0);
    }
    if (tetherAtoms.empty()) {
      throw RbtError(_WHERE_,
		     "SMARTS query should have at least one match in the reference SD file");
    }
    if (tetherAtoms.size() > 1) {
      //throw RbtError(_WHERE_,
      cout << "More than one SMARTS match found with the reference SD file"
	   << endl;
      cout << "Will only tether to the first matching path..." << endl;
    }
    
    //DM 1 Jul 2002 - calculate principal axes and COM of tethered atoms
    //for the reference structure. Use to prealign each ligand with the reference.
    RbtAtomList tetheredAtomList(*(tetherAtoms.begin()));
    RbtPrincipalAxes refAxes = Rbt::GetPrincipalAxes(tetheredAtomList);
    
    //Prepare the SD file sink for saving the docked conformations 
    //for each ligand
    if (bOutput) {
      //writing down the reference sd file for the query
      ostrstream ost;
      print_atoms(tetheredAtomList, ost);
      RbtVariant vTetherAtoms(ost.str());
      RbtMolecularFileSinkPtr spRefMdlFileSink(new RbtMdlFileSink(strRunName+"_reference.sd", 
								  spReferenceMdl));
      spReferenceMdl->SetDataValue("TETHERED ATOMS", vTetherAtoms);
      spRefMdlFileSink->Render();
    }
    
    //DM 18 May 1999
    //Variants describing the library version, exe version, 
    //parameter file, and current directory
    //Will be stored in the ligand SD files
    RbtVariant vLib(Rbt::GetProduct()+" ("+Rbt::GetVersion()
		    +", Build"+Rbt::GetBuild()+")");
    RbtVariant vExe(strExeName+EXEVERSION);
    RbtVariant vDir(Rbt::GetCurrentDirectory());
    ///////////////////////////////////
      //MAIN LOOP OVER LIGAND RECORDS
      ///////////////////////////////////
    //DM 20 Apr 1999 - add explicit bPosIonise and bNegIonise flags 
    //to MdlFileSource constructor
    RbtMolecularFileSinkPtr spMdlFileSink
                    (new RbtMdlFileSink(strRunName+".sd",RbtModelPtr()));
    RbtMolecularFileSourcePtr spMdlFileSource
         (new RbtMdlFileSource(strLigandMdlFile,bPosIonise,bNegIonise,bImplH));
    for (RbtInt nRec=1; spMdlFileSource->FileStatusOK(); 
                        spMdlFileSource->NextRecord(), nRec++) 
    {
      cout.setf(ios_base::left,ios_base::adjustfield);
      cout << endl
           << "**************************************************" << endl
           << "RECORD #" << nRec << endl;
      RbtError molStatus = spMdlFileSource->Status();
      if (!molStatus.isOK()) 
      {
        cout << endl << molStatus << endl
             << "************************************************" << endl;
        continue;
      }
      
      //DM 26 Jul 1999 - only read the largest segment 
      //(guaranteed to be called H)
      spMdlFileSource->SetSegmentFilterMap(Rbt::ConvertStringToSegmentMap("H"));
      
      if (spMdlFileSource->isDataFieldPresent("Name"))
        cout << "NAME:   " << spMdlFileSource->GetDataValue("Name") << endl;
      if (spMdlFileSource->isDataFieldPresent("REG_Number"))
        cout << "REG_Num:" << spMdlFileSource->GetDataValue("REG_Number") 
             << endl;
      
      //Create the ligand model
      RbtModelPtr spLigand(new RbtModel(spMdlFileSource));
      RbtString strMolName = spLigand->GetName();
      
      tetherAtoms = DT::QueryModel(spLigand, strQuery, strSmiles);
      cout << "SMILES: " << strSmiles << endl;
      cout << "SMARTS: " << strQuery << endl;
      if (tetherAtoms.empty()) {
	cout << "No SMARTS match, structure not output" << endl;
	continue;  // ligand filtered out
      }
      RbtAtomList ligAtomList = spLigand->GetAtomList();
      for (RbtAtomListListIter alli = tetherAtoms.begin(); alli != tetherAtoms.end() ; alli++) {
	cout << "Path: ";
	for (RbtAtomListIter ali = alli->begin() ; ali != alli->end() ; ali++) {
	  cout << (*ali)->GetAtomName() << "\t";
	}
	cout << endl;
	/////////////////////////////////
	  //DM 1 Jul 2002 - prealign each ligand with the reference tether
	  //Calculate quat needed to overlay principal axes of tethered fragment with reference tether
	  RbtPrincipalAxes prAxes = Rbt::GetPrincipalAxes(*alli);
	  RbtQuat q = Rbt::GetQuatFromAlignAxes(prAxes, refAxes);
	  //Apply the transformation to the whole ligand
	  //1. Translate COM of tethered fragment to the origin
	  //2. Perform the rotation
	  //3. Translate COM of tethered fragment to reference tether COM
	  std::for_each(ligAtomList.begin(),ligAtomList.end(),Rbt::TranslateAtom(-prAxes.com));
	  std::for_each(ligAtomList.begin(),ligAtomList.end(),Rbt::RotateAtomUsingQuat(q));
	  std::for_each(ligAtomList.begin(),ligAtomList.end(),Rbt::TranslateAtom(refAxes.com));
	  ///////////////////////////////////
	    
	    ostrstream ost;
	    print_atoms(*alli, ost);
	    RbtVariant vTetherAtoms(ost.str());
	    //DM 18 May 1999 - store run info in model data
	    //Clear any previous Rbt.* data fields
	    spLigand->ClearAllDataFields("Rbt.");
	    spLigand->ClearAllDataFields("TETHERED ATOMS");
	    spLigand->SetDataValue("Rbt.Library",vLib);
	    spLigand->SetDataValue("Rbt.Executable",vExe);
	    spLigand->SetDataValue("Rbt.Current_Directory",vDir);
	    spLigand->SetDataValue("TETHERED ATOMS", vTetherAtoms);
	    spMdlFileSink->SetModel(spLigand);
	    spMdlFileSink->Render();
      }
    }
    //END OF MAIN LOOP OVER LIGAND RECORDS
    ////////////////////////////////////////////////////
    cout << endl << "END OF RUN" << endl;
  }
  catch (RbtError& e) 
    {
      cout << e << endl;
    }
  catch (...) 
    {
      cout << "Unknown exception" << endl;
    }
  
  _RBTOBJECTCOUNTER_DUMP_(cout)
    
    return 0;
}

void print_atoms(RbtAtomList &atoms, ostrstream &ost)
{
  ost.clear();
  for (RbtInt iter = 0 ; iter < atoms.size() ; iter++)
    {
      ost << atoms[iter]->GetAtomId();
      //    if ((((iter + 1) % 35) == 0) || (iter == (atoms.size() - 1)))
      if (iter == (atoms.size() - 1))
	ost << "\n";
      else if (((iter + 1) % 35) == 0)
      ost << ",\n";
      else
	ost << ",";
    }
  ost << ends;
}
