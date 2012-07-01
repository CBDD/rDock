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

//Calculates vdW grids for use by RbtVdwGridSF scoring function class

#include <fstream>
#include <iomanip>
#include <cstring>
#include "RbtBiMolWorkSpace.h"
#include "RbtParameterFileSource.h"
#include "RbtPRMFactory.h"
#include "RbtSFFactory.h"
#include "RbtRealGrid.h"
#include "RbtTriposAtomType.h"

const RbtString EXEVERSION = " ($Id: //depot/dev/client3/rdock/2006.1/src/exe/rbcalcgrid.cxx#3 $)";
const RbtString _ROOT_SF = "SCORE";

//Creates list of probe models
//NOTE: MUST BE IN ORDER OF ASCENDING RADII
RbtModelList CreateProbes() {
  RbtModelList probes;
  RbtTriposAtomTypeList atomTypes;
  atomTypes.push_back(RbtTriposAtomType::UNDEFINED);
  atomTypes.push_back(RbtTriposAtomType::Br);
  atomTypes.push_back(RbtTriposAtomType::C_cat);  
  atomTypes.push_back(RbtTriposAtomType::C_1);  
  atomTypes.push_back(RbtTriposAtomType::C_1_H1);  
  atomTypes.push_back(RbtTriposAtomType::C_2);  
  atomTypes.push_back(RbtTriposAtomType::C_2_H1);  
  atomTypes.push_back(RbtTriposAtomType::C_2_H2);  
  atomTypes.push_back(RbtTriposAtomType::C_3);  
  atomTypes.push_back(RbtTriposAtomType::C_3_H1);  
  atomTypes.push_back(RbtTriposAtomType::C_3_H2);  
  atomTypes.push_back(RbtTriposAtomType::C_3_H3);  
  atomTypes.push_back(RbtTriposAtomType::C_ar);  
  atomTypes.push_back(RbtTriposAtomType::C_ar_H1);  
  atomTypes.push_back(RbtTriposAtomType::Cl);  
  atomTypes.push_back(RbtTriposAtomType::F);  
  atomTypes.push_back(RbtTriposAtomType::H);  
  atomTypes.push_back(RbtTriposAtomType::H_P);  
  atomTypes.push_back(RbtTriposAtomType::I);  
  atomTypes.push_back(RbtTriposAtomType::N_1);  
  atomTypes.push_back(RbtTriposAtomType::N_2);  
  atomTypes.push_back(RbtTriposAtomType::N_3);  
  atomTypes.push_back(RbtTriposAtomType::N_4);  
  atomTypes.push_back(RbtTriposAtomType::N_am);  
  atomTypes.push_back(RbtTriposAtomType::N_ar);  
  atomTypes.push_back(RbtTriposAtomType::N_pl3);  
  atomTypes.push_back(RbtTriposAtomType::O_2);  
  atomTypes.push_back(RbtTriposAtomType::O_3);  
  atomTypes.push_back(RbtTriposAtomType::O_co2); 
  atomTypes.push_back(RbtTriposAtomType::P_3);  
  atomTypes.push_back(RbtTriposAtomType::S_2);  
  atomTypes.push_back(RbtTriposAtomType::S_3);  
  atomTypes.push_back(RbtTriposAtomType::S_o);  
  atomTypes.push_back(RbtTriposAtomType::S_o2);  
  for (RbtTriposAtomTypeListConstIter iter = atomTypes.begin(); iter != atomTypes.end(); iter++) {
    RbtAtomList atomList;
    RbtBondList bondList;
    RbtAtomPtr spAtom(new RbtAtom(1));
    spAtom->SetTriposType(*iter);
    spAtom->SetAtomicMass(12.0);//Mass is irrelevant as we are not rotating models around COM
    atomList.push_back(spAtom);
    probes.push_back(new RbtModel(atomList,bondList));
  }
  return probes;  	
}

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
  RbtString strSuffix(".grd");
  RbtString strReceptorPrmFile;//Receptor param file
  RbtString strSFFile("calcgrid_attr.prm");//Scoring function file
  RbtDouble gs(0.5);//grid step
  RbtDouble border(1.0);//grid border around docking site
  
  //Brief help message
  if (argc == 1) {
    cout << endl << "rbcalcgrid - calculates vdw grids for each atom type" << endl;
    cout << endl << "Usage:\trbcalcgrid -o<OutputRoot> -r<ReceptorPrmFile> -p<SFPrmFile> [-g<GridStep>]" << endl;
    cout << endl << "Options:\t-o<OutputSuffix> - suffix for grid (.grd IS required)" << endl;
    cout << "\t\t-r<ReceptorPrmFile> - receptor param file (contains active site params)" << endl;
    cout << "\t\t-p<SFPrmFile> - scoring function param file (either calcgrid_attr.prm or calcgrid_repul.prm)" << endl;
    cout << "\t\t-g<GridStep> - grid step (default=0.5A)" << endl;
    cout << "\t\t-b<Border> - grid border around docking site (default=1.0A)" << endl;
    return 1;
  }

  //Check command line arguments
  cout << endl << "Command line args:" << endl;
  for (RbtInt iarg = 1; iarg < argc; iarg++) {
    cout << argv[iarg];
    RbtString strArg(argv[iarg]);
    if (strArg.find("-o")==0)
      strSuffix = strArg.substr(2);
    else if (strArg.find("-r")==0)
      strReceptorPrmFile = strArg.substr(2);
    else if (strArg.find("-p")==0)
      strSFFile = strArg.substr(2);
    else if (strArg.find("-g")==0) {
      RbtString strGridStep = strArg.substr(2);
      gs = atof(strGridStep.c_str());
    }
    else if (strArg.find("-b")==0) {
      RbtString strBorder = strArg.substr(2);
      border = atof(strBorder.c_str());
    }
    else {
      cout << " ** INVALID ARGUMENT" << endl;
      return 1;
    }
    cout << endl;
  }

  cout << endl;

  try {
    //Create a bimolecular workspace
    RbtBiMolWorkSpacePtr spWS(new RbtBiMolWorkSpace());
    //Set the workspace name to the root of the receptor .prm filename
    RbtStringList componentList = Rbt::ConvertDelimitedStringToList(strReceptorPrmFile,".");
    RbtString wsName = componentList.front();
    spWS->SetName(wsName);
    
    //Read the receptor parameter file
    RbtParameterFileSourcePtr spRecepPrmSource(new RbtParameterFileSource(Rbt::GetRbtFileName("data/receptors",strReceptorPrmFile)));
    cout << endl << "RECEPTOR:" << endl << spRecepPrmSource->GetFileName() << endl << spRecepPrmSource->GetTitle() << endl;

    //Read the scoring function file
    RbtParameterFileSourcePtr spSFSource(new RbtParameterFileSource(Rbt::GetRbtFileName("data/sf",strSFFile)));
    RbtSFFactoryPtr spSFFactory(new RbtSFFactory());//Factory class for scoring functions
    RbtSFAggPtr spSF(spSFFactory->CreateAggFromFile(spSFSource,_ROOT_SF));//Root SF aggregate

    //Register the scoring function with the workspace
    //Dump details to cout
    spWS->SetSF(spSF);
    cout << endl << "SCORING FUNCTION DETAILS:" << endl << *spSF << endl;
    
    //Create the receptor model from the file names in the receptor parameter file
    spRecepPrmSource->SetSection();
    RbtPRMFactory prmFactory(spRecepPrmSource);
    RbtModelPtr spReceptor = prmFactory.CreateReceptor();
    //Trap multiple receptor conformations here: this SF does not support them yet
    RbtBool bEnsemble = (spReceptor->GetNumSavedCoords() > 1);
    if (bEnsemble) {
      RbtString message("rbcalcgrid does not support multiple receptor conformations yet");
      throw RbtInvalidRequest(_WHERE_,message);
    }

    //Read docking site from file and register with workspace
    RbtString strASFile = spWS->GetName()+".as";
    RbtString strInputFile = Rbt::GetRbtFileName("data/grids",strASFile);
    //DM 26 Sep 2000 - ios_base::binary is invalid with IRIX CC
#if defined(__sgi) && !defined(__GNUC__)
    ifstream istr(strInputFile.c_str(),ios_base::in);
#else
    ifstream istr(strInputFile.c_str(),ios_base::in|ios_base::binary);
#endif
    RbtDockingSitePtr spDS(new RbtDockingSite(istr));
    istr.close();
    spWS->SetDockingSite(spDS);

    //Register receptor with workspace
    spWS->SetReceptor(spReceptor);
    cout << endl << "DOCKING SITE" << endl << (*spDS) << endl;
    
    //Create a grid covering the docking site, plus user-defined border
    RbtCoord minCoord = spDS->GetMinCoord()-border;
    RbtCoord maxCoord = spDS->GetMaxCoord()+border;
    RbtVector recepExtent = maxCoord-minCoord;
    RbtVector gridStep(gs,gs,gs);
    RbtUInt nX = int(recepExtent.x/gridStep.x)+1;
    RbtUInt nY = int(recepExtent.y/gridStep.y)+1;
    RbtUInt nZ = int(recepExtent.z/gridStep.z)+1;
    cout << "Constructing grid of size " << nX << " x " << nY << " x " << nZ << endl;
    RbtRealGridPtr spGrid(new RbtRealGrid(minCoord,gridStep,nX,nY,nZ));
    float* gridData = spGrid->GetGridData();
    
    //Create probes
    RbtModelList probes = CreateProbes();
    
    //Open output file
    RbtString strOutputFile(spWS->GetName()+strSuffix);
#if defined(__sgi) && !defined(__GNUC__)
    ofstream ostr(strOutputFile.c_str(),ios_base::out|ios_base::trunc);
#else
    ofstream ostr(strOutputFile.c_str(),ios_base::out|ios_base::binary|ios_base::trunc);
#endif
    //Write header string (RbtVdwGridSF)
    const char* const header = "RbtVdwGridSF";
    RbtInt length = strlen(header);
    Rbt::WriteWithThrow(ostr, (const char*) &length, sizeof(length));
    Rbt::WriteWithThrow(ostr, header, length);
    //Write number of grids
    RbtInt nGrids = probes.size();
    Rbt::WriteWithThrow(ostr, (const char*) &nGrids, sizeof(nGrids));

    //Store regular pointers to avoid smart pointer dereferencing overheads
    RbtRealGrid* pGrid(spGrid);
    RbtSFAgg* pSF(spSF);
    RbtTriposAtomType triposType;
    //Main loop over each probe model
    for (RbtModelListConstIter mIter = probes.begin(); mIter != probes.end(); mIter++) {
      RbtModelPtr spLigand(*mIter);
      RbtAtom* pAtom = spLigand->GetAtomList().front();
      RbtTriposAtomType::eType atomType = pAtom->GetTriposType();
      RbtString strType = triposType.Type2Str(atomType);
      cout << "Atom type=" << strType << endl;
      //Register ligand with workspace
      spWS->SetLigand(spLigand);
      pGrid->SetAllValues(0.0);
      //Loop over all grid coords and calculate the score at each position		
      for (RbtUInt i = 0; i < spGrid->GetN(); i++ ) {
	pAtom->SetCoords(pGrid->GetCoord(i));
	gridData[i] = pSF->Score();
      }
      //Write the atom type string to the grid file, before the grid itself
      const char* const szType = strType.c_str();
      RbtInt l = strlen(szType);
      Rbt::WriteWithThrow(ostr, (const char*) &l, sizeof(l));
      Rbt::WriteWithThrow(ostr, szType, l);	  
      pGrid->Write(ostr);
    }
    ostr.close();		
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
