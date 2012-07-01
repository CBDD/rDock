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

//Main docking application
#include <iomanip>
using std::setw;

#include <popt.h>		// for command-line parsing
#include <errno.h>

#include "RbtBiMolWorkSpace.h"
#include "RbtMdlFileSource.h"
#include "RbtMdlFileSink.h"
#include "RbtCrdFileSink.h"
#include "RbtParameterFileSource.h"
#include "RbtPRMFactory.h"
#include "RbtSFFactory.h"
#include "RbtTransformFactory.h"
#include "RbtRand.h"
#include "RbtModelError.h"
#include "RbtDockingError.h"
#include "RbtLigandError.h"
#include "RbtFilter.h"
#include "RbtSFRequest.h"
#include "RbtFileError.h"

const RbtString EXEVERSION = " ($Id: //depot/dev/client3/rdock/2006.1/src/exe/rbdock.cxx#4 $)";
//Section name in docking prm file containing scoring function definition
const RbtString _ROOT_SF = "SCORE";
const RbtString _RESTRAINT_SF = "RESTR";
const RbtString _ROOT_TRANSFORM = "DOCK";

void PrintUsage(void)
{
  cout << endl << "Usage:" << endl;
  cout << "rbdock -i <sdFile> -o <outputRoot> -r <recepPrmFile> -p <protoPrmFile> [-n <nRuns>] [-ap] [-an] [-allH]" << endl;
  cout << "       [-t <targetScore|targetFilterFile>] [-c] [-T <traceLevel>] [-s <rndSeed>]" << endl;
  cout << endl << "Options:\t-i <sdFile> - input ligand SD file" << endl;
  cout << "\t\t-o <outputRoot> - root name for output file(s)" << endl;
  cout << "\t\t-r <recepPrmFile> - receptor parameter file " << endl;
  cout << "\t\t-p <protoPrmFile> - docking protocol parameter file" << endl;
  cout << "\t\t-n <nRuns> - number of runs/ligand (default=1)" << endl;
  cout << "\t\t-ap - protonate all neutral amines, guanidines, imidazoles (default=disabled)" << endl;
  cout << "\t\t-an - deprotonate all carboxylic, sulphur and phosphorous acid groups (default=disabled)" << endl;
  cout << "\t\t-allH - read all hydrogens present (default=polar hydrogens only)" << endl;
  cout << "\t\t-t - score threshold OR filter file name" << endl;
  cout << "\t\t-c - continue if score threshold is met (use with -t <targetScore>, default=terminate ligand)" << endl;
  cout << "\t\t-T <traceLevel> - controls output level for debugging (0 = minimal, >0 = more verbose)" << endl;
  cout << "\t\t-s <rndSeed> - random number seed (default=from sys clock)" << endl;
}

/////////////////////////////////////////////////////////////////////
// MAIN PROGRAM STARTS HERE
/////////////////////////////////////////////////////////////////////

int main(int argc, const char* argv[])
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
	RbtString	strLigandMdlFile;
	RbtBool		bOutput(false);
	RbtString	strRunName;
	RbtString	strReceptorPrmFile;		//Receptor param file
	RbtString	strParamFile;			//Docking run param file
	RbtString       strFilterFile; // Filter file
	RbtInt		nDockingRuns(0);//Init to zero, so can detect later whether user explictly typed -n

	//Params for target score
	RbtBool		bTarget(false);
	RbtBool		bStop(true);			//DM 25 May 2001 - if true, stop once target is met
	RbtBool         bDockingRuns(false); // is argument -n present?
	RbtDouble	dTargetScore(0.0);
	RbtBool         bFilter(false);

	RbtBool		bPosIonise(false);
	RbtBool		bNegIonise(false);
	RbtBool		bImplH(true);			//if true, read only polar hydrogens from SD file, else read all H's present
	RbtBool		bSeed(false);			//Random number seed (default = from system clock)
	RbtInt		nSeed(0);
	RbtBool         bTrace(false);
	RbtInt		iTrace(0);//Trace level, for debugging

	// variables for popt command-line parsing
	char 			c;						// for argument parsing
	poptContext		optCon;					// ditto

	char 			*inputFile=NULL;		// will be 'strLigandMdlFile'
	char 			*outputFile=NULL;		// will be 'strRunName'
	char 			*receptorFile=NULL;		// will be 'strReceptorPrmFile'
	char 			*protocolFile=NULL;		// will be 'strParamFile'
	char 			*strTargetScr=NULL;		// will be 'dTargetScore' 
	struct poptOption optionsTable[] = {	// command line options
		{"input",		'i',POPT_ARG_STRING|POPT_ARGFLAG_ONEDASH,&inputFile,   'i',"input file"},
		{"output",		'o',POPT_ARG_STRING|POPT_ARGFLAG_ONEDASH,&outputFile,  'o',"output file"},
		{"receptor",	'r',POPT_ARG_STRING|POPT_ARGFLAG_ONEDASH,&receptorFile,'r',"receptor file"},
		{"protocol",	'p',POPT_ARG_STRING|POPT_ARGFLAG_ONEDASH,&protocolFile,'p',"protocol file"},
		{"runs",		'n',POPT_ARG_INT   |POPT_ARGFLAG_ONEDASH,&nDockingRuns,'n',"number of runs"},
		{"trace",		'T',POPT_ARG_INT   |POPT_ARGFLAG_ONEDASH,&iTrace,'T',"trace level for debugging"},
		{"seed",		's',POPT_ARG_INT   |POPT_ARGFLAG_ONEDASH,&nSeed,       's',"random seed"},
		{"ap",			'P',POPT_ARG_NONE  |POPT_ARGFLAG_ONEDASH,0,            'P',"protonate groups"},
		{"an",			'D',POPT_ARG_NONE  |POPT_ARGFLAG_ONEDASH,0,            'D',"DEprotonate groups"},
		{"allH",        'H',POPT_ARG_NONE  |POPT_ARGFLAG_ONEDASH,0,            'H',"read all Hs"},
		{"target",      't',POPT_ARG_STRING|POPT_ARGFLAG_ONEDASH,&strTargetScr,'t',"target score"},
		{"cont",        'C',POPT_ARG_NONE  |POPT_ARGFLAG_ONEDASH,0,            'C',"continue even if target met"},
		POPT_AUTOHELP
		{NULL,0,0,NULL,0}
	};

	optCon	= poptGetContext(NULL, argc, argv, optionsTable, 0);
	poptSetOtherOptionHelp(optCon, "-r<receptor.prm> -p<protocol.prm> -i<infile> -o<outfile> [options]");
	//Brief help message
	if (argc < 2) {
		PrintUsage();
		return 1;
	}
	RbtDouble val(0.0);
	while ( (c=poptGetNextOpt(optCon))>= 0) {
		switch(c) {
			case 'P':	// protonate
				bPosIonise = true;
				break;
			case 'D':	// deprotonate
				bNegIonise = true;
				break;
			case 'H':	// all-H model
				bImplH = false;
				break;
			case 'C':
				bStop = false;
				break;
			case 't':
			  // If str can be translated to an integer, I assume is a
			  // threshold. Otherwise, I assume is the filter file name
			  char *error;
			  errno = 0;
			  val = strtod(strTargetScr, &error);
			  if (!errno && !*error)  // Is it a number?
			    {
			      dTargetScore = val;
			      bTarget = true;
			    }
			  else // Assume is the filter file name
			    {
			      bFilter = true;
			      strFilterFile = strTargetScr;
			    }
			  break;
			case 's':
				bSeed = true;
				break;
			case 'T':
				bTrace = true;
				break;
			default:
				break;
		}
	}
	cout << endl;
	poptFreeContext(optCon);
	
	// print out arguments
	// input ligand file, receptor and parameter is compulsory
	cout << endl << "Command line args:" << endl;
	if(!inputFile || !receptorFile || !protocolFile) { // if any of them is missing
		poptPrintUsage(optCon, stderr, 0);
		exit(1);
	} else {
		strLigandMdlFile	= inputFile;
		strReceptorPrmFile	= receptorFile;
		strParamFile		= protocolFile;
		cout << " -i " << strLigandMdlFile		<< endl;
		cout << " -r " << strReceptorPrmFile	<< endl;
		cout << " -p " << strParamFile			<< endl;
	}
	// output is not that important but good to have
	if(outputFile) {
		strRunName	= outputFile;
		bOutput = true;
		cout << " -o " << strRunName << endl;
	} else {
		cout << "WARNING: output file name is missing."<< endl;
	}
	// docking runs
	if(nDockingRuns >= 1) {//User typed -n explicitly, so set bDockingRuns to true
	  bDockingRuns = true;
	  cout << " -n " << nDockingRuns << endl;
	} else {
	  nDockingRuns = 1;//User didn't type -n explicitly, so fall back to the default of n=1
	  cout << " -n " << nDockingRuns << " (default) " <<endl;
	}
	if(bSeed)		// random seed (if provided)
		cout << " -s " << nSeed << endl;
	if(bTrace)		// random seed (if provided)
		cout << " -T " << iTrace << endl;
	if(bPosIonise)	// protonate 
		cout << " -ap " << endl;
	if(bNegIonise)	// deprotonate
		cout << " -an " << endl;
	if(!bImplH)		// all-H
		cout << " -allH " << endl;
	if(!bStop)		// stop after target
		cout << " -cont " << endl;
	if(bTarget)
		cout << " -t " << dTargetScore << endl;

  //BGD 26 Feb 2003 - Create filters to simulate old rbdock
  //behaviour
  ostrstream strFilter;
  if (!bFilter)
  {
    if (bTarget) // -t<TS>
    {
      if (!bDockingRuns) // -t<TS> only
      {
        strFilter << "0 1 - SCORE.INTER " << dTargetScore << endl;
      }
      else  // -t<TS> -n<N> need to check if -cont present
            // for all other cases it doesn't matter
        if (!bStop) // -t<TS> -n<N> -cont
        {
          strFilter << "1 if - SCORE.NRUNS " << (nDockingRuns - 1)
                    << " 0.0 -1.0,\n1 - SCORE.INTER " << dTargetScore << endl;
        }
        else // -t<TS> -n<N>
        {
          strFilter << "1 if - " << dTargetScore << " SCORE.INTER 0.0 " 
                    << "if - SCORE.NRUNS " << (nDockingRuns - 1) 
                    << " 0.0 -1.0,\n1 - SCORE.INTER " << dTargetScore << endl;
        }
    } // no target score, no filter
    else if (bDockingRuns) // -n<N>
    {
      strFilter << "1 if - SCORE.NRUNS " << (nDockingRuns - 1) 
                << " 0.0 -1.0,\n0";
    }
    else // no -t no -n
    {
      strFilter << "0 0\n";
    }
  }

  //DM 20 Apr 1999 - set the auto-ionise flags
  if (bPosIonise)
    cout << "Automatically protonating positive ionisable groups (amines, imidazoles, guanidines)" << endl;
  if (bNegIonise)
    cout << "Automatically deprotonating negative ionisable groups (carboxylic acids, phosphates, sulphates, sulphonates)" << endl;
  if (bImplH)
    cout << "Reading polar hydrogens only from ligand SD file" << endl;
  else
    cout << "Reading all hydrogens from ligand SD file" << endl;

  if (bTarget) {
    cout << endl << "Lower target intermolecular score = " << dTargetScore << endl;
  }

  try {
    //Create a bimolecular workspace
    RbtBiMolWorkSpacePtr spWS(new RbtBiMolWorkSpace());
    //Set the workspace name to the root of the receptor .prm filename
    RbtStringList componentList = Rbt::ConvertDelimitedStringToList(strReceptorPrmFile,".");
    RbtString wsName = componentList.front();
    spWS->SetName(wsName);
    
    //Read the docking protocol parameter file
    RbtParameterFileSourcePtr spParamSource(new RbtParameterFileSource(Rbt::GetRbtFileName("data/scripts",strParamFile)));
    //Read the receptor parameter file
    RbtParameterFileSourcePtr spRecepPrmSource(new RbtParameterFileSource(Rbt::GetRbtFileName("data/receptors",strReceptorPrmFile)));
    cout << endl << "DOCKING PROTOCOL:" << endl << spParamSource->GetFileName() << endl << spParamSource->GetTitle() << endl;
    cout << endl << "RECEPTOR:" << endl << spRecepPrmSource->GetFileName() << endl << spRecepPrmSource->GetTitle() << endl;
    
    //Create the scoring function from the SCORE section of the docking protocol prm file
    //Format is:
    //SECTION SCORE
    //    INTER    RbtInterSF.prm
    //    INTRA RbtIntraSF.prm
    //END_SECTION
    //
    //Notes:
    //Section name must be SCORE. This is also the name of the root SF aggregate
    //An aggregate is created for each parameter in the section.
    //Parameter name becomes the name of the subaggregate (e.g. SCORE.INTER)
    //Parameter value is the file name for the subaggregate definition
    //Default directory is $RBT_ROOT/data/sf
    RbtSFFactoryPtr spSFFactory(new RbtSFFactory());//Factory class for scoring functions
    RbtSFAggPtr spSF(new RbtSFAgg(_ROOT_SF));//Root SF aggregate
    spParamSource->SetSection(_ROOT_SF);
    RbtStringList sfList(spParamSource->GetParameterList());
    //Loop over all parameters in the SCORE section
    for (RbtStringListConstIter sfIter = sfList.begin(); sfIter != sfList.end(); sfIter++) {
      //sfFile = file name for scoring function subaggregate
      RbtString sfFile(Rbt::GetRbtFileName("data/sf",spParamSource->GetParameterValueAsString(*sfIter)));
      RbtParameterFileSourcePtr spSFSource(new RbtParameterFileSource(sfFile));
      //Create and add the subaggregate
      spSF->Add(spSFFactory->CreateAggFromFile(spSFSource,*sfIter));
    }
    
    //Add the RESTRAINT subaggregate scoring function from any SF definitions in the receptor prm file
    spSF->Add(spSFFactory->CreateAggFromFile(spRecepPrmSource,_RESTRAINT_SF));
    
    //Create the docking transform aggregate from the transform definitions in the docking prm file
    RbtTransformFactoryPtr spTransformFactory(new RbtTransformFactory());
    spParamSource->SetSection();
    RbtTransformAggPtr spTransform(spTransformFactory->CreateAggFromFile(spParamSource,_ROOT_TRANSFORM));
    
    //Override the TRACE levels for the scoring function and transform
    //Dump details to cout
    //Register the scoring function and the transform with the workspace
    if (bTrace) {
      RbtRequestPtr spTraceReq(new RbtSFSetParamRequest("TRACE",iTrace));
      spSF->HandleRequest(spTraceReq);
      spTransform->HandleRequest(spTraceReq);
    }
    if (iTrace > 0) {
      cout << endl << "SCORING FUNCTION DETAILS:" << endl << *spSF << endl;
      cout << endl << "SEARCH DETAILS:" << endl << *spTransform << endl;
    }
    spWS->SetSF(spSF);
    spWS->SetTransform(spTransform);
    
    //DM 18 May 1999
    //Variants describing the library version, exe version, parameter file, and current directory
    //Will be stored in the ligand SD files
    RbtVariant vLib(Rbt::GetProduct()+" ("+Rbt::GetVersion()+", Build"+Rbt::GetBuild()+")");
    RbtVariant vExe(strExeName+EXEVERSION);
    RbtVariant vRecep(spRecepPrmSource->GetFileName());
    RbtVariant vPrm(spParamSource->GetFileName());
    RbtVariant vDir(Rbt::GetCurrentDirectory());
    
    spRecepPrmSource->SetSection();
   //Read docking site from file and register with workspace
    RbtString strASFile = spWS->GetName()+".as";
    RbtString strInputFile = Rbt::GetRbtFileName("data/grids",strASFile);
    //DM 26 Sep 2000 - ios_base::binary is invalid with IRIX CC
#if defined(__sgi) && !defined(__GNUC__)
    ifstream istr(strInputFile.c_str(),ios_base::in);
#else
    ifstream istr(strInputFile.c_str(),ios_base::in|ios_base::binary);
#endif
    //DM 14 June 2006 - bug fix to one of the longest standing rDock issues
    //(the cryptic "Error reading from input stream" message, if cavity file was missing) 
    if (!istr) {
      RbtString message = "Cavity file (" + strASFile + ") not found in current directory or $RBT_HOME";
      message += " - run rbcavity first"; 
      throw RbtFileReadError(_WHERE_,message);
    } 
    RbtDockingSitePtr spDS(new RbtDockingSite(istr));
    istr.close();
    spWS->SetDockingSite(spDS);
    cout << endl << "DOCKING SITE" << endl << (*spDS) << endl;

    RbtPRMFactory prmFactory(spRecepPrmSource, spDS);
    prmFactory.SetTrace(iTrace);
    //Create the receptor model from the file names in the receptor parameter file
    RbtModelPtr spReceptor = prmFactory.CreateReceptor();
    spWS->SetReceptor(spReceptor);

    //Register any solvent
    RbtModelList solventList = prmFactory.CreateSolvent();
 	spWS->SetSolvent(solventList);
 	if (spWS->hasSolvent()) {
 		RbtInt nSolvent = spWS->GetSolvent().size();
 		cout << endl << nSolvent << " solvent molecules registered" << endl;
 	}
 	else {
 		cout << endl << "No solvent" << endl;
 	}
 	
    //Prepare the SD file sink for saving the docked conformations for each ligand
    //DM 3 Dec 1999 - replaced ostrstream with RbtString in determining SD file name
    if (bOutput) {
      RbtMolecularFileSinkPtr spMdlFileSink(new RbtMdlFileSink(strRunName+".sd",RbtModelPtr()));
      spWS->SetSink(spMdlFileSink);
    }
    
    //Seed the random number generator
    RbtRand& theRand = Rbt::GetRbtRand();//ref to random number generator
    if (bSeed) {
      theRand.Seed(nSeed);
    }
    
    //Create the filter object for controlling early termination of protocol
    RbtFilterPtr spfilter;
    if (bFilter) {
      spfilter = new RbtFilter(strFilterFile);
      if (bDockingRuns) {
	spfilter->SetMaxNRuns(nDockingRuns);
      }
    }
    else {
      spfilter = new RbtFilter(strFilter.str(), true);
    }
    if (bTrace) {
      RbtRequestPtr spTraceReq(new RbtSFSetParamRequest("TRACE",iTrace));
      spfilter->HandleRequest(spTraceReq);
    }

    //Register the Filter with the workspace
    spWS->SetFilter(spfilter);

    //MAIN LOOP OVER LIGAND RECORDS
    //DM 20 Apr 1999 - add explicit bPosIonise and bNegIonise flags to MdlFileSource constructor
    RbtMolecularFileSourcePtr spMdlFileSource(new RbtMdlFileSource(strLigandMdlFile,bPosIonise,bNegIonise,bImplH));
    for (RbtInt nRec=1; spMdlFileSource->FileStatusOK(); spMdlFileSource->NextRecord(), nRec++) {
      cout.setf(ios_base::left,ios_base::adjustfield);
      cout << endl
       << "**************************************************" << endl
       << "RECORD #" << nRec << endl;
      RbtError molStatus = spMdlFileSource->Status();
      if (!molStatus.isOK()) {
        cout << endl << molStatus << endl
             << "************************************************" << endl;
        continue;
      }
      
      //DM 26 Jul 1999 - only read the largest segment (guaranteed to be called H)
      //BGD 07 Oct 2002 - catching errors created by the ligands,
      //so rbdock continues with the next one, instead of
      //completely stopping
      try
      {
        spMdlFileSource->SetSegmentFilterMap
                                  (Rbt::ConvertStringToSegmentMap("H"));
      
        if (spMdlFileSource->isDataFieldPresent("Name"))
          cout << "NAME:   " << spMdlFileSource->GetDataValue("Name") << endl;
        if (spMdlFileSource->isDataFieldPresent("REG_Number"))
          cout << "REG_Num:" << spMdlFileSource->GetDataValue("REG_Number") 
               << endl;
        cout << setw(30) << "RANDOM_NUMBER_SEED:" << theRand.GetSeed() << endl;
      
        //Create and register the ligand model
        RbtModelPtr spLigand = prmFactory.CreateLigand(spMdlFileSource);
        RbtString strMolName = spLigand->GetName();
        spWS->SetLigand(spLigand);
        //Update any model coords from embedded chromosomes in the ligand file
        spWS->UpdateModelCoordsFromChromRecords(spMdlFileSource, iTrace);
      
        //DM 18 May 1999 - store run info in model data
        //Clear any previous Rbt.* data fields
        spLigand->ClearAllDataFields("Rbt.");
        spLigand->SetDataValue("Rbt.Library",vLib);
        spLigand->SetDataValue("Rbt.Executable",vExe);
        spLigand->SetDataValue("Rbt.Receptor",vRecep);
        spLigand->SetDataValue("Rbt.Parameter_File",vPrm);
        spLigand->SetDataValue("Rbt.Current_Directory",vDir);
      
        //DM 10 Dec 1999 - if in target mode, loop until target score is reached
        RbtBool bTargetMet = false;
        
        ////////////////////////////////////////////////////
        //MAIN LOOP OVER EACH SIMULATED ANNEALING RUN
        //Create a history file sink, just in case it's needed by any 
        //of the transforms
        RbtInt iRun = 1;
	// need to check this here. The termination 
	// filter is only run once at least
	// one docking run has been done.
        if (nDockingRuns < 1) 	
          bTargetMet = true;
        while (!bTargetMet) {
	  //Catching errors with this specific run
          try {
	    if (bOutput) {
              ostrstream histr;
              histr << strRunName << "_" << strMolName << nRec << "_his_" 
                    << iRun << ".sd" << ends;
              RbtMolecularFileSinkPtr spHistoryFileSink
		(new RbtMdlFileSink(histr.str(),spLigand));
              delete histr.str();
              spWS->SetHistorySink(spHistoryFileSink);
            }
            spWS->Run();//Dock!
	    RbtBool bterm = spfilter->Terminate();
	    RbtBool bwrite = spfilter->Write();
	    if (bterm)
	      bTargetMet = true;
	    if (bOutput && bwrite) {
	      spWS->Save();
	    }
	    iRun++;
          }
          catch (RbtDockingError& e) {
	    cout << e << endl;
	  }
        }
	//END OF MAIN LOOP OVER EACH SIMULATED ANNEALING RUN
	////////////////////////////////////////////////////
      } 
      //END OF TRY
      catch (RbtLigandError& e) {
	cout << e << endl;
      }
    }
    //END OF MAIN LOOP OVER LIGAND RECORDS
    ////////////////////////////////////////////////////
    cout << endl << "END OF RUN" << endl;
//    if (bOutput && flexRec) {
//      RbtMolecularFileSinkPtr spRecepSink(new RbtCrdFileSink(strRunName+".crd",spReceptor));
//      spRecepSink->Render();
//    }
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
