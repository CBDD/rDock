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

// Main docking application
#include <errno.h>

//#include <cxxopts.hpp>
#include <iomanip>

#include "RbtArgParser.h"
#include "RbtBiMolWorkSpace.h"
#include "RbtCrdFileSink.h"
#include "RbtDockingError.h"
#include "RbtFileError.h"
#include "RbtFilter.h"
#include "RbtLigandError.h"
#include "RbtMdlFileSink.h"
#include "RbtMdlFileSource.h"
#include "RbtModelError.h"
#include "RbtPRMFactory.h"
#include "RbtParameterFileSource.h"
#include "RbtPlatformCompatibility.h"
#include "RbtRand.h"
#include "RbtSFFactory.h"
#include "RbtSFRequest.h"
#include "RbtTransformFactory.h"
#include "RbtVersion.h"

const RbtString EXEVERSION = RBT_VERSION;
// Section name in docking prm file containing scoring function definition
const RbtString _ROOT_SF = "SCORE";
const RbtString _RESTRAINT_SF = "RESTR";
const RbtString _ROOT_TRANSFORM = "DOCK";

RbtArgParser::RbtArgParser get_arguments_parser() {
    using std::string;

    RbtArgParser::RbtArgParser parser("rbdock", "Docking application");
    parser.add<string>("i,input", "input ligand SD file (mandatory)");
    parser.add<string>("o,output-root", "root name for output file(s)", "");
    parser.add<string>("r,receptor", "receptor parameter file (mandatory)");
    parser.add<string>("p,protocol", "docking protocol parameter file", "dock.prm");
    parser.add<int>("n,runs", "number of runs/ligand", "1");
    parser.add_flag("P,ap", "protonate all neutral amines, guanidines, imidazoles");
    parser.add_flag("D,an", "deprotonate all carboxylic, sulphur and phosphorous acid groups");
    parser.add_flag("H,allH", "read all hydrogens present. If disabled, read polar hydrogens only");
    parser.add<string>("t,target", "score threshold OR filter file name");
    parser.add_flag("C,cont", "if enabled, continue if score threshold is met (use with -t <targetScore>)");
    parser.add<int>("T,trace", "controls output level for debugging (0 = minimal, >0 = more verbose)", "0");
    parser.add<int>("s,seed", "random number seed (default=from sys clock)");

    return parser;
}

struct RBDockConfig {
    RbtString strLigandMdlFile;
    RbtBool bOutput = false;
    RbtString strRunName;
    RbtString strReceptorPrmFile;
    RbtString strParamFile;
    RbtString strFilterFile;
    RbtInt nDockingRuns = 1;
    RbtBool bTarget = false;
    RbtBool bContinue = false;
    RbtBool bDockingRuns = false;
    RbtDouble dTargetScore;
    RbtBool bFilter = false;
    RbtBool bPosIonise = false;
    RbtBool bNegIonise = false;
    RbtBool bAllH = false;
    RbtBool bSeed = false;
    RbtInt nSeed = 0;
    RbtBool bTrace = false;
    RbtInt iTrace = 0;

    friend std::ostream &operator<<(std::ostream &os, const RBDockConfig &config);

    void validate() {
        using RbtArgParser::ValidationError;
        if (strLigandMdlFile.empty()) throw ValidationError("input ligand SD file is mandatory");
        if (strReceptorPrmFile.empty()) throw ValidationError("receptor parameter file is mandatory");
        if (nDockingRuns < 1) throw ValidationError("number of runs must be greater than 0");
    }
};

std::ostream &operator<<(std::ostream &os, const RBDockConfig &config) {
    os << "Command line args:" << std::endl;
    os << " -i " << config.strLigandMdlFile << std::endl;
    os << " -r " << config.strReceptorPrmFile << std::endl;
    os << " -p " << config.strParamFile << std::endl;
    if (config.bOutput) {
        os << " -o " << config.strRunName << std::endl;
    } else {
        os << "WARNING: output file name is missing." << std::endl;
    }
    auto def_aux = config.bDockingRuns ? " (default) " : "";
    os << " -n " << config.nDockingRuns << def_aux << std::endl;
    if (config.bSeed) os << " -s " << config.nSeed << std::endl;
    if (config.bTrace) os << " -T " << config.iTrace << std::endl;
    if (config.bPosIonise) os << " -ap " << std::endl;
    if (config.bNegIonise) os << " -an " << std::endl;
    if (!config.bAllH) os << " -allH " << std::endl;
    if (config.bTarget) os << " -t " << config.dTargetScore << std::endl;
    if (config.bContinue) os << " -cont " << std::endl;

    return os;
}

RBDockConfig parse_args(int argc, const char *argv[]) {
    auto parser = get_arguments_parser();
    try {
        auto parser_result = RbtArgParser::RbtParseResult(parser.parse(argc, argv));
        RBDockConfig config;
        parser_result["input"] >> config.strLigandMdlFile;
        parser_result["receptor"] >> config.strReceptorPrmFile;
        parser_result["protocol"] >> config.strParamFile;
        config.bOutput = parser_result["output-root"].is_present();
        parser_result["output-root"] >> config.strRunName;
        parser_result["runs"] >> config.nDockingRuns;
        parser_result["ap"] >> config.bPosIonise;
        parser_result["an"] >> config.bNegIonise;
        parser_result["allH"] >> config.bAllH;
        config.bSeed = parser_result["seed"].is_present();
        parser_result["seed"] >> config.nSeed;
        config.bTrace = parser_result["trace"].is_present();
        parser_result["trace"] >> config.iTrace;
        if (parser_result["target"].is_present()) {
            std::string target;
            parser_result["target"] >> target;
            if (target.find(".prm") != string::npos) {
                config.bFilter = true;
                config.strFilterFile = target;
            } else {
                config.bTarget = true;
                config.dTargetScore = std::stod(target);
            }
        }
        parser_result["cont"] >> config.bContinue;

        config.validate();
        return config;
    } catch (RbtArgParser::ParsingError &e) {
        std::cerr << "Error parsing options: " << e.what() << std::endl;
    } catch (RbtArgParser::ValidationError &e) {
        std::cerr << "Invalid configuration: " << e.what() << std::endl;
    }
    std::cerr << parser.help() << std::endl;
    exit(1);
}

/////////////////////////////////////////////////////////////////////
// MAIN PROGRAM STARTS HERE
/////////////////////////////////////////////////////////////////////

std::string get_filter_string(const RBDockConfig &config) {
    ostringstream strFilter;
    if (!config.bFilter) {
        if (config.bTarget) {            // -t<TS>
            if (!config.bDockingRuns) {  // -t<TS> only
                strFilter << "0 1 - SCORE.INTER " << config.dTargetScore << endl;
            } else  // -t<TS> -n<N> need to check if -cont present
                    // for all other cases it doesn't matter
                if (config.bContinue) {  // -t<TS> -n<N> -cont
                    strFilter << "1 if - SCORE.NRUNS " << (config.nDockingRuns - 1) << " 0.0 -1.0,\n1 - SCORE.INTER "
                              << config.dTargetScore << endl;
                } else {  // -t<TS> -n<N>
                    strFilter << "1 if - " << config.dTargetScore << " SCORE.INTER 0.0 "
                              << "if - SCORE.NRUNS " << (config.nDockingRuns - 1) << " 0.0 -1.0,\n1 - SCORE.INTER "
                              << config.dTargetScore << endl;
                }
        }                                // no target score, no filter
        else if (config.bDockingRuns) {  // -n<N>
            strFilter << "1 if - SCORE.NRUNS " << (config.nDockingRuns - 1) << " 0.0 -1.0,\n0";
        } else {  // no -t no -n
            strFilter << "0 0\n";
        }
    }
    return strFilter.str();
}

void RBDock(const RBDockConfig &config, const RbtString &strExeName) {
    // Create a bimolecular workspace
    RbtBiMolWorkSpacePtr spWS(new RbtBiMolWorkSpace());
    // Set the workspace name to the root of the receptor .prm filename
    RbtStringList componentList = Rbt::ConvertDelimitedStringToList(config.strReceptorPrmFile, ".");
    RbtString wsName = componentList.front();
    spWS->SetName(wsName);

    // Read the docking protocol parameter file
    RbtParameterFileSourcePtr spParamSource(
        new RbtParameterFileSource(Rbt::GetRbtFileName("data/scripts", config.strParamFile))
    );
    // Read the receptor parameter file
    RbtParameterFileSourcePtr spRecepPrmSource(
        new RbtParameterFileSource(Rbt::GetRbtFileName("data/receptors", config.strReceptorPrmFile))
    );
    cout << endl
         << "DOCKING PROTOCOL:" << endl
         << spParamSource->GetFileName() << endl
         << spParamSource->GetTitle() << endl;
    cout << endl
         << "RECEPTOR:" << endl
         << spRecepPrmSource->GetFileName() << endl
         << spRecepPrmSource->GetTitle() << endl;

    // Create the scoring function from the SCORE section of the docking protocol prm file
    // Format is:
    // SECTION SCORE
    //     INTER    RbtInterSF.prm
    //     INTRA RbtIntraSF.prm
    // END_SECTION
    //
    // Notes:
    // Section name must be SCORE. This is also the name of the root SF aggregate
    // An aggregate is created for each parameter in the section.
    // Parameter name becomes the name of the subaggregate (e.g. SCORE.INTER)
    // Parameter value is the file name for the subaggregate definition
    // Default directory is $RBT_ROOT/data/sf
    RbtSFFactoryPtr spSFFactory(new RbtSFFactory());  // Factory class for scoring functions
    RbtSFAggPtr spSF(new RbtSFAgg(_ROOT_SF));         // Root SF aggregate
    spParamSource->SetSection(_ROOT_SF);
    RbtStringList sfList(spParamSource->GetParameterList());
    // Loop over all parameters in the SCORE section
    for (RbtStringListConstIter sfIter = sfList.begin(); sfIter != sfList.end(); sfIter++) {
        // sfFile = file name for scoring function subaggregate
        RbtString sfFile(Rbt::GetRbtFileName("data/sf", spParamSource->GetParameterValueAsString(*sfIter)));
        RbtParameterFileSourcePtr spSFSource(new RbtParameterFileSource(sfFile));
        // Create and add the subaggregate
        spSF->Add(spSFFactory->CreateAggFromFile(spSFSource, *sfIter));
    }

    // Add the RESTRAINT subaggregate scoring function from any SF definitions in the receptor prm file
    spSF->Add(spSFFactory->CreateAggFromFile(spRecepPrmSource, _RESTRAINT_SF));

    // Create the docking transform aggregate from the transform definitions in the docking prm file
    RbtTransformFactoryPtr spTransformFactory(new RbtTransformFactory());
    spParamSource->SetSection();
    RbtTransformAggPtr spTransform(spTransformFactory->CreateAggFromFile(spParamSource, _ROOT_TRANSFORM));

    // Override the TRACE levels for the scoring function and transform
    // Dump details to cout
    // Register the scoring function and the transform with the workspace
    if (config.bTrace) {
        RbtRequestPtr spTraceReq(new RbtSFSetParamRequest("TRACE", config.iTrace));
        spSF->HandleRequest(spTraceReq);
        spTransform->HandleRequest(spTraceReq);
    }
    if (config.iTrace > 0) {
        cout << endl << "SCORING FUNCTION DETAILS:" << endl << *spSF << endl;
        cout << endl << "SEARCH DETAILS:" << endl << *spTransform << endl;
    }
    spWS->SetSF(spSF);
    spWS->SetTransform(spTransform);

    // DM 18 May 1999
    // Variants describing the library version, exe version, parameter file, and current directory
    // Will be stored in the ligand SD files
    RbtVariant vLib(Rbt::GetProduct() + " (" + Rbt::GetVersion() + ", Build" + Rbt::GetBuild() + ")");
    RbtVariant vExe(strExeName + " - " + EXEVERSION);
    RbtVariant vRecep(spRecepPrmSource->GetFileName());
    RbtVariant vPrm(spParamSource->GetFileName());
    RbtVariant vDir(Rbt::GetCurrentDirectory());

    spRecepPrmSource->SetSection();
    // Read docking site from file and register with workspace
    RbtString strASFile = spWS->GetName() + ".as";
    RbtString strInputFile = Rbt::GetRbtFileName("data/grids", strASFile);
    // DM 26 Sep 2000 - ios_base::binary is invalid with IRIX CC
    ifstream istr(strInputFile.c_str(), Rbt::inputMode);
    // DM 14 June 2006 - bug fix to one of the longest standing rDock issues
    //(the cryptic "Error reading from input stream" message, if cavity file was missing)
    if (!istr) {
        RbtString message = "Cavity file (" + strASFile + ") not found in current directory or $RBT_HOME";
        message += " - run rbcavity first";
        throw RbtFileReadError(_WHERE_, message);
    }
    RbtDockingSitePtr spDS(new RbtDockingSite(istr));
    istr.close();
    spWS->SetDockingSite(spDS);
    cout << endl << "DOCKING SITE" << endl << (*spDS) << endl;

    // Prepare the SD file sink for saving the docked conformations for each ligand
    // DM 3 Dec 1999 - replaced ostringstream with RbtString in determining SD file name
    // SRC 2014 moved here this block to allow WRITE_ERROR TRUE
    if (config.bOutput) {
        RbtMolecularFileSinkPtr spMdlFileSink(new RbtMdlFileSink(config.strRunName + ".sd", RbtModelPtr()));
        spWS->SetSink(spMdlFileSink);
    }

    RbtPRMFactory prmFactory(spRecepPrmSource, spDS);
    prmFactory.SetTrace(config.iTrace);
    // Create the receptor model from the file names in the receptor parameter file
    RbtModelPtr spReceptor = prmFactory.CreateReceptor();
    spWS->SetReceptor(spReceptor);

    // Register any solvent
    RbtModelList solventList = prmFactory.CreateSolvent();
    spWS->SetSolvent(solventList);
    if (spWS->hasSolvent()) {
        RbtInt nSolvent = spWS->GetSolvent().size();
        cout << endl << nSolvent << " solvent molecules registered" << endl;
    } else {
        cout << endl << "No solvent" << endl;
    }

    // SRC 2014 removed sector bOutput from here to some blocks above, for WRITEERRORS TRUE

    // Seed the random number generator
    RbtRand &theRand = Rbt::GetRbtRand();  // ref to random number generator
    if (config.bSeed) {
        theRand.Seed(config.nSeed);
    }

    // Create the filter object for controlling early termination of protocol
    RbtFilterPtr spfilter;
    if (config.bFilter) {
        spfilter = new RbtFilter(config.strFilterFile);
        if (config.bDockingRuns) {
            spfilter->SetMaxNRuns(config.nDockingRuns);
        }
    } else {
        spfilter = new RbtFilter(get_filter_string(config), true);
    }
    if (config.bTrace) {
        RbtRequestPtr spTraceReq(new RbtSFSetParamRequest("TRACE", config.iTrace));
        spfilter->HandleRequest(spTraceReq);
    }

    // Register the Filter with the workspace
    spWS->SetFilter(spfilter);

    // MAIN LOOP OVER LIGAND RECORDS
    // DM 20 Apr 1999 - add explicit bPosIonise and bNegIonise flags to MdlFileSource constructor
    RbtMolecularFileSourcePtr spMdlFileSource(
        new RbtMdlFileSource(config.strLigandMdlFile, config.bPosIonise, config.bNegIonise, !config.bAllH)
    );
    for (RbtInt nRec = 1; spMdlFileSource->FileStatusOK(); spMdlFileSource->NextRecord(), nRec++) {
        cout.setf(ios_base::left, ios_base::adjustfield);
        cout << endl << "**************************************************" << endl << "RECORD #" << nRec << endl;
        RbtError molStatus = spMdlFileSource->Status();
        if (!molStatus.isOK()) {
            cout << endl << molStatus << endl << "************************************************" << endl;
            continue;
        }

        // DM 26 Jul 1999 - only read the largest segment (guaranteed to be called H)
        // BGD 07 Oct 2002 - catching errors created by the ligands,
        // so rbdock continues with the next one, instead of
        // completely stopping
        try {
            spMdlFileSource->SetSegmentFilterMap(Rbt::ConvertStringToSegmentMap("H"));

            if (spMdlFileSource->isDataFieldPresent("Name"))
                cout << "NAME:   " << spMdlFileSource->GetDataValue("Name") << endl;
            if (spMdlFileSource->isDataFieldPresent("REG_Number"))
                cout << "REG_Num:" << spMdlFileSource->GetDataValue("REG_Number") << endl;
            cout << std::setw(30) << "RANDOM_NUMBER_SEED:" << theRand.GetSeed() << endl;

            // Create and register the ligand model
            RbtModelPtr spLigand = prmFactory.CreateLigand(spMdlFileSource);
            RbtString strMolName = spLigand->GetName();
            spWS->SetLigand(spLigand);
            // Update any model coords from embedded chromosomes in the ligand file
            spWS->UpdateModelCoordsFromChromRecords(spMdlFileSource, config.iTrace);

            // DM 18 May 1999 - store run info in model data
            // Clear any previous Rbt.* data fields
            spLigand->ClearAllDataFields("Rbt.");
            spLigand->SetDataValue("Rbt.Library", vLib);
            spLigand->SetDataValue("Rbt.Executable", vExe);
            spLigand->SetDataValue("Rbt.Receptor", vRecep);
            spLigand->SetDataValue("Rbt.Parameter_File", vPrm);
            spLigand->SetDataValue("Rbt.Current_Directory", vDir);

            // DM 10 Dec 1999 - if in target mode, loop until target score is reached
            RbtBool bTargetMet = false;

            ////////////////////////////////////////////////////
            // MAIN LOOP OVER EACH SIMULATED ANNEALING RUN
            // Create a history file sink, just in case it's needed by any
            // of the transforms
            RbtInt iRun = 1;
            // need to check this here. The termination
            // filter is only run once at least
            // one docking run has been done.
            if (config.nDockingRuns < 1) bTargetMet = true;
            while (!bTargetMet) {
                // Catching errors with this specific run
                try {
                    if (config.bOutput) {
                        ostringstream histr;
                        histr << config.strRunName << "_" << strMolName << nRec << "_his_" << iRun << ".sd";
                        RbtMolecularFileSinkPtr spHistoryFileSink(new RbtMdlFileSink(histr.str(), spLigand));
                        spWS->SetHistorySink(spHistoryFileSink);
                    }
                    spWS->Run();  // Dock!
                    RbtBool bterm = spfilter->Terminate();
                    RbtBool bwrite = spfilter->Write();
                    if (bterm) bTargetMet = true;
                    if (config.bOutput && bwrite) {
                        spWS->Save();
                    }
                    iRun++;
                } catch (RbtDockingError &e) {
                    cout << e << endl;
                }
            }
            // END OF MAIN LOOP OVER EACH SIMULATED ANNEALING RUN
            ////////////////////////////////////////////////////
        }
        // END OF TRY
        catch (RbtLigandError &e) {
            cout << e << endl;
        }
    }
    // END OF MAIN LOOP OVER LIGAND RECORDS
    ////////////////////////////////////////////////////
    cout << endl << "END OF RUN" << endl;
    //    if (bOutput && flexRec) {
    //      RbtMolecularFileSinkPtr spRecepSink(new RbtCrdFileSink(strRunName+".crd",spReceptor));
    //      spRecepSink->Render();
    //    }
}

int main(int argc, const char *argv[]) {
    cout.setf(ios_base::left, ios_base::adjustfield);

    // Strip off the path to the executable, leaving just the file name
    RbtString strExeName(argv[0]);
    RbtString::size_type i = strExeName.rfind("/");
    if (i != RbtString::npos) strExeName.erase(0, i + 1);

    // Print a standard header
    Rbt::PrintStdHeader(cout, strExeName + " - " + EXEVERSION);

    // Parse command line arguments
    RBDockConfig config = parse_args(argc, argv);
    std::cout << "Command line args:\n";
    std::cout << config << std::endl;

    // BGD 26 Feb 2003 - Create filters to simulate old rbdock
    // behaviour

    // DM 20 Apr 1999 - set the auto-ionise flags
    if (config.bPosIonise)
        cout << "Automatically protonating positive ionisable groups (amines, imidazoles, guanidines)" << endl;
    if (config.bNegIonise)
        cout << "Automatically deprotonating negative ionisable groups (carboxylic acids, phosphates, sulphates, "
                "sulphonates)"
             << endl;
    if (!config.bAllH)
        cout << "Reading polar hydrogens only from ligand SD file" << endl;
    else
        cout << "Reading all hydrogens from ligand SD file" << endl;

    if (config.bTarget) {
        cout << endl << "Lower target intermolecular score = " << config.dTargetScore << endl;
    }

    try {
        RBDock(config, strExeName);
    } catch (RbtError &e) {
        cout << e << endl;
        return 1;
    } catch (...) {
        cout << "Unknown exception" << endl;
        return 1;
    }

    _RBTOBJECTCOUNTER_DUMP_(cout)

    return 0;
}
