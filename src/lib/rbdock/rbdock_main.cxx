#include "rbdock/rbdock_main.h"

#include <algorithm>
#include <iomanip>

#include "rbdock/rbdock_argparser.h"
#include "rbdock/rbdock_config.h"

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

// Section name in docking prm file containing scoring function definition
const RbtString _ROOT_SF = "SCORE";
const RbtString _RESTRAINT_SF = "RESTR";
const RbtString _ROOT_TRANSFORM = "DOCK";

std::string RBDock::get_filter_string(const RBDock::RBDockConfig &config) {
    std::ostringstream strFilter;
    if (!config.bFilter) {
        if (config.bTarget) {            // -t<TS>
            if (!config.bDockingRuns) {  // -t<TS> only
                strFilter << "0 1 - SCORE.INTER " << config.dTargetScore << std::endl;
            } else  // -t<TS> -n<N> need to check if -cont present
                    // for all other cases it doesn't matter
                if (config.bContinue) {  // -t<TS> -n<N> -cont
                    strFilter << "1 if - SCORE.NRUNS " << (config.nDockingRuns - 1) << " 0.0 -1.0,\n1 - SCORE.INTER "
                              << config.dTargetScore << std::endl;
                } else {  // -t<TS> -n<N>
                    strFilter << "1 if - " << config.dTargetScore << " SCORE.INTER 0.0 "
                              << "if - SCORE.NRUNS " << (config.nDockingRuns - 1) << " 0.0 -1.0,\n1 - SCORE.INTER "
                              << config.dTargetScore << std::endl;
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

void RBDock::RBDock(const RBDock::RBDockConfig &config, const RbtString &strExeName) {
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
    RbtVariant vExe(strExeName + " - " + RBT_VERSION);
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