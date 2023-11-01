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

// Standalone executable for generating docking site .as files for rbdock

#include <string>
#include <iomanip>
using std::setw;

#include <popt.h>  // for popt command-line parsing

#include <algorithm>

#include "RbtBiMolWorkSpace.h"
#include "RbtCrdFileSink.h"
#include "RbtDockingSite.h"
#include "RbtPRMFactory.h"
#include "RbtParameterFileSource.h"
#include "RbtPsfFileSink.h"
#include "RbtSiteMapperFactory.h"

const RbtString EXEVERSION = " ($Id: //depot/dev/client3/rdock/2021.1/src/exe/rbcavity.cxx#3 $)";

void PrintUsage(void) {
    cout << "rbcavity - calculate docking cavities" << endl;
    cout << "Usage:\trbcavity -r<ReceptorPrmFile> [-was -ras -d -l<dist> -b<border>" << endl;
    cout << "Options:" << endl;
    cout << "\t\t-r<PrmFile> - receptor param file (contains active site params)" << endl;
    cout << "\t\t-was [-W]   - write docking cavities (plus distance grid) to .as file" << endl;
    cout << "\t\t-ras [-R]   - read docking cavities (plus distance grid) from .as file" << endl;
    cout << "\t\t-d          - dump InsightII grids for each cavity for visualisation" << endl;
    cout << "\t\t-v          - dump target PSF/CRD files for rDock Viewer" << endl;
    cout << "\t\t-l<dist>    - list receptor atoms with <dist> A of any cavity" << endl;
    cout << "\t\t-s          - print SITE descriptors (counts of exposed atoms)" << endl;
    cout << "\t\t-b<border>  - set the border around the cavities for the distance grid (default=8A)" << endl;
    cout << "\t\t-m          - write active site into a MOE grid" << endl;
}

struct RBCavityConfig {
    std::string prmFile;
    bool readAS;
    bool writeAS;
    bool dump;
    bool viewer;
    bool list;
    bool site;
    bool moeGrid;
    double border;
    double dist;

};

RBCavityConfig parseCommandLine(int argc, const char *argv[]) {
    RBCavityConfig config;
    config.readAS = false;
    config.writeAS = false;
    config.dump = false;
    config.viewer = false;
    config.list = false;
    config.site = false;
    config.moeGrid = false;
    config.border = 8.0;
    config.dist = 5.0;

    poptContext optCon;       // ditto
    char *prmFile = NULL;     // will be strReceptorPrmFile
    char *listDist = NULL;    // will be 'dist'
    char *borderDist = NULL;  // will be 'border'
    struct poptOption optionsTable[] = {
        // command line options
        {"receptor", 'r', POPT_ARG_STRING | POPT_ARGFLAG_ONEDASH, &prmFile, 0, "receptor file"},
        {"was", 'W', POPT_ARG_NONE | POPT_ARGFLAG_ONEDASH, 0, 'W', "write active site"},
        {"ras", 'R', POPT_ARG_NONE | POPT_ARGFLAG_ONEDASH, 0, 'R', "read active site"},
        {"dump-insight", 'd', POPT_ARG_NONE | POPT_ARGFLAG_ONEDASH, 0, 'd', "dump InsightII/PyMol grids"},
        //{"dump-moe",    'm',POPT_ARG_NONE  |POPT_ARGFLAG_ONEDASH,0,          'm',"dump MOE grids"}, //not working
        // right now so commenting it
        {"viewer", 'v', POPT_ARG_NONE | POPT_ARGFLAG_ONEDASH, 0, 'v', "dump Viewer PSF/CRD files"},
        {"list", 'l', POPT_ARG_STRING | POPT_ARGFLAG_ONEDASH, &listDist, 'l', "list receptor atoms within <dist>"},
        {"site", 's', POPT_ARG_NONE | POPT_ARGFLAG_ONEDASH, 0, 's', "print site descriptors"},
        {"border",
         'b',
         POPT_ARG_STRING | POPT_ARGFLAG_ONEDASH,
         &borderDist,
         'b',
         "set the border around the cavities"},
        POPT_AUTOHELP{NULL, 0, 0, NULL, 0}};
    optCon = poptGetContext(NULL, argc, argv, optionsTable, 0);
    poptSetOtherOptionHelp(optCon, "-r<receptor.prm> [options]");
    if (argc < 2) {
        poptPrintUsage(optCon, stderr, 0);
        exit(1);
    }
    char c;
    while ((c = poptGetNextOpt(optCon)) >= 0) {
        switch (c) {
            case 'b':
                config.border = atof(borderDist);
                break;
            case 'R':  // also for -ras
                config.readAS = true;
                break;
            case 'W':
                config.writeAS = true;
                break;
            case 'd':
                config.dump = true;
                break;
            case 'v':
                config.viewer = true;
                break;
            case 'l':
                config.list = true;
                config.dist = atof(listDist);
                break;
            case 'm':
                config.moeGrid = true;
                break;
            case 's':
                config.site = true;
                break;
            default:
                cout << "WARNING: unknown argument: " << c << endl;
                ;
                break;
        }
    }
    cout << endl;
    poptFreeContext(optCon);
    if (prmFile) {
        config.prmFile = prmFile;
    } else {
        cout << "Missing receptor parameter file name" << endl;
        PrintUsage();
        exit(1);
    }
    return config;
}

/////////////////////////////////////////////////////////////////////
// MAIN PROGRAM STARTS HERE
/////////////////////////////////////////////////////////////////////

void showConfig(const RBCavityConfig &config) {
    cout << "Command line arguments:" << endl;
    cout << "-r " << config.prmFile << endl;
    if (config.list) cout << "-l " << config.dist << endl;
    if (config.border) cout << "-b " << config.border << endl;
    if (config.writeAS) cout << "-was" << endl;
    if (config.readAS) cout << "-ras" << endl;
    if (config.moeGrid) cout << "-m" << endl;
    if (config.dump) cout << "-d" << endl;
    if (config.site) cout << "-s" << endl;
    if (config.viewer) cout << "-v" << endl;
}

ios_base::openmode addBinaryMode(ios_base::openmode mode) {
#if defined(__sgi) && !defined(__GNUC__)
    return mode;
#else
    return mode | ios_base::binary;
#endif
}

void rbcavity(const RBCavityConfig &config) {
    // Create a bimolecular workspace
    RbtBiMolWorkSpacePtr spWS(new RbtBiMolWorkSpace());
    // Set the workspace name to the root of the receptor .prm filename
    RbtStringList componentList = Rbt::ConvertDelimitedStringToList(config.prmFile, ".");
    RbtString wsName = componentList.front();
    spWS->SetName(wsName);

    // Read the protocol parameter file
    RbtParameterFileSourcePtr spRecepPrmSource(
        new RbtParameterFileSource(Rbt::GetRbtFileName("data/receptors", config.prmFile))
    );

    // Create the receptor model from the file names in the parameter file
    spRecepPrmSource->SetSection();
    RbtPRMFactory prmFactory(spRecepPrmSource);
    RbtModelPtr spReceptor = prmFactory.CreateReceptor();

    RbtDockingSitePtr spDockSite;
    RbtString strASFile = wsName + ".as";

    // Either read the docking site from the .as file
    if (config.readAS) {
        RbtString strInputFile = Rbt::GetRbtFileName("data/grids", strASFile);
        auto mode = addBinaryMode(ios_base::in);
        ifstream istr(strInputFile.c_str(), mode);
        spDockSite = RbtDockingSitePtr(new RbtDockingSite(istr));
        istr.close();
    }
    // Or map the site using the prescribed mapping algorithm
    else {
        RbtSiteMapperFactoryPtr spMapperFactory(new RbtSiteMapperFactory());
        RbtSiteMapperPtr spMapper = spMapperFactory->CreateFromFile(spRecepPrmSource, "MAPPER");
        spMapper->Register(spWS);
        spWS->SetReceptor(spReceptor);
        cout << *spMapper << endl;

        RbtInt nRI = spReceptor->GetNumSavedCoords() - 1;
        if (nRI == 0) {
            spDockSite = RbtDockingSitePtr(new RbtDockingSite((*spMapper)(), config.border));
        } else {
            RbtCavityList allCavs;
            for (RbtInt i = 1; i <= nRI; i++) {
                spReceptor->RevertCoords(i);
                RbtCavityList cavs((*spMapper)());
                std::copy(cavs.begin(), cavs.end(), std::back_inserter(allCavs));
            }
            spDockSite = RbtDockingSitePtr(new RbtDockingSite(allCavs, config.border));
        }
    }

    cout << endl << "DOCKING SITE" << endl << (*spDockSite) << endl;
    if (config.writeAS) {
        auto mode = addBinaryMode(ios_base::out | ios_base::trunc);
        ofstream ostr(strASFile.c_str(), mode);
        spDockSite->Write(ostr);
        ostr.close();
    }

    // Write PSF/CRD files to keep the rDock Viewer happy (it doesn't read MOL2 files yet)
    if (config.viewer) {
        RbtMolecularFileSinkPtr spRecepSink = new RbtPsfFileSink(wsName + "_for_viewer.psf", spReceptor);
        cout << "Writing PSF file: " << spRecepSink->GetFileName() << endl;
        spRecepSink->Render();
        spRecepSink = new RbtCrdFileSink(wsName + "_for_viewer.crd", spReceptor);
        cout << "Writing CRD file: " << spRecepSink->GetFileName() << endl;
        spRecepSink->Render();
    }

    // Write an ASCII InsightII grid file for each defined cavity
    if (config.dump) {
        RbtCavityList cavList = spDockSite->GetCavityList();
        for (RbtInt i = 0; i < cavList.size(); i++) {
            ostringstream filename;
            filename << wsName << "_cav" << i + 1 << ".grd" << ends;
            ofstream dumpFile(filename.str());
            if (dumpFile) {
                cavList[i]->GetGrid()->PrintInsightGrid(dumpFile);
                dumpFile.close();
            }
        }
    }
    // writing active site into MOE grid
    if (config.moeGrid) {
        cout << "MOE grid feature not yet implemented, sorry." << endl;
    }
    // List all receptor atoms within given distance of any cavity
    if (config.list) {
        RbtRealGridPtr spGrid = spDockSite->GetGrid();
        RbtAtomList atomList = spDockSite->GetAtomList(spReceptor->GetAtomList(), 0.0, config.dist);
        cout << atomList.size() << " receptor atoms within " << config.dist << " A of any cavity" << endl;
        cout << endl << "DISTANCE,ATOM" << endl;
        for (RbtAtomListConstIter iter = atomList.begin(); iter != atomList.end(); iter++) {
            cout << spGrid->GetSmoothedValue((*iter)->GetCoords()) << "\t" << **iter << endl;
        }
        cout << endl;
    }

    // DM 15 Jul 2002 - print out SITE descriptors
    // Use a crude measure of solvent accessibility - count #atoms within 4A of each atom
    // Use an empirical threshold to determine if atom is exposed or not
    if (config.site) {
        RbtDouble cavDist = 4.0;  // Use a fixed definition of cavity atoms - all those within 4A of docking volume
        RbtDouble neighbR = 4.0;  // Sphere radius for counting nearest neighbours
        RbtDouble threshold = 15;  // Definition of solvent exposed: neighbours < threshold
        // RbtRealGridPtr spGrid = spDockSite->GetGrid();
        RbtAtomList recepAtomList = spReceptor->GetAtomList();
        RbtAtomList cavAtomList = spDockSite->GetAtomList(recepAtomList, cavDist);
        RbtAtomList exposedAtomList;  // The list of exposed cavity atoms
        cout << endl << "SOLVENT EXPOSED CAVITY ATOMS" << endl;
        cout << "1) Consider atoms within " << cavDist << "A of docking site" << endl;
        cout << "2) Determine #neighbours within " << neighbR << "A of each atom" << endl;
        cout << "3) If #neighbours < " << threshold << " => exposed" << endl;
        cout << "4) Calculate SITE* descriptors over exposed cavity atoms only" << endl;
        cout << endl << "ATOM NAME\t#NEIGHBOURS" << endl;

        // Get the list of solvent exposed cavity atoms
        for (RbtAtomListConstIter iter = cavAtomList.begin(); iter != cavAtomList.end(); iter++) {
            RbtInt nNeighb =
                Rbt::GetNumAtoms(recepAtomList, Rbt::isAtomInsideSphere((*iter)->GetCoords(), neighbR));
            nNeighb--;
            if (nNeighb < threshold) {
                cout << (*iter)->GetFullAtomName() << "\t" << nNeighb << endl;
                exposedAtomList.push_back(*iter);
            }
        }

        // Total +ve and -ve charges
        RbtDouble posChg(0.0);
        RbtDouble negChg(0.0);
        for (RbtAtomListConstIter iter = exposedAtomList.begin(); iter != exposedAtomList.end(); iter++) {
            RbtDouble chg = (*iter)->GetGroupCharge();
            if (chg > 0.0) {
                posChg += chg;
            } else if (chg < 0.0) {
                negChg += chg;
            }
        }

        // Atom type counts
        Rbt::isHybridState_eq bIsArom(RbtAtom::AROM);
        RbtInt nAtoms = exposedAtomList.size();
        RbtInt nLipoC = Rbt::GetNumAtoms(exposedAtomList, Rbt::isAtomLipophilic());
        RbtInt nArom = Rbt::GetNumAtoms(exposedAtomList, bIsArom);
        RbtInt nNHBD = Rbt::GetNumAtoms(exposedAtomList, Rbt::isAtomHBondDonor());
        RbtInt nMetal = Rbt::GetNumAtoms(exposedAtomList, Rbt::isAtomMetal());
        RbtInt nGuan = Rbt::GetNumAtoms(exposedAtomList, Rbt::isAtomGuanidiniumCarbon());
        RbtInt nNHBA = Rbt::GetNumAtoms(exposedAtomList, Rbt::isAtomHBondAcceptor());

        // Cavity volume
        cout << endl << wsName << ",SITE_VOL," << spDockSite->GetVolume() << endl;
        // Atom type counts
        cout << wsName << ",SITE_NATOMS," << nAtoms << endl;
        cout << wsName << ",SITE_NLIPOC," << nLipoC << endl;
        cout << wsName << ",SITE_NAROMATOMS," << nArom << endl;
        cout << wsName << ",SITE_NHBD," << nNHBD << endl;
        cout << wsName << ",SITE_NMETAL," << nMetal << endl;
        cout << wsName << ",SITE_NGUAN," << nGuan << endl;
        cout << wsName << ",SITE_NHBA," << nNHBA << endl;
        // Atom type percentages
        cout << wsName << ",SITE_PERC_LIPOC," << 100.0 * nLipoC / nAtoms << endl;
        cout << wsName << ",SITE_PERC_AROMATOMS," << 100.0 * nArom / nAtoms << endl;
        cout << wsName << ",SITE_PERC_HBD," << 100.0 * nNHBD / nAtoms << endl;
        cout << wsName << ",SITE_PERC_METAL," << 100.0 * nMetal / nAtoms << endl;
        cout << wsName << ",SITE_PERC_GUAN," << 100.0 * nGuan / nAtoms << endl;
        cout << wsName << ",SITE_PERC_HBA," << 100.0 * nNHBA / nAtoms << endl;
        // Charges
        cout << wsName << ",SITE_POS_CHG," << posChg << endl;
        cout << wsName << ",SITE_NEG_CHG," << negChg << endl;
        cout << wsName << ",SITE_TOT_CHG," << posChg + negChg << endl;
    }
}

void printHeader(const char *argv[]) {
    // Strip off the path to the executable, leaving just the file name
    RbtString strExeName(argv[0]);
    RbtString::size_type i = strExeName.rfind("/");
    if (i != RbtString::npos) strExeName.erase(0, i + 1);

    // Print a standard header
    Rbt::PrintStdHeader(cout, strExeName + EXEVERSION);
}

int main(int argc, const char *argv[]) {
    auto config = parseCommandLine(argc, argv);

    printHeader(argv);

    showConfig(config);

    cout.setf(ios_base::left, ios_base::adjustfield);

    try {
        rbcavity(config);
    } catch (RbtError &e) {
        cout << e << endl;
    } catch (...) {
        cout << "Unknown exception" << endl;
    }

    _RBTOBJECTCOUNTER_DUMP_(cout)

    return 0;
}
