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
#include <set>
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

// this should go to a common header file, as it is used in several places across the code
ios_base::openmode addBinaryMode(ios_base::openmode mode) {
#if defined(__sgi) && !defined(__GNUC__)
    return mode;
#else
    return mode | ios_base::binary;
#endif
}

// until we can use c++20, we need to implement this ourselves
bool contains(std::set<char> & charSet, char c) { return charSet.find(c) != charSet.end(); }


struct RBCavityConfig {
    std::string prmFile;
    bool readAS {false};
    bool writeAS {false};
    bool dump {false};
    bool viewer {false};
    bool list {false};
    bool site {false};
    bool moeGrid {false};
    double border {8.0};
    double dist {5.0};
};

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

class RBCavityCommandLineParser {
 protected:
    char *prmFile;
    char *listDist;
    char *borderDist;
    std::vector<poptOption> optionsTable;

    // this is quite ugly, but sets in c++ didn't have a contains method until c++20
    // this function will be removed once we can use c++20...
 public:

    RBCavityCommandLineParser();
    RBCavityConfig parse(int argc, const char * argv[]);

 protected:
    std::vector<poptOption> createOptionsTable();
    std::set<char> getKnownFlags();
    std::set<char> getFlags(poptContext & optCon);
    void warnUnknownFlags(std::set<char> & flags);

};

/////////////////////////////////////////////////////////////////////
// MAIN PROGRAM STARTS HERE
/////////////////////////////////////////////////////////////////////

class SITEDescriptors {
 protected:
    RbtDouble cavDist = 4.0;  // Use a fixed definition of cavity atoms - all those within 4A of docking volume
    RbtDouble neighbR = 4.0;  // Sphere radius for counting nearest neighbours
    RbtDouble threshold = 15;  // Definition of solvent exposed: neighbours < threshold
    RbtString wsName;
    RbtDockingSitePtr spDockSite;
    RbtModelPtr spReceptor;

    // properties
    RbtDouble posChg;
    RbtDouble negChg;
    RbtInt nAtoms;
    RbtInt nLipoC;
    RbtInt nArom;
    RbtInt nNHBD;
    RbtInt nMetal;
    RbtInt nGuan;
    RbtInt nNHBA;

 public:

    SITEDescriptors(RbtString wsName, RbtDockingSitePtr spDockSite, RbtModelPtr spReceptor):
        cavDist(4.0),
        neighbR(4.0),
        threshold(15),
        wsName(wsName),
        spDockSite(spDockSite),
        spReceptor(spReceptor),
        posChg(0.0),
        negChg(0.0),
        nAtoms(0),
        nLipoC(0),
        nArom(0),
        nNHBD(0),
        nMetal(0),
        nGuan(0),
        nNHBA(0)
    {}

    void show(){
        RbtAtomList exposedAtomList = getExposedAtomList();

        calculateCharges(exposedAtomList);
        calculateAtomCounts(exposedAtomList);

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

 protected:
    RbtAtomList getExposedAtomList() {
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
        return exposedAtomList;
    }

    void calculateCharges(const RbtAtomList & exposedAtomList) {
        // Total +ve and -ve charges
        posChg = 0.0;
        negChg = 0.0;
        for (RbtAtomListConstIter iter = exposedAtomList.begin(); iter != exposedAtomList.end(); iter++) {
            RbtDouble chg = (*iter)->GetGroupCharge();
            if (chg > 0.0) {
                posChg += chg;
            } else if (chg < 0.0) {
                negChg += chg;
            }
        }
    }

    void calculateAtomCounts(const RbtAtomList & exposedAtomList) {
        // Atom type counts
        Rbt::isHybridState_eq bIsArom(RbtAtom::AROM);

        nAtoms = exposedAtomList.size();
        nLipoC = Rbt::GetNumAtoms(exposedAtomList, Rbt::isAtomLipophilic());
        nArom = Rbt::GetNumAtoms(exposedAtomList, bIsArom);
        nNHBD = Rbt::GetNumAtoms(exposedAtomList, Rbt::isAtomHBondDonor());
        nMetal = Rbt::GetNumAtoms(exposedAtomList, Rbt::isAtomMetal());
        nGuan = Rbt::GetNumAtoms(exposedAtomList, Rbt::isAtomGuanidiniumCarbon());
        nNHBA = Rbt::GetNumAtoms(exposedAtomList, Rbt::isAtomHBondAcceptor());
    }
};

class RBCavity {
    RBCavityConfig config;
    RbtDockingSitePtr spDockSite;
    RbtBiMolWorkSpacePtr spWS;
    RbtString wsName;
    RbtString strASFile;
    RbtParameterFileSourcePtr spRecepPrmSource;
    RbtModelPtr spReceptor;

 public:
    RBCavity(RBCavityConfig config):
        config(config),
        spWS(new RbtBiMolWorkSpace()),
        wsName(getWsName()),
        strASFile(wsName + ".as"),
        spRecepPrmSource(new RbtParameterFileSource(Rbt::GetRbtFileName("data/receptors", config.prmFile)))
    {
        spWS->SetName(getWsName());

    }

    RbtString getWsName() {
        RbtStringList componentList = Rbt::ConvertDelimitedStringToList(config.prmFile, ".");
        return componentList.front();
    }

    void loadReceptor() {
        spRecepPrmSource->SetSection();
        RbtPRMFactory prmFactory(spRecepPrmSource);
        spReceptor = prmFactory.CreateReceptor();
    }

    void loadDockingSite() {
        if(config.readAS) {
            loadDockingSiteFromAS();
        } else {
            loadDockingSiteFromMapper();
        }
    }

    void loadDockingSiteFromAS() {
        RbtString strInputFile = Rbt::GetRbtFileName("data/grids", strASFile);
        auto mode = addBinaryMode(ios_base::in);
        ifstream istr(strInputFile.c_str(), mode);
        spDockSite = RbtDockingSitePtr(new RbtDockingSite(istr));
        istr.close();
    }

    void loadDockingSiteFromMapper() {
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

    void setup(){
        loadReceptor();
        loadDockingSite();
    }

    void showDockingSite() {
        cout << endl << "DOCKING SITE" << endl << (*spDockSite) << endl;
    }

    void writeAS() {
        auto mode = addBinaryMode(ios_base::out | ios_base::trunc);
        ofstream ostr(strASFile.c_str(), mode);
        spDockSite->Write(ostr);
        ostr.close();
    }

    void writeMoeGrid() {
        cout << "MOE grid feature not yet implemented, sorry." << endl;
    }

    void writePSF() {
        RbtMolecularFileSinkPtr spRecepSink = new RbtPsfFileSink(wsName + "_for_viewer.psf", spReceptor);
        cout << "Writing PSF file: " << spRecepSink->GetFileName() << endl;
        spRecepSink->Render();
    }

    void writeCRD() {
        RbtMolecularFileSinkPtr spRecepSink = new RbtCrdFileSink(wsName + "_for_viewer.crd", spReceptor);
        cout << "Writing CRD file: " << spRecepSink->GetFileName() << endl;
        spRecepSink->Render();
    }

    void writeInisghtGrids() {
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

    void listReceptorAtoms() {
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
    void showSITEDescriptors(){
        SITEDescriptors(wsName, spDockSite, spReceptor).show();
    }


    void run() {
        setup();

        showDockingSite();
        if (config.writeAS) writeAS();

        // Write PSF/CRD files to keep the rDock Viewer happy (it doesn't read MOL2 files yet)
        if (config.viewer) {
            writePSF();
            writeCRD();
        }

        // Write an ASCII InsightII grid file for each defined cavity
        if (config.dump) writeInisghtGrids();
        // writing active site into MOE grid
        if (config.moeGrid) writeMoeGrid();
        // List all receptor atoms within given distance of any cavity
        if (config.list) listReceptorAtoms();
        if (config.site) showSITEDescriptors();

    }

};

void printHeader(const char *argv[]) {
    // Strip off the path to the executable, leaving just the file name
    RbtString strExeName(argv[0]);
    RbtString::size_type i = strExeName.rfind("/");
    if (i != RbtString::npos) strExeName.erase(0, i + 1);

    // Print a standard header
    Rbt::PrintStdHeader(cout, strExeName + EXEVERSION);
}

int main(int argc, const char *argv[]) {
    printHeader(argv);
    int exitCode = 0;
    try {
        auto config = RBCavityCommandLineParser().parse(argc, argv);
        showConfig(config);
        cout.setf(ios_base::left, ios_base::adjustfield);
        RBCavity(config).run();
    } catch (RbtError &e) {
        cout << e << endl;
        exitCode = 1;
    } catch (...) {
        cout << "Unknown exception" << endl;
        exitCode = 1;
    }

    _RBTOBJECTCOUNTER_DUMP_(cout)

    return exitCode;
}


std::vector<poptOption> RBCavityCommandLineParser::createOptionsTable() {
    return {
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
        POPT_AUTOHELP{NULL, 0, 0, NULL, 0}
    };
}

RBCavityCommandLineParser::RBCavityCommandLineParser():
    prmFile(nullptr),
    listDist(nullptr),
    borderDist(nullptr),
    optionsTable(createOptionsTable())
{}

std::set<char> RBCavityCommandLineParser::getFlags(poptContext & optCon) {
    std::set<char> flags;
    char c;
    while ((c = poptGetNextOpt(optCon)) >= 0){
        flags.emplace(c);
        cerr << "FLAG: " << c << endl;
        if (borderDist != nullptr) {
        cerr << string(borderDist) << endl;
        } else {
            cerr << "NULL" << endl;
        }
    }
    auto a = flags.find('r');
    return flags;
}

std::set<char> RBCavityCommandLineParser::getKnownFlags() {
    std::set<char> knownFlags;
    for (auto & option : optionsTable) {
        if (option.shortName) knownFlags.emplace(option.shortName);
    }
    return knownFlags;
}

RBCavityConfig RBCavityCommandLineParser::parse(int argc, const char * argv[]) {
    RBCavityConfig config{};
    poptContext optCon = poptGetContext(NULL, argc, argv, optionsTable.data(), 0);
    poptSetOtherOptionHelp(optCon, "-r<receptor.prm> [options]");
    if (argc < 2) {
        poptPrintUsage(optCon, stderr, 0);
        throw RbtError(_WHERE_, "Not enough parameters");
    }

    // extract configuration
    auto flags = getFlags(optCon);
    config.readAS = contains(flags, 'R');
    config.writeAS = contains(flags, 'W');
    config.dump = contains(flags, 'd');
    config.viewer = contains(flags, 'v');
    config.list = contains(flags, 'l');
    config.site = contains(flags, 's');
    config.moeGrid = contains(flags, 'm');
    showConfig(config);
    if (config.list) config.dist = atof(listDist);
    if (contains(flags, 'b')) config.border = atof(borderDist);
    if (!prmFile) {
        PrintUsage();
        throw RbtError(_WHERE_, "Missing receptor parameter file name");
    }
    config.prmFile = prmFile;

    
    // cleanup
    poptFreeContext(optCon);

    // warn the user about possible typos
    warnUnknownFlags(flags);

    return config;
}

void RBCavityCommandLineParser::warnUnknownFlags(std::set<char> & flags) {
    auto knownFlags = getKnownFlags();
    for (auto & flag : flags) {
        if (!contains(knownFlags, flag)) {
            cout << "WARNING: unknown argument: " << flag << endl;
        }
    }
}
