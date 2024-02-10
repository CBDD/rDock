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

#include <algorithm>
#include <cxxopts.hpp>
#include <iomanip>

#include "RbtArgParser.h"
#include "RbtBiMolWorkSpace.h"
#include "RbtCrdFileSink.h"
#include "RbtDockingSite.h"
#include "RbtPRMFactory.h"
#include "RbtParameterFileSource.h"
#include "RbtPlatformCompatibility.h"
#include "RbtPsfFileSink.h"
#include "RbtSiteMapperFactory.h"
#include "RbtVersion.h"

const RbtString EXEVERSION = RBT_VERSION;

cxxopts::options get_options_parser() {
    using RbtArgParser::add_flag;
    using RbtArgParser::add_float;
    using RbtArgParser::add_string;

    cxxopts::options options("rbcavity", "calculate docking cavities");
    add_string(options, "r,receptor", "receptor param file (contains active site params)");
    add_float(options, "l,list", "list receptor atoms within a distance in angstrom of any cavity", "5.0f");
    add_float(options, "b,border", "set the border (in angstrom) around the cavities for the distance grid", "8.0f");
    add_flag(options, "W,was", "write docking cavities (plus distance grid) to .as file");
    add_flag(options, "R,ras", "read docking cavities (plus distance grid) from .as file");
    add_flag(options, "d,dump-insight", "dump InsightII/PyMOL grids for each cavity for visualisation");
    add_flag(options, "v,viewer", "dump target PSF/CRD files for rDock Viewer");
    add_flag(options, "s,site", "print SITE descriptors (counts of exposed atoms)");
    add_flag(options, "m", "write active site into a MOE grid");
    return options;
}

struct RBCavityConfig {
    RbtString strReceptorPrmFile;
    RbtBool bReadAS = false;   // If true, read Active Site from file
    RbtBool bWriteAS = false;  // If true, write Active Site to file
    RbtBool bDump = false;     // If true, dump cavity grids in Insight format
    RbtBool bViewer = false;   // If true, dump PSF/CRD files for rDock Viewer
    RbtBool bList = false;     // If true, list atoms within 'distance' of cavity
    RbtBool bBorder = false;   // If true, set the border around the cavities for the distance grid
    RbtBool bSite = false;     // If true, print out "SITE" descriptors  (counts of exposed atoms)
    RbtBool bMOEgrid = false;  // If true, create a MOE grid file for AS visualisation
    RbtDouble border = 8.0;    // Border to allow around cavities for distance grid
    RbtDouble dist = 5.0;      // Distance to cavity for atom listing

    class ValidationError: public std::runtime_error {
     public:
        using std::runtime_error::runtime_error;
    };

    friend std::ostream &operator<<(std::ostream &os, const RBCavityConfig &config);

    void validate() const {
        if (strReceptorPrmFile.empty()) throw ValidationError("Missing receptor parameter file name");
        if (bList && dist <= 0) throw ValidationError("Invalid distance to cavity. must be a positive number");
        if (bBorder && border <= 0) throw ValidationError("Invalid border distance. must be a positive number");
    }
};

std::ostream &operator<<(std::ostream &os, const RBCavityConfig &config) {
    os << "Command line arguments:" << endl;
    os << "-r " << config.strReceptorPrmFile << endl;
    if (config.bList) os << "-l " << config.dist << endl;
    if (config.bBorder) os << "-b " << config.border << endl;
    if (config.bWriteAS) os << "-was" << endl;
    if (config.bReadAS) os << "-ras" << endl;
    if (config.bMOEgrid) os << "-m" << endl;
    if (config.bDump) os << "-d" << endl;
    if (config.bSite) os << "-s" << endl;
    if (config.bViewer) os << "-v" << endl;
    return os;
}

// all this is just for retrocompatibility with the original rbcavity
// we need to replace the arguments with the long version
const RbtArgParser::ArgsSubstitutions ARGS_SUBSTITUTIONS = {
    {"-was", "--was"},
    {"-ras", "--ras"},
    {"-receptor", "--receptor"},
    {"-dump-insight", "--dump-insight"},
    {"-viewer", "--viewer"},
    {"-list", "--list"},
    {"-site", "--site"},
    {"-border", "--border"},
};

RBCavityConfig parse_args(int argc, const char *argv[]) {
    auto options = get_options_parser();
    auto args = RbtArgParser::preprocessArgs(argc, argv, ARGS_SUBSTITUTIONS);
    try {
        auto result = options.parse(argc, args.data());
        RBCavityConfig config;
        config.strReceptorPrmFile = result["receptor"].as<std::string>();
        config.bReadAS = result["ras"].as<bool>();
        config.bWriteAS = result["was"].as<bool>();
        config.bDump = result["dump-insight"].as<bool>();
        config.bViewer = result["viewer"].as<bool>();
        config.bMOEgrid = result["m"].as<bool>();
        config.bBorder = result.count("border") > 0;
        if (config.bBorder) config.border = result["border"].as<float>();
        config.bList = result.count("list") > 0;
        if (config.bList) config.dist = result["list"].as<float>();
        config.bSite = result["site"].as<bool>();

        config.validate();
        return config;
    } catch (cxxopts::parse_error &e) {
        std::cerr << "Error parsing options: " << e.what() << std::endl;
    } catch (RBCavityConfig::ValidationError &e) {
        std::cerr << "Invalid configuration: " << e.what() << std::endl;
    }
    // if we reach this point, something went wrong. Print the help and exit
    std::cerr << options.help() << std::endl;
    exit(1);
}

void RBCavity(const RBCavityConfig & config) {
    // Create a bimolecular workspace
    RbtBiMolWorkSpacePtr spWS(new RbtBiMolWorkSpace());
    // Set the workspace name to the root of the receptor .prm filename
    RbtStringList componentList = Rbt::ConvertDelimitedStringToList(config.strReceptorPrmFile, ".");
    RbtString wsName = componentList.front();
    spWS->SetName(wsName);

    // Read the protocol parameter file
    RbtParameterFileSourcePtr spRecepPrmSource(
        new RbtParameterFileSource(Rbt::GetRbtFileName("data/receptors", config.strReceptorPrmFile))
    );

    // Create the receptor model from the file names in the parameter file
    spRecepPrmSource->SetSection();
    RbtPRMFactory prmFactory(spRecepPrmSource);
    RbtModelPtr spReceptor = prmFactory.CreateReceptor();

    RbtDockingSitePtr spDockSite;
    RbtString strASFile = wsName + ".as";

    // Either read the docking site from the .as file
    if (config.bReadAS) {
        RbtString strInputFile = Rbt::GetRbtFileName("data/grids", strASFile);
        ifstream istr(strInputFile.c_str(), Rbt::inputMode);
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

    if (config.bWriteAS) {
        ofstream ostr(strASFile.c_str(), Rbt::outputMode);
        spDockSite->Write(ostr);
        ostr.close();
    }

    // Write PSF/CRD files to keep the rDock Viewer happy (it doesn't read MOL2 files yet)
    if (config.bViewer) {
        RbtMolecularFileSinkPtr spRecepSink = new RbtPsfFileSink(wsName + "_for_viewer.psf", spReceptor);
        cout << "Writing PSF file: " << spRecepSink->GetFileName() << endl;
        spRecepSink->Render();
        spRecepSink = new RbtCrdFileSink(wsName + "_for_viewer.crd", spReceptor);
        cout << "Writing CRD file: " << spRecepSink->GetFileName() << endl;
        spRecepSink->Render();
    }

    // Write an ASCII InsightII grid file for each defined cavity
    if (config.bDump) {
        RbtCavityList cavList = spDockSite->GetCavityList();
        for (RbtUInt i = 0; i < cavList.size(); i++) {
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
    if (config.bMOEgrid) {
        cout << "MOE grid feature not yet implemented, sorry." << endl;
    }
    // List all receptor atoms within given distance of any cavity
    if (config.bList) {
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
    if (config.bSite) {
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

/////////////////////////////////////////////////////////////////////
// MAIN PROGRAM STARTS HERE
/////////////////////////////////////////////////////////////////////

int main(int argc, const char *argv[]) {
    // Strip off the path to the executable, leaving just the file name
    RbtString exeFullPath(argv[0]);
    RbtString strExeName = exeFullPath.substr(exeFullPath.find_last_of("/\\") + 1);
    Rbt::PrintStdHeader(cout, strExeName + " - " + EXEVERSION);

    RBCavityConfig config = parse_args(argc, argv);
    // writing command line arguments
    std::cout << config << std::endl;

    cout.setf(ios_base::left, ios_base::adjustfield);

    try {
        RBCavity(config);
    } catch (RbtError &e) {
        cout << e << endl;
    } catch (...) {
        cout << "Unknown exception" << endl;
    }

    _RBTOBJECTCOUNTER_DUMP_(cout)

    return 0;
}
