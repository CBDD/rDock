#ifndef _RBDOCK_CONFIG_H_
#define _RBDOCK_CONFIG_H_

#include <iostream>
#include <optional>
#include <sstream>
#include <string>

#include "RbtValidationError.h"

namespace RBDock {

struct RBDockConfig {
    // mandatory parameters
    std::string strLigandMdlFile;
    std::string strReceptorPrmFile;
    std::string strParamFile;

    // optional parameters
    std::optional<std::string> strOutputPrefix;
    std::optional<std::string> strFilterFile;
    std::optional<int> nDockingRuns;
    std::optional<double> dTargetScore;
    std::optional<int> nSeed;
    std::optional<int> iTrace;

    // flags
    bool bContinue = false;
    bool bPosIonise = false;
    bool bNegIonise = false;
    bool bAllH = false;

    friend std::ostream &operator<<(std::ostream &os, const RBDockConfig &config) {
        os << "Command line args:" << std::endl;
        os << " -i " << config.strLigandMdlFile << std::endl;
        os << " -r " << config.strReceptorPrmFile << std::endl;
        os << " -p " << config.strParamFile << std::endl;
        if (config.strOutputPrefix.has_value()) {
            os << " -o " << config.strOutputPrefix.value() << std::endl;
        } else {
            os << "WARNING: output file name is missing." << std::endl;
        }
        auto def_aux = config.nDockingRuns.has_value() ? " (default) " : "";
        os << " -n " << config.nDockingRuns.value_or(1) << def_aux << std::endl;
        if (config.nSeed.has_value()) os << " -s " << config.nSeed.value() << std::endl;
        if (config.iTrace.has_value()) os << " -T " << config.iTrace.value() << std::endl;
        if (config.bPosIonise) os << " --ap " << std::endl;
        if (config.bNegIonise) os << " --an " << std::endl;
        if (!config.bAllH) os << " --allH " << std::endl;
        if (config.dTargetScore.has_value()) os << " -t " << config.dTargetScore.value() << std::endl;
        if (config.bContinue) os << " --cont " << std::endl;
        if (config.bPosIonise)
            os << "Automatically protonating positive ionisable groups (amines, imidazoles, guanidines)" << std::endl;
        if (config.bNegIonise)
            os << "Automatically deprotonating negative ionisable groups (carboxylic acids, phosphates, sulphates, "
                  "sulphonates)"
               << std::endl;
        if (!config.bAllH)
            os << "Reading polar hydrogens only from ligand SD file" << std::endl;
        else
            os << "Reading all hydrogens from ligand SD file" << std::endl;

        if (config.dTargetScore.has_value()) {
            os << std::endl << "Lower target intermolecular score = " << config.dTargetScore.value() << std::endl;
        }
        if (config.strFilterFile.has_value()) os << " -f " << config.strFilterFile.value() << std::endl;

        return os;
    }

    void validate() {
        if (strLigandMdlFile.empty()) throw Rbt::ValidationError("input ligand SD file is mandatory");
        if (strReceptorPrmFile.empty()) throw Rbt::ValidationError("receptor parameter file is mandatory");
        if (nDockingRuns.has_value() && nDockingRuns < 1)
            throw Rbt::ValidationError("number of runs must be greater than 0");
    }

    bool operator==(const RBDockConfig &rhs) const {
        return strLigandMdlFile == rhs.strLigandMdlFile && strOutputPrefix == rhs.strOutputPrefix
               && strReceptorPrmFile == rhs.strReceptorPrmFile && strParamFile == rhs.strParamFile
               && strFilterFile == rhs.strFilterFile && nDockingRuns == rhs.nDockingRuns && bContinue == rhs.bContinue
               && dTargetScore == rhs.dTargetScore && bPosIonise == rhs.bPosIonise && bNegIonise == rhs.bNegIonise
               && bAllH == rhs.bAllH && nSeed == rhs.nSeed && iTrace == rhs.iTrace;
    }

    std::string get_filter_string() const {
        std::ostringstream strAuxFilter;
        if (strFilterFile.has_value()) return "";
        if (dTargetScore.has_value()) {       // -t<TS>
            if (!nDockingRuns.has_value()) {  // -t<TS> only
                strAuxFilter << "0 1 - SCORE.INTER " << dTargetScore.value() << std::endl;
            } else  // -t<TS> -n<N> need to check if -cont present
                    // for all other cases it doesn't matter
                if (bContinue) {  // -t<TS> -n<N> -cont
                    strAuxFilter << "1 if - SCORE.NRUNS " << (nDockingRuns.value() - 1)
                                << " 0.0 -1.0,\n1 - SCORE.INTER " << dTargetScore.value() << std::endl;
                } else {  // -t<TS> -n<N>
                    strAuxFilter << "1 if - " << dTargetScore.value() << " SCORE.INTER 0.0 "
                                << "if - SCORE.NRUNS " << (nDockingRuns.value() - 1)
                                << " 0.0 -1.0,\n1 - SCORE.INTER " << dTargetScore.value() << std::endl;
                }
        }                                     // no target score, no filter
        else if (nDockingRuns.has_value()) {  // -n<N>
            strAuxFilter << "1 if - SCORE.NRUNS " << (nDockingRuns.value() - 1) << " 0.0 -1.0,\n0";
        } else {  // no -t no -n
            strAuxFilter << "0 0\n";
        }
    return strAuxFilter.str();
    }
};

}  // namespace RBDock

#endif  // _RBDOCK_CONFIG_H_
