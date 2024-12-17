#ifndef _RBDOCK_CONFIG_H_
#define _RBDOCK_CONFIG_H_

#include <iostream>
#include <string>

#include "RbtValidationError.h"

namespace RBDock {

struct RBDockConfig {
    std::string strLigandMdlFile;
    bool bOutput = false;
    std::string strRunName;
    std::string strReceptorPrmFile;
    std::string strParamFile;
    std::string strFilterFile;
    int nDockingRuns = 1;
    bool bTarget = false;
    bool bContinue = false;
    bool bDockingRuns = false;
    double dTargetScore;
    bool bFilter = false;
    bool bPosIonise = false;
    bool bNegIonise = false;
    bool bAllH = false;
    bool bSeed = false;
    int nSeed = 0;
    bool bTrace = false;
    int iTrace = 0;

    friend std::ostream &operator<<(std::ostream &os, const RBDockConfig &config) {
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

        if (config.bTarget) {
            os << std::endl << "Lower target intermolecular score = " << config.dTargetScore << std::endl;
        }
        return os;
    }

    void validate() {
        if (strLigandMdlFile.empty()) throw Rbt::ValidationError("input ligand SD file is mandatory");
        if (strReceptorPrmFile.empty()) throw Rbt::ValidationError("receptor parameter file is mandatory");
        if (nDockingRuns < 1) throw Rbt::ValidationError("number of runs must be greater than 0");
    }

    bool operator==(const RBDockConfig &rhs) const {
        return strLigandMdlFile == rhs.strLigandMdlFile && bOutput == rhs.bOutput && strRunName == rhs.strRunName
               && strReceptorPrmFile == rhs.strReceptorPrmFile && strParamFile == rhs.strParamFile
               && strFilterFile == rhs.strFilterFile && nDockingRuns == rhs.nDockingRuns && bTarget == rhs.bTarget
               && bContinue == rhs.bContinue && bDockingRuns == rhs.bDockingRuns && dTargetScore == rhs.dTargetScore
               && bFilter == rhs.bFilter && bPosIonise == rhs.bPosIonise && bNegIonise == rhs.bNegIonise
               && bAllH == rhs.bAllH && bSeed == rhs.bSeed && nSeed == rhs.nSeed && bTrace == rhs.bTrace
               && iTrace == rhs.iTrace;
    }
};

}  // namespace RBDock

#endif  // _RBDOCK_CONFIG_H_
