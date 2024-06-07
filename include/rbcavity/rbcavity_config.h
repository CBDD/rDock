#ifndef _RBCAVITY_CONFIG_H_
#define _RBCAVITY_CONFIG_H_

#include <iostream>
#include <string>

#include "RbtValidationError.h"

namespace RBCavity {

struct RBCavityConfig {
    std::string strReceptorPrmFile;
    bool bReadAS = false;   // If true, read Active Site from file
    bool bWriteAS = false;  // If true, write Active Site to file
    bool bDump = false;     // If true, dump cavity grids in Insight format
    bool bViewer = false;   // If true, dump PSF/CRD files for rDock Viewer
    bool bList = false;     // If true, list atoms within 'distance' of cavity
    bool bBorder = false;   // If true, set the border around the cavities for the distance grid
    bool bSite = false;     // If true, print out "SITE" descriptors  (counts of exposed atoms)
    bool bMOEgrid = false;  // If true, create a MOE grid file for AS visualisation
    float border = 8.0;     // Border to allow around cavities for distance grid
    float dist = 5.0;       // Distance to cavity for atom listing

    friend std::ostream &operator<<(std::ostream &os, const RBCavityConfig &config) {
        os << "Command line arguments:" << std::endl;
        os << "-r " << config.strReceptorPrmFile << std::endl;
        if (config.bList) os << "-l " << config.dist << std::endl;
        if (config.bBorder) os << "-b " << config.border << std::endl;
        if (config.bWriteAS) os << "--was" << std::endl;
        if (config.bReadAS) os << "--ras" << std::endl;
        if (config.bMOEgrid) os << "-m" << std::endl;
        if (config.bDump) os << "-d" << std::endl;
        if (config.bSite) os << "-s" << std::endl;
        if (config.bViewer) os << "-v" << std::endl;
        return os;
    }

    inline void validate() const {
        using Rbt::ValidationError;
        if (strReceptorPrmFile.empty()) throw ValidationError("Missing receptor parameter file name");
        if (bList && dist <= 0) throw ValidationError("Invalid distance to cavity. must be a positive number");
        if (bBorder && border <= 0) throw ValidationError("Invalid border distance. must be a positive number");
    }

    friend bool operator==(const RBCavity::RBCavityConfig & c1, const RBCavity::RBCavityConfig & c2) {
        return (
            c1.strReceptorPrmFile == c2.strReceptorPrmFile
            && c1.bReadAS == c2.bReadAS
            && c1.bWriteAS == c2.bWriteAS
            && c1.bDump == c2.bDump
            && c1.bViewer == c2.bViewer
            && c1.bList == c2.bList
            && c1.bBorder == c2.bBorder
            && c1.bSite == c2.bSite
            && c1.bMOEgrid == c2.bMOEgrid
            && c1.border == c2.border
            && c1.dist == c2.dist
        );
    }
};

}  // namespace RBCavity

#endif  // _RBCAVITY_CONFIG_H_