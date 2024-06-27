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

#include <iomanip>
#include <iostream>
#include <string>

#include "rbcavity/rbcavity_argparser.h"
#include "rbcavity/rbcavity_config.h"
#include "rbcavity/rbcavity_main.h"

#include "Rbt.h"
#include "RbtDebug.h"
#include "RbtVersion.h"

const std::string EXEVERSION = RBT_VERSION;

using std::cerr;
using std::cout;
using std::endl;
/////////////////////////////////////////////////////////////////////
// MAIN PROGRAM STARTS HERE
/////////////////////////////////////////////////////////////////////

int main(int argc, const char *argv[]) {
    // Strip off the path to the executable, leaving just the file name
    std::string exeFullPath(argv[0]);
    std::string strExeName = exeFullPath.substr(exeFullPath.find_last_of("/\\") + 1);
    Rbt::PrintStdHeader(cout, strExeName + " - " + EXEVERSION);

    try {
        RBCavity::RBCavityConfig config = RBCavity::parse_args(argc, argv);
        // writing command line arguments
        cout << config << endl;
        cout.setf(std::ios_base::left, std::ios_base::adjustfield);
        RBCavity::RBCavity(config);
    } catch (RbtError &e) {
        cerr << e << endl;
        return 1;
    } catch (std::exception &e) {
        cerr << "Unknown exception" << endl;
        cerr << typeid(e).name() << endl;
        cerr << e.what() << endl;
        return 1;
    }
    _RBTOBJECTCOUNTER_DUMP_(cout)
    return 0;
}
