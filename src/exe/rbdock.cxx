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
#include "rbdock/rbdock_argparser.h"
#include "rbdock/rbdock_config.h"
#include "rbdock/rbdock_main.h"

#include "Rbt.h"
#include "RbtDebug.h"
#include "RbtVersion.h"

const RbtString EXEVERSION = RBT_VERSION;

/////////////////////////////////////////////////////////////////////
// MAIN PROGRAM STARTS HERE
/////////////////////////////////////////////////////////////////////

int main(int argc, const char *argv[]) {
    // Strip off the path to the executable, leaving just the file name
    std::string exeFullPath(argv[0]);
    std::string strExeName = exeFullPath.substr(exeFullPath.find_last_of("/\\") + 1);
    Rbt::PrintStdHeader(std::cout, strExeName + " - " + EXEVERSION);

    try {
        RBDock::RBDockConfig config = RBDock::parse_args(argc, argv);
        std::cout << config << std::endl;
        std::cout.setf(std::ios_base::left, std::ios_base::adjustfield);
        RBDock::RBDock(config, strExeName);
    } catch (RbtError &e) {
        std::cerr << e << std::endl;
        return 1;
    } catch (std::exception &e) {
        std::cerr << "Unknown exception" << std::endl;
        std::cerr << typeid(e).name() << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    }

    _RBTOBJECTCOUNTER_DUMP_(std::cout)
    return 0;
}
