#include "rbdock/rbdock_argparser.h"

#include <iostream>
#include <string>

#include "rbdock/rbdock_config.h"

#include "RbtArgParser.h"
#include "RbtError.h"

RbtArgParser::RbtArgParser RBDock::get_options_parser() {
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

RBDock::RBDockConfig RBDock::parse_args(int argc, const char *argv[]) {
    auto parser = get_options_parser();
    try {
        auto parser_result = RbtArgParser::RbtParseResult(parser.parse(argc, argv));
        RBDock::RBDockConfig config;
        parser_result["input"] >> config.strLigandMdlFile;
        parser_result["receptor"] >> config.strReceptorPrmFile;
        if (parser_result["protocol"].is_present()) {
            parser_result["protocol"] >> config.strParamFile;
        }
        config.bOutput = parser_result["output-root"].is_present();
        parser_result["output-root"] >> config.strRunName;
        parser_result["runs"] >> config.nDockingRuns;
        parser_result["ap"] >> config.bPosIonise;
        parser_result["an"] >> config.bNegIonise;
        parser_result["allH"] >> config.bAllH;
        config.bSeed = parser_result["seed"].is_present();
        if (config.bSeed) parser_result["seed"] >> config.nSeed;
        config.bTrace = parser_result["trace"].is_present();
        if (config.bTrace) parser_result["trace"] >> config.iTrace;
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