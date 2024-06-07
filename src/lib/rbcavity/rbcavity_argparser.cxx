#include "rbcavity/rbcavity_argparser.h"

#include <iostream>
#include <string>

#include "rbcavity/rbcavity_config.h"

#include "RbtArgParser.h"
#include "RbtError.h"

RbtArgParser::RbtArgParser RBCavity::get_options_parser() {
    using std::string;

    RbtArgParser::RbtArgParser parser("rbcavity", "calculate docking cavities");
    parser.add<string>("r,receptor", "receptor param file (contains active site params)");
    parser.add<float>("l,list", "list receptor atoms within a distance in angstrom of any cavity", "5.0f");
    parser.add<float>("b,border", "set the border (in angstrom) around the cavities for the distance grid", "8.0f");
    parser.add_flag("W,was", "write docking cavities (plus distance grid) to .as file");
    parser.add_flag("R,ras", "read docking cavities (plus distance grid) from .as file");
    parser.add_flag("d,dump-insight", "dump InsightII/PyMOL grids for each cavity for visualisation");
    parser.add_flag("v,viewer", "dump target PSF/CRD files for rDock Viewer");
    parser.add_flag("s,site", "print SITE descriptors (counts of exposed atoms)");
    parser.add_flag("m", "write active site into a MOE grid");
    return parser;
}

RBCavity::RBCavityConfig RBCavity::parse_args(int argc, const char *argv[]) {
    auto parser = RBCavity::get_options_parser();
    try {
        auto arguments = RbtArgParser::RbtParseResult(parser.parse(argc, argv));
        RBCavity::RBCavityConfig config{};
        arguments["receptor"] >> config.strReceptorPrmFile;
        arguments["ras"] >> config.bReadAS;
        arguments["was"] >> config.bWriteAS;
        arguments["dump-insight"] >> config.bDump;
        arguments["viewer"] >> config.bViewer;
        arguments["m"] >> config.bMOEgrid;
        config.bBorder = arguments["border"].is_present();
        if (config.bBorder) arguments["border"] >> config.border;
        config.bList = arguments["list"].is_present();
        if (config.bList) arguments["list"] >> config.dist;
        arguments["site"] >> config.bSite;

        config.validate();
        return config;
    } catch (RbtArgParser::ParsingError &e) {
        std::cerr << "Error parsing options: " << e.what() << std::endl;
    } catch (RbtArgParser::ValidationError &e) {
        std::cerr << "Invalid configuration: " << e.what() << std::endl;
    }
    // if we reach this point, something went wrong. Print the help and exit
    std::cerr << parser.help() << std::endl;
    throw RbtError(_WHERE_, "Log Expressions only have 1 operand");
}
