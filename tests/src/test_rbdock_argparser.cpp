#include <catch2/catch_amalgamated.hpp>
#include <tuple>
#include <vector>

#include "rbdock/rbdock_argparser.h"
#include "rbdock/rbdock_config.h"
#include "test_utils.hpp"

using RBDock::RBDockConfig;

typedef std::tuple<std::vector<std::string>, RBDockConfig> TestConfig;

TEST_CASE("rbdock", "[argparser]") {
    std::vector<TestConfig> configs{
        // base config
        {{"-r", "receptor.prm", "-i", "ligand.sdf"},
         RBDockConfig{.strLigandMdlFile = "ligand.sdf", .strReceptorPrmFile = "receptor.prm"}},
        // config with output
        {{"-r", "receptor.prm", "-i", "ligand.sdf", "-o", "output"},

         RBDockConfig{
             .strLigandMdlFile = "ligand.sdf",
             .strReceptorPrmFile = "receptor.prm",
             .strOutputPrefix = "output",
         }},
        // config with protocol
        {{"-r", "receptor.prm", "-i", "ligand.sdf", "-p", "protocol.prm"},

         RBDockConfig{
             .strLigandMdlFile = "ligand.sdf", .strReceptorPrmFile = "receptor.prm", .strParamFile = "protocol.prm"}},
        // config with filter
        {{"-r", "receptor.prm", "-i", "ligand.sdf", "-t", "filter.txt"},

         RBDockConfig{
             .strLigandMdlFile = "ligand.sdf",
             .strReceptorPrmFile = "receptor.prm",
             .strFilterFile = "filter.txt",
         }},
        // config with target score
        {{"-r", "receptor.prm", "-i", "ligand.sdf", "-t", "0.5"},

         RBDockConfig{.strLigandMdlFile = "ligand.sdf", .strReceptorPrmFile = "receptor.prm", .dTargetScore = 0.5}},
        // config with continue
        {{"-r", "receptor.prm", "-i", "ligand.sdf", "-C"},
         RBDockConfig{.strLigandMdlFile = "ligand.sdf", .strReceptorPrmFile = "receptor.prm", .bContinue = true}},
        // config with docking runs
        {{"-r", "receptor.prm", "-i", "ligand.sdf", "-n", "5"},
         RBDockConfig{.strLigandMdlFile = "ligand.sdf", .strReceptorPrmFile = "receptor.prm", .nDockingRuns = 5}},
        // config with protonate
        {{"-r", "receptor.prm", "-i", "ligand.sdf", "-P"},
         RBDockConfig{.strLigandMdlFile = "ligand.sdf", .strReceptorPrmFile = "receptor.prm", .bPosIonise = true}},
        {{"-r", "receptor.prm", "-i", "ligand.sdf", "-ap"},
         RBDockConfig{.strLigandMdlFile = "ligand.sdf", .strReceptorPrmFile = "receptor.prm", .bPosIonise = true}},
        // config with deprotonate
        {{"-r", "receptor.prm", "-i", "ligand.sdf", "-D"},
         RBDockConfig{.strLigandMdlFile = "ligand.sdf", .strReceptorPrmFile = "receptor.prm", .bNegIonise = true}},
        {{"-r", "receptor.prm", "-i", "ligand.sdf", "-an"},

         RBDockConfig{.strLigandMdlFile = "ligand.sdf", .strReceptorPrmFile = "receptor.prm", .bNegIonise = true}},
        // config with all hydrogens
        {{"-r", "receptor.prm", "-i", "ligand.sdf", "-H"},

         RBDockConfig{.strLigandMdlFile = "ligand.sdf", .strReceptorPrmFile = "receptor.prm", .bAllH = true}},
        {{"-r", "receptor.prm", "-i", "ligand.sdf", "-allH"},

         RBDockConfig{.strLigandMdlFile = "ligand.sdf", .strReceptorPrmFile = "receptor.prm", .bAllH = true}},
        // config with seed
        {{"-r", "receptor.prm", "-i", "ligand.sdf", "-s", "123"},

         RBDockConfig{.strLigandMdlFile = "ligand.sdf", .strReceptorPrmFile = "receptor.prm", .nSeed = 123}},
        // config with trace
        {{"-r", "receptor.prm", "-i", "ligand.sdf", "-T", "2"},

         RBDockConfig{.strLigandMdlFile = "ligand.sdf", .strReceptorPrmFile = "receptor.prm", .iTrace = 2}},
        // config with filter
        {{"-r", "receptor.prm", "-i", "ligand.sdf", "-f", "filter.txt"},
         RBDockConfig{
             .strLigandMdlFile = "ligand.sdf", .strReceptorPrmFile = "receptor.prm", .strFilterFile = "filter.txt"}},

        // config with all options
        {
            {"-r", "receptor.prm", "-i", "ligand.sdf", "-o", "output", "-p", "protocol.prm", "-t", "0.5",
             "-C", "-n",           "5",  "-P",         "-D", "-H",     "-s", "123",          "-T", "2"},

            RBDockConfig{
                .strLigandMdlFile = "ligand.sdf",
                .strReceptorPrmFile = "receptor.prm",
                .strParamFile = "protocol.prm",
                .strOutputPrefix = "output",
                .nDockingRuns = 5,
                .dTargetScore = 0.5,
                .nSeed = 123,
                .iTrace = 2,
                .bContinue = true,
                .bPosIonise = true,
                .bNegIonise = true,
                .bAllH = true,
            },
        },
    };
    NullBuffer null_buffer;
    CoutRedirect redirect_cout(&null_buffer);
    CerrRedirect redirect_cerr(&null_buffer);

    auto [test_input, expected_config] = GENERATE_COPY(from_range(configs));

    auto inputs = std::vector<const char *>{"rbdock"};
    for (auto &input: test_input) inputs.push_back(input.c_str());
    std::cerr << "inputs: ";
    for (auto &input: inputs) std::cerr << input << " ";
    std::cerr << std::endl;
    auto result = RBDock::parse_args(inputs.size(), inputs.data());
    REQUIRE(result == expected_config);
}
