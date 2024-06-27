#include <catch2/catch_amalgamated.hpp>

#include "rbcavity/rbcavity_argparser.h"
#include "rbcavity/rbcavity_config.h"

#include "test_utils.hpp"

#include <vector>
#include <tuple>


struct TestConfig{
    std::vector<std::string> input;
    RBCavity::RBCavityConfig expected_config;
};


TEST_CASE("rbcavity", "[argparser]") {
    auto config_base = RBCavity::RBCavityConfig{.strReceptorPrmFile = "receptor.prm"};
    auto config_with_list = RBCavity::RBCavityConfig{.strReceptorPrmFile = "receptor.prm", .bList = true, .dist = 3.0f};
    auto config_with_border = RBCavity::RBCavityConfig{.strReceptorPrmFile = "receptor.prm", .bBorder = true, .border = 5.0f};
    auto config_with_read_as = RBCavity::RBCavityConfig{.strReceptorPrmFile = "receptor.prm", .bReadAS = true};
    auto config_with_write_as = RBCavity::RBCavityConfig{.strReceptorPrmFile = "receptor.prm", .bWriteAS = true};
    auto config_with_dump_insight = RBCavity::RBCavityConfig{.strReceptorPrmFile = "receptor.prm", .bDump = true};
    auto config_with_viewer = RBCavity::RBCavityConfig{.strReceptorPrmFile = "receptor.prm", .bViewer = true};
    auto config_with_site = RBCavity::RBCavityConfig{.strReceptorPrmFile = "receptor.prm", .bSite = true};
    auto config_with_moe_grid = RBCavity::RBCavityConfig{.strReceptorPrmFile = "receptor.prm", .bMOEgrid = true};

    std::vector<std::string> all_settings{"-r", "receptor.prm", "-l", "3", "-b", "5", "-R", "-W", "-d", "-v", "-s", "-m"};
    auto config_all_settings = RBCavity::RBCavityConfig{
        .strReceptorPrmFile = "receptor.prm",
        .bReadAS = true,
        .bWriteAS = true,
        .bDump = true,
        .bViewer = true,
        .bList = true,
        .bBorder = true,
        .bSite = true,
        .bMOEgrid = true,
        .border = 5.0f,
        .dist = 3.0f,
    };

    std::vector<TestConfig> configs {
        TestConfig{.input = {"-r", "receptor.prm"}, .expected_config = config_base},
        TestConfig{.input = {"-rreceptor.prm"}, .expected_config = config_base},
        TestConfig{.input = {"-receptor=receptor.prm"}, .expected_config = config_base},
        TestConfig{.input = {"--receptor", "receptor.prm"}, .expected_config = config_base},
        TestConfig{.input = {"--receptor=receptor.prm"}, .expected_config = config_base},
        TestConfig{.input = {"-r", "receptor.prm", "-l", "3"}, .expected_config = config_with_list},
        TestConfig{.input = {"-r", "receptor.prm", "-l3"}, .expected_config = config_with_list},
        TestConfig{.input = {"-r", "receptor.prm", "-list=3"}, .expected_config = config_with_list},
        TestConfig{.input = {"-r", "receptor.prm", "--list", "3"}, .expected_config = config_with_list},
        TestConfig{.input = {"-r", "receptor.prm", "--list=3"}, .expected_config = config_with_list},
        TestConfig{.input = {"-r", "receptor.prm", "-b", "5"}, .expected_config = config_with_border},
        TestConfig{.input = {"-r", "receptor.prm", "-b5"}, .expected_config = config_with_border},
        TestConfig{.input = {"-r", "receptor.prm", "-border=5"}, .expected_config = config_with_border},
        TestConfig{.input = {"-r", "receptor.prm", "--border", "5"}, .expected_config = config_with_border},
        TestConfig{.input = {"-r", "receptor.prm", "--border=5"}, .expected_config = config_with_border},
        TestConfig{.input = {"-r", "receptor.prm", "-R"}, .expected_config = config_with_read_as},
        TestConfig{.input = {"-r", "receptor.prm", "-ras"}, .expected_config = config_with_read_as},
        TestConfig{.input = {"-r", "receptor.prm", "--ras"}, .expected_config = config_with_read_as},
        TestConfig{.input = {"-r", "receptor.prm", "-W"}, .expected_config = config_with_write_as},
        TestConfig{.input = {"-r", "receptor.prm", "-was"}, .expected_config = config_with_write_as},
        TestConfig{.input = {"-r", "receptor.prm", "--was"}, .expected_config = config_with_write_as},
        TestConfig{.input = {"-r", "receptor.prm", "-d"}, .expected_config = config_with_dump_insight},
        TestConfig{.input = {"-r", "receptor.prm", "-dump-insight"}, .expected_config = config_with_dump_insight},
        TestConfig{.input = {"-r", "receptor.prm", "--dump-insight"}, .expected_config = config_with_dump_insight},
        TestConfig{.input = {"-r", "receptor.prm", "-v"}, .expected_config = config_with_viewer},
        TestConfig{.input = {"-r", "receptor.prm", "-viewer"}, .expected_config = config_with_viewer},
        TestConfig{.input = {"-r", "receptor.prm", "--viewer"}, .expected_config = config_with_viewer},
        TestConfig{.input = {"-r", "receptor.prm", "-s"}, .expected_config = config_with_site},
        TestConfig{.input = {"-r", "receptor.prm", "-site"}, .expected_config = config_with_site},
        TestConfig{.input = {"-r", "receptor.prm", "--site"}, .expected_config = config_with_site},
        TestConfig{.input = {"-r", "receptor.prm", "-m"}, .expected_config = config_with_moe_grid},
        TestConfig{.input = {"-r", "receptor.prm", "-dump-moe"}, .expected_config = config_with_moe_grid},
        TestConfig{.input = {"-r", "receptor.prm", "--dump-moe"}, .expected_config = config_with_moe_grid},
        TestConfig{.input = all_settings, .expected_config = config_all_settings}
    };
    NullBuffer null_buffer;
    CoutRedirect redirect_cout(&null_buffer);
    CerrRedirect redirect_cerr(&null_buffer);
    
    for (auto test_config: configs) {
        auto inputs = std::vector<const char *>{"rbcavity"};
        for (auto & input: test_config.input) inputs.push_back(input.c_str());
        auto result = RBCavity::parse_args(inputs.size(), inputs.data());
        REQUIRE(result == test_config.expected_config);
    }
}
