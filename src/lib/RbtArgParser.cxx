#include "RbtArgParser.h"

#include <cxxopts.hpp>
#include <iostream>
#include <string>
#include <vector>

namespace RbtArgParser {

RbtArgParser::RbtArgParser(const std::string &program, const std::string &description):
    parser(program, description),
    substitutions({}) {}

void RbtArgParser::add_flag(const std::string &opts, const std::string &dest, const char *default_value) {
    add<bool>(opts, dest, default_value);
}

cxxopts::parse_result RbtArgParser::parse(int argc, const char *argv[]) {
    auto fixed_args = preprocess_args(argc, argv);
    auto fixed_char_ptr_args = cast_args(fixed_args);
    return parser.parse(fixed_char_ptr_args.size(), fixed_char_ptr_args.data());
}

// all the functions below are a retrocompatibility layer to replace single dash long
// arguments with double dash long arguments. This is a deprecated feature and will be
// removed in the future.
std::string RbtArgParser::fix_value(const char *arg) {
    std::string arg_str(arg);
    for (auto &substitution: substitutions) {
        if (arg_str.find(substitution) == 0) {
            warn_deprecated_argument(arg);
            arg_str.replace(0, substitution.size(), "-" + substitution);
            return arg_str;
        }
    }
    return arg;
}

void RbtArgParser::warn_deprecated_argument(const char *arg) {
    std::cerr << "WARNING: Single dash long argument detected: '" << arg << "'."
              << " This is a deprecated feature and will be removed in the future."
              << " Please use either double dashed arguments or single dashed abbreviations instead." << std::endl;
}

std::vector<std::string> RbtArgParser::preprocess_args(int argc, const char *argv[]) {
    std::vector<std::string> args;
    for (int i = 0; i < argc; i++) {
        args.push_back(fix_value(argv[i]));
    }
    return args;
}

std::vector<const char *> RbtArgParser::cast_args(const std::vector<std::string> &args) const {
    std::vector<const char *> casted_args;
    for (const auto &arg: args) {
        casted_args.push_back(arg.c_str());
    }
    return casted_args;
}

// populate the substitution table for single dash long arguments automatically
// whenever a long option is added
void RbtArgParser::set_substitution(const std::string &opts) {
    auto comma_position = opts.find(",");
    std::string long_option_name;
    if (comma_position != std::string::npos) {
        long_option_name = opts.substr(comma_position + 1);
    } else if (opts.length() > 1) {
        long_option_name = opts;
    }
    if (!long_option_name.empty()) {
        substitutions.push_back({"-" + long_option_name});
    }
}

}  // namespace RbtArgParser
