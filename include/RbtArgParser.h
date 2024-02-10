#ifndef _RBTARGPARSER_H_
#define _RBTARGPARSER_H_

#include <cxxopts.hpp>
#include <iostream>
#include <string>
#include <vector>

namespace RbtArgParser {

typedef std::vector<std::pair<std::string, const char *>> ArgsSubstitutions;

// convenience function to add a boolean flag to the options
// greatly improves readability
void add_flag(
    cxxopts::options &options, const std::string &opts, const std::string &dest, const char *default_value = "false"
) {
    options.add_options()(opts, dest, cxxopts::value<bool>()->default_value(default_value));
}

template <typename T>
void add(
    cxxopts::options &options, const std::string &opts, const std::string &dest, const char *default_value = nullptr
) {
    if (default_value == nullptr) {
        options.add_options()(opts, dest, cxxopts::value<T>());
    } else {
        options.add_options()(opts, dest, cxxopts::value<T>()->default_value(default_value));
    }
}

// retrocompatibility function to replace single dash long arguments with double dash long arguments
const char *replace_value(const char *arg, ArgsSubstitutions substitutions) {
    for (auto &substitution: substitutions) {
        if (std::string(arg) == substitution.first) {
            std::cerr << "WARNING: Single dash long argument detected: '" << arg << "'."
                      << " This is a deprecated feature and will be removed in the future."
                      << " Please use double dash long arguments instead." << std::endl;
            return substitution.second;
        }
    }
    return arg;
}

std::vector<const char *> preprocessArgs(int argc, const char *argv[], const ArgsSubstitutions &substitutions) {
    std::vector<const char *> args;
    for (int i = 0; i < argc; i++) {
        args.push_back(replace_value(argv[i], substitutions));
    }
    return args;
}

}  // namespace RbtArgParser
#endif  //_RBTARGPARSER_H_
