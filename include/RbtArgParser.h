#ifndef _RBTARGPARSER_H_
#define _RBTARGPARSER_H_

#include <cxxopts.hpp>
#include <string>
#include <vector>

namespace RbtArgParser {

using parsing_error = cxxopts::exceptions::parsing;
typedef std::vector<std::string> ArgsSubstitutions;

class RbtArgParser {
 public:
    cxxopts::Options parser;
    ArgsSubstitutions substitutions;

 public:
    RbtArgParser(const std::string &program, const std::string &description = "");

    cxxopts::ParseResult parse(int argc, const char *argv[]);
    inline std::string help() const { return parser.help(); }

    // convenience functions to improve readability when adding options
    // add a boolean flag to the options
    inline void add_flag(const std::string &opts, const std::string &dest, const char *default_value = "false") {
        add<bool>(opts, dest, default_value);
    }

    // add whatever type of option
    template <typename T>
    void add(const std::string &opts, const std::string &dest, const char *default_value = nullptr) {
        add_substitution(opts);
        if (default_value == nullptr) {
            parser.add_options()(opts, dest, cxxopts::value<T>());
        } else {
            parser.add_options()(opts, dest, cxxopts::value<T>()->default_value(default_value));
        }
    }

    // all the functions below are a retrocompatibility layer to replace single dash long
    // arguments with double dash long arguments.
    // This is a deprecated feature and will be removed in the future, and so are these functions.
    std::string fix_value(const char *arg);
    void warn_deprecated_argument(const char *arg);
    std::vector<std::string> preprocess_args(int argc, const char *argv[]);
    std::vector<const char *> cast_args(const std::vector<std::string> &args) const;
    void add_substitution(const std::string &opts);
};

class ValidationError: public std::runtime_error {
 public:
    using std::runtime_error::runtime_error;
};

}  // namespace RbtArgParser
#endif  //_RBTARGPARSER_H_