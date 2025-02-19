#ifndef _RBTARGPARSER_H_
#define _RBTARGPARSER_H_

#include <cxxopts.hpp>
#include <optional>
#include <string>
#include <vector>

namespace RbtArgParser {

using ParsingError = cxxopts::exceptions::parsing;
typedef std::vector<std::string> ArgsSubstitutions;

class RbtArgParser {
 public:
    cxxopts::Options parser;
    ArgsSubstitutions substitutions;

 public:
    RbtArgParser(const std::string &program, const std::string &description = "");

    cxxopts::ParseResult parse(int argc, const char *argv[]);
    std::string help() const { return parser.help(); }

    // convenience functions to improve readability when adding options
    // add a boolean flag to the options
    RbtArgParser &add_flag(const std::string &opts, const std::string &dest, const char *default_value = "false") {
        return add<bool>(opts, dest, default_value);
    }

    // add whatever type of option
    template <typename T>
    RbtArgParser &add(const std::string &opts, const std::string &dest, const char *default_value = nullptr) {
        add_substitution(opts);
        if (default_value == nullptr) {
            parser.add_options()(opts, dest, cxxopts::value<T>());
        } else {
            parser.add_options()(opts, dest, cxxopts::value<T>()->default_value(default_value));
        }
        return *this;
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

class RbtOptionValue {
 public:
    cxxopts::OptionValue value;

    inline RbtOptionValue(cxxopts::OptionValue value): value(value) {}

    template <typename T>
    RbtOptionValue &operator>>(T &v) {
        v = value.as<T>();
        return *this;
    }

    template <typename T>
    RbtOptionValue &operator>>(std::optional<T> &v) {
        v = value.as<T>();
        return *this;
    }

    inline bool is_present() { return value.count() > 0; }
};

class RbtParseResult {
 public:
    cxxopts::ParseResult result;

    RbtParseResult(cxxopts::ParseResult result): result(result) {}

    inline RbtOptionValue operator[](const std::string &key) { return RbtOptionValue(result[key]); }
};

class ValidationError: public std::runtime_error {
 public:
    using std::runtime_error::runtime_error;
};

}  // namespace RbtArgParser
#endif  //_RBTARGPARSER_H_
