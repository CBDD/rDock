#ifndef _RBT_RBDOCK_ARGPARSER_H_
#define _RBT_RBDOCK_ARGPARSER_H_

#include "rbdock/rbdock_config.h"

#include "RbtArgParser.h"

namespace RBDock {

RbtArgParser::RbtArgParser get_options_parser();

RBDock::RBDockConfig parse_args(int argc, const char *argv[]);

}  // namespace RBDock

#endif  // _RBT_RBDOCK_ARGPARSER_H_