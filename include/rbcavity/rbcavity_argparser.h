#ifndef _RBT_RBCAVITY_ARGPARSER_H_
#define _RBT_RBCAVITY_ARGPARSER_H_

#include "rbcavity/rbcavity_config.h"

#include "RbtArgParser.h"

namespace RBCavity {

RbtArgParser::RbtArgParser get_options_parser();

RBCavity::RBCavityConfig parse_args(int argc, const char *argv[]);

}  // namespace RBCavity

#endif  // _RBT_RBCAVITY_ARGPARSER_H_