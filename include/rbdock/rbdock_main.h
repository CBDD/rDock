#ifndef _RBT_RBDOCK_MAIN_H_
#define _RBT_RBDOCK_MAIN_H_

#include <string>

#include "rbdock/rbdock_config.h"

namespace RBDock {
void RBDock(const RBDockConfig &config, const std::string &strExeName);
}  // namespace RBDock

#endif  // _RBT_RBDOCK_MAIN_H_