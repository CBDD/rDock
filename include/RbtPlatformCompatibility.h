#ifndef _RBT_PLATFORM_COMPATIBILITY_H_
#define _RBT_PLATFORM_COMPATIBILITY_H_
#include <ios>
namespace Rbt {

#if defined(__sgi) && !defined(__GNUC__)
const auto inputMode = std::ios_base::in;
#else
const auto inputMode = std::ios_base::in | std::ios_base::binary;
#endif

#if defined(__sgi) && !defined(__GNUC__)
const auto outputMode = std::ios_base::out | std::ios_base::trunc;
#else
const auto outputMode = ios_base::out | ios_base::binary | ios_base::trunc;
#endif

}  // namespace Rbt
#endif  // _RBT_PLATFORM_COMPATIBILITY_H_