
#ifndef _RBT_BINARY_IO_H_
#define _RBT_BINARY_IO_H_

#include <iostream>

#include "RbtError.h"
#include "RbtFileError.h"

namespace Rbt {

template <typename T>
inline void bin_write(std::ostream& ostr, const T& data) {
    try {
        ostr.write(reinterpret_cast<const char*>(&data), sizeof(T));
    } catch (std::ios_base::failure& e) {
        throw RbtFileWriteError(_WHERE_, e.what());
    }
}

template <typename T>
inline void bin_write(std::ostream& ostr, const T* data, size_t n) {
    try {
        ostr.write(reinterpret_cast<const char*>(data), n * sizeof(T));
    } catch (std::ios_base::failure& e) {
        throw RbtFileWriteError(_WHERE_, e.what());
    }
}

template <typename T>
inline void bin_read(std::istream& istr, T& data) {
    try {
        istr.read(reinterpret_cast<char*>(&data), sizeof(T));
    } catch (std::ios_base::failure& e) {
        throw RbtFileReadError(_WHERE_, e.what());
    }
}

template <typename T>
inline void bin_read(std::istream& istr, T* data, size_t n) {
    try {
        istr.read(reinterpret_cast<char*>(data), n * sizeof(T));

    } catch (std::ios_base::failure& e) {
        throw RbtFileReadError(_WHERE_, e.what());
    }
}

}  // namespace Rbt

#endif  // _RBT_BINARY_IO_H_
