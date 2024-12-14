
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

template <>
inline void bin_write(std::ostream& ostr, const std::string& data) {
    // casting size_t to int, and storing the null character at the end
    // for retro-compatibility with old code only. This is not the best way to do it
    // and will be modified in the near future
    bin_write(ostr, (int)(data.size() + 1));
    bin_write(ostr, data.data(), data.size() + 1);
}

template <>
inline void bin_read(std::istream& istr, std::string& data) {
    // similar to the bin_write(std::ostream& ostr, const std::string& data) function
    // this is a temporary solution for retro-compatibility with old code
    // the size_t is casted to int, and the null character added at the end must be removed
    int size;
    bin_read(istr, size);
    data.resize(size);
    bin_read(istr, data.data(), size);
    data.pop_back();
}

}  // namespace Rbt

#endif  // _RBT_BINARY_IO_H_
