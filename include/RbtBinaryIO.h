
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

//speciallization for std::string

template <>
inline void bin_write(std::ostream& ostr, const std::string& data) {
    bin_write(ostr, data.size());
    bin_write(ostr, data.data(), data.size());
}

template <>
inline void bin_read(std::istream& istr, std::string& data) {
    size_t size;
    bin_read(istr, size);
    data.resize(size);
    bin_read(istr, &data[0], size);
}

}  // namespace Rbt

#endif  // _RBT_BINARY_IO_H_
