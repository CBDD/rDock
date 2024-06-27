#ifndef _RBT_TEST_UTILS_H_
#define _RBT_TEST_UTILS_H_

#include <iostream>
#include <streambuf>
#include <sstream>

class CoutRedirect {
public:
    inline CoutRedirect(std::streambuf* new_buffer): old(std::cout.rdbuf(new_buffer)) {}
    inline ~CoutRedirect() {std::cout.rdbuf(old);}

private:
    std::streambuf* old;
};


class CerrRedirect {
public:
    inline CerrRedirect(std::streambuf* new_buffer): old(std::cerr.rdbuf(new_buffer)) {}
    inline ~CerrRedirect() {std::cerr.rdbuf(old);}

private:
    std::streambuf* old;
};

class NullBuffer : public std::streambuf {
public:
    inline int overflow(int c) override { return c; }
};

#endif  // _RBT_TEST_UTILS_H_
