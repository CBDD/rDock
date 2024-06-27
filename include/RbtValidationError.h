#ifndef _RBT_VALIDATION_ERROR_H_
#define _RBT_VALIDATION_ERROR_H_

#include <exception>
#include <string>

namespace Rbt {

class ValidationError: public std::exception {
 public:
    explicit ValidationError(const std::string message): message(message) {}
    const char* what() const noexcept override { return message.c_str(); }

 protected:
    std::string message;
};

}  // namespace Rbt

#endif  // _RBT_VALIDATION_ERROR_H_