#pragma once

#include <memory>
#include <string>
#include <vector>

#include "grammar_typed.h"

namespace macrodr {

namespace interface {

struct IObject {
    virtual std::string kind() const = 0;
    virtual std::string describe() const = 0;

    virtual ~IObject() = default;
};
}  // namespace interface
}  // namespace macrodr
