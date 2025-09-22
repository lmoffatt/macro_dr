#pragma once

#include <macrodr/dsl/grammar_typed.h>

#include <memory>
#include <string>
#include <vector>

namespace macrodr::interface {

struct IObject {
    // [[nodiscard]] virtual std::string kind() const = 0;
    // virtual std::string describe() const = 0;

    virtual ~IObject() = default;
};
}  // namespace macrodr::interface