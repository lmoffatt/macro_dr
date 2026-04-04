#pragma once

#include <indexed.h>
#include <maybe_error.h>

#include <string>
#include <vector>

namespace macrodr::cmd {

Maybe_error<var::Axis> axis(std::string name, std::vector<std::string> labels);

template <class T>
Maybe_error<var::Indexed<T>> indexed_by(var::Axis axis, std::vector<T> values) {
    auto valid_axis = axis.validate();
    if (!valid_axis) {
        return valid_axis.error();
    }

    var::IndexSpace space{{axis}};
    auto valid_space = space.validate();
    if (!valid_space) {
        return valid_space.error();
    }

    var::Indexed<T> indexed(std::move(space), std::move(values));
    auto valid_indexed = indexed.validate();
    if (!valid_indexed) {
        return valid_indexed.error();
    }
    return indexed;
}

}  // namespace macrodr::cmd
