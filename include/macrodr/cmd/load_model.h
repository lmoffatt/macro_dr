#pragma once

#include <string>

#include "maybe_error.h"
#include "parameters.h"

// Forward declarations to avoid include cycles with IModel.h -> qmodel.h -> CLI headers.
namespace macrodr::interface {
template <typename... ParamValues>
struct IModel;
}  // namespace macrodr::interface

namespace var {
template <class, class>
class Derivative;  // forward decl only
class Parameters_transformed;  // forward decl only
}  // namespace var

namespace macrodr::cmd {

using model_handle = Maybe_unique<macrodr::interface::IModel<var::Parameters_values>>;

using dmodel_handle = Maybe_unique<macrodr::interface::IModel<
    var::Parameters_values, var::Derivative<var::Parameters_values, var::Parameters_transformed>>>;

model_handle load_model(const std::string& model_name);

dmodel_handle load_dmodel(const std::string& model_name);

}  // namespace macrodr::cmd
