#pragma once

#include <macrodr/interface/IModel.h>

#include <string>

#include "maybe_error.h"
#include "parameters.h"

namespace macrodr::cmd {

using model_handle =
    Maybe_unique<typename macrodr::interface::IModel<typename var::Parameters_values>>;

using dmodel_handle = Maybe_unique<typename macrodr::interface::IModel<
    typename var::Parameters_values,
    var::Derivative<var::Parameters_values, var::Parameters_transformed>>>;

model_handle load_model(const std::string& model_name);

dmodel_handle load_dmodel(const std::string& model_name);

}  // namespace macrodr::cmd
