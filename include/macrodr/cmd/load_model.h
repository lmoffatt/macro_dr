#pragma once

#include <macrodr/interface/IModel.h>

#include <string>

#include "maybe_error.h"
#include "parameters.h"

namespace macrodr::cmd {
#if 0
using model_handle =
    Maybe_unique<typename macrodr::interface::IModel<typename var::Parameters_values>>;

model_handle load_model(const std::string& model_name);
#endif
}  // namespace macrodr::cmd
