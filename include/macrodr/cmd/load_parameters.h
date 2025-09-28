#pragma once

#include <memory>
#include <string>
#include <utility>

#include <macrodr/interface/IModel.h>

#include "maybe_error.h"
#include "parameters.h"
#include "qmodel.h"
#include "patch_model.h"  // reuse ModelPtr alias

namespace macrodr::cmd {

// Load parameter transformations from a file (filename, separator) using the
// provided model to validate names/order. Returns the full Parameters_Transformations
// object so callers can derive standard values or perform further ops.
Maybe_error<var::Parameters_Transformations> load_parameters(
    const ModelPtr& model, const std::pair<std::string, std::string>& parameters_file_and_sep);

// Convenience: legacy alias returning transformations rather than raw values (safer across steps)
Maybe_error<var::Parameters_Transformations> load_parameter_values(
    const ModelPtr& model, const std::pair<std::string, std::string>& parameters_file_and_sep);

}  // namespace macrodr::cmd
