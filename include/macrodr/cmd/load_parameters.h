#pragma once

#include <macrodr/interface/IModel.h>

#include <string>
#include <utility>

#include "maybe_error.h"
#include "parameters.h"
#include "patch_model.h"  // reuse ModelPtr alias

namespace macrodr::cmd {

// Load parameter transformations from a file (filename, separator) using the
// provided model to validate names/order. Returns the full Parameters_Transformations
// object so callers can derive standard values or perform further ops.
Maybe_error<var::Parameters_Transformations> load_parameters(const ModelPtr& model,
                                                             const std::string& parameters_file);

// Convenience: legacy alias returning transformations rather than raw values (safer across steps)
var::Parameters_values get_standard_parameter_values(var::Parameters_Transformations const& tr);

}  // namespace macrodr::cmd
