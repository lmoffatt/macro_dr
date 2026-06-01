#pragma once

#include <macrodr/interface/IModel.h>

#include <cstddef>
#include <string>
#include <utility>

#include "maybe_error.h"
#include "parameters.h"
#include "parameters_distribution.h"
#include "patch_model.h"  // reuse ModelPtr alias

namespace macrodr::cmd {

// Load parameter transformations from a file (filename, separator) using the
// provided model to validate names/order. Returns the full Parameters_Transformations
// object so callers can derive standard values or perform further ops.
Maybe_error<var::Parameters_Transformations> load_parameters(const ModelPtr& model,
                                                             const std::string& parameters_file);


//model_name,i_par,parameter_name,parameter_transformation,parameter_value
Maybe_error<var::Parameters_Transformations>
create_parameters(const ModelPtr& model,
    std::vector<std::tuple<std::string, std::string, double>>  parameters_info);


// In-line construction of a Parameters_Normal_Distribution prior. Counterpart
// to load_prior (CSV-based); use this from .macroir to build a diagonal
// Gaussian prior without an external file. Each tuple is
// (parameter_name, transformation, value, transformed_variance) — the prior
// mean in transformed space is derived as tr(value).
Maybe_error<var::Parameters_Normal_Distribution>
create_prior(const ModelPtr& model,
    std::vector<std::tuple<std::string, std::string, double, double>> prior_info);



// Convenience: legacy alias returning transformations rather than raw values (safer across steps)
var::Parameters_values get_standard_parameter_values(var::Parameters_Transformations const& tr);

var::Parameters_transformed get_standard_parameter_transformed_values(var::Parameters_Transformations const& tr);

}  // namespace macrodr::cmd
