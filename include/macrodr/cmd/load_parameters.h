#pragma once

#include <macrodr/interface/IModel.h>

#include <cstddef>
#include <string>
#include <utility>

#include "indexed.h"
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

// Two-step perturbation builder for the per-sample FD diagnostic (kept OUTSIDE
// the diagnostic command so it stays agnostic).
//
// Step 1 — by_parameter_coordinate: the raw, axis-aligned displacements Δ (NOT
// yet applied to θ), one Parameters_transformed per (parameter coordinate,
// relative step) that is zero except slot i = h (i.e. h·eᵢ). Indexed over
// `parameter_coordinate` × `h_rel`. The `by_…` prefix marks a function that
// RETURNS an Indexed keyed by the axis it DERIVES (vs the `…_by` suffix family
// where you pass the axis). Plain args + Indexed return → native-indexed
// dispatch (no lifting).
Maybe_error<var::Indexed<var::Parameters_transformed>> by_parameter_coordinate(
    var::Parameters_transformed const& par, std::vector<double> h_rels);

// Step 2 — apply_relative_perturbation: θ′ᵢ = θᵢ + Δᵢ·max(|θᵢ|,1) (identical to
// calculate_mnumerical_fisher_information). Returns a PLAIN Parameters_
// transformed, so the DSL lifts it over Δ's axes (and θ's, if Indexed),
// combining them.
Maybe_error<var::Parameters_transformed> apply_relative_perturbation(
    var::Parameters_transformed const& par, var::Parameters_transformed const& delta);

}  // namespace macrodr::cmd
