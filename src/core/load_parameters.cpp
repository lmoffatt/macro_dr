// src/core/load_parameters.cpp

#include <macrodr/cmd/load_parameters.h>
#include <macrodr/cmd/indexed_construction.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include "parameters.h"
#include "parameters_distribution.h"

namespace macrodr::cmd {

Maybe_error<var::Parameters_Transformations> load_parameters(const ModelPtr& model,
                                                             const std::string& parameters_file) {
    if (!model) {
        return error_message("model is null");
    }
    return var::load_Parameters(parameters_file, ",", model->model_name(), model->names());
}

//model_name,i_par,parameter_name,parameter_transformation,parameter_value
Maybe_error<var::Parameters_Transformations>
create_parameters(const ModelPtr& model,
    std::vector<std::tuple<std::string, std::string, double>>  parameters_info){
    if (!model) {
        return error_message("model is null");
        }
    return var::create_parameters(model->model_name(),model->names(), std::move(parameters_info));
}

Maybe_error<var::Parameters_Normal_Distribution>
create_prior(const ModelPtr& model,
    std::vector<std::tuple<std::string, std::string, double, double>> prior_info) {
    if (!model) {
        return error_message("model is null");
    }
    return var::create_prior(model->model_name(), model->names(), std::move(prior_info));
}

var::Parameters_values get_standard_parameter_values(var::Parameters_Transformations const& tr) {
    return var::Parameters_values(tr, tr.standard_values());
}

var::Parameters_transformed get_standard_parameter_transformed_values(var::Parameters_Transformations const& tr) {
    return tr.standard_parameter_transformed();
}

// Step 1 of the per-sample FD diagnostic perturbation, built OUTSIDE the
// diagnostic so the command stays agnostic. The `by_…` prefix marks a function
// that RETURNS an Indexed value, keyed by the axis it derives (here, the
// parameter coordinates) — counterpart to the `…_by` suffix family where you
// pass the axis in.
//
// Returns the raw, axis-aligned displacements Δ (NOT yet applied to θ): for each
// (parameter coordinate i, relative step h) a Parameters_transformed that is
// zero everywhere except slot i = h — i.e. h·eᵢ. The result is Indexed over two
// axes — `parameter_coordinate` (one coord per parameter, labeled by name) ×
// `h_rel` (one coord per value in h_rels, e.g. [+h,−h] or several scales). The
// `max(|θᵢ|,1)` scaling is NOT here — apply_relative_perturbation adds it. h_rels
// is a PLAIN vector and the function returns Indexed → it dispatches as a
// native-indexed (exact) function, no lifting. Agnostic: reads the parameter
// count and names from `par`.
//
// Flat-index order is "axis 0 fastest" (indexed.h flat_index), and the index
// space is {parameter_coordinate, h_rel}, so `parameter_coordinate` varies
// fastest → the values vector is filled with the parameter loop INNER.
Maybe_error<var::Indexed<var::Parameters_transformed>> by_parameter_coordinate(
    var::Parameters_transformed const& par, std::vector<double> h_rels) {
    const std::size_t p = par.size();
    const std::vector<std::string>& names = par.parameters().names();
    if (names.size() != p) {
        return error_message("by_parameter_coordinate: parameter-name count ", names.size(),
                             " differs from parameter count ", p);
    }
    if (h_rels.empty()) {
        return error_message("by_parameter_coordinate: h_rels is empty");
    }

    auto ax_param = axis("parameter_coordinate", names);
    if (!ax_param) {
        return ax_param.error();
    }

    std::vector<std::string> h_labels;
    h_labels.reserve(h_rels.size());
    for (double h : h_rels) {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%g", h);
        h_labels.emplace_back(buf);
    }
    auto ax_h = axis("h_rel", h_labels);
    if (!ax_h) {
        return ax_h.error();
    }

    var::IndexSpace space{{ax_param.value(), ax_h.value()}};  // axis 0 = parameter_coordinate (fast)
    auto valid_space = space.validate();
    if (!valid_space) {
        return valid_space.error();
    }

    std::vector<var::Parameters_transformed> deltas;
    deltas.reserve(p * h_rels.size());
    for (double h : h_rels) {              // outer = h_rel (slow axis)
        for (std::size_t i = 0; i < p; ++i) {  // inner = parameter_coordinate (fast axis)
            auto delta = par();            // copy the Matrix<double> to get the p×1 shape
            for (std::size_t k = 0; k < p; ++k) {
                delta[k] = 0.0;
            }
            delta[i] = h;                  // raw relative step along this parameter coordinate
            deltas.push_back(par.create(std::move(delta)));
        }
    }

    var::Indexed<var::Parameters_transformed> indexed(std::move(space), std::move(deltas));
    auto valid = indexed.validate();
    if (!valid) {
        return valid.error();
    }
    return indexed;
}

// Step 2: apply a relative perturbation. Given θ (`par`) and a raw displacement
// Δ (from by_parameter_coordinate), returns θ′ with θ′ᵢ = θᵢ + Δᵢ·max(|θᵢ|, 1)
// — exactly the perturbation calculate_mnumerical_fisher_information finite-
// differences (hᵢ = h_rel·max(|θᵢ|,1)). Returns a PLAIN Parameters_transformed,
// so the DSL lifts this over Δ's axes (parameter_coordinate × h_rel) — and over
// θ's axes if θ is Indexed — combining them. The `max(|θᵢ|,1)` scaling lives
// here (not in Δ) so it uses each θ's own values and stays identical to the
// numerical-Fisher logic.
Maybe_error<var::Parameters_transformed> apply_relative_perturbation(
    var::Parameters_transformed const& par, var::Parameters_transformed const& delta) {
    const std::size_t n = par.size();
    if (delta.size() != n) {
        return error_message("apply_relative_perturbation: perturbation size ", delta.size(),
                             " differs from parameter count ", n);
    }
    auto values = par();  // copy the underlying Matrix<double>
    for (std::size_t i = 0; i < n; ++i) {
        values[i] = par[i] + delta[i] * std::max(std::abs(par[i]), 1.0);
    }
    return par.create(std::move(values));
}


}  // namespace macrodr::cmd
