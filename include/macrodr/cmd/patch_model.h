#pragma once

#include <memory>
#include <string>
#include <utility>

#include <macrodr/interface/IModel.h>

#include "maybe_error.h"
#include "parameters.h"
#include "qmodel.h"

namespace macrodr::cmd {

// Convenient aliases for readability
using ModelPtr = std::unique_ptr<macrodr::interface::IModel<var::Parameters_values>>;
using PatchModel = macrodr::Transfer_Op_to<var::Parameters_values, macrodr::Patch_Model>;
using QxEig = macrodr::Transfer_Op_to<PatchModel, macrodr::Qx_eig>;
using PMean = macrodr::Transfer_Op_to<PatchModel, macrodr::P_mean>;

// Build a Patch_Model by applying parameter values to a loaded model.
// `parameters_file_and_sep` is the pair returned by load_Parameter(filename, sep)
Maybe_error<PatchModel> patch_model(const ModelPtr& model,
                                   const std::pair<std::string, std::string>& parameters_file_and_sep);

// Expose a couple of low-level qmodel helpers
Maybe_error<QxEig> calc_eigen(const PatchModel& pm, double atp_concentration);

PMean calc_peq(const QxEig& qx_eig, const PatchModel& pm);

}  // namespace macrodr::cmd

