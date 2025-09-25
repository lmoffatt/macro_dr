// src/core/patch_model.cpp

#include <macrodr/cmd/patch_model.h>

#include <string>
#include <utility>

#include "parameters.h"
#include "qmodel.h"

namespace macrodr::cmd {

Maybe_error<PatchModel> patch_model(
    const ModelPtr& model, const std::pair<std::string, std::string>& parameters_file_and_sep) {
    if (!model) {
        return error_message("model is null");
    }
    // Load parameter values using the legacy helper
    auto maybe_params =
        var::load_Parameters(parameters_file_and_sep.first, parameters_file_and_sep.second,
                             model->model_name(), model->names());
    if (!maybe_params) {
        return maybe_params.error();
    }

    // Use the model’s natural operator() to produce a Patch_Model
    auto pvalues = maybe_params.value().standard_parameter();  // var::Parameters_values
    auto maybe_pm = (*model)(pvalues);
    if (!maybe_pm) {
        return maybe_pm.error();
    }
    return maybe_pm.value();
}

Maybe_error<QxEig> calc_eigen(const PatchModel& pm, double atp_concentration) {
    auto m = Macro_DMR{};
    return m.calc_eigen(pm, ::macrodr::ATP_concentration(atp_concentration));
}

PMean calc_peq(const QxEig& qx_eig, const PatchModel& pm) {
    auto m = Macro_DMR{};
    return m.calc_Peq_(qx_eig, pm);
}

}  // namespace macrodr::cmd
