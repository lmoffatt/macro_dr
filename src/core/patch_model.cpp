// src/core/patch_model.cpp

#include <experiment.h>
#include <macrodr/cmd/patch_model.h>

#include <string>
#include <utility>

#include "parameters.h"
#include "qmodel.h"

namespace macrodr::cmd {

static Maybe_error<void> validate_patch_model_dims(const PatchModel& m) {
    // Basic dimensional consistency checks to avoid undefined behavior downstream
    const auto nst = get<N_St>(m)();
    const auto& Q0m = get<Q0>(m)();
    const auto& Qam = get<Qa>(m)();
    const auto& gm = get<g>(m)();
    const auto& Pini = get<P_initial>(m)();

    auto fail = [&](const std::string& msg) { return error_message(std::string{"patch_model invalid: "} + msg); };

    if (nst == 0 || nst > 1000000) {
        return fail("N_St out of bounds: " + std::to_string(nst));
    }
    if (Q0m.nrows() != nst || Q0m.ncols() != nst) {
        return fail("Q0 dims " + std::to_string(Q0m.nrows()) + "x" + std::to_string(Q0m.ncols()) +
                    " != N_St=" + std::to_string(nst));
    }
    if (Qam.nrows() != nst || Qam.ncols() != nst) {
        return fail("Qa dims " + std::to_string(Qam.nrows()) + "x" + std::to_string(Qam.ncols()) +
                    " != N_St=" + std::to_string(nst));
    }
    if (!((gm.nrows() == nst && gm.ncols() == 1) || (gm.nrows() == 1 && gm.ncols() == nst))) {
        return fail("g dims " + std::to_string(gm.nrows()) + "x" + std::to_string(gm.ncols()) +
                    " incompatible with N_St=" + std::to_string(nst));
    }
    if (!((Pini.nrows() == 1 && Pini.ncols() == nst) || (Pini.nrows() == nst && Pini.ncols() == 1))) {
        return fail("P_initial dims " + std::to_string(Pini.nrows()) + "x" + std::to_string(Pini.ncols()) +
                    " incompatible with N_St=" + std::to_string(nst));
    }
    return {};
}

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

Maybe_error<PatchModel> patch_model(const ModelPtr& model, const var::Parameters_values& values) {
    if (!model) {
        return error_message("model is null");
    }
    auto maybe_pm = (*model)(values);
    if (!maybe_pm) {
        return maybe_pm.error();
    }
    // Validate dimensions early so failures show up at step_2 instead of crashing later
    if (auto chk = validate_patch_model_dims(maybe_pm.value()); !chk) {
        return chk.error();
    }
    return maybe_pm.value();
}

Maybe_error<PatchModel> patch_model(const ModelPtr& model, const var::Parameters_Transformations& tr) {
    if (!model) {
        return error_message("model is null");
    }
    auto pvalues = tr.standard_parameter();
    auto maybe_pm = (*model)(pvalues);
    if (!maybe_pm) {
        return maybe_pm.error();
    }
    if (auto chk = validate_patch_model_dims(maybe_pm.value()); !chk) {
        return chk.error();
    }
    return maybe_pm.value();
}

Maybe_error<macrodr::Patch_State> path_state(const PatchModel& pm, double initial_atp) {
    if (auto chk = validate_patch_model_dims(pm); !chk) {
        return chk.error();
    }
    auto m = Macro_DMR{};
    // predictions<0> returns Patch_State (no evolution stored)
    auto maybe_init = m.init<return_predictions<0>>(pm,
                                                   initial_ATP_concentration(ATP_concentration(initial_atp)));
    if (!maybe_init) {
        return maybe_init.error();
    }
    if constexpr (std::is_same_v<decltype(maybe_init.value()), Patch_State>) {
        return maybe_init.value();
    } else {
        // Transfer_Op_to may wrap Patch_State for type propagation; normalize via construction
        return static_cast<Patch_State>(maybe_init.value());
    }
}

Maybe_error<macrodr::Patch_State> path_state(const ModelPtr& model,
                                             const var::Parameters_values& values,
                                             double initial_atp) {
    auto maybe_pm = patch_model(model, values);
    if (!maybe_pm) {
        return maybe_pm.error();
    }
    return path_state(maybe_pm.value(), initial_atp);
}

Maybe_error<QxEig> calc_eigen(const PatchModel& pm, double atp_concentration) {
    if (auto chk = validate_patch_model_dims(pm); !chk) {
        return chk.error();
    }
    auto m = Macro_DMR{};
    return m.calc_eigen(pm, ::macrodr::ATP_concentration(atp_concentration));
}

PMean calc_peq(const QxEig& qx_eig, const PatchModel& pm) {
    auto m = Macro_DMR{};
    return m.calc_Peq_(qx_eig, pm);
}

}  // namespace macrodr::cmd
