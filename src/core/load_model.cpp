// src/core/load_model.cpp

#include <macrodr/cmd/load_model.h>
#include <macrodr/interface/IModel.h>

#include <type_traits>
#include <utility>
#include <variant>

#include "allosteric_models.h"
#include "models_MoffattHume_allosteric.h"
#include "models_MoffattHume_linear.h"
#include "parameters.h"

namespace macrodr::cmd {

model_handle load_model(const std::string& model_name) {
    auto Maybe_model = macrodr::get_model(model_name);
    if (!Maybe_model) {
        return Maybe_model.error();
    }

    // Move out the variant of owning pointers to concrete model types
    auto model = Maybe_model.value();

    // Wrap the concrete model into a polymorphic interface erased to IModel<var::Parameters_values>
    return std::visit(
        [](auto&& model0ptr) -> model_handle {
            using ModelType = std::decay_t<decltype(*model0ptr)>;
            auto iface =
                macrodr::interface::make_model_interface<ModelType, var::Parameters_values>(
                    *model0ptr);
            return iface;
        },
        model);
}

dmodel_handle load_dmodel(const std::string& model_name) {
    auto Maybe_model = macrodr::get_model(model_name);
    if (!Maybe_model) {
        return Maybe_model.error();
    }

    // Move out the variant of owning pointers to concrete model types
    auto model = Maybe_model.value();

    // Wrap the concrete model into a polymorphic interface erased to IModel<var::Parameters_values>
    return std::visit(
        [](auto&& model0ptr) -> dmodel_handle {
            using ModelType = std::decay_t<decltype(*model0ptr)>;
            auto iface = macrodr::interface::make_model_interface<
                ModelType, var::Parameters_values,
                var::Derivative<var::Parameters_values, var::Parameters_transformed>>(*model0ptr);
            return iface;
        },
        model);
}

}  // namespace macrodr::cmd
