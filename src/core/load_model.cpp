#include <macrodr/cmd/load_model.h>

#include "allosteric_models.h"
#include "models_MoffattHume_allosteric.h"
#include "parameters.h"

namespace macrodr::cmd {
#if 0    
template <typename Scheme>
auto make_model(Scheme scheme) {
    return std::make_unique < macrodr::interface::ConcreteModel < Scheme,
           var::Parameters_values >>> (std::move(scheme));
}  // namespace macrodr::interface

model_handle load_model(const std::string& model_name) {
    auto Maybe_model = macrodr::get_model(model_name);
    if (!Maybe_model) {
        return Maybe_model.error();
    }
    auto model = Maybe_model.value();

    return std::visit(
        [](auto model0ptr) {
            using ModelType = std::decay_t<decltype(*model0ptr)>;
            return macrodr::interface::make_model_interface<ModelType, var::Parameters_values>(
                *model0ptr);
        },
        model);
}
#endif
}  // namespace macrodr::cmd
