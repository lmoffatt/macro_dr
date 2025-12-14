// src/core/load_parameters.cpp

#include <macrodr/cmd/load_parameters.h>

#include <string>
#include <utility>

#include "parameters.h"

namespace macrodr::cmd {

Maybe_error<var::Parameters_Transformations> load_parameters(const ModelPtr& model,
                                                             const std::string& parameters_file) {
    if (!model) {
        return error_message("model is null");
    }
    return var::load_Parameters(parameters_file, ",", model->model_name(), model->names());
}

var::Parameters_values get_standard_parameter_values(var::Parameters_Transformations const& tr) {
    return var::Parameters_values(tr, tr.standard_values());
}

var::Parameters_transformed get_standard_parameter_transformed_values(var::Parameters_Transformations const& tr) {
    return tr.standard_parameter_transformed();
}


}  // namespace macrodr::cmd
