// src/core/load_parameters.cpp

#include <macrodr/cmd/load_parameters.h>

#include <string>
#include <utility>

#include "parameters.h"

namespace macrodr::cmd {

Maybe_error<var::Parameters_Transformations> load_parameters(
    const ModelPtr& model, const std::pair<std::string, std::string>& parameters_file_and_sep) {
    if (!model) {
        return error_message("model is null");
    }
    return var::load_Parameters(parameters_file_and_sep.first, parameters_file_and_sep.second,
                                model->model_name(), model->names());
}

Maybe_error<var::Parameters_Transformations> load_parameter_values(
    const ModelPtr& model, const std::pair<std::string, std::string>& parameters_file_and_sep) {
    return load_parameters(model, parameters_file_and_sep);
}

}  // namespace macrodr::cmd
