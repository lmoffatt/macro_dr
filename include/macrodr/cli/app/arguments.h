#pragma once

#include <string>
#include <vector>

namespace macrodr::cli::app {

using CliArguments = std::vector<std::string>;

CliArguments collect_arguments(int argc, char** argv);

}  // namespace macrodr::cli::app

