#ifndef COMMAND_MANAGER_H
#define COMMAND_MANAGER_H

#include <macrodr/dsl/lexer_typed.h>
#include <macrodr/dsl/Lexer.h>
#include <macrodr/dsl/Compiler.h>

namespace macrodr::cli {

dsl::Compiler<dsl::Lexer> make_compiler_new();

}  // namespace macrodr::cli

#endif  // COMMAND_MANAGER_H
