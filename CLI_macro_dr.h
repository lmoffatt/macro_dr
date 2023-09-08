#ifndef CLI_MACRO_DR_H
#define CLI_MACRO_DR_H


#include "lexer_typed.h"

namespace dcli {

template <class... Args, class F, class... String>
  requires(sizeof...(Args) == sizeof...(String)) &&
          (std::is_constructible_v<std::string, String> && ...)
auto to_typed_function(F t_f, String &&...args) {
    
    return new function_compiler<Lexer, Compiler, F, Args...>(
      t_f, new typed_identifier<Lexer, Compiler, Args>(
                 to_Identifier<Lexer>(std::string(args)).value())...);
}



} // namespace dcli

#endif // CLI_MACRO_DR_H
