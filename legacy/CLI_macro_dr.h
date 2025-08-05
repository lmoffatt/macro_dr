#ifndef CLI_MACRO_DR_H
#define CLI_MACRO_DR_H

#include <string>
#include <type_traits>

//#include "lexer_typed.h"
#include <macrodr/dsl/lexer_untyped.h>
namespace macrodr {

namespace dsl {

template <class... Args, class F, class... String>
    requires((sizeof...(Args) == sizeof...(String)) &&
             (std::is_constructible_v<std::string, String> && ...) &&
             (std::is_void_v<std::invoke_result_t<F, Args...>> ||
              std::is_object_v<std::invoke_result_t<F, Args...>>))
auto to_typed_function(F t_f, String&&... args) {
    return new function_compiler<Lexer, Compiler, F, Args...>(
        t_f,
        field_compiler<Lexer, Compiler, Args>(to_Identifier<Lexer>(std::string(args)).value())...);
}


template <class T, class P>
    requires( std::is_same_v<Maybe_error<T>, std::invoke_result_t<P,T>>)
auto to_typed_predicate(P t_f, std::string arg) {
    return new predicate_compiler<Lexer, Compiler, P, T>(
        t_f,
        field_compiler<Lexer, Compiler,T>(to_Identifier<Lexer>(arg).value()));
}




}  // namespace dsl
}  // namespace macrodr
#endif  // CLI_MACRO_DR_H
