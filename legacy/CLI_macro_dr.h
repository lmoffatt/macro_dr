#pragma once
#include <macrodr/dsl/grammar_Identifier.h>
#include <macrodr/dsl/grammar_untyped.h>
#include <macrodr/dsl/lexer_typed.h>
#include <macrodr/dsl/lexer_untyped.h>

#include <memory>
#include <string>
#include <type_traits>

namespace macrodr::dsl {

template <class... Args, class F, class... String>
    requires((sizeof...(Args) == sizeof...(String)) &&
             (std::is_constructible_v<std::string, String> && ...) &&
             (std::is_void_v<std::invoke_result_t<F, Args...>> ||
              std::is_object_v<std::invoke_result_t<F, Args...>>))
auto to_typed_function(F t_f, String&&... args) {
    return std::make_unique<function_compiler<Lexer, Compiler, F, Args...>>(
        t_f, field_compiler<Lexer, Compiler, detail::function_argument_storage_t<Args>>(
                 to_Identifier<Lexer>(std::forward<String>(args)).value())...);
}

template <class T, class NameString, class LabelsString, class ValuesString>
    requires(std::is_constructible_v<std::string, NameString> &&
             std::is_constructible_v<std::string, LabelsString> &&
             std::is_constructible_v<std::string, ValuesString>)
auto to_typed_indexed_constructor(NameString&& name_arg, LabelsString&& labels_arg,
                                  ValuesString&& values_arg) {
    return std::make_unique<indexed_constructor_compiler<Lexer, Compiler, T>>(
        to_Identifier<Lexer>(std::forward<NameString>(name_arg)).value(),
        to_Identifier<Lexer>(std::forward<LabelsString>(labels_arg)).value(),
        to_Identifier<Lexer>(std::forward<ValuesString>(values_arg)).value());
}

template <class T, class P>
    requires(std::is_same_v<Maybe_error<T>, std::invoke_result_t<P, T>>)
auto to_typed_predicate(P t_f, std::string arg) {
    return std::make_unique<predicate_compiler<Lexer, Compiler, P, T>>(
        t_f, field_compiler<Lexer, Compiler, T>(to_Identifier<Lexer>(std::move(arg)).value()));
}

}  // namespace macrodr::dsl
