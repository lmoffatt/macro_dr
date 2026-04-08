#include <maybe_error.h>
#include "json_spec.h"

namespace macrodr::dsl {
template <class Lexer>
class Identifier;

template <class Lexer, class Compiler> class Environment;
template <class Lexer, class Compiler> class base_typed_statement;
template <class Lexer, class Compiler> class base_typed_expression;
template <class Lexer, class Compiler> class base_typed_argument;
template <class Lexer, class Compiler> class base_Identifier_compiler;
template <class Lexer, class Compiler> class base_function_compiler;

template <class Lexer, class Compiler, class T> class typed_expression;
template <class Lexer, class Compiler, class T> class typed_argument;
template <class Lexer, class Compiler, class T> class typed_literal;
template <class Lexer, class Compiler, class T> class Identifier_compiler;
template <class Lexer, class Compiler, class T> class parameter_compiler;
template <class Lexer, class Compiler>
class Environment;

template <class Lexer, class Compiler>
class base_typed_argument;


template <class Lexer, class Compiler, class T>
class expression_argument;

template <class Lexer, class Compiler, class T>
class pointwise_argument;

template <class Lexer, class Compiler, class T>
class borrowed_expression_argument;

template <class Lexer, class Compiler, class T>
class pointwise_borrowed_argument;

template <class Lexer, class Compiler, class T>
class typed_vector_construction;

template <class Lexer, class Compiler, class T>
class typed_lifted_vector_construction;

template <class Lexer, class Compiler, class T>
class typed_set_construction;

template <class Lexer, class Compiler, class T>
class typed_lifted_set_construction;


template <class Lexer, class Compiler, class... T>
class typed_tuple_construction;

template <class Lexer, class Compiler, class... T>
class typed_lifted_tuple_construction;


template <class Lexer, class Compiler, class T>
class typed_identifier;

template <class Lexer, class Compiler, class T>
class typed_identifier_ref_const;

template <class Lexer, class Compiler, class T>
class typed_identifier_ref;



template <class Lexer, class Compiler>
class untyped_argument_list;

template <class Lexer, class Compiler, class T, class S>
    requires(std::convertible_to<S, T>)
class typed_conversion;

template <class Lexer, class Compiler, class F, class... Args>
    requires(std::is_object_v<std::invoke_result_t<F, Args...>> ||
             std::is_void_v<std::invoke_result_t<F, Args...>>)
class typed_scalar_function_evaluation;

template <class Lexer, class Compiler, class F, class... Args>
    requires(std::is_object_v<std::invoke_result_t<F, Args...>> ||
             std::is_void_v<std::invoke_result_t<F, Args...>>)
class typed_lifted_function_evaluation;

template <class Lexer, class Compiler, class F, class... Args>
    requires(std::is_object_v<std::invoke_result_t<F, Args...>> ||
             std::is_void_v<std::invoke_result_t<F, Args...>>)
class typed_exact_indexed_function_evaluation;

template <class Lexer, class Compiler>
class untyped_program;

template <class Lexer, class Compiler>
class typed_program;

template <class Lexer, class Compiler>
class untyped_expression;

template <class Lexer, class Compiler>
class untyped_statement;

template <class Lexer, class Compiler>
class untyped_numeric_literal;

template <class Lexer, class Compiler>
class untyped_string_literal;

template <class Lexer, class Compiler>
class untyped_identifier;

template <class Lexer, class Compiler>
class untyped_assignment;

template <class Lexer, class Compiler>
class untyped_literal;

template <class Lexer, class Compiler>
class untyped_function_evaluation;

template <class Lexer, class Compiler>
class untyped_vector_construction;

template <class Lexer, class Compiler>
class untyped_tuple_construction;

template <class Lexer, class Compiler, class P, class T>
    requires(std::is_same_v<Maybe_error<T>, std::invoke_result_t<P, T>>)
class typed_predicate_evaluation;


template <class Lexer, class Compiler,  class... Ts>
 class tuple_compiler;

 template <class Lexer, class Compiler, class T>
Maybe_error<void> load_literal_from_json_helper(const typename json_spec<Lexer>::Json& value,
                                                const std::string& path,
                                                typename json_spec<Lexer>::TagPolicy policy,
                                                const Identifier<Lexer>& id,
                                                Environment<Lexer, Compiler>& env);


}