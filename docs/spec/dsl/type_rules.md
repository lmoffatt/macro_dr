
---

# MacroIR Type Deduction Rules (C++-style, grounded)

**Status:** SPEC

**Derived from:** `compile_expression`, `compile_statement`, `to_typed_function`

**Last updated:** 2025-07-17

---

```cpp
typed_expression_ptr compile_expression(const untyped_expression& expr,
                                        Environment<Lexer,Compiler>& cm) {
    if (expr.is_numeric_literal()) {
        double value = std::stod(expr.str());
        return new typed_literal<Lexer,Compiler,double>(value);
    }

    if (expr.is_string_literal()) {
        return new typed_literal<Lexer,Compiler,std::string>(expr.str());
    }

    if (expr.is_identifier()) {
        auto resolved = cm.get_Identifier(expr.id());
        if (!resolved) return resolved.error();
        return resolved.value();  // already a typed_expression<...>
    }

    if (expr.is_function_call()) {
        auto fn = cm.get_function(expr.fid());
        if (!fn) return fn.error();
        return fn.value()->compile_function_evaluation(cm, expr.args());
    }

    if (expr.is_argument_list()) {
        // Returns typed_argument_list<...>
        return expr.compile_argument_list(cm);
    }

    if (expr.is_vector_literal()) {
        // Parsed as untyped_vector_construction<Lexer,Compiler>
        auto& node = expr.as_vector_construction();
        // The element type T is determined by the target function parameter or
        // assignment context (e.g. std::vector<double>, std::vector<Model>, ...)
        using T = /* context-dependent element type */;
        return vector_compiler<Lexer,Compiler,T>{}
            .compile_vector_construction(cm, node.args());
    }

    if (expr.is_tuple_literal()) {
        // Parsed as untyped_tuple_construction<Lexer,Compiler>
        auto& node = expr.as_tuple_construction();
        // The tuple type is std::tuple<Ts...>, inferred from the context
        using tuple_t = /* context-dependent tuple type */;
        using compiler_t = tuple_compiler<Lexer,Compiler,Ts...>;
        return compiler_t{}
            .compile_tuple_construction(cm, node);
    }

    return error_message("Unhandled expression type");
}
````
---

## Composite Type Rules (Vectors and Tuples)

### Storage Conversion

For any parameter or element type `U`, the internal storage type used during
compilation is:

```cpp
using storage_t = detail::function_argument_storage_t<U>;
````

with the following cases:

* `U` or `const U` → `storage_t = std::remove_cvref_t<U>`
* `U&` or `const U&` → `storage_t = std::reference_wrapper<U>` /
  `std::reference_wrapper<const U>`
* `const IModel<ParamValues...>&` → `storage_t = std::unique_ptr<IModel<ParamValues...>>`

`field_compiler` (for function parameters) and `element_compiler` (for
vector/tuple elements) always compile to `typed_expression<storage_t>` and
later adapt to `U` at runtime.

### Vector Type Rule

* If every element `ei` in `[e1, ..., en]` has type `T` in the target context,
  then the literal has type `std::vector<T>`.
* The compiled node is `typed_vector_construction<Lexer,Compiler,T>`, whose
  `run(env)` returns `std::vector<T>`.

### Tuple Type Rule

* If elements `{e1, ..., en}` have types `T1, ..., Tn` in the target context,
  then the literal has type `std::tuple<T1, ..., Tn>`.
* The compiled node is `typed_tuple_construction<Lexer,Compiler,T1,...,Tn>`,
  whose `run(env)` returns `std::tuple<T1,...,Tn>`.

### Reference Elements

When `U` is a reference type (`T&` or `const T&`):

* The storage type is `std::reference_wrapper<T>` / `std::reference_wrapper<const T>`.
* Only **identifier** expressions are accepted in that position.
* Literals, function calls, and nested composites are rejected by
  `element_compiler` / `field_compiler` with a clear diagnostic.

````

---

### ✅ Highlights From the Actual Code

* `compile_expression` dispatches on untyped node kind and delegates to either `typed_literal`, `typed_identifier`, or `typed_function_evaluation`.
* Function types are resolved from the compiler via `get_function(...)` and implemented using `to_typed_function<F, Args...>()`.
* Runtime evaluation calls `run(env)` inside `typed_expression<T>` or `run_expression(env)` to return wrapped result literals.
* Assignments store evaluated results in the environment using `env.insert(...)`.

---

Would you like the companion runtime rule block (`run_rules.md`) formatted next?
