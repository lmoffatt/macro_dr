
---

# MacroIR Type Deduction Rules (C++-style, grounded)

**Status:** SPEC

**Derived from:** `compile_expression`, `compile_statement`, `to_typed_function`

**Last updated:** 2025-07-17

---

```cpp
typed_expression_ptr compile_expression(const untyped_expression& expr, Environment<Lexer,Compiler>& cm) {
    if (expr.is_numeric_literal()) {
        double value = std::stod(expr.str());
        return new typed_literal<double>(value);  // type: double
    }

    if (expr.is_string_literal()) {
        return new typed_literal<std::string>(expr.str());  // type: string
    }

    if (expr.is_identifier()) {
        auto resolved = cm.get_Identifier(expr.id());
        if (!resolved) return resolved.error();
        return resolved.value();
    }

    if (expr.is_function_call()) {
        auto fn = cm.get_function(expr.fid());
        if (!fn) return fn.error();
        return fn.value()->compile_function_evaluation(cm, expr.args());
    }

    if (expr.is_argument_list()) {
        return expr.compile_argument_list(cm);  // typed_argument_list<...>
    }

    return error_message("Unhandled expression type");
}
```

---

### ✅ Highlights From the Actual Code

* `compile_expression` dispatches on untyped node kind and delegates to either `typed_literal`, `typed_identifier`, or `typed_function_evaluation`.
* Function types are resolved from the compiler via `get_function(...)` and implemented using `to_typed_function<F, Args...>()`.
* Runtime evaluation calls `run(env)` inside `typed_expression<T>` or `run_expression(env)` to return wrapped result literals.
* Assignments store evaluated results in the environment using `env.insert(...)`.

---

Would you like the companion runtime rule block (`run_rules.md`) formatted next?
