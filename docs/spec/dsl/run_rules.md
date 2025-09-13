
---

## MacroIR Runtime Execution Rules (C++-style, grounded)

**Status:** SPEC

**Derived from:** `typed_expression::run`, `typed_assignment::run_statement`, environment dispatch

**Last updated:** 2025-07-17

---
### 🔹 Program Execution


```cpp
Environment run(const typed_program& program) {
    Environment env;
    for (const auto& stmt : program.statements()) {
        auto result = stmt->run_statement(env);
        if (!result) throw RuntimeError(result.error());
    }
    return env;
}
```

### 🔹 Assignment Execution

```cpp
bool typed_assigment::run_statement(Environment& env) const {
    auto result = m_expr->run_expression(env);  // evaluates RHS
    if (!result) return result.error();
    env.insert(m_id, result.value());           // bind result to name
    return true;
}
```

### 🔹 Literal Evaluation

```cpp
template <typename T>
T typed_literal<T>::run(const Environment&) const {
    return m_value;  // constant result
}
```

### 🔹 Identifier Evaluation

```cpp
template <typename T>
T typed_identifier<T>::run(const Environment& env) const {
    auto bound = env.get(m_id);
    auto expr = dynamic_cast<typed_expression<T>*>(bound.value());
    return expr->run(env);  // re-run stored expression
}
```

### 🔹 Function Call Evaluation

```cpp
template <typename F, typename... Args>
auto typed_function_evaluation<F, Args...>::run(const Environment& env) const {
    return std::apply(
        [this, &env](auto&... arg_ptrs) {
            return m_f(arg_ptrs->run(env)...);
        }, m_args
    );
}
```

---

### ✅ Highlights From the Actual Code

* Each compiled statement is a `typed_statement` with a `run_statement()` method.
* The runtime `Environment` holds all bound values and intermediate results.
* Function calls evaluate arguments recursively and then invoke registered functions.
* All `typed_expression<T>` subclasses implement `T run(const Environment&) const`.
* Assignments are executed by evaluating RHS and binding the result in `env`.

---

           
