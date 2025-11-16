
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

### 🔹 Vector Literal Evaluation

```cpp
template <class Lexer, class Compiler, class T>
Maybe_error<std::vector<T>>
typed_vector_construction<Lexer,Compiler,T>::run(const Environment<Lexer,Compiler>& env) const {
    using storage_t = detail::function_argument_storage_t<T>;

    std::vector<T> out;
    out.reserve(m_args.size());
    std::string err;

    for (std::size_t i = 0; i < m_args.size(); ++i) {
        // Each arg: typed_expression<storage_t>
        auto maybe_elem = m_args[i]->run(env);  // Maybe_error<storage_t>
        if (!maybe_elem) {
            err += std::to_string(i) + ": " + maybe_elem.error()();
            continue;
        }
        // Same adaptation logic as function arguments
        out.emplace_back(adapt<T>(maybe_elem.value()));
    }

    if (!err.empty()) {
        return error_message(err);
    }
    return out;
}
````

### 🔹 Tuple Literal Evaluation

```cpp
template <class Lexer, class Compiler, class... Ts>
Maybe_error<std::tuple<Ts...>>
typed_tuple_construction<Lexer,Compiler,Ts...>::run(const Environment<Lexer,Compiler>& env) const {
    using storage_t = detail::function_argument_storage_t;

    // Evaluate each element → Maybe_error<storage_t<Tk>>
    auto maybe_storage = std::apply(
        [&](auto&... exprs) {
            return std::tuple(exprs->run(env)...);
        }, m_args);

    bool ok = true;
    std::string msg;

    std::apply([&](auto&... me) {
        (([&]{
            if (!me.valid()) {
                ok = false;
                msg += me.error()();
            }
        }()), ...);
    }, maybe_storage);

    if (!ok) {
        return error_message(msg);
    }

    // Build the exposed tuple<Ts...>, adapting each storage element
    auto builder = [&](auto&... me) {
        return std::tuple<Ts...>( adapt<Ts>(me.value())... );
    };

    return std::apply(builder, maybe_storage);
}
```

### 🔹 Reference Semantics in Composites

* Elements compiled as `std::reference_wrapper<T>` or `std::reference_wrapper<const T>`
  behave like function arguments:

  * At runtime they **read** from or **write** to values stored in the `Environment`.
  * For vectors/tuples of references, individual elements remain aliases to
    the underlying environment literals.
* Elements compiled as owning types (`T`, `std::unique_ptr<T>`, etc.) behave
  as value copies or ownership transfers, depending on the underlying storage
  type and `adapt<>()` rule.

```

---


---

### ✅ Highlights From the Actual Code

* Each compiled statement is a `typed_statement` with a `run_statement()` method.
* The runtime `Environment` holds all bound values and intermediate results.
* Function calls evaluate arguments recursively and then invoke registered functions.
* All `typed_expression<T>` subclasses implement `T run(const Environment&) const`.
* Assignments are executed by evaluating RHS and binding the result in `env`.

---

           
