# `field_compiler` and Argument Adaptation

This note compares the current DSL behavior for:

- value-like storage: `T`
- reference-like storage: `std::reference_wrapper<T>`
- const-reference-like storage: `std::reference_wrapper<const T>`

It focuses on the function-argument path.

## Naming clarification

There is no standalone class named `argument_compiler` in the current code.

For function calls, the effective argument pipeline is:

1. `detail::function_argument_storage_t<Arg>`
2. `field_compiler<Lexer, Compiler, storage_t<Arg>>`
3. `typed_expression<..., storage_t<Arg>>`
4. `function_compiler::adapt_argument<Arg>(storage)`
5. actual C++ function call

So the "argument compiler" role is split between:

- compile time: `field_compiler`
- call-time adaptation: `function_compiler::adapt_argument(...)`

Code anchors:

- [lexer_typed.h](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/lexer_typed.h)
- [grammar_typed.h](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/grammar_typed.h)

## Storage rule first

The storage choice happens before `field_compiler` is instantiated:

```cpp
template <class T>
struct function_argument_storage<const T&> {
    using type = std::reference_wrapper<const T>;
};

template <class T>
struct function_argument_storage<T&> {
    using type = std::reference_wrapper<T>;
};
```

So:

- parameter `T` uses `field_compiler<T>`
- parameter `const T&` uses `field_compiler<std::reference_wrapper<const T>>`
- parameter `T&` uses `field_compiler<std::reference_wrapper<T>>`

This is why the comparison is really between different `field_compiler`
specializations.

## Side by side

| Parameter shape | Storage shape | Compiler used | What input is accepted | Runtime result |
| --- | --- | --- | --- | --- |
| `T` | `T` | `field_compiler<T>` | identifiers, decodable literals, function calls, and composites when `T` is composite | a value of type `T` |
| `const T&` | `std::reference_wrapper<const T>` | `field_compiler<std::reference_wrapper<const T>>` | identifier expressions only in practice | a borrowed `const T&` via `storage.get()` |
| `T&` | `std::reference_wrapper<T>` | `field_compiler<std::reference_wrapper<T>>` | identifier expressions only in practice | a borrowed mutable `T&` via `storage.get()` |

The key difference is:

- value arguments compile a value-producing typed expression
- reference arguments compile an identifier-reference typed expression

## 1. `field_compiler<T>`

Code:

- [lexer_typed.h:267](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/lexer_typed.h#L267)

### Accepted inputs

`field_compiler<T>` accepts:

- identifiers
- literals, if `literal_decodable<T>::value`
- function calls whose compiled result is exactly `typed_expression<..., T>`
- vectors if `T` is `std::vector<...>`
- sets if `T` is `std::set<...>`
- tuples if `T` is `std::tuple<...>`

### Identifier path

For an identifier:

1. it calls `cm.get_Identifier(...)`
2. it obtains a precompiled typed expression from the environment
3. it checks that the expression is actually `typed_expression<..., T>`
4. if so, it releases and returns it

This means value arguments can bind to:

- literals already stored in the environment
- variables produced by previous function calls
- composite expressions already compiled to the right type

### Literal path

If literals are decodable for `T`, it builds:

```cpp
typed_literal<Lexer, Compiler, T>
```

If not, it rejects literals with a type-specific diagnostic.

### Function-call path

For a function call:

1. it compiles the nested call
2. checks that the nested result is `typed_expression<..., T>`
3. returns that expression if types match

### Statement wrapper path

If the DSL argument is written as an assignment-like node, `field_compiler<T>`
also checks that the explicit argument name matches the expected field id.

That check is specific to function arguments and is not part of vector/tuple
element compilation.

## 2. `field_compiler<std::reference_wrapper<const T>>`

Code:

- [lexer_typed.h:433](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/lexer_typed.h#L433)
- [grammar_typed.h:388](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/grammar_typed.h#L388)

### Accepted inputs

Conceptually this specialization accepts only identifier-backed references.

Its real behavior is:

- identifiers: accepted if they name an environment value stored as
  `typed_literal<..., T>`
- literals: rejected
- function calls: rejected
- vectors/tuples: not a meaningful supported path here

So this is not a value compiler. It is a binder to an existing environment slot.

### Identifier path

For an identifier:

1. it calls `cm.get(...)`, not `cm.get_Identifier(...)`
2. it checks the stored runtime value
3. that value must currently be a `typed_literal<..., T>`
4. if so, it builds `typed_identifier_ref_const<..., T>(id)`

This is stricter than the value case.

It does not accept "anything that can compile to `T`".
It accepts "an environment variable currently stored as a literal value of `T`".

### Runtime path

At execution time, `typed_identifier_ref_const`:

1. calls `env.get(id)`
2. checks again that the stored value is `typed_literal<..., T>`
3. returns `std::cref(literal->value_ref())`

So the actual runtime storage is still the environment value.
The DSL node only remembers the identifier.

## 3. `field_compiler<std::reference_wrapper<T>>`

Code:

- [lexer_typed.h:518](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/lexer_typed.h#L518)
- [grammar_typed.h:417](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/grammar_typed.h#L417)

This is the mutable analogue of the previous case.

### Accepted inputs

Same practical rule:

- identifiers to existing environment variables of literal type `T`
- no literals
- no function calls
- no meaningful composite path

### Runtime path

`typed_identifier_ref`:

1. loads the environment entry
2. checks for `typed_literal<..., T>`
3. `const_cast`s that literal node
4. returns `std::ref(literal->value_ref())`

So mutable references work only because the underlying stored object is a
mutable literal node in the environment.

## 4. Side-by-side behavior

### Compile-time difference

`field_compiler<T>`:

- "compile an expression that yields a value of type `T`"

`field_compiler<std::reference_wrapper<const T>>`:

- "bind this argument to an already-stored environment variable of type `T`"

`field_compiler<std::reference_wrapper<T>>`:

- same as above, but mutable

### What each one queries

`field_compiler<T>` uses:

- `cm.get_Identifier(...)`

This follows the identifier compiler path and accepts any typed expression of
the right resulting type.

Reference-wrapper specializations use:

- `cm.get(...)`

This follows the runtime variable store and insists on an actual stored value.

### What node each one produces

`field_compiler<T>` may produce:

- `typed_literal<..., T>`
- a released `typed_expression<..., T>` from an identifier
- a compiled nested function expression
- a composite typed node

Reference-wrapper specializations produce:

- `typed_identifier_ref_const<..., T>`
- `typed_identifier_ref<..., T>`

They do not produce independent value expressions.

## 5. Argument adaptation inside `function_compiler`

Code:

- [lexer_typed.h:1027](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/lexer_typed.h#L1027)

After each argument expression is run, `function_compiler::adapt_argument` maps
storage to the actual parameter type.

### Value case

If `Param = T` and storage is `T`, the final branch is used:

```cpp
return static_cast<Param>(std::move(storage));
```

So the callee receives a value.

### Const-reference case

If `Param = const T&` and storage is `std::reference_wrapper<const T>`, the
first branch is used:

```cpp
return static_cast<Param>(storage.get());
```

So the callee receives a borrowed reference.

### Mutable-reference case

If `Param = T&` and storage is `std::reference_wrapper<T>`, the same first
branch is used:

```cpp
return static_cast<Param>(storage.get());
```

Again, the callee receives a borrowed reference, but mutable.

## 6. The same pattern for composite elements

For vectors, sets, and tuples, the analogous compiler is `element_compiler`,
not `field_compiler`.

Code:

- [lexer_typed.h:598](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/lexer_typed.h#L598)
- [lexer_typed.h:750](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/lexer_typed.h#L750)
- [lexer_typed.h:825](/home/lmoffatt/Code/macro_dr/macro_dr/include/macrodr/dsl/lexer_typed.h#L825)

The structure is intentionally parallel:

- `element_compiler<T>` mirrors `field_compiler<T>`
- `element_compiler<std::reference_wrapper<const T>>` mirrors the const-ref
  `field_compiler`
- `element_compiler<std::reference_wrapper<T>>` mirrors the mutable-ref
  `field_compiler`

The main difference is positional semantics:

- function arguments may check argument-name congruence
- vector/tuple/set elements are positional and do not carry argument ids

## 7. Main conclusion

The current DSL has two very different compilation modes:

### Value mode

- compile expressions
- allow literals and nested calls
- produce a value expression
- move or copy that value into the callee

### Reference mode

- bind to an existing environment variable
- reject literals and computed temporaries
- produce an identifier-reference expression
- hand the callee a borrowed reference through `reference_wrapper`

So `std::reference_wrapper<T>` is not "just another storage type".
It selects a different compilation regime.
