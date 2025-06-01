
# Compile-Time Semantic Values for Template Parameters in C++: Minimalist and Type-Safe

## Motivation

When writing high-level generic C++ code, especially scientific or library code, you often want to control behavior at compile time with parameters—recursion, algorithm variants, etc.
The standard way is to use plain `bool`, `int`, or `enum` template parameters:

```cpp
template <bool recursive, int averaging>
void my_func() { ... }
```

But as soon as you have more than two or three such parameters, problems appear:

* **Semantic confusion:** It’s easy to forget what the third `true` or `2` means.
* **Accidental misuse:** Mixing up the order of arguments or types is not caught by the compiler.
* **Readability:** Call sites are hard to read and review.

What if you could **encode intent and meaning directly in the type**, so your template parameters are self-documenting, type-checked, and remain fully constexpr?

---

## The Solution: Semantic Constexpr Values as Template Parameters

This pattern introduces a *semantic type* for each logical template parameter.
Each semantic value is:

* A unique type (so you can't mix them up)
* A constexpr value (so `if constexpr` works as usual)
* Instantiated in a single, minimal line per semantic parameter (no extra tag structs, no macros)
* Used as a template parameter; it can also be used as a runtime variable, although it will lose its constexpr properties in that context.

---

### 1. The Minimal Pattern: `semantic_constexpr_value`

```cpp

template <class Id, class T>
class semantic_constexpr_value {
public:
    T value;
    static constexpr bool is_variable = true;
    constexpr semantic_constexpr_value(T t_x) : value{t_x} {}
    constexpr semantic_constexpr_value() = default;
    constexpr auto& operator[](semantic_constexpr_value<Id>) const { return *this; }
};
```

This class acts as a semantic wrapper: the *type* is unique to each parameter, the *value* can be any compile-time constant.

---

### 2. Defining Semantic Compile-Time Parameters

For each parameter, define a semantic value type **in one line**:

```cpp
class uses_recursive_approximation
    : public semantic_constexpr_value<uses_recursive_approximation, bool> {};

class uses_averaging_approximation
    : public semantic_constexpr_value<uses_averaging_approximation, int> {};
```

---

### 3. Using Semantic Values in Templates

You can now write generic code that is:

* **Self-documenting** (names appear at call site)
* **Type-safe** (compiler enforces parameter meaning)
* **Constexpr-friendly** (can use `if constexpr (recursive.value)` directly)

Example:

```cpp
template <uses_recursive_approximation recursive, uses_averaging_approximation averaging>
void my_algorithm() {
    if constexpr (recursive.value) {
        // Recursive case
    }
    if constexpr (averaging.value == 2) {
        // Averaging==2 case
    }
}
```

Usage:

```cpp
my_algorithm<uses_recursive_approximation(true), uses_averaging_approximation(2)>();
```

---

### 4. Semantic Types and the C++ Community

This is a special case of “**semantic types**” or “**strong types**,” used to make code more robust and self-explanatory.
By encoding meaning into the *type* of the template parameter, you make your code safer and easier to understand, even as the template interface grows.
The key is that each parameter is *distinguished* at the type level, not just by position or value.

---

### 5. Why Not Use Enums or Structs?

* `enum class` is good for enumerations, but not for generic parameters (e.g., arbitrary `int` or `bool` values).
* Using plain `struct` types requires extra boilerplate and does not scale as elegantly or provide direct `.value` access.
* This pattern gives you the *minimal* syntax and full type safety.

---

### 6. Benefits Recap

* **Self-documenting template calls**: You see `uses_recursive_approximation(true)` at the call site, not just a `true`.
* **No mixups:** Compiler catches any misordering or misnaming.
* **Zero runtime cost:** Everything happens at compile time.
* **No boilerplate:** One line per semantic parameter type.

---

### 7. Full Minimal Example

```cpp
// Base pattern
template <class Id, class T>
class semantic_constexpr_value {
public:
    T value;
    static constexpr bool is_variable = true;
    constexpr semantic_constexpr_value(T t_x) : value{t_x} {}
    constexpr semantic_constexpr_value() = default;
    constexpr auto& operator[](semantic_constexpr_value<Id>) const { return *this; }
};

// Semantic parameter types (one line each)
class uses_recursive_approximation
    : public semantic_constexpr_value<uses_recursive_approximation, bool> {};
class uses_averaging_approximation
    : public semantic_constexpr_value<uses_averaging_approximation, int> {};

// Usage in template
template <uses_recursive_approximation recursive, uses_averaging_approximation averaging>
void my_algorithm() {
    if constexpr (recursive.value) {
        // Recursive logic
    }
    if constexpr (averaging.value == 2) {
        // Averaging=2 logic
    }
}

// Instantiation:
int main() {
    my_algorithm<uses_recursive_approximation(true), uses_averaging_approximation(2)>();
}
```

---

## 8. About “Semantic Types”

The idea of semantic types (sometimes called “strong typedefs” or “named types”) is a way to encode *meaning* into the type system, not just the value.
In this pattern, the *type* enforces the semantic role of each template parameter, so your code remains robust even as complexity grows.

For more, see:

* [Fluent C++: Strong Types for Strong Interfaces](https://www.fluentcpp.com/2016/12/08/strong-types-for-strong-interfaces/)
* [Modernescpp: Strongly Typed Typedefs](https://www.modernescpp.com/index.php/strongly-typed-typedefs/)

---

## Summary

* Use `semantic_constexpr_value` to declare template parameters with semantic meaning.
* Pass those semantic values as template parameters for self-documenting, type-safe, and constexpr-friendly generic code.
* This approach is ultra-minimal, has no runtime penalty, and avoids template confusion and accidental mixups—*especially as your codebase grows*.

---

