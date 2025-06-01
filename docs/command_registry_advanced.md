# MacroDR Command Registry — Advanced Developer Guide

---

## Table of Contents

- [1. Motivation and Architecture](#1-motivation-and-architecture)
- [2. Core Concepts](#2-core-concepts)
- [3. Key Classes and Files](#3-key-classes-and-files)
- [4. Registration Workflow (Under the Hood)](#4-registration-workflow-under-the-hood)
    - [4.1. Registration Flow Diagram](#41-registration-flow-diagram)
- [5. Execution Workflow (From Script to C++)](#5-execution-workflow-from-script-to-c)
    - [5.1. Execution Flow Diagram](#51-execution-flow-diagram)
- [6. Naming Conventions and Design Choices](#6-naming-conventions-and-design-choices)
- [7. Extending and Refactoring](#7-extending-and-refactoring)
- [8. Troubleshooting & FAQ](#8-troubleshooting--faq)
- [9. Performance and Template Expansion](#9-performance-and-template-expansion)
- [10. Future Directions and Migration Notes](#10-future-directions-and-migration-notes)

---

## 1. Motivation and Architecture

MacroDR exposes a custom, type-safe, and extensible DSL for scientific modeling. The Command Registry enables users to call C++ functions by name from scripts, mapping arguments and types automatically. This design ensures:
- Flexible scripting and automation.
- Strong type guarantees.
- Ease of extension for new commands or models.

**Key design features:**
- Each command is registered with a unique name and argument list.
- Argument types are enforced at registration and at script execution.
- Errors are reported with meaningful messages.

---

## 2. Core Concepts

- **Command:** A C++ function made accessible via the DSL (e.g., `get_Experiment_file`).
- **Handler:** The actual C++ function or functor.
- **Argument Binding:** Mapping DSL script arguments to C++ types and names.
- **Type Checking:** Automatic validation at both registration and runtime.

---

## 3. Key Classes and Files

- **Registration/Dispatch:**
  - `Compiler`: Holds the registry of commands (`m_func`).
  - `push_function`: Adds a command to the registry.
- **Function Wrappers:**
  - `base_function_compiler`, `function_compiler`: Templates that handle type erasure and dispatch.
  - `to_typed_function`: The factory that generates wrappers for function registration.
- **Arguments:**
  - `field_compiler`: Maps argument names/types.
- **Relevant Files:**
  - `lexer_typed.h`, `grammar_typed.h`, `CLI_function_table.h`, `CLI_macro_dr.h`, `CLI_likelihood.h`, `CLI_thermo_evidence_dts.h`, `main.cpp`

---

## 4. Registration Workflow (Under the Hood)

1. **Define your handler:**  
   ```cpp
   std::string get_Experiment_file(std::string filename, double frequency_of_sampling, double initial_ATP);
````

2. **Register it:**

   ```cpp
   cm.push_function(
       "get_Experiment_file",
       dcli::to_typed_function<std::string, double, double>(
           &get_Experiment_file,
           "filename", "frequency_of_sampling", "initial_ATP"
       )
   );
   ```
3. **What happens internally:**

   * `to_typed_function` generates a `function_compiler` template instance for the exact argument signature.
   * The function is stored in a registry as a `base_function_compiler` pointer (type-erased).
   * Argument names are stored for runtime mapping and documentation.

### 4.1. Registration Flow Diagram

```text
+----------------------------+
| Define C++ Handler         |
|----------------------------|
| get_Experiment_file(...)   |
+----------------------------+
           |
           v
+-----------------------------+
| Register with Registry      |
|-----------------------------|
| cm.push_function( ... )     |
+-----------------------------+
           |
           v
+-------------------------------------+
| to_typed_function Factory           |
|-------------------------------------|
| Creates function_compiler instance  |
| with type info and arg names        |
+-------------------------------------+
           |
           v
+----------------------------+
| Registry (Compiler.m_func) |
|----------------------------|
| Stores as base_function_   |
| compiler pointer           |
+----------------------------+
```

---

## 5. Execution Workflow (From Script to C++)

1. **User writes a script:**

   ```dsl
   exp = get_Experiment_file("experiment1.csv", 10.0, 2.5)
   ```
2. **Script is parsed:**

   * Lexer and grammar modules parse the command and its arguments.
3. **Lookup and dispatch:**

   * The command registry finds the registered command by name.
   * Arguments are type-checked and mapped by position/name.
   * The underlying C++ function is invoked with the converted arguments.
4. **Return value (or error) is handled** and can be used in further script expressions.

### 5.1. Execution Flow Diagram

```text
+------------------+
|  DSL Script      |
|------------------|
| exp =            |
|   get_Experiment_file("exp.csv", 10, 2.5)  |
+------------------+
           |
           v
+---------------------+
|   Lexer/Parser      |
|---------------------|
| Tokenizes and       |
| parses the script   |
+---------------------+
           |
           v
+-------------------------------+
| Command Registry (Compiler)   |
|-------------------------------|
| Looks up "get_Experiment_file"|
| and argument names/types      |
+-------------------------------+
           |
           v
+-----------------------------+
| Function Wrapper            |
| (function_compiler,         |
|  type erasure)              |
|-----------------------------|
| Checks types, converts args |
+-----------------------------+
           |
           v
+--------------------------+
| C++ Handler Function     |
|--------------------------|
| get_Experiment_file(...) |
+--------------------------+
           |
           v
+--------------------------+
|   Return value/result    |
|--------------------------|
| Used in next DSL step    |
+--------------------------+
```

---

## 6. Naming Conventions and Design Choices

* **Command names:** Should be descriptive and unique.
* **Argument names:** Use domain-relevant, unambiguous terms.
* **Handler signatures:** Prefer simple types (primitives or well-documented structs).
* **Template usage:** Template expansion is linear per unique argument signature.

---

## 7. Extending and Refactoring

* **Adding a new command:**

  * Define the handler.
  * Register with `push_function`.
  * Test in scripts.

* **Refactoring tips:**

  * Identify overlapping commands and unify where possible.
  * Remove redundant or deprecated commands.
  * Group commands by topic for maintainability.

* **Advanced:**

  * You can add support for new argument types by extending `field_compiler` logic.
  * Complex types are possible, but require careful mapping and documentation.

---

## 8. Troubleshooting & FAQ

* **Argument count/type mismatch:**

  * Check both registration and script for alignment.

* **Unclear errors:**

  * Look at the error message details and at argument order/names.

* **Adding new types:**

  * Extend `field_compiler` and document new usage.

* **Common pitfalls:**

  * Registering lambdas with unique signatures leads to more template expansion.
  * Inconsistent naming between registration and script.
  * Missing argument documentation.

---

## 9. Performance and Template Expansion

* **Each unique command signature expands a new template instance.**
* **Growth is linear** (not exponential): hundreds of commands are feasible.
* **Debug builds** are larger due to symbol retention; release builds are compact.
* **For minimal build times:**

  * Reuse signatures when possible.
  * Avoid unnecessary overloads.

---

## 10. Future Directions and Migration Notes

* **If you outgrow the custom registry:**

  * pybind11, sol2 (Lua), Chaiscript, and others offer robust alternatives for scripting and command registration.
  * Migration is possible: most commands are already modular.
  * Consider community libraries if you need Python integration or dynamic scripting.

---

## See Also

* [Quick Start Guide](#)
* Example command registrations in `CLI_function_table.h`
* Example scripts in `tests/` or `examples/`
* For questions or help, see `README.md` or contact maintainers.

---

