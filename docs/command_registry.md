# MacroDR Command Registry — Quick Start Guide

## What is this?

MacroDR’s Command Registry is a system for exposing C++ functions as *commands* in a custom Domain-Specific Language (DSL). It allows users and collaborators to script and automate complex model operations, data analysis, and simulation — all in a way that is both flexible and type-safe.

---

## How does it work?

1. **C++ functions (“handlers”) are registered as commands** with the registry, using clear string names and argument names.
2. **User scripts (in the MacroDR DSL)** call these commands by name, passing arguments as needed.
3. **The registry parses, checks types, and dispatches the function call** at runtime, returning results or errors.

---

## A. Example: Registering a Command

Suppose you have the following function in C++:
```cpp
std::string get_Experiment_file(std::string filename, double frequency_of_sampling, double initial_ATP);
````

You can register it as a DSL command like this:

```cpp
cm.push_function(
    "get_Experiment_file",
    dcli::to_typed_function<std::string, double, double>(
        &get_Experiment_file,
        "filename", "frequency_of_sampling", "initial_ATP"
    )
);
```

Now the DSL user can write in their script:

```
get_Experiment_file("experiment1.csv", 10.0, 2.5)
```

---

## B. What does a script look like?

A typical MacroDR DSL script might look like:

```
exp = get_Experiment_file("experiment1.csv", 10.0, 2.5)
obs = get_Observations("obs1.csv")
nparams = get_num_parameters("my_model")
sim = simulate_trajectory("my_model", exp, 100.0)
lik = compute_likelihood(obs, "my_model", sim)
evi = compute_thermo_evidence("my_model", sim, 1.0)
```

---

## C. How do I add a new command?

1. **Define the C++ function** (use only supported types as arguments!).
2. **Register it with `push_function`**, specifying argument names as strings, matching their order.
3. **Test in a script.** If the argument count or types do not match, the system will return an error with details.

---

## D. Tips and Caveats

* **Command and argument names should be clear** — the DSL is user-facing!
* If two commands have very similar signatures, use descriptive names to avoid confusion.
* If you get type errors, check both the registration and the script for argument mismatches.
* It is possible to register functions with complex argument types, but prefer primitive or well-documented structures for clarity.
* To remove redundancies or refactor, look for commands with overlapping functionality and unify their signatures if possible.

---

## E. Where is everything?

* Core registration logic: `lexer_typed.h`, `grammar_typed.h`, `CLI_function_table.h`
* Function handlers: `CLI_macro_dr.h`, `CLI_likelihood.h`, `CLI_thermo_evidence_dts.h`, etc.
* See also: `main.cpp` for real usage patterns.

---

## F. If you get stuck…

* Check the error message (they usually tell you the argument name/type problem).
* See example registrations in `CLI_function_table.h` or other `CLI_*.h` files.
* Ask for help or look at this guide’s advanced section!

---

*For more detail, see the [Advanced Developer Guide](#).*




