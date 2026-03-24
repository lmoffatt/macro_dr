

# MacroIR Codebase Map

This document provides an overview of the key modules, major classes/templates, and their relationships in the MacroIR scientific computing repository. Utilities, math libraries, and modeling logic are mapped separately for clarity.

---

## Top-Level Structure

```
.
├── main.cpp
├── *.h (core headers)
├── scripts/                  # Data/parameter scripts (no code)
├── slurm*/                   # Batch scripts for HPC jobs
├── test_macro_dr/            # Tests and validation
├── Testing/                  # Temporary and log files
├── models*/                  # Model definitions and parameters
├── experiments/              # Experimental data
├── docs/                     # Documentation
├── catch2/                   # Third-party testing library
├── clusters/                 # Cluster execution scripts
├── command_lines/            # Example CLI input scripts
├── macro_dr/                 # (possibly legacy binary/data)
├── multi_task/               # Multi-model execution scripts
├── parameters.h/.h           # Model parameter handling
└── ... (see below)
```

---

## Core Modules

### 1. **Modeling & Inference Logic**

* **qmodel.h**
  *Core classes and templates for model representation and simulation.*
* **allosteric\_models.h, models\_MoffattHume\_linear.h, models\_Ag\_log\_baseline/, models\_posterior/, models\_different\_inact/**
  *Definitions for specific biophysical models, kinetic schemes, and posterior distributions.*
* **mcmc.h, parallel\_tempering.h, parallel\_tempering\_fraction.h, parallel\_tempering\_linear\_regression.h**
  *MCMC, tempering, and Bayesian inference engines.*
* **bayesian\_linear\_regression.h**
  *Bayesian regression model implementation.*

### 2. **CLI and Scripting**

* **CLI\_macro\_dr.h, CLI\_macro\_dr\_base.h, CLI\_thermo\_evidence*.h, CLI\_function\_table.h, CLI\_regular\_types.h, CLI\_grammar.h, CLI\_base.h*\*
  *Command-line interfaces for model specification, evidence calculation, and type handling.*
* **command\_lines/**
  *Input scripts for running models through the CLI.*
* **scripts/**
  *Parameter sets, experiment inputs, and simulation configs.*

### 3. **Math and Scientific Computing**

* **multivariate\_normal\_distribution.h, matrix\_random.h, matrix\_derivative.h, matrix.h, exponential\_matrix.h**
  *Linear algebra, random matrix, and distribution tools.*
* **lapack\_headers.h, gsl\_integrate.h, lgamma.h, gsl\_\*.h**
  *Bindings for GSL/LAPACK special functions and integration.*
* **distributions.h, random\_samplers.h, function\_memoization.h, function\_measure\_verification\_and\_optimization.h**
  *Random sampling, distribution tools, and function optimization.*
* **type\_algebra.h, variables.h, variables\_derivative.h, parameters\_distribution.h, parameters\_derivative.h**
  *Template metaprogramming for model parameters and algebra.*

### 4. **Utilities**

* **maybe\_error.h, general\_algorithm\_on\_containers.h, general\_output\_operator.h, derivative\_test.h, derivative\_operator.h, indexed.h, fold.h, continuation.h**
  *Generic utilities for error handling, algorithm abstraction, output, and functional programming.*
* **lexer\_untyped.h, lexer\_typed.h, grammar\_Identifier.h, grammar\_typed.h, grammar\_untyped.h**
  *Custom lexers and grammar definitions (possibly for parsing DSLs or model input).*

### 5. **Testing**

* **test\_macro\_dr/**
  *C++ tests, validation scripts, and example outputs (see `main.cpp`, `tst_*.cpp`, and `examples/`).*
* **Testing/**
  *Temporary test logs and files.*
* **catch2/**
  *Third-party unit testing framework (Catch2).*

---

## Major Classes & Templates

> **Note:** For full class lists, see the individual headers. Here are key abstractions inferred from file names and conventions:

* `QModel` (in `qmodel.h`): Main model abstraction (likely templated).
* `MCMC` and `ParallelTempering` classes: Inference engine for sampling.
* `AllostericModel` (in `allosteric_models.h`): Biophysical model class.
* `BayesianLinearRegression`: Bayesian linear regression logic.
* Utility templates for error handling (`MaybeError`), functional (`Fold`, `Continuation`), algebra (`TypeAlgebra`), and variable management.

---

## Relationships

* **Model classes** (e.g., `QModel`, `AllostericModel`) depend on utilities (e.g., random samplers, algebra templates) and are parameterized by data from the `models*` folders.
* **Inference engines** (MCMC, parallel tempering) are built around these model classes and interact with parameter and data scripts.
* **CLI modules** serve as a bridge between command-line usage and internal modeling/inference logic.
* **Testing** imports core modules and runs end-to-end validation with sample data/parameters.

---

## Special Purpose/Legacy

* **macrodr\_h.zip, macro\_dr/**: Legacy binaries or bundled headers.
* **old\_runs/, Testing/Temporary/, .qtc\_clangd/**: Output, logs, IDE/project config, and cache.
* \**slurm*/ and clusters/\*\*: Scripts for running jobs on clusters.

---

## Math Libraries & Scientific Utilities

* **LAPACK/GSL wrappers:**
  `lapack_headers.h`, `gsl_integrate.h`, `gsl_*.h`
  Provide numerical integration, linear algebra, and special functions.
* **Random sampling and distributions:**
  `random_samplers.h`, `distributions.h`, `multivariate_normal_distribution.h`
* **Matrix calculus and algebra:**
  `matrix_derivative.h`, `exponential_matrix.h`, `type_algebra.h`

---

## Modeling Logic

* All model/parameter definitions are under:

  * `models/`, `models_Ag/`, `models_posterior/`, etc.
  * Each file represents a kinetic scheme, parameter set, or prior/posterior.
* Templates and class logic for combining models with data and inference engines are in headers like `qmodel.h`, `allosteric_models.h`, and `mcmc.h`.

---

## Testing & Validation

* All C++ test logic in `test_macro_dr/`
* Output logs, example runs, and benchmarks in `test_macro_dr/examples/` and `Testing/`
* Parameterized test data in `scripts/` and `models*/`

---

## Documentation

* Markdown docs and command registry under `docs/`
* Includes blog posts, template instantiation reports, and CLI documentation.

---

## Build System

* `CMakeLists.txt` in root and subdirs.
* `.clang-format` and `.clang-format-llvm` for style.

---

## Quick Reference Table

| Module/Folder            | Purpose                              |
| ------------------------ | ------------------------------------ |
| `qmodel.h`               | Core model logic and templates       |
| `allosteric_models.h`    | Biophysical model classes            |
| `mcmc.h`                 | MCMC inference logic                 |
| `parallel_tempering*`    | Parallel tempering (sampling)        |
| `CLI_*`                  | Command-line interfaces              |
| `models*/`               | Model definitions/parameters         |
| `matrix_*.h`             | Matrix and linear algebra utilities  |
| `random_samplers.h`      | Random number/distribution utilities |
| `function_memoization.h` | Function caching/memoization         |
| `test_macro_dr/`         | C++ testing and example validation   |
| `slurm*/`                | Batch/cluster job scripts            |
| `catch2/`                | Unit test framework                  |
| `docs/`                  | Documentation and blog               |

---

## How to Extend/Explore Further

* For class diagrams or dependency graphs, consider using [Doxygen](https://www.doxygen.nl/) or [clangd-index](https://clangd.llvm.org/).
* For a test coverage or missing test map, explore `test_macro_dr/` and `Testing/`.

---

**Feel free to request a finer breakdown of any submodule or a graphical map (class/flow diagram) as needed.**

---

*Generated automatically from file tree and filename heuristics. For module-level details, see headers and docstrings.*

