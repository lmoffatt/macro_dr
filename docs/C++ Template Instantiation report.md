# C++ Template Instantiation Analysis Report

## Executive Summary

This report summarizes the investigation into template instantiation counts in the MacroDR project, focusing on the risks of exponential template expansion and the impact on compilation time and binary size. The analysis demonstrates that the current template usage, while extensive, does not exhibit exponential growth, and is mostly driven by multilinear combinations of models, flags, and STL adaptors. Recommendations and observations are provided for future maintainability and performance.

---

## 1. Context & Motivation

Template metaprogramming is a central design technique in MacroDR, allowing for flexible and efficient code generation for many model/flag combinations at compile time. However, C++ templates can cause compilation slowdowns or excessive binary bloat if combinatorial or exponential instantiation patterns emerge. This report aims to:

* Quantify template instantiation in the project.
* Detect pathological cases (exponential explosion).
* Guide future refactors, e.g., towards type-erasure if needed.

---

## 2. Methodology

1. Template instantiation counts were collected (e.g., via `-fdump-ipa-all`, `clang -Xclang -ast-dump` or equivalent tools).
2. Results were analyzed and aggregated by template/function/class name, including instantiation context and example parameters.
3. CSV analysis was performed using Python/Pandas for summary statistics and histograms.

---

## 3. Key Findings (Summary)

* **No exponential explosion detected:** The highest number of instantiations for project-specific templates is around 1,500 (for core algorithmic routines). STL/lambda helpers (e.g., `operator()`, `__invoke_impl`) have >10,000, but this is standard and not problematic.
* **Template combinatorics are multilinear, not exponential:** Instantiations reflect the number of model classes × flag combinations × STL adaptors/lambdas.
* **STL and lambda machinery dominate the count:** Most instantiations are due to `std::variant`, lambdas, `operator()`, and STL adapters.
* **Main project-specific templates (with >100 instances):**

  * `Likelihood_Model`, `thermo_evidence`, `step_stretch_thermo_mcmc`, `log_Likelihood`, etc.
* **Models involved:** The majority are from combinations of `Allost1`, `Model0`, lambdas, and their variant-generated classes.

---

## 4. Detailed Results

### Top 20 Templates by Instantiations

| Template                           | Instantiations |
| ---------------------------------- | -------------- |
| operator()                         | 23,092         |
| \_\_invoke\_impl                   | 18,858         |
| remove\_reference                  | 17,655         |
| \_\_invoke                         | 17,376         |
| forward                            | 16,541         |
| \_Variadic\_union                  | 12,640         |
| \_\_gen\_vtable\_impl              | 10,026         |
| \_\_element\_by\_index\_or\_cookie | 8,544          |
| invoke                             | 6,927          |
| \_\_get\_n                         | 6,576          |
| report                             | 6,324          |
| \~                                 | 6,042          |
| variant                            | 5,122          |
| report\_model\_all                 | 3,988          |
| f                                  | 3,960          |
| get                                | 3,830          |
| fold                               | 3,591          |
| log\_Likelihood                    | 3,465          |
| \_\_get                            | 3,396          |
| \_Tuple\_impl                      | 3,105          |

*(See CSV for full breakdown)*

#### Instantiation Histogram

* Most templates are instantiated only a few times (2-4).
* Only 100 templates exceed 100 instantiations.
* Only STL/lambda helpers reach >10,000.

#### Main Models/Flags in Combinations

Most frequent: `Allost1`, `Model0`, and lambdas as variant types.

---

## 5. Conclusions & Recommendations

* **No evidence of uncontrolled combinatorial explosion.** The template usage is heavy but multilinear and under control.
* **Compile times (7–10 min) are reasonable for this template complexity.**
* **No urgent need to refactor for type-erasure or reduce template use, unless faster iteration or smaller binaries are required.**
* If future code growth increases the number of models or flags significantly, consider:

  * Refactoring interface code with type-erasure for commonly passed objects (e.g., `std::function`, pImpl, or similar).
  * Replacing some deep template layers with runtime polymorphism for algorithm dispatch, if needed.

---

## 6. Appendices

### A. Example Template Instantiation Record

```
Template: Likelihood_Model
Instantiations: 676
Context: function
Example: <macrodr::uses_adaptive_aproximation{...}, macrodr::uses_recursive_aproximation{...}, macrodr::uses_averaging_aproximation{...}, macrodr::uses_variance_aproximation{...}, macrodr::uses_variance_correction_aproximation{...}, macrodr::Model_Patch<macrodr::Allost1>::Model<...>>
```

### B. Analysis Script (Python/Pandas)

```python
import pandas as pd
# Place results.csv in the same directory as this script
# ... [Script snippet as in previous messages] ...
```

---

## 7. References

* C++ Best Practices: [CppCoreGuidelines](https://github.com/isocpp/CppCoreGuidelines)
* [Type Erasure in Modern C++](https://www.modernescpp.com/index.php/type-erasure-the-ultimate-guide/)
* [Abseil type erasure patterns](https://abseil.io/tips/152)

---
* Generated with ChatGPT, reviewed by Luciano, 2025-05-31.*

