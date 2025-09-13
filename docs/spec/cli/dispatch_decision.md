# MacroDR Dispatch‑Mechanism Decision  
_Why we deliberately keep **`dynamic_cast`** for `typed_expression` down‑casts and postpone double‑dispatch visitors._

---

## 1  Context & Goal  
MacroDR’s DSL compiler produces a heterogeneous tree of `typed_expression<Lexer, Compiler, T>` nodes (`T = int, double, Filename …`).  Algorithms occasionally need the concrete `T`.  The question is **how to down‑cast** from `base_typed_expression` to the desired concrete node.

Our goal is:
1. **Rapid type‑set growth** – new `T`s appear frequently during research.  
2. **Minimal maintenance friction** – avoid touching central headers whenever a new type is added.  
3. **Safe failure** – mis‑casts must yield clean diagnostic paths, not UB.  
4. **Acceptable runtime cost** – but numerical kernels dominate, so dispatch micro‑cost is secondary.

---

## 2  Options Compared  
| ID | Mechanism | Runtime path | Add new `T` requires… | Cross‑project subset reuse | Failure semantics |
|----|-----------|--------------|-----------------------|----------------------------|-------------------|
| **A** | **`dynamic_cast` + RTTI** | 1 virtual call + RTTI lookup | _Nothing_ – new header only | ✅ **Easy** (base never changes) | `nullptr → Maybe_error` |
| **B** | Double‑dispatch visitor | 1 virtual + overload | **Edit `base_typed_expression`** _and_ every visitor | �� Requires recompiling a new base for each subset – hard to share | compile‑time exhaustiveness |
| **C** | Manual enum/string tag + `reinterpret_cast` | virtual + int/str compare | Add enum value & plumbing | �� Same central registry issue | UB if tag drifts |

*Visitor cost* is **code‑maintenance quadratic** (algorithms × types).  Runtime for all three is ~O(1) per node and negligible compared to simulations.

---

## 3  Key Observations
* **Growth friction**: Option B/C introduce a _"central registry"_ that must know **every concrete `T` in advance**.  If a partner project wants only a subset, they must fork the base class or keep unused v‑table slots – a maintenance “death‑kiss.”
* **Dynamic‐cast safety**: RTTI is compiled‑in; a bad cast returns `nullptr` that we already wrap in `Maybe_error`, meeting requirement 3 without UB.
* **Performance**: On modern compilers a `dynamic_cast` of pointer to polymorphic base is a single cached lookup; profile shows it dwarfed by statistical kernels.

---

## 4  Decision  
> **Stay with Option A (`dynamic_cast`) until a concrete need emerges.**

* **Zero friction** when new `typed_expression<T>` is added – no central edits.  
* **Cross‑project reuse**: different MacroDR plugins can include only the `T`s they need, sharing the same base class.  
* **Runtime & binary size**: RTTI tables are marginal next to Eigen/GSL/LAPACK.  
* **Safe diagnostics** already integrated via `Maybe_error`.

---

## 5  Situations That Could Re‑open the Question  
Re‑evaluate only if **all** of the following align:
1. **RTTI must be disabled** for platform or size reasons _and_  
2. The set of `typed_expression` types stabilises (few additions per year) _and_  
3. We introduce ≥ N (>10) **new independent algorithms** that each traverse the AST, making visitor boilerplate worth its cost.

Until those triggers occur, **avoid redesigning dispatch** – the maintenance overhead outweighs speculative benefits.

_Last updated: 2025‑07‑24_

