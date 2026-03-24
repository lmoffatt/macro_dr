# MacroIR Proof-Oriented Design Manifesto

---

## 0. Introduction

This document formalizes the design philosophy that has emerged during the conceptual refactoring of the MacroIR scientific modeling codebase.

We propose a new paradigm: **Proof-Oriented Design (POD)** — a method for building scientific software that encodes mathematical invariants directly into the code structure, enabling correctness by construction.

---

## 1. The Problem

Traditional scientific software (and MacroIR's initial codebase) suffers from:

* Highly coupled components
* Weak enforcement of mathematical correctness
* Manual and scattered invariant checks
* Heavy reliance on testing without formal guarantees

Our goal is to refactor MacroIR to achieve:

* Semantic correctness
* Maintainability
* Extensibility
* Verifiability

---

## 2. Foundational Insights

### 2.1 Dependencies vs Concerns

* **Dependencies** are logical and objective (who calls whom, includes, data flow).
* **Concerns** are semantic and subjective (which module owns which part of the model).
* Good design requires both clear dependency graphs and clean separation of concerns.

### 2.2 Invariants as Core Semantics

* Invariants define scientific correctness:

  * Stochastic matrices: rows sum to zero
  * Probabilities: sum to one
  * Transition rates: non-negative
* Invariants drive code structure.

### 2.3 Types as Propositions (Curry–Howard Correspondence)

* Type == Proposition
* Valid object == Proof of invariant
* Type construction == Proof of correctness

### 2.4 Functions as Proof Transformers

* Functions transform proofs of correctness.
* Function signatures encode input/output invariants.
* Testing becomes empirical sampling of correctness.

### 2.5 Testing as Partial Proof

* Tests provide existential witnesses for correctness.
* Complete proofs would eliminate testing; testing remains useful when full proofs are impractical.

---

## 3. The Long-Term Ideal: Homotopy Type Theory (HoTT)

* HoTT unifies types, proofs, and computation.
* In HoTT:

  * Types represent mathematical spaces.
  * Equalities are paths.
  * Proofs are first-class citizens.
* HoTT provides the ultimate foundation for scientific software correctness.

---

## 4. C++ as the Host Language: Strengths and Limits

### 4.1 Strengths

* Templates can simulate dependent types.
* Compile-time computations possible via constexpr and metaprogramming.
* Can partially encode invariants at compile-time.

### 4.2 Limitations

* No native dependent types.
* No first-class proofs.
* No semantic enforcement of invariants.
* Unsafe operations remain possible.

---

## 5. Two Enforcement Strategies

### 5.1 External Checker

* Build a static analysis tool:

  * Registers HoTT-safe types and functions.
  * Checks categorical closure.
  * Forbids escape from safe subcategory.
* Analogous to borrow checkers or refinement type systems.

### 5.2 Internal Monad Discipline

* Embed invariants into the type system via monadic structure:

```cpp
template<typename T>
struct HoTT {
    T value;
    // only safe operations allowed
};
```

* Functions become morphisms: `T -> HoTT<U>`
* Composition guarantees safety.
* Compiler enforces safety locally.

---

## 6. The Category Theory Perspective

* We define a safe subcategory $\mathcal{C} \subset \mathcal{CPP}$
* Objects: HoTT-safe types
* Morphisms: Invariant-preserving functions
* Composition: Chained safe operations
* Identity: Safe identity morphisms
* Monad: Reflective embedding functor into $\mathcal{C}$

---

## 7. Proof-Oriented Design Principles

1. **Types encode invariants**
2. **Construction equals proof**
3. **Functions transform proofs**
4. **Testing approximates missing proofs**
5. **Monads preserve safety**
6. **Precompilers enforce global structure**
7. **Software mirrors algebraic models**

---

## 8. Scientific Software as Executable Mathematics

By following Proof-Oriented Design, MacroIR aims to:

* Translate scientific theories directly into type-safe code.
* Ensure correctness through construction rather than testing.
* Build models whose correctness is encoded and checkable.
* Provide both scientific rigor and software robustness.

---

## 9. Future Work

* Prototype HoTT-monadic core in C++
* Build external semantic precompiler
* Explore embedding HoTT via DSLs or proof assistants
* Move toward fully verified scientific modeling frameworks

---

## 10. Summary Diagram

```text
Scientific Model  →  Invariants  →  Types  →  Proofs  →  Code

C++  →  HoTT Subcategory  →  Monad Discipline / External Checker  →  Safety
```

---

## 11. Conclusion

Proof-Oriented Design is not merely a coding style; it is a new architecture for trustworthy scientific software. It blends:

* Category theory
* Type theory
* Scientific modeling
* Software architecture

MacroIR stands as a pioneering attempt to embody these principles in real, high-integrity scientific computation.
