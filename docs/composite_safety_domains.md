# Composite Safety Domains: The Scalability Layer

---

## 0. Introduction

This document extends the Contagiously Safe Design Pattern and the Safety Explosion Limit by introducing a pragmatic solution to managing type explosion: **Composite Safety Domains (CSDs)**.

By coarsening safety granularity, we preserve contagious correctness while controlling the combinatorial growth of type instantiations.

---

## 1. The Fundamental Insight

> **Collapse multiple fine-grained safety properties into composite safety domains that enforce multiple invariants simultaneously as atomic correctness zones.**

This trades fine-grained diagnosis for scalable enforcement.

---

## 2. The Safety Granularity Coarsening Principle

- Instead of encoding every safety property as an independent type layer:
  - HoTT safety
  - Units consistency
  - Memory safety
  - Numerical stability
  - Side-effect safety

- We aggregate them into bundled safety classes:

```cpp
template<typename T>
struct ScientificSafe {
    T value;
    // Composite enforcement of all critical safety dimensions
};
```

---

## 3. The Key Tradeoff

| Gain | Loss |
| ---- | ---- |
| Controlled type growth | Coarser safety granularity |
| Simpler contagion rules | Loss of per-property failure diagnostics |
| Linear instantiation growth | Reduced orthogonality |
| Debuggable type system | Less introspection on individual invariants |

---

## 4. Type Explosion Control

### Type Growth Comparison

| Strategy | Type Instantiations |
| -------- | ------------------- |
| Fine-grained nesting | \(M \times 2^N\) |
| Composite domains | \(M \times K\) (where \(K\) is small) |

This dramatically reduces the combinatorial explosion while preserving correctness propagation.

---

## 5. Contagion Discipline Remains Intact

- The core principle of contagious safety still applies:
  - Once inside a composite safety domain, all downstream operations must remain inside.
  - Violations become type errors.
  - Composition closure remains enforced.

---

## 6. The Meta-Architecture

```text
+---------------------------------+
| ScientificSafe<T>               |
|  - Algebraic Safety (HoTT)      |
|  - Memory Safety                |
|  - Units Consistency            |
|  - Numerical Stability          |
|  - Side-effect control          |
+---------------------------------+

Contagion applies at ScientificSafe<T> level.
```

---

## 7. Practical Language Parallels

| Language | Analogy |
| -------- | ------- |
| Rust | `unsafe {}` vs fully safe contexts |
| Haskell | Monadic zones: `IO`, `STM`, `Safe`, etc. |
| Ada/SPARK | Modular contract zones |
| F* / Dafny | Verification modes per module |

---

## 8. Meta-Theorem Extension

> **Contagiously Safe Systems remain practically scalable when safety dimensions are bundled into composite correctness domains that maintain closure and enforce contagion at coarser granularity.**

---

## 9. Scientific Software Implications

- This design allows:
  - Full safety propagation
  - Compile-time tractability
  - Zero-cost abstraction
  - Simplified cognitive model for developers
  - Easy incremental adoption in large codebases

- The domain-specific composite bundles can reflect scientific modeling boundaries naturally.

---

## 10. The Ultimate Safety Design Stack

```text
Proof-Oriented Design
    ↓
HoTT Contagion Principle
    ↓
Contagiously Safe Design Pattern
    ↓
Safety Explosion Limit
    ↓
Composite Safety Domains (Scalable Enforcement)
```

---

## 11. Next Steps

- Define canonical safety bundles for MacroIR:
  - `ScientificSafe`
  - `PhysicalSafe`
  - `NumericalSafe`
- Build first minimal prototypes in C++20 using concepts.
- Investigate hybrid designs with selective finer-grained subdomains.
- Evaluate migration pathways for legacy code integration.

---

## 12. Conclusion

Composite Safety Domains offer a tractable, enforceable, and scalable framework for building proof-oriented, correctness-preserving scientific software.

They balance mathematical rigor with practical engineering constraints, allowing high-integrity scientific computation to be both formally grounded and computationally feasible in real-world production environments.

