# The HoTT Contagion Principle

---

## 0. Introduction

This document extends the Proof-Oriented Design Manifesto by introducing a key organizing principle that naturally arises from our architectural philosophy: **The HoTT Contagion Principle**.

---

## 1. The Core Idea

> **Once you operate on a HoTT-safe object, any operation downstream must remain HoTT-safe, or fail to compile.**

This principle makes invariant safety *infectious* throughout the codebase, ensuring that once correctness is established at any point, it propagates forward automatically unless explicitly broken.

---

## 2. Formal Description

- We define a safe subcategory \(\mathcal{C} \subset \mathcal{CPP}\) inside C++.
- Objects: HoTT-safe types (invariant-correct data structures).
- Morphisms: Functions that preserve invariants.
- Once a computation begins inside \(\mathcal{C}\), it must remain inside \(\mathcal{C}\).

If a function attempts to apply an unsafe operation on a HoTT-safe object, the compiler rejects it as a type error.

---

## 3. Monadic Formulation

We achieve contagion inside C++ via a **HoTT Monad**:

```cpp
template<typename T>
struct HoTT {
    T value; // Proven to satisfy invariants
};

// Invariant-preserving functions:
template<typename T, typename U>
HoTT<U> bind(const HoTT<T>& x, function<T, HoTT<U>> f) {
    return f(x.value);
}
```

- All safe functions have signature `T -> HoTT<U>`.
- Composition is forced to respect the monadic structure.
- The type system enforces that once inside HoTT, you stay inside HoTT.

---

## 4. Category-Theoretic View

- The HoTT Monad creates a **reflective subcategory** inside the host language.
- Composition of HoTT-safe morphisms is closed.
- Non-HoTT functions are not composable with HoTT-safe computations without explicit unsafe conversions.
- This embeds a safe world inside an unsafe host.

---

## 5. Why This Works

- The monadic structure enforces local safety through types.
- The compiler becomes an active proof assistant.
- Once a type is known to satisfy its invariants, downstream correctness is automatic.
- Violations become type errors, not runtime bugs.

---

## 6. The Contagion Mechanism

| Step                                | Effect                              |
| ----------------------------------- | ----------------------------------- |
| Start with a HoTT-safe object       | `HoTT<T>`                           |
| Apply invariant-preserving function | `T -> HoTT<U>`                      |
| Result remains HoTT-safe            | `HoTT<U>`                           |
| Attempt to apply unsafe function    | Compilation error                   |
| Escape from HoTT                    | Requires explicit, unsafe operation |

---

## 7. Practical Benefits for Scientific Code

- Guarantees mathematical model consistency.
- Eliminates entire classes of runtime errors.
- Reduces need for scattered manual checks.
- Shifts correctness from testing to type enforcement.
- Aligns code structure directly with the algebra of the scientific domain.

---

## 8. Functional Language Analogies

| Language   | Contagion Mechanism              |
| ---------- | -------------------------------- |
| Haskell    | IO monad, STM monad              |
| Rust       | `Result<T, E>`, ownership system |
| Idris/Agda | Dependent type proofs            |
| Scala      | Tagless final encoding           |

---

## 9. Summary Statement

> **HoTT Contagion allows correctness to propagate automatically by design.**
>
> Once inside HoTT safety, all computations remain safe unless explicitly broken.
>
> This forms the foundation for building mathematically trustworthy scientific software.

---

## 10. Next Steps

- Formalize the HoTT-safe type system for MacroIR.
- Enumerate approved invariant-safe types and morphisms.
- Implement the monadic structure inside C++.
- Begin enforcing contagion discipline throughout MacroIR's refactor.

---

## 11. Closing Remark

The HoTT Contagion Principle turns correctness into an architectural invariant of the entire system, allowing scientific software to approach the rigor of formal mathematics while remaining pragmatically executable inside existing host languages.

