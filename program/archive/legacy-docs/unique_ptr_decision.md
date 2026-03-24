# MacroDR Memory‑Ownership Decision

*Why the codebase ****continues**** to use **`std::unique_ptr`** for IR nodes today and what the ****future migration path**** looks like.*

---

## 1  Current Contract

- **All factory routines** (`Compiler::get_Identifier`, literal builders, etc.) return
  ```cpp
  Maybe_error<std::unique_ptr<base_typed_expression<Lexer,Compiler>>>
  ```
  The caller gains **exclusive ownership** and is free to move, clone, or destroy.
- Internal passes pass nodes further by **moving** the `unique_ptr` or by borrowing a raw pointer.

This interface has shipped for months and underpins every working CLI command.

---

## 2  Why we *stay* with `unique_ptr` **for now**

| Requirement                                  | `unique_ptr` status                                                 |
| -------------------------------------------- | ------------------------------------------------------------------- |
| **Zero friction when new node types appear** | *Met* – no central edits needed.                                    |
| **Fail‑fast lifetime errors**                | *Met* – double‑free impossible; use‑after‑move caught by compiler.  |
| **Performance inside tight numeric loops**   | *Met* – no atomic ref‑counts.                                       |
| **Cross‑project subset reuse**               | *Met* – base class never changes, plugins add only their own types. |
| **RTTI support for **``                      | *Met* – status quo.                                                 |

Changing `get_Identifier` now would ripple through every call‑site without delivering an immediate performance win – cloning cost on identifiers is not a measured bottleneck.

---

## 3  Observed Pain‑Points (not urgent)

1. **Identifier duplication** – every lookup returns a fresh copy, wasting a small amount of heap.
2. **Branch complexity** – callers sometimes juggle `unique_ptr` vs borrowed pointer semantics when literals are built in place.

These are inconveniences, *not* correctness or performance blockers today.

---

## 4  Recommended Future Improvement (NodePtr borrow)

When the AST stabilises and refactors calm down, switch only the *lookup* path (`Compiler::get_Identifier`, potential literal‑arena helper) to:

```cpp
using NodePtr = base_typed_expression<Lexer,Compiler> const*;
Maybe_error<NodePtr> get_Identifier(StringView name) const;  // borrow, no copy
```

- **Borrowed lifetime** – the compiler’s symbol table (or a small literal arena) owns the nodes.
- **Zero allocations per lookup**, no ref‑count, still no central registry.
- **Migration path** – callers clone only when they truly need ownership:
  ```cpp
  auto raw = get_Identifier(name).value();
  auto own = std::make_unique<Node>(*raw);  // explicit copy if needed
  ```

This change is **postponed** until profiling proves the duplication cost noticeable or the codebase stops heavy churn.

---

## 5  Re‑evaluation Triggers

Revisit the interface when **all** of these are true:

1. DSL node types stabilise (few new `typed_expression<T>` per quarter).
2. Profiling shows ≥10 % wall‑time or ≥5 % heap in identifier cloning.
3. We introduce an immutable node pool (literal & identifier arena) inside `Compiler`.

Until then, keep the current `unique_ptr` contract to maximise development velocity and minimise merge friction.

*Last updated: 2025‑07‑24*

