# Tech debt

Open architectural items where the codebase works correctly but a future
refactor would yield material improvements (compile time, binary size,
maintainability, etc.).

Each note explains: current architecture, why it's suboptimal, the cost of
keeping it, and a "trigger" condition for prioritizing the refactor.

## Items

- [`constexpr_flags_analysis.md`](constexpr_flags_analysis.md) — algorithm
  flags as compile-time `constexpr_Var` template parameters vs runtime
  fields. Variant explosion, compile-time cost, binary size.
