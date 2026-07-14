# MacroIR Theory Docs

Curated baseline-facing MacroIR theory is promoted here incrementally.

Current promoted material includes:

- `Macro_IR/`
  - core MacroIR specification, derivations, tutorials, supplement drafts, and
    implementation-facing documents for the implemented MacroIR line
- `Likelihood_Information_Distortion/`
  - core Likelihood Information Distortion theory documents and the current
    implementation note for subspace handling

This area should contain stable conceptual and mathematical baseline material,
plus implementation-facing documents for methods that are actually realized in
the code. It should not contain superseded draft predecessors or generated
theory artifacts.

## Status vs the current eLife paper (2026-07)

The current paper is **likelihood-only, two-state, non-stationary** (see `papers/macroir-elife-2025/02_decision_log.md`). For agents working on it:

- **Canonical for the current paper:** `Macro_IR/`, `Likelihood_Information_Distortion/`, `Gaussian_Fisher_Distortion_Family.md`.
- **Cut from the current paper, kept for future components:** `Posterior_Information_Distortion/` (the posterior framework was dropped; the likelihood-side evidence correction stays as motivation only), `Macro_IRT/`, `Macro_MRT/`. Do not treat these as the current method.
