# MacroIR Theory Docs

Curated baseline-facing MacroIR theory is promoted here incrementally.

Current promoted material includes:

- `Macro_IR/`
  - core MacroIR specification, derivations, tutorials, supplement drafts, and
    implementation-facing documents for the implemented MacroIR line
- `Macro_IRT/`
  - MacroIR + state-dependent (open-channel) noise via Laplace/Taylor in
    boundary-state space; the "T" extension of the IR family
- `Macro_MRT/`
  - MacroR with state-dependent noise via Laplace/Taylor in macroscopic
    occupancy space; the instantaneous-current ("MR") counterpart of MacroIRT,
    used as the controlled baseline against which interval-averaging gains
    in MacroIRT must be measured
- `Information_Distortion_Matrix/`
  - core Information Distortion Matrix theory documents and the current
    implementation note for subspace handling

This area should contain stable conceptual and mathematical baseline material,
plus implementation-facing documents for methods that are actually realized in
the code. It should not contain superseded draft predecessors or generated
theory artifacts.
