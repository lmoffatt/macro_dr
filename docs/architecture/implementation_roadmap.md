# MacroDR Implementation Roadmap

This roadmap turns the architecture decisions into incremental, low‑risk changes. It keeps the current behavior intact while moving code into stable boundaries.

Goals
- Achieve the functionality we had before start moving things around (this is done in branch rr_scheme, also located in folder ~/Code/macro_dr_rr_scheme/macro_dr)
- One first change should be that we have Abstract bases for Models and algorithms, we hope to reduce the compiling time in this way. 
- Keep legacy compiling; wrap it via Core until migrated.
- Add minimal contract tests and a CLI smoke test.

Principles
- Make small, reversible changes per step.
- Prefer adapters at boundaries over deep refactors.
- Do not change behavior unless explicitly stated.

Milestones
- M0: Docs finalized (this roadmap + ADR‑001/002).
- M1: Registry consolidation; CLI simplified to aggregator.
- M2: Priority CMD/Core surfaces: load_experiment, simulate.
- M3: Lingua franca pilot for 2 entities + boundary checks.
- M4: Broader CMD migration (likelihood, thermo/dts config) from legacy.
- M5: Cleanup, performance passes, and deprecations.

---

Phase 1 — Registry Consolidation (M1)
Deliverables
- `include/macrodr/cli/command_manager.h` is the only place that aggregates registries.
- Domain registries live in `include/macrodr/cmd/*` (small headers only).
- `src/cli/main.cpp` no longer contains ad‑hoc `push_function` calls.

Tasks
1. Extract all registrations from `src/cli/main.cpp` into per‑domain builders:
   - utilities/io/experiment/model/simulation/likelihood/dts (start with those already present in legacy builders).
2. Align `make_compiler_new` declaration/definition:
   - Keep it in `namespace macrodr::cli` and update/remove `include/macrodr/cmd/init_commands.h` accordingly.
3. Ensure build compiles with registrations coming solely from CMD headers.

Acceptance
- CLI builds and runs existing scripts unchanged.
- No direct `legacy/*` includes in `src/cli/` (only via public headers).

---

Phase 2 — CMD Surfaces + Core Implementations (M2)
Deliverables
- Stable CMD headers for priority commands.
- Core `.cpp` files implement them by orchestrating legacy code.

Targets (priority)
- `cmd::load_experiment` (already present) — verify header/API and complete implementation.
- `cmd::simulate` — verify header/API and complete implementation.

Tasks
1. Define argument types, pre/post stubs in CMD headers (no heavy logic).
2. Implement orchestration in `src/core/*.cpp` using legacy as needed.
3. Add CLI smoke test: a small `.macroir` script that calls load+simulate and produces output.

Acceptance
- Behavior identical to current runs for sample scripts.
- Minimal pre/post stubs compiled (may be no‑ops initially).

---

Phase 3 — Lingua Franca Pilot (M3)
Deliverables
- `include/macrodr/data/` and `include/macrodr/spec/` header‑only modules.
- Encoders for 2 entities (e.g., `Recording`, `Experiment`) and a minimal Schema for each.
- Boundary checks (pre/post) enabled in the CMD surface for these two flows.

Tasks
1. Add headers: `data/value.h`, `data/schema.h`, `data/encode.h`, `data/vectorize.h` (prototypes only).
2. Add `spec/spec.h` with `Spec<T>` and simple predicate helpers.
3. Implement `to_value/from_value` for `Recording` and `Experiment`.
4. Validate in CMD: preconditions on inputs; postconditions on outputs.
5. Optional: a DSL predicate to assert basic invariants in a test script.

Acceptance
- Conversions used only at boundaries (no use in inner loops).
- CLI smoke test extended to assert a basic postcondition.

---

Phase 4 — Broader CMD Migration (M4)
Deliverables
- Migrate additional commands out of legacy registries:
  - Likelihood configuration (`set_Likelihood_algorithm`).
  - Thermo/dts config (`set_ThermoAlgorithm_dts`).
  - Other high‑value commands referenced in `legacy/CLI_function_table.h`, `legacy/CLI_macro_dr_base.h`, `legacy/CLI_likelihood.h`.

Tasks
1. For each command group, define CMD header with argument schema and pre/post stubs.
2. Implement Core orchestration that wraps legacy functionality.
3. Add or extend smoke tests to exercise the new surfaces.

Acceptance
- Existing scripts continue to run with no behavior change.
- Registrations originate from CMD only.

---

Phase 5 — Cleanup & Performance (M5)
Deliverables
- Optional vectorization for parameter‑like entities (`to_vector/from_vector` + metadata).
- Clear separation in CMake targets for data/spec (INTERFACE) and core.
- Deprecation notes for direct legacy inclusions.

Tasks
1. Introduce `macrodr_data` and `macrodr_spec` interface targets; link CLI/Core against them.
2. Add `VectorSpec` for parameters; keep conversions local to boundaries.
3. Audit includes in CLI/CMD for accidental legacy dependencies.
4. Document deprecations and the preferred paths.

Acceptance
- No `legacy/*` includes in CLI/CMD; only Core links legacy.
- Optional vector path available for hot numeric flows.

---

Testing Plan (ongoing)
- CLI smoke: runs a minimal `.macroir` script that loads data and simulates; verifies outputs exist and are parseable.
- Contract tests: 1–2 predicates per priority command (e.g., sampling interval > 0; time series lengths match expectations).
- Unit tests: small Core units around IO adapters and parameter loading.

Risk & Mitigation
- Conversion overhead: keep Value at boundaries only; vectorize hot paths.
- Schema drift: versioned fields and clear ownership in CMD/Core.
- Compile‑time cost: keep new headers minimal; reuse existing legacy where possible.

Backout Plan
- Each phase is independently revertible.
- Legacy pathways remain intact; CMD/Core can be pointed back to legacy registries if needed.

Ownership & Sign‑off
- Architecture docs: maintainer (you) approves before each phase.
- Phase leads: assign per domain (experiment/simulation/likelihood/dts) as needed.

