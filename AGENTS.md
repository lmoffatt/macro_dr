# AGENTS.md — Working With This Repo

Purpose: give coding agents a fast, accurate orientation to the codebase so you can find the right files quickly, follow local conventions, and avoid common pitfalls.

Scope: applies to the whole repository.

## Top‑Level Layout

- `include/macrodr/`
  - `dsl/` — DSL lexer/grammar/compiler + function wrappers. Key: `lexer_typed.h` defines the registry machinery (e.g., `function_compiler`).
  - `cli/` — CLI aggregation only (no domain logic). Key: `cli/command_manager.h` builds a `dsl::Compiler` by pushing/merging commands.
  - `cmd/` — Public command surfaces (headers) callable from scripts; light orchestration and argument types. Examples: `cmd/load_experiment.h`, `cmd/simulate.h`.
  - `interface/` — Stable interfaces (`IModel.h`, `IObject.h`).
  - `io/` — File IO and related interfaces.
- `src/`
  - `cli/` — Entrypoint (`main.cpp`) and small helpers; should not contain registration logic long‑term.
  - `core/` — Implementations behind `cmd/*` (e.g., `load_experiment.cpp`, `simulate.cpp`). These may wrap legacy components.
- `legacy/` — Heavy, legacy headers: models, inference engines, math, utilities, and older CLI registries (e.g., `CLI_*`, `qmodel.h`, `parallel_tempering.h`). Leave in place; wrap via `src/core` or `include/macrodr/cmd`.
- `projects/p2x2/ops/` — Operational scripts for batch/cluster runs.
  - `slurm/` — sbatch wrappers; submit jobs and continuations.
  - `multi_task/M_scheme_N_tasks.sh` — constructs macro_dr invocations per local task using env arrays.
  - `MacroIR/simulation.macroir` — example DSL script(s) used by ops.
- `docs/` — Architecture, ADRs, design notes. New: `docs/architecture/modules.md`, `docs/adr/ADR-001*`, `docs/adr/ADR-002*`. CLI overview lives in `docs/cli.md`.
- `docs/cli.md` — CLI usage and migration notes.
- `docs/testing.md` — How to build and run unit/integration tests.
- `tests/` — Unit/CLI tests (if present); keep new tests focused and fast.
- `CMakeLists.txt`, `CMakePresets.json` — Build config (targets: `macrodr_core`, `macrodr_cli`).

## Key Concepts & Where To Edit

- DSL registry
  - Add a new script‑callable function by registering it into a `dsl::Compiler` via `to_typed_function(...)` (factory currently in `legacy/CLI_macro_dr.h`).
  - Prefer building registries in small domain builders under `include/macrodr/cmd/` and aggregating them in `include/macrodr/cli/command_manager.h`.

- Commands (CMD)
  - Public surface that scripts call; headers live in `include/macrodr/cmd/` and should declare clean argument types and (eventually) pre/postconditions.
  - Implementations live in `src/core/` and may call into `legacy/` code while we migrate.

- CLI Entrypoint
  - `src/cli/main.cpp` constructs a `dsl::Compiler` (currently merges several registries, some from legacy). Goal: keep it thin; move registrations into CMD builders.

## Ops / HPC Workflow (p2x2)

Example: `projects/p2x2/ops/slurm/run_32_scheme_10_DR_nyquist.sh`
- Sets SLURM resources and experiment parameters (schemes, experiment, likelihood, scouts, betas, runtime, continuation count).
- Submits an initial job and a chain of continuations via `slurm/M_scheme_N_tasks.sh`.
- `slurm/M_scheme_N_tasks.sh`
  - Sources cluster profile `${PATH_MACRO}/macro_dr/clusters/${CLUSTER}.sh` (external to this repo), sets threads env, then `srun`s the multi‑task script.
- `multi_task/M_scheme_N_tasks.sh`
  - Selects per‑task variables using `${SLURM_LOCALID}`.
  - Builds an ID and invokes the MacroDR binary with a sequence of script files and inline arguments, e.g.:
    - `${PATH_MACRO}/${PATH_MACRO_DRX}/macro_dr` \
      `${PATH_MACRO}/macro_dr/${SCHEME_DIR}/${SCHEME}.txt` \
      `${PATH_MACRO}/macro_dr/scripts/${EXPERIMENT}.txt` \
      `${PATH_MACRO}/macro_dr/scripts/simulation.txt` \
      `"--runIdName= \"${IDNAME}\""` \
      `"--num_scouts_per_ensemble = get_number(n=${N_SCOUTS})"` \
      `"--max_iter_equilibrium = get_number(n=${MAX_ITER})"` \
      `${PATH_MACRO}/macro_dr/scripts/${LIKELIHOOD}.txt` \
      `${PATH_MACRO}/macro_dr/scripts/beta_${N_BETA}.txt` \
      `${PATH_MACRO}/macro_dr/scripts/evidence_${EVIDENCE_ALGORITHM}_(data|continuation).txt`

Notes:
- These ops scripts rely on environment variables: `PATH_MACRO`, `PATH_MACRO_DRX`, `CLUSTER`, `PARTITION`, etc. Some referenced files (e.g., `macro_dr/scripts/*`) may live outside this repo in deployment environments.

## Build & Run

- Presets: use `CMakePresets.json` to configure and build; `compile_commands.json` symlink points to `build/gcc-debug/`.
- Targets:
  - `macrodr_core` — static library with core logic and legacy headers included.
  - `macrodr_cli` — CLI binary; links `macrodr_core`.
- Typical run (local): pass one or more `.macroir` or `.txt` script files and/or inline commands prefixed with `--`.

## Conventions & Do/Don’t

- Do
  - Keep CLI thin; put registrations in CMD builders.
  - Wrap legacy via Core; don’t move heavy code without need.
  - Use `Maybe_error<>` for error propagation; don’t throw.
  - Follow `.clang-format` and `.clang-tidy` (config in repo).

- Don’t
  - Include `legacy/*` directly from `src/cli/`.
  - Add unrelated refactors in the same change as a feature/fix.
  - Break existing ops scripts without providing a compatibility layer.

## Quick Pointers

- Registering commands: see `include/macrodr/cli/command_manager.h` and `legacy/CLI_macro_dr.h` (factory helpers).
- Implementations: `src/core/load_experiment.cpp`, `src/core/simulate.cpp` are examples that wrap legacy code.
- DSL machinery: `include/macrodr/dsl/lexer_typed.h` (search for `class function_compiler`).
- Docs: `docs/architecture/modules.md`, `docs/adr/ADR-001*`, `docs/adr/ADR-002*`.

## Pitfalls / Historical Quirks

- Namespaces: `make_compiler_new` exists under `macrodr::cli`; avoid duplicates under `macrodr::cmd`.
- Command registrations historically lived in `src/cli/main.cpp`; the plan is to centralize them under CMD.
- Legacy headers are large; add includes judiciously to keep compile times down.

## If You’re Lost

- Search tips: prefer ripgrep (`rg -n "symbol"`) if available; fallback `grep -RIn` across the repo.
- Start from the command surface you need under `include/macrodr/cmd/`, then jump to its implementation in `src/core/` and referenced legacy.
