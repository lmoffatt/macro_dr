# MacroDR Architecture & Module Organization

Purpose: Decouple modules, stabilize public APIs, and enable correctness-by-construction while preserving high-performance inner loops.

Core Decision: Introduce a single, implementation-agnostic “lingua franca” (canonical data model) used at module boundaries, with an optional numeric projection for hot paths.

## Principles

- Single lingua: one canonical data model for boundaries; internal representations are free.
- Boundary contracts: commands and IO validate pre/postconditions and schemas at edges.
- One-way dependencies: higher layers depend on lower layers; no cycles.
- Fast paths: inner loops use typed/vector APIs; adapt to lingua franca only at edges.

## Modules

DSL
- Responsibilities: parser, typed/untyped AST, function registry/wrappers. No domain logic.
- Public APIs: `dsl::Compiler`, typed/predicate wrappers, identifier/types utilities.
- Depends on: STL only.

CLI
- Responsibilities: entrypoint, script ingestion, registry aggregation only.
- Public APIs: `make_compiler()` aggregator.
- Depends on: DSL, CMD.

CMD (Command Surface)
- Responsibilities: definitive command set (args, pre/post), light orchestration, registration with DSL.
- Public APIs: per-domain headers exposing commands + registry builders.
- Depends on: Core, Domain Entities, IO, DSL (registration only).

Core (Implementations)
- Responsibilities: final implementations behind CMD; orchestrates Models, Inference, IO, Math.
- Public APIs: none (internal to library); invoked by CMD.
- Depends on: Models, Inference, IO, Probability, Math, Utils, Domain Entities.

Domain Entities
- Responsibilities: canonical data/value types independent of specific models (e.g., `Experiment`, `Recording`, `ParameterSet`, `LikelihoodSpec`, `SimulationAlgo`).
- Public APIs: strongly typed objects + builders; no IO or sampling.
- Depends on: Utils.

IO
- Responsibilities: serialization/deserialization (CSV/JSON/etc.), adapters to/from Domain Entities; no business logic.
- Public APIs: load/save functions; DTO mappers.
- Depends on: Domain Entities, Utils.

Models
- Responsibilities: biophysical models, names/parameters, deterministic forward simulations; model-level invariants.
- Public APIs: `get_model`, traits (`names()`, `model_name()`).
- Depends on: Math, Utils.

Inference
- Responsibilities: priors/likelihood assembly, posteriors, samplers/schedulers (MCMC/PT/dts).
- Public APIs: configuration/builders for likelihood/prior; sampling and evaluation entrypoints.
- Depends on: Models, Probability, Math, Utils, Domain Entities.

Probability
- Responsibilities: distribution primitives, RNG, density/log-likelihood components (not samplers).
- Public APIs: distributions, random samplers, density helpers.
- Depends on: Math, Utils.

Math
- Responsibilities: linear algebra, calculus/AD, numerical integration; system libraries.
- Public APIs: matrix ops, derivatives, special functions, BLAS/LAPACK/GSL wrappers.
- Depends on: system libs, Utils.

Utils
- Responsibilities: error/result (`Maybe_error`), semantics/variables, metaprogramming (fold/continuation/indexed), memoization, general helpers.
- Depends on: STL only.

Tests
- Responsibilities: HoTT/contract tests (postconditions/invariants), standard unit/integration tests, CLI smoke scripts.
- Depends on: public surfaces (CLI/CMD/Core/Domain Entities).

## Lingua Franca (Canonical Data Model)

- Data Value: JSON-like value tree (null, bool, number, string, array, object), trivially JSON-encodable. Sits below all modules; DSL is a client, not the owner.
- Schema: minimal schema documenting fields, types, units, bounds, required/optional.
- Contracts: pre/postcondition predicates checked at boundaries; integrates with tests and optional DSL predicates.
- Vector View (optional): standard numeric projection (dimension + labels metadata) for performance-critical flows.

## Allowed Dependencies

- CLI → DSL, CMD
- CMD → Core, Domain Entities, IO, DSL (registration only)
- Core → Models, Inference, IO, Probability, Math, Utils, Domain Entities
- Inference → Models, Probability, Math, Utils, Domain Entities
- Models → Math, Utils
- Probability → Math, Utils
- IO → Domain Entities, Utils
- DSL → data (adapters only)
- Utils, Domain Entities → STL only

## Folder/Namespace Layout

- `include/macrodr/dsl/` (`macrodr::dsl`): lexer, grammar, compiler, function wrappers
- `include/macrodr/cli/` (`macrodr::cli`): registry aggregator
- `include/macrodr/cmd/` (`macrodr::cmd`): command headers + registry builders by domain
- `include/macrodr/data/` (`macrodr::data`): lingua franca (`value.h`, `schema.h`, `encode.h`, `vectorize.h`, optional `json_adapter.h`)
- `include/macrodr/spec/` (`macrodr::spec`): contracts (`spec.h` for pre/post/invariants)
- `include/macrodr/interface/`: stable public interfaces
- `include/macrodr/io/`: IO interfaces
- `include/macrodr/core/`: promoted common helpers/types (as needed)
- `src/cli/`: `main.cpp` and small CLI helpers
- `src/core/`: implementations behind CMD
- `legacy/`: existing heavy subsystems (kept compiling, wrapped by Core/CMD)
- `docs/`: architecture docs and ADRs

## Related Architecture Notes

- `docs/architecture/macro_state_types.md` — rationale for `Macro_State`, `dMacro_State`, and `ddMacro_State`.
- `docs/architecture/moment_statistics.md` — moment/variance tracking helpers and usage notes.

## Interactions (Examples)

Load Experiment
- CLI script → CMD `load_experiment` (pre checks) → Core → IO parsers → Domain Entities (schema validated) → return.

Simulate
- CLI → CMD `simulate` → Core orchestrates model lookup + parameter loading + RNG → uses Models/Inference + Math → IO persists outputs → return.

## Bypass Rules

- Use typed/vector APIs inside Core↔Models↔Inference.
- Use Data Value + Schema at CMD and IO boundaries.
- Use zero-copy buffers for large numeric arrays; avoid copying large data into value trees.

## Migration Plan (Phased)

1) Add `include/macrodr/data/` and `include/macrodr/spec/` (headers only).
2) Centralize command registration under CMD; remove ad-hoc registrations from `src/cli/main.cpp`.
3) Port high-value commands to CMD/Core surfaces: `load_experiment`, `simulate`.
4) Wrap legacy in Core as needed; keep behavior identical.
5) Add CLI smoke tests and initial contract tests (postconditions) for key commands.
6) Document with ADRs: Modules & Dependencies; Lingua Franca & Contracts.

## Open Decisions

- Numerics: keep int/double distinct vs normalize to double.
- Units: embed in Schema vs enforce via invariants in Spec.
- Validation policy: always on vs debug/configurable per command.
- Vectorization scope: which entities must be vectorizable (parameters, states).
- Naming: “Domain Entities” vs “Domain Model”.
