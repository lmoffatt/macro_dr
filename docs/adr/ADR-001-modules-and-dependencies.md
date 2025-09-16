# ADR-001: Modules & Dependencies

Status: Draft
Date: 2025-09-13

## Context

MacroDR evolved with significant functionality embedded in legacy headers and CLI registration code. We need a clear module structure with one-way dependencies, a stable public surface, and room for high-performance internals.

## Decision

Adopt the following modules and dependency rules:

- DSL: parser/AST/registry; no domain logic.
- CLI: entrypoint and registry aggregation; no domain logic.
- CMD: definitive command surface (arguments + pre/postconditions) with light orchestration; registers with DSL.
- Core: final implementations behind CMD; orchestrates Models/Inference/IO/Math.
- Domain Entities: canonical data/value types; no IO or algorithms.
- IO: serialization/deserialization; adapters to/from Domain Entities.
- Models: biophysical models and parameterization; deterministic forward simulation.
- Inference: likelihood/prior assembly, posteriors, samplers/schedulers.
- Probability: distributions/RNG/density primitives (not samplers).
- Math: linear algebra, calculus/AD, integration; system libs.
- Utils: error/result, semantics/variables, metaprogramming, memoization.
- Tests: HoTT/contract tests + unit/integration + CLI smoke.

Allowed dependencies (one-way):
- CLI → DSL, CMD
- CMD → Core, Domain Entities, IO, DSL (registration only)
- Core → Models, Inference, IO, Probability, Math, Utils, Domain Entities
- Inference → Models, Probability, Math, Utils, Domain Entities
- Models → Math, Utils
- Probability → Math, Utils
- IO → Domain Entities, Utils
- DSL → data (adapters only)
- Utils, Domain Entities → STL only

## Rationale

- Decoupling: keeps boundaries stable while allowing internal migrations.
- Correctness: enables pre/postcondition checks at the surface (HoTT path) without slowing inner loops.
- Performance: hot paths remain typed and vectorized; boundary conversions are explicit.
- Evolution: domains can add capabilities without cross-cutting edits.

## Consequences

- CLI simplifies to aggregation; CMD becomes the authoritative surface.
- Legacy stays compiling; Core wraps it as needed.
- Tests can target the surface, leaving internals free to change.

## Notes

- See ADR-002 for the lingua franca used at boundaries.

