# ADR-002: Lingua Franca & Contracts

Status: Draft
Date: 2025-09-13

## Context

Core must serve two independent clients: CMD (command surface/DSL) and IO (persistence/interchange). We need a shared data representation to communicate structure and intent without coupling modules or sacrificing performance in inner loops.

## Decision

- Canonical data model: a JSON-like Value + Schema used at module boundaries.
  - Value: null, bool, number, string, array, object; trivially JSON-encodable.
  - Schema: fields, types, units, bounds, required/optional.
- Contracts: pre/postcondition predicates checked at boundaries; integrate with tests and optional DSL predicates.
- Vector view (optional): standardized numeric projection (dimension + labels metadata) as a documented facet of the same schema for hot numeric paths.
- Placement: the data/contract layer sits below all modules; DSL is a client, not the provider.

## Rationale

- Decoupling: boundaries speak a single format; internals remain typed and optimized.
- Inspectability: easy to log, diff, and serialize; good for debugging and tests.
- Correctness: contracts enforce invariants at the surface (HoTT path) without infecting internals.
- Evolution: schemas are versionable; unknown fields can be preserved; easy adapter swaps.

## Consequences

- Boundary-only conversions: keep Value at edges; avoid using it in inner loops.
- Fast-path escape: allow direct typed/vector APIs where performance demands it.
- Adapter discipline: `to_value/from_value` and `to_vector/from_vector` are free functions near types; swapping encodings later is localized.

## Bypass Rules

- Use typed/vector APIs inside Core↔Models↔Inference.
- Reserve Value + Schema for CMD and IO boundaries.
- Use zero-copy buffers or borrowed views for large arrays.

## Open Questions

- Numerics: keep int/double distinct vs normalize to double.
- Units: embed in Schema vs enforce via invariants.
- Validation policy: always on vs debug/configurable per command.
- Vectorization scope: which entities must be vectorizable (parameters, states)?
- Naming: “Domain Entities” vs “Domain Model”.

