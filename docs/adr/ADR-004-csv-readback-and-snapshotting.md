# ADR-004 — CSV Read-Back and Intermediate Result Snapshotting

Status: Deferred. Design recorded for future revival; no implementation beyond a
small foundational utility ([fingerprint.h](../../include/macrodr/io/fingerprint.h))
that is not yet used by callers.

## Context

MacroDR scripts produce artifacts via `write_csv` at several points in a pipeline
(e.g. [figure_2.macroir](../../projects/eLife_2025/ops/local/figure_2.macroir)):

- `simulations` (upstream) — simulate 1024 traces from a model + params.
- `dlikelihood_predictions` (intermediate) — derivative-aware likelihood state,
  input to downstream diagnostics.
- `likelihood_analysis` (terminal) — bootstrap + cross-correlation diagnostics
  for plots.

A recurring question: should we be able to **read these CSVs back** into the
corresponding C++ objects, to enable either

1. **Caching** — skip expensive recomputation on reruns (the expensive step is
   `calc_dlikelihood_predictions`), or
2. **Snapshotting** — resume a pipeline after a late-stage failure without
   redoing earlier stages.

Today there is no read_csv. Scripts rerun end-to-end.

## The design we worked out before deferring

We explored the full shape of a read/cache layer before deciding not to build it.
The design is captured here so it can be revived without re-deriving it.

### CSV schema

Already in place (see [write_csv_common.h](../../include/macrodr/cmd/detail/write_csv_common.h)).
Long format, one value per row:

```
scope, simulation_index, sample_index, segment_index, sub_index,
n_step, step_start, step_end, step_middle, agonist, patch_current,
component_path, value_row, value_col,
probit, calculus, statistic, quantile_level,
param_index, param_col, param_name,
[axis columns…],
value
```

`scope` ∈ {`simulation`, `state`, `evolution`, `summary`}. `calculus` ∈
{`primitive`, `derivative`}. Self-descriptive enough that any templated structure
built from primitives + derivatives + `Evolution` sub-records is writable by the
same emission machinery.

### Five primitives + one wrapper

| Op | Reads sidecar? | Computes? | Writes? |
|---|---|---|---|
| `calc_*(inputs…) → Target` | — | yes | no |
| `write_csv(target, ctx…, path)` | — | no | CSV (+ sidecar) |
| `read_csv<Target>(path, metadata…)` | **strict** | no | no |
| `read_csv_unchecked<Target>(path, metadata…)` | **skip** | no | no |
| `write_if_absent(target, ctx…, path)` | checks | no | only if missing/partial |
| `load_or_calc(path, inputs…)` *(wrapper)* | strict | on miss | no |

### Reader signature

Factory form, not out-param:

```cpp
template <class Target>
Maybe_error<Target> read_csv(std::string path, /* non-reconstructable metadata… */);
```

Write arguments split on read into two roles:

- **Reconstructable from rows** (recording values, derivative entries, experiment
  scaffolding) → comes from the CSV.
- **Non-reconstructable metadata** (template target, `Parameters_transformed`
  domain) → passed as const input.

`Indexed<T>` and `vector<T>` are dispatched at runtime from the coord columns —
one reader per leaf type, not 4× per-wrapper overloads.

### Fingerprint vs naming discipline

Two ways to detect a stale cache:

1. **Content fingerprint.** Hash the inputs that produced the CSV at write time,
   record the hash in a sidecar JSON, re-hash current inputs on read, compare.
2. **Naming discipline.** User is responsible for giving different outputs
   different filenames. Same name ⇒ same inputs, by author contract.

We chose **naming discipline** (see "Decision" below). Fingerprinting would have
required one of:

- In-object fingerprint (Option C): extend `Simulated_Recording<SimTag>` and
  `dMacro_State<Vars…>` with a fingerprint member, plumbed through every
  constructor and serialization site. Heavy type surgery because `dMacro_State`
  is a `Vector_Space` of variadic `please_include` tags.
- DSL-level fingerprint (Option E): add `Environment::fingerprint_for(id)` and a
  fingerprint function per value-producing command (~15 commands for figure_2).
  No type changes, but still substantial DSL surgery in
  [grammar_typed.h](../../include/macrodr/dsl/grammar_typed.h).

Even with naming discipline, the sidecar JSON is still useful as a
**transaction-commit marker**: `write_csv` writes `.csv.tmp`, writes the sidecar,
then renames `.csv.tmp` → `.csv`. Presence of matching `.csv` + `.meta.json` ⇒
the writer ran to completion and the schema is known. Prevents `load_or_calc`
from serving a truncated CSV left behind by a crashed run.

### Which artifacts get a reader

- **Upstream/intermediate artifacts** (feeds further C++ computation): reader
  justified. Here: `simulations`, `dlikelihood_predictions`.
- **Terminal artifacts** (plot input, human diagnostics): no reader. Write-only.
  Here: `likelihood_analysis`.

## The alternatives and why CSV read-back is the wrong tool for debugging

The user scenario that *seemed* to justify this work was "reproduce errors from
overnight runs without re-running expensive stages". CSV round-trip does not
fit that scenario well:

- CSV is a **lossy human-readable view**. `dMacro_State_Ev_gradient_all`
  contains `Parameters_transformed` with transform metadata
  (log10/standard/etc., bounds, names) that the CSV only captures by name. A
  reconstructed object is an approximation, not the same object the failing
  downstream stage would have received.
- Round-tripping templated types requires either pre-declaring the exact template
  instantiation at the call site or growing runtime reflection. Either is real
  engineering.
- The *full-fidelity* snapshotting path the scenario actually calls for is
  **extending `--env-save` / [environment_io.h](../../include/macrodr/io/json/environment_io.h)**
  to cover `Simulated_Recording` and `dMacro_State` via `to_json`/`from_json`.
  Reuses existing infrastructure. Schema drift is localised to one pair of
  functions per type instead of a CSV parser per type.

## Decision

**Path A — do nothing.**

For figure_2 scale (~8 hour overnight compute, late-stage bugs rare), rerunning
from scratch is the right answer. Running end-to-end also maximises confidence
that no silent serialization bug is affecting results.

Implications:

- No CSV readers.
- No cache wrappers.
- `write_csv` keeps its current shape; no sidecar emission from callers today.
- For debugging: users make a small-scale version of the script (fewer channels,
  fewer simulations) that reproduces the bug in seconds.

## What's already in tree

- [include/macrodr/io/fingerprint.h](../../include/macrodr/io/fingerprint.h) —
  `Fingerprint`, `Hasher`, and `Sidecar` read/write built on the existing
  [minijson.h](../../include/macrodr/io/json/minijson.h). Header-only, no users.
- [tests/io/test_fingerprint.cpp](../../tests/io/test_fingerprint.cpp) — hash
  stability, sidecar round-trip, hex round-trip. Passes on `[fingerprint]`
  filter.

These were built before the decision to defer. They do no harm sitting unused;
they are the foundation a future revival would start from. No other changes were
made to `write_csv`/`read_csv` surfaces.

## If you come back to this

**Revival trigger — speed caching:** compute time of an intermediate step grows
past the point where "rerun overnight" is acceptable (~> 1 day), or iteration
on downstream code accelerates to the point where multiple reruns per day
becomes the norm.

**Revival trigger — snapshotting:** the specific scenario "I have an expensive
intermediate result on disk and want to feed it to a C++ function that was
updated after the compute ran." At that point the right answer is probably
extending `env-save` (full-fidelity JSON), not CSV readers.

**Steps in the order they should be tackled:**

1. Pilot the reader on `Simulated_Recording<please_include<>>` — simplest case,
   CSV is fully self-describing (every row has timing, agonist, patch_current).
   Skip the sidecar entirely in v1 (unchecked mode by default) to keep the
   design pressure on the reader, not the cache semantics.
2. Add atomic-write + sidecar emission to `write_csv`. Define the sidecar
   `scheme` field ("simulate:v1", "dlikelihood_predictions:v1"). Strict read
   mode checks it; unchecked mode skips.
3. Build a minimal `load_or_calc(path, compute_fn)` wrapper at C++ level. DSL
   exposure can wait.
4. Tackle `dMacro_State_Ev_gradient_all` read only after step 1 has proven the
   reader schema. The derivative columns share the same row structure as
   simulation values; the hard part is reconstructing the template shape and
   `Parameters_transformed` metadata (which must be passed in as a const input).
5. Decide then whether to also do fingerprinting, or stay with naming
   discipline. Naming discipline is the 80/20 answer until a real stale-cache
   incident happens.

**What not to do:**

- Do not put fingerprints on the computed types themselves. The surgery is
  disproportionate for a property that only matters at the write/read boundary.
- Do not build a reader for terminal artifacts ("completeness" is not a
  justification — maintenance cost outweighs it).
- Do not introduce a DSL metalanguage to auto-generate cached variants until
  you have at least three cross-cutting concerns needing the same pattern. Until
  then, hand-written thin wrappers are cheaper than the metalanguage investment.

## Consequences

- Any CSV already on disk stays usable by R/plotting code unchanged.
- Scripts that fail late cost a full rerun. Mitigation is out-of-band: build a
  small-scale variant of the script when iterating.
- The deferred design is recoverable: resuming from this ADR plus
  [fingerprint.h](../../include/macrodr/io/fingerprint.h) lands the first
  reader in days, not weeks. The conceptual work is not lost.
