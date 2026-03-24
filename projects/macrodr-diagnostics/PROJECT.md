# MacroDR Diagnostics Project

## Identity

- Project id: `macrodr-diagnostics`
- Program: `macroir`
- Type: code / diagnostics / developer-experience
- Status: active
- Publication relevance: irrelevant
- Primary output type: improved baseline diagnostics and error reporting

## Purpose

Improve MacroDR diagnostics so failures are easier to understand, reproduce,
and fix.

This includes better reporting for parsing, evaluation, execution, and
project-level workflows such as the current `eLife_2025` figure scripts.

## Scope

### In scope

- runtime and parser-facing error reporting improvements
- clearer file, command, and script context in failures
- better diagnostics for failing `.macroir` workflows
- making error output more actionable for humans and agents
- project-level acceptance using current failing cases such as `figure_1`

### Out of scope

- broad semantic redesign of the DSL itself
- unrelated algorithmic changes
- paper writing or figure interpretation

## Key problem

Current diagnostics are not strong enough for active workflows. When a script
fails, the output often does not make clear enough:

- where the failure occurred
- which command or expression failed
- which file or evaluation context triggered it
- what the likely next debugging step is

This slows down both project work and baseline development.

## Conceptual core

This is a baseline-improvement project.

`projects/eLife_2025/` owns the triggering use case, but diagnostics
improvements should be reusable across the whole program.

Part of this work depends on clarifying or extending DSL-facing structures, so
it is closely related to `projects/macrodr-dsl/`.

## Deliverables

### D1. Diagnostics project home

Tracked project definition for diagnostics work.

### D2. Current failure inventory

A list of representative failing cases and missing diagnostic information.

### D3. Improved diagnostic output

Concrete improvements to MacroDR’s reporting for parser/runtime failures.

### D4. Acceptance on real workflows

Verification that at least one current failing workflow, such as `figure_1`,
becomes easier to diagnose.

## Success criteria

This project is successful if:

1. diagnostics work has an explicit project identity
2. representative failures are documented
3. new diagnostics expose more useful context than the current baseline
4. at least one real project workflow benefits directly

## Current status

### Completed

- need for a dedicated diagnostics project identified from active `eLife_2025`
  work

### In progress

- creating the project home and manifest

### Not yet done

- current failure inventory
- concrete code changes
- acceptance on current failing workflows

## Validation and checks

This project should be considered validated when:

- the project docs exist
- representative failure cases are documented
- improved diagnostics are implemented
- one real workflow shows better failure reporting

## Exit criteria

This project is established when all of the following exist:

1. `README.md`
2. `PROJECT.md`
3. `PROJECT.yaml`
4. registration in the program manifest
5. a first failure-inventory step recorded

## Dependencies

### Upstream dependencies

- `eLife_2025`
- `macrodr-dsl`

### Downstream dependencies

- better debugging of current execution projects
- future baseline release readiness

## Summary sentence

This project makes MacroDR failures easier to understand by improving the
baseline diagnostic and error-reporting layer.
