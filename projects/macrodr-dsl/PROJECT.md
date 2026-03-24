# MacroDR DSL Project

## Identity

- Project id: `macrodr-dsl`
- Program: `macroir`
- Type: code / language / interface
- Status: active
- Publication relevance: optional
- Primary output type: improved DSL behavior and DSL-facing baseline structure

## Purpose

Develop the MacroDR DSL as an explicit baseline project.

This includes the parts of the language surface that affect:

- script usability
- parser behavior
- command/expression structure
- diagnostics and evaluation context

## Scope

### In scope

- bounded improvements to the DSL needed by active workflows
- DSL-facing structure that improves diagnostics
- language/interface cleanup needed to support reliable scripting
- clarification of how DSL constructs are surfaced in the baseline

### Out of scope

- full replacement of the current language
- broad runtime refactors unrelated to the DSL
- project-specific figure generation logic

## Key problem

Diagnostics quality depends in part on DSL structure. If the language surface
does not carry enough context, error reporting remains weak even if the runtime
tries to improve.

The repository also needs an explicit project home for DSL work so it can be
tracked as a bounded baseline effort rather than as incidental fixes.

## Conceptual core

This project owns bounded DSL work at the baseline level.

It is distinct from `projects/macrodr-diagnostics/`, which owns the reporting
and acceptance side, but the two are tightly coupled.

## Deliverables

### D1. DSL project home

Tracked project definition for DSL work.

### D2. DSL improvement inventory

A shortlist of DSL problems or limitations relevant to current workflows.

### D3. DSL reconstruction document

A current DSL description grounded in the existing docs and the current code
entry points.

### D4. DSL documentation delta

A clear record of what changed, what is stale, and what old docs remain
useful.

### D5. First DSL improvement set

A first bounded set of DSL improvements justified by current project needs.

### D6. Diagnostics-enabling DSL changes

Specific DSL-side changes that make better diagnostics possible or easier.

## Success criteria

This project is successful if:

1. DSL work has an explicit project identity
2. current DSL pain points are documented
3. the current DSL is reconstructed from docs and code
4. DSL drift from older docs is explicit
5. a first bounded DSL improvement set is identified
6. DSL work is clearly connected to diagnostics and scripting quality

## Current status

### Completed

- need for a dedicated DSL project identified from diagnostics and scripting work
- first DSL inventory file created
- current DSL reconstructed from docs and implementation anchors
- DSL documentation drift recorded explicitly

### In progress

- creating the project home and manifest
- seeding the DSL pain-point inventory from active project failures

### Not yet done

- DSL improvement inventory beyond the first seeded case
- first bounded DSL improvement set
- explicit diagnostics-enabling DSL changes

## Validation and checks

This project should be considered validated when:

- the project docs exist
- the DSL improvement inventory exists
- the reconstructed DSL document exists
- the DSL delta document exists
- at least one bounded DSL improvement set is defined
- the DSL-to-diagnostics linkage is explicit

## Exit criteria

This project is established when all of the following exist:

1. `README.md`
2. `PROJECT.md`
3. `PROJECT.yaml`
4. registration in the program manifest
5. a first DSL-improvement inventory step recorded
6. a reconstructed DSL document
7. a DSL changes document

## Dependencies

### Upstream dependencies

- current code baseline

### Downstream dependencies

- `macrodr-diagnostics`
- future script usability improvements

## Summary sentence

This project gives bounded DSL work an explicit home so language improvements
can be tracked as baseline development rather than incidental fixes.
