# Project - eLife 2025 Execution

## Identity

- Project id: `eLife_2025`
- Program: `macroir`
- Type: execution / validation / paper-support
- Status: active
- Publication relevance: direct
- Primary output type: validation runs, datasets, figures, and reproducible paper support

## Purpose

Execute the validation and figure-generation workflow that supports the current
MacroIR validation paper.

This project owns the operational path from selected models and experiments to:

- reproducible runs
- analyzed outputs
- figure-ready datasets
- manuscript-support artifacts

## Scope

### In scope

- local and scripted execution of MacroIR figure and validation workflows
- management of model and experiment inputs
- generation of intermediate and final datasets for figures
- run provenance and stored outputs
- figure-generation notebooks and scripts
- paper-support execution work for the active eLife paper pack

### Out of scope

- program-level roadmap governance
- canonical paper planning documents
- the full theory-to-implementation mapping for the whole program
- unrelated execution projects outside the current paper effort

## Key problem

The current MacroIR validation path already has scripts, data, notebooks, and
many recorded runs, but it lacks an explicit project identity and status.

Without that project definition, the roadmap and the paper pack refer to work
that is already happening here, but ownership remains implicit.

## Conceptual core

This project is the execution-side counterpart of the active MacroIR paper
effort.

It does not replace the paper pack under `papers/`, and it does not replace the
program roadmap. Instead, it owns the concrete operational workflow that
produces the validation outputs used by the paper.

## Deliverables

### D1. Validation run corpus

Recorded runs with sufficient provenance to understand how outputs were
generated.

### D2. Figure datasets

Stable datasets used by the current figure notebooks and manuscript drafts.

### D3. Figure-generation workflow

Executable scripts or notebooks that transform run outputs into paper-facing
figures.

### D4. Paper-support execution record

A documented bridge between:

- model/experiment inputs
- runs
- figure datasets
- manuscript-facing outputs

## Success criteria

This project is successful if:

1. the validation and figure workflows are explicitly owned here
2. the main figure-generation scripts are runnable from documented inputs
3. runs and outputs are stored with enough provenance to be interpretable
4. the project clearly supports the active eLife paper without absorbing the
   paper pack itself

## Current status

### Completed

- execution assets already exist under `ops/`, `figures/`, `runs/`, and `data/`
- many recorded runs already exist
- model inputs and figure notebooks are present
- figure-support datasets already exist in the project tree
- the project has now been explicitly identified as the owner of the current
  validation-paper execution path
- root-level figure datasets were normalized into `figures/data/`
- root-level model inputs were normalized into `data/models/`
- debug and scratch residue were moved into `archive/debug/`

### In progress

- consolidating the current execution path into an explicit project definition
- clarifying which roadmap milestones are operationally owned here

### Not yet done

- explicit project-local runbook
- explicit project-local roadmap or workboard, if needed
- project-local mapping from runs to manuscript figures

## Validation and checks

This project should be considered validated when:

- the project docs exist
- the main execution scripts are identified
- run records exist
- figure notebooks and data exist
- the project can be named explicitly as the owner of the operational paper path

## Exit criteria

This project is complete when all of the following exist:

1. explicit project identity documents
2. a documented execution path from inputs to figure outputs
3. clear linkage to the active paper pack
4. enough provenance to rerun or interpret the core figure workflows

## Dependencies

### Upstream dependencies

- `macro-dr-reorganization`
- `macrodr-program-roadmap`
- `papers/macroir-elife-2025/`

### Downstream dependencies

- active MacroIR validation paper
- future figure regeneration
- future validation reruns or expansions

## Risks

- mixing paper planning and execution artifacts in one place
- unclear linkage between runs and manuscript outputs
- a large volume of historical runs without explicit curation

## Next actions

1. identify the main figure and validation workflows currently owned here
2. map the roadmap execution milestones to this project
3. add a project-local runbook or workboard if the workflow needs tighter tracking
4. document the linkage from runs and figure datasets to manuscript outputs
