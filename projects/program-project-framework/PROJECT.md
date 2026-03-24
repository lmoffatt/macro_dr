# Project 1 - Program/Project Framework

## Identity

- Project id: `program-project-framework`
- Program: `macroir`
- Type: governance / architecture / documentation
- Status: drafted
- Publication relevance: optional
- Primary output type: framework, documentation, schemas

## Purpose

Define a practical, stable, and extensible framework for organizing the
MacroIR program and its repository.

This project establishes the vocabulary, structure, and validation model needed
to distinguish:

- the long-lived program
- the current baseline status of theory + code + canonical documentation
- bounded projects
- papers and releases
- runs, results, and artifacts

Its goal is to reduce ambiguity, make progress assessable, and allow the
repository to grow coherently.

## Scope

### In scope

- Define the core vocabulary:
  - `program`
  - `baseline`
  - `project`
  - `paper`
  - `release`
  - `run`
  - `artifact`
  - `archive`
- Define the role of the repository as a container for a program, not
  necessarily a single project.
- Define the relationship between:
  - program-level organization
  - theory
  - code
  - execution projects
  - papers
- Define minimal machine-readable manifests:
  - `PROGRAM.yaml`
  - `PROJECT.yaml`
- Define project and program stages, checks, and readiness concepts.
- Define a pragmatic folder organization for `macro_dr`.
- Define rules for when a `README.md` is enough and when a folder should be
  split into multiple documents.

### Out of scope

- Full implementation of the repository reorganization.
- Moving files or editing the existing repo layout.
- Designing a rigid ontology for every future use case.
- Building automation beyond defining what it should validate.

## Key problem

The current repository mixes several kinds of work:

- evolving theory
- evolving code
- execution-specific projects
- paper production
- archival material

Without a stable vocabulary and a clear organizational theory, these layers are
easy to confuse. Exploratory theory can be mistaken for canonical software
truth, project artifacts can be mistaken for program baseline, and papers can
be mistaken for the source of truth.

## Conceptual core

The central idea is:

> MacroIR is a program, not just a project.

Within that program:

- the current accepted state of theory + code + canonical docs is the baseline
- bounded efforts that test, extend, validate, or apply the baseline are
  projects
- publications are papers, which are one possible class of deliverable, but
  not the only one

This framework should support organic growth: start small, remain usable when
incomplete, and allow new features, workstreams, and checks to be added later
without redefining everything.

## Deliverables

### D1. Vocabulary document

A short canonical document defining:

- program
- baseline
- project
- paper
- release
- run
- artifact
- archive

This document is the single source of truth for the meaning of these terms.

### D2. Repository organization model

A documented proposal for organizing the repository into major layers, likely
along these lines:

- `program/`
- `theory/`
- `code/`
- `projects/`
- `papers/`
- `releases/` 

### D3. Program/project lifecycle model

A draft of:

- project stages
- program health/coherence dimensions
- validation and readiness concepts

### D4. `PROJECT.yaml` draft schema

A minimal machine-readable schema describing one bounded project, including:

- identity
- purpose
- scope
- deliverables
- status/stage
- checks
- success criteria
- publication relevance
- exit criteria

The schema terms must be consistent with D1 and should not introduce a
parallel vocabulary.

### D5. `PROGRAM.yaml` draft schema

A minimal machine-readable schema describing the MacroIR program baseline,
including:

- identity
- purpose
- scope
- baseline
- workstreams
- active projects
- checks
- readiness / coherence concepts

The schema terms must be consistent with D1 and should not introduce a
parallel vocabulary.

### D6. Document splitting rule

A simple rule of thumb for:

- when `README.md` is enough
- when to split into additional documents
- when a concept should become a folder with an index `README.md`

## Success criteria

This project is successful if:

1. The vocabulary is clear enough to use consistently across repo discussions.
2. `PROGRAM.yaml` and `PROJECT.yaml` are useful and justified by checks, not
   merely descriptive.
3. The distinction between baseline and project is operationally clear.
4. The resulting structure is simple enough to use immediately, not just
   theoretically elegant.
5. The framework supports non-paper outcomes, including code improvements,
   validations, benchmarks, technical assessments, demonstrations, reusable
   modules, and negative results.
6. The framework is extensible without needing a redesign when the program
   evolves.

## Non-paper deliverables explicitly allowed

Projects in this framework may end in any of the following:

- paper
- code
- benchmark
- validation report
- technical assessment
- demonstration
- reusable module
- negative result
- decision to stop or archive a line of work

## Current status

### Completed

- Initial distinction between program and project explored.
- MacroIR identified as a program, not merely a project.
- Need for a distinct term for the current evolving accepted state identified.
- `baseline` currently preferred for that role.
- Need for both `PROGRAM.yaml` and `PROJECT.yaml` established.
- Need for checks at both project and program level established.
- Need to separate theory, code, projects, papers, and archival material
  identified as a core architectural principle.

### In progress

- Finalizing vocabulary.
- Finalizing what belongs in `PROGRAM.yaml` vs `PROJECT.yaml`.
- Finalizing the practical repository layout.

### Not yet done

- Written canonical vocabulary document.
- Written schemas.
- Written stage/check model.
- Applied framework to actual `macro_dr` reorganization.

## Validation and checks

This project should be considered validated when the following checks pass.

### Conceptual checks

- The vocabulary distinguishes all major entities without overlap.
- `program`, `baseline`, and `project` are not interchangeable.
- Papers are correctly modeled as optional deliverables rather than the sole
  target of projects.

### Structural checks

- The proposed directory model can classify all current major materials in
  `macro_dr`.
- The model gives each artifact class a clear home.
- The model distinguishes canonical material from exploratory and archival
  material.

### Operational checks

- `PROJECT.yaml` can describe at least this project and one execution project
  without awkwardness.
- `PROGRAM.yaml` can describe the MacroIR program in a way that includes
  meaningful checks.
- The framework is simple enough to start using before the repository is fully
  reorganized.

### Extensibility checks

- New projects can be added without changing the vocabulary.
- New checks can be added without breaking the model.
- New workstreams can be introduced without collapsing the distinction between
  program and project.

## Exit criteria

This project is complete when all of the following exist in draft form:

1. A canonical vocabulary document.
2. A first draft `PROJECT.yaml` template.
3. A first draft `PROGRAM.yaml` template.
4. A short lifecycle/check model for both project and program.
5. A practical proposed top-level directory structure for `macro_dr`.
6. Enough decision stability to begin Project 2: repository reorganization.

## Dependencies

### Upstream dependencies

- None strict; this project is foundational.

### Downstream dependencies

- `macro-dr-reorganization`
- future project status automation
- future program/project dashboards
- future documentation governance work

## Risks

- Overdesign: creating a framework too abstract to use.
- Vocabulary inflation: too many terms for the same concept.
- Bureaucracy: requiring too much metadata too early.
- Premature formalization: trying to define all future use cases now.
- Confusing repository structure with conceptual structure.

## Guiding constraints

- Be pragmatic before elegant.
- Prefer a small stable kernel over a complete theory.
- Every YAML manifest must justify itself through checks or coordination value.
- The framework must help real work, not merely describe it.
- The system must support organic growth.

## Next actions

1. Write the canonical vocabulary document.
2. Draft `PROJECT.yaml` and `PROJECT_DEFINITIONS.md`. 
3. Draft `PROGRAM.yaml` and `PROGRAM_DEFINITIONS.md`.
4. Define project stages and program health/coherence dimensions.
5. Write the proposed top-level directory structure for `macro_dr`.
6. Start Project 2 using this framework, then refine the framework from real
   friction.

## Summary sentence

This project defines the governance architecture of the MacroIR program: its
vocabulary, baseline concept, project model, checks, and repository
organization rules.
