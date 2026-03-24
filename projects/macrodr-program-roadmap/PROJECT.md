# Project 3 - MacroDR Program Roadmap

## Identity

- Project id: `macrodr-program-roadmap`
- Program: `macroir`
- Type: governance / planning / roadmap
- Status: active
- Publication relevance: optional
- Primary output type: canonical program roadmap

## Purpose

Define the medium-term roadmap of the MacroIR / MacroDR program.

This project turns the current planning material into a canonical, current, and
program-level roadmap that can guide priorities across theory, code,
validation, execution projects, releases, and papers.

## Scope

### In scope

- identify the main active workstreams of the program
- define priority milestones and sequencing
- identify cross-project dependencies
- identify near-term and medium-term goals
- define roadmap-facing success criteria
- consolidate the current roadmap material into a more canonical form
- document the mapping between theory documents and implemented algorithms
- classify roadmap-relevant documents into canonical vs archive

### Out of scope

- detailed implementation planning for every downstream project
- executing the roadmap itself
- rewriting scientific theory or code as part of roadmap drafting
- replacing local project plans with one central monolith

## Key problem

The repository now has a clearer structure, but the forward program plan is
still dispersed across notes, existing projects, and older roadmap material.

Without a canonical roadmap, it is harder to decide:

- what comes next
- which workstreams are primary
- which projects are prerequisites for others
- what should count as medium-term progress for the program

## Conceptual core

This project assumes:

- MacroIR / MacroDR is a program with multiple simultaneous workstreams
- projects are bounded efforts inside that program
- a roadmap should coordinate, not replace, those projects
- the roadmap should remain pragmatic and revisable rather than pretend to fix
  the entire future in advance

## Deliverables

### D1. Canonical program roadmap

A roadmap document, likely promoted into `program/ROADMAP.md`, describing the
next program phase.

### D2. Milestone structure

A milestone model covering near-term and medium-term goals across:

- theory
- code
- validation
- execution projects
- papers
- releases

### D3. Dependency map

A compact account of which program efforts depend on which others.

### D4. Active-project prioritization

A shortlist of active and next projects and their relative importance.

### D5. Theory-to-implementation mapping

A compact mapping between the main MacroIR and validation-theory documents and
the implemented algorithm families they support.

### D6. Canonical-vs-archive classification for roadmap sources

A classification of the main roadmap-relevant theory and manuscript materials
into:

- canonical
- active draft
- archive / historical support

## Success criteria

This project is successful if:

1. the program has a current roadmap that is concrete enough to guide work
2. the roadmap distinguishes near-term priorities from later ambitions
3. cross-project dependencies are explicit
4. theory, code, validation, papers, and execution projects are all reflected
5. the roadmap is concise enough to remain usable

## Current status

### Completed

- the repository now has a program-oriented structure
- the earlier roadmap note has been preserved at
  `program/notes/Project Roadmap MacroIR.md`
- the need for a canonical roadmap project has been identified
- a first canonical roadmap draft has been written at `program/ROADMAP.md`
- the current roadmap focus has been narrowed to the MacroIR validation path
- the roadmap source-of-truth inputs have been identified as:
  - `papers/macroir-elife-2025/docs/manuscript-drafts/`
  - `theory/macroir/docs/Macro_IR/`
  - `theory/macroir/docs/Information_Distortion_Matrix/`

### In progress

- refining the roadmap milestones and dependencies
- reconciling older roadmap notes with the manuscript and curated theory docs
- defining the theory-to-implementation mapping for the validation path
- classifying roadmap-relevant source documents into canonical vs archive

### Not yet done

- refining milestone groupings into concrete projects and deliverables
- deciding which future directions should become separate projects
- writing the first explicit theory-to-implementation map
- writing the first explicit canonical-vs-archive classification for roadmap sources

## Validation and checks

This project should be considered validated when:

- a roadmap document exists
- it covers the main program workstreams
- it identifies milestones and dependencies
- it is consistent with the current repository organization
- it identifies the main theory-to-implementation mapping for the current
  validation path
- it distinguishes canonical source documents from archive or historical residue

## Exit criteria

This project is complete when all of the following exist:

1. a canonical roadmap document
2. a milestone structure
3. an explicit dependency summary
4. a current prioritization of active and next projects
5. a theory-to-implementation mapping for the current validation path
6. a canonical-vs-archive classification of the roadmap source documents

## Dependencies

### Upstream dependencies

- `program-project-framework`
- `macro-dr-reorganization`

### Downstream dependencies

- future release planning
- project prioritization
- validation planning
- paper sequencing

## Risks

- producing a roadmap that is too broad to guide decisions
- producing a roadmap that is too detailed to remain useful
- confusing roadmap milestones with fixed commitments

## Next actions

1. refine the validation-example milestone into concrete execution targets
2. identify the active projects needed to complete the validation paper
3. document the mapping between theory/manuscript documents and implemented algorithms
4. classify the main roadmap source documents into canonical vs archive
5. decide which future directions belong inside the current roadmap horizon
6. update the roadmap after the next project/prioritization pass
