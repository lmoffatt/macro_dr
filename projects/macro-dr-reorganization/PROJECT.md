# Project 2 - MacroDR Reorganization

## Identity

- Project id: `macro-dr-reorganization`
- Program: `macroir`
- Type: governance / architecture / repository migration
- Status: active
- Publication relevance: irrelevant
- Primary output type: reorganization plan and repository structure

## Purpose

Reorganize the `macro_dr` repository so it reflects the MacroIR program
framework defined in Project 1.

This project translates the abstract framework into a concrete repository
structure that clearly separates:

- program-level identity and governance
- baseline theory
- baseline code
- bounded projects
- papers
- archive material

Its goal is to make the repository easier to navigate, easier to evolve, and
less likely to confuse exploratory material with canonical baseline material.

## Scope

### In scope

- Inventory current major materials in the repository.
- Classify them by role:
  - program
  - theory
  - code
  - project
  - paper
  - release
  - archive
- Define the target top-level structure for `macro_dr`.
- Define migration rules from the current layout to the target layout.
- Identify which existing documents are:
  - canonical
  - active draft
  - duplicated
  - superseded
  - archival
- Seed the first program/project documents needed for the new structure.

### Out of scope

- Full cleanup of every legacy or generated file in one step.
- Rewriting large bodies of theory or manuscript content.
- Changing scientific content during the organizational migration.
- Tooling for automatic migration beyond what is needed immediately.

## Key problem

The current repository already contains:

- software implementation
- architecture and spec documentation
- theory drafts
- execution projects
- paper planning material
- archival notes and transcripts

These materials coexist in ways that make it hard to distinguish current
baseline truth from historical or exploratory content.

## Conceptual core

This project assumes the framework from Project 1:

- MacroIR / MacroDR is a program.
- The current accepted code + theory + canonical docs form the baseline.
- Projects are bounded efforts that act on the baseline.
- Papers are publication-facing deliverables derived from projects.
- Archive material should remain available without being mistaken for baseline.

The reorganization should be incremental, not a one-step big-bang rewrite.

## Deliverables

### D1. Repository classification document

A document that classifies current major directories and document groups by the
framework vocabulary.

### D2. Target directory model for `macro_dr`

A concrete repository structure specifying where program, theory, code,
projects, papers, and archive materials should live.

### D3. Migration plan

A phased plan for moving from the current structure to the target structure.

### D4. Canonical-vs-archive mapping

A list or document identifying which current documents are:

- canonical
- active draft
- duplicate
- former version
- archive

### D5. Initial moved/created scaffolding

The first structural scaffolding needed for the new model, including program
and project definitions.

## Success criteria

This project is successful if:

1. Every major material class in the repo has a clear destination.
2. Active and archival materials are clearly distinguished.
3. Code-facing documentation is separated from exploratory theory.
4. Paper material is separated from execution projects.
5. The migration can be performed incrementally without losing traceability.
6. The resulting structure is simple enough for daily use by you and by agents.

## Current status

### Completed

- Project 1 defined the initial vocabulary, lifecycle model, and organization
  principles.
- The repository has already been surveyed at a high level.
- Repository classification document written.
- Repo-specific target structure written.
- Migration plan written.
- Canonical-vs-archive mapping written.
- First top-level migration scaffolding created:
  - `README.md`
  - `program/`
  - `theory/`
  - `papers/`
- First conservative promotion step applied:
  - scientific-software theory copied into `theory/scientific-software/docs/`
  - first MacroIR theory families classified into `theory/macroir/docs/`,
    `theory/macroir/notes/`, and `theory/macroir/archive/`
  - active eLife paper pack copied into `papers/macroir-elife-2025/`
  - historical behind-the-paper material copied into `papers/archive/`
- Former `docs/eLife 2025/` material migrated into the active paper-local
  split:
  - manuscript draft sources in `papers/macroir-elife-2025/docs/manuscript-drafts/`
  - compiled/history residue in `papers/macroir-elife-2025/archive/version-inicial/`
- Former `docs/theoretical results/` material migrated into `theory/macroir/`
- Former `docs/audios/` material migrated into `program/source-notes/`
- Stale legacy docs archived under `program/archive/legacy-docs/`

### In progress

- Refining the first promoted contents for `program/`, `theory/`, and `papers/`.
- Recording accepted reorganization decisions and pending curation choices.
- Further classifying the promoted MacroIR theory families.
- Cleaning the remaining overloaded `docs/` surface.
- Reclassifying source-note and archive material at program level.
- Closing out the remaining support-folder and archive decisions.

### Not yet done

- Promoted the first canonical program documents into `program/`.
- Reduced the remaining mixed role of `docs/theoretical results/` and
  `docs/papers/` at the source side.
- Remaining overload in `docs/` is now concentrated in archive/history material
  and support/reference folders rather than theory/paper source trees.
- Reduced the mixed role of `docs/` substantially.

## Validation and checks

This project should be considered validated when the following checks pass.

### Structural checks

- The current repository can be classified without major leftover ambiguity.
- The target structure covers all major material types.
- The migration plan does not require destructive cleanup as a first step.

### Operational checks

- The new structure remains navigable for code work.
- Existing build and test workflows are not conceptually broken by the new
  model.
- Existing project directories can be mapped into the framework.

### Governance checks

- Program, project, and baseline vocabulary are used consistently.
- The reorganization does not collapse project, paper, and archive into one
  layer again.

## Exit criteria

This project is complete when all of the following exist:

1. A repository classification document.
2. A target directory model specific to `macro_dr`.
3. A phased migration plan.
4. A canonical-vs-archive mapping for the main documentation areas.
5. Enough of the repo has been reorganized that the framework is real, not only
   theoretical.

## Dependencies

### Upstream dependencies

- `program-project-framework`

### Downstream dependencies

- future release criteria
- future program health automation
- future repository cleanup and archival projects

## Risks

- Over-migrating too early and losing orientation.
- Treating reorganization as content editing.
- Creating a structure that looks elegant but is not convenient in daily work.
- Mixing generated outputs with curated deliverables.

## Guiding constraints

- Prefer incremental moves over a big-bang rewrite.
- Preserve traceability.
- Distinguish canonical, active draft, and archive.
- Do not break the repository for current work.
- Keep the structure usable by humans and agents.

## Next actions

1. Confirm the remaining support-folder decisions for `docs/`.
2. Decide whether any extra Communications Biology material still needs to be
   promoted into `papers/archive/comms-biol/`.
3. Close out the machine-readable program manifest transition around
   `program/PROGRAM.yaml`.
4. Mark the project as delivered once those final decisions are made.

## First applied moves

The project has now applied an initial low-risk structural step:

1. created a root `README.md`
2. created `program/` scaffolding
3. created `theory/` scaffolding
4. created `papers/` scaffolding
5. kept the code tree in place to avoid unnecessary disruption

This makes the target organization visible without yet moving large bodies of
content.

## First conservative promotions

The project has now applied a first conservative copy-based promotion step:

1. copied the accepted scientific-software baseline documents into
   `theory/scientific-software/docs/`
2. classified and copied the first MacroIR theory families into:
   - `theory/macroir/docs/`
   - `theory/macroir/notes/`
   - `theory/macroir/archive/`
3. copied the active eLife paper pack into `papers/macroir-elife-2025/`
4. copied published behind-the-paper material into `papers/archive/`
5. migrated `docs/eLife 2025/` into the active paper-local `docs/` plus
   `archive/` split
6. migrated the remaining `docs/theoretical results/` families into
   `theory/macroir/notes/`

The early promotions were copy-based for traceability. The later promotions of
`docs/eLife 2025/` and `docs/theoretical results/` were completed as actual
moves once the structure had stabilized.

## Recorded decisions

Project-level accepted decisions and pending curation choices should be
recorded in:

- `DECISIONS.md`
- `REMAINING.md`

## Summary sentence

This project converts the Program/Project Framework into a concrete repository
organization for `macro_dr`.
