# Reorganization Decisions

This document records decisions that materially constrain or guide Project 2.

## Accepted decisions

### D1. Target model accepted

The target model in `TARGET.md` is accepted with the following pragmatic
interpretation:

- the code baseline remains physically in:
  - `include/`
  - `src/`
  - `tests/`
  - `tools/`
  - `legacy/`
- the explicit new top-level destinations are:
  - `program/`
  - `theory/`
  - `papers/`

Rationale:

- this preserves build and include-path stability
- it avoids a disruptive `code/` relocation
- it still makes the program structure explicit

### D2. Migration mode

The migration mode for Project 2 is:

- `conservative`

Meaning:

- create destinations first
- classify before moving
- promote a small number of documents at a time
- avoid large folder moves until the mapping is stable

Rationale:

- the repository is active
- the code tree should not be disrupted prematurely
- mixed documentation areas need curation before bulk moves

### D3. Scope of `docs/`

`docs/` should remain the home of canonical code-facing and
architecture/program-facing documentation.

This includes:

- architecture
- ADRs
- specs
- canonical guides such as CLI and testing
- documents that explain how the code relates to or realizes the theory

Theory and paper material should gradually leave `docs/` and be promoted into:

- `theory/`
- `papers/`

Rationale:

- `docs/` is currently overloaded
- theory and paper material need their own identity and archive boundaries

### D4. Scope of `theory/`

`theory/` should include both:

- MacroIR-specific mathematical and algorithmic theory
- scientific software theory that functions as part of the conceptual
  foundation of the program

These should be distinguished internally rather than collapsed into one
undifferentiated body of documents.

Working implication:

- material that defines the conceptual foundation of the program belongs in
  `theory/`
- material that explains current architecture, implementation behavior,
  specifications, or guides belongs in `docs/`

Rationale:

- in this repository, scientific software theory is not merely external
  philosophy; part of it functions as baseline conceptual infrastructure
- separating foundational theory from code-facing documentation reduces the
  overload currently present in `docs/`

### D5. Internal structure of `theory/`

`theory/` should be split into two first-level theory families:

- `theory/macroir/`
- `theory/scientific-software/`

Each family should have its own internal lifecycle structure:

- `docs/`
- `notes/`
- `archive/`

Target shape:

```text
theory/
  macroir/
    README.md
    docs/
    notes/
    archive/
  scientific-software/
    README.md
    docs/
    notes/
    archive/
```

Rationale:

- MacroIR theory and scientific software theory are both part of the program's
  conceptual foundation
- they should remain clearly distinguished
- each family needs its own curated baseline, active notes, and archive

### D6. Information Distortion Matrix is theory

The Information Distortion Matrix line of work belongs to `theory/`, not to
`docs/`.

Working implication:

- theoretical IDM documents should be promoted under `theory/`
- only implementation-facing documents that explain how the code realizes or
  diagnoses IDM-related quantities belong in `docs/`

Rationale:

- importance does not determine whether something belongs in `docs/`
- IDM is part of the mathematical and conceptual framework of the program
- `docs/` should remain the software-facing interpretive layer rather than a
  container for all important writing

### D7. Source notes are program-level material

Raw transcripts and similar exploratory source material should live in a
dedicated program-level source-notes area, not in `docs/`.

Working implication:

- audio transcripts and related free-form note captures belong in
  `program/source-notes/`
- they remain there until analyzed and extracted
- only after analysis should they be promoted into structured program, theory,
  project, or paper documentation, or moved into archive

Rationale:

- these materials are not code-facing documentation
- they are not yet canonical theory or project documents
- they remain valuable as source material while analysis is pending

### D8. Cross-cutting legacy docs archive

Cross-cutting legacy repository documents should be archived under
`program/archive/legacy-docs/`.

Rationale:

- these documents are historically important
- they should no longer appear as current canonical docs
- they do not belong to one paper, one project, or one theory family

### D9. Active paper-local docs vs paper archive

Active paper manuscript sources should live in paper-local `docs/`, not inside
paper-local `archive/`.

Working implication:

- active manuscript draft sources belong under `papers/<name>/docs/`
- planning documents remain at the top of the paper folder
- compiled outputs, generated residue, and superseded paper-local history
  belong under `papers/<name>/archive/`

Rationale:

- manuscript drafts are active working documents, not archive material
- papers need the same internal distinction already established for theory:
  active material vs historical residue

## Pending curation decisions

### P0. New projects and papers to create

Before promoting large amounts of material, decide which new top-level project
and paper folders should exist in the reorganized repository.

This matters because some current material in `docs/` should not simply be
moved by topic. It should first be assigned to the correct owning project or
paper.

Questions to resolve:

- which current efforts deserve explicit project folders
- which publication efforts deserve explicit paper folders
- which current mixed areas are better treated as archive rather than as new
  active projects or papers

Current resolved additions:

#### New projects accepted

- `projects/scientific-software/`
  - rationale: the scientific-software line of work has its own conceptual
    center, document family, and likely future internal split into several
    independent ideas or subprojects

Likely candidates still to evaluate:


#### New projects

- a theory-curation project
- a release-definition or release-readiness project
- a documentation-governance or docs-curation project
- whether Proof-Oriented Design / algebraic CLI should remain one internal line
  of work inside `projects/scientific-software/` or later become an explicit
  subproject within that area

- project-specific cleanup under existing efforts such as `projects/eLife_2025/`
  or `projects/p2x2/`

#### New papers

- promotion of `docs/papers/macroir-elife-2025/` into a top-level paper folder 
- promotion of `docs/behind the paper/` it is a genuine
  publication-facing workstream (published https://communities.springernature.com/manage/posts/300507) 

  

Working rule:

- create a new `projects/<name>/` folder when the effort has its own goal,
  deliverables, checks, and lifecycle
- create a new `papers/<name>/` folder when the material is publication-facing
  and has its own narrative, sources, drafts, or archive

When a first explicit list is chosen, record it here.

### P1. First theory documents to promote

Record here the first theory documents that are considered baseline-facing
enough to move into either theory family.

When a first promotion set is chosen, record it here explicitly.

#### Assigned to `theory/macroir/`

- the theory families currently under `docs/theoretical results/`, including:
  - `Macro_IR/`
  - `Information_Distortion_Matrix/`
  - `Adaptive_MacroR/`
  - `Gmean_ij_gvarij/`
  - `Macro_Taylor/`
  - `Macro_TaylorIR/`
  - `MacroLogit/`

Working rule:

- these families belong conceptually to `theory/macroir/`
- they still require internal classification into:
  - `theory/macroir/docs/`
  - `theory/macroir/notes/`
  - `theory/macroir/archive/`
- draft predecessors, `archive_*` files, and generated theory artifacts should
  not all be promoted into `theory/macroir/docs/`

#### Assigned to `theory/scientific-software/docs/`
#### Assigned to `program/notes/`

- `macro_dr_vision_and_realization.md`
  - treated as preliminary program vision rather than baseline MacroIR theory

#### Assigned to `theory/scientific-software/notes/`

- `docs/Proof_Oriented_Design_Manifesto.md`
- `docs/general_scientific_software_design_theory.md`
- `docs/provenance_types_hott.md`
- `docs/algebraic_CLI.md`
- `docs/contagiously_safe_design.md`
- `docs/composite_safety_domains.md`
- `docs/hott_contagion_principle.md`
- `docs/runtime_vs_semantic_safety_narrative.md`
- `docs/safety_explosion_limit.md`

Working rule:

- these are currently treated as vision, proposal, or not-yet-implemented
  conceptual material
- they should remain in `notes/` until a smaller accepted subset justifies
  promotion into `theory/scientific-software/docs/`




### P2. First paper pack to promote

Record here which active paper pack should be promoted first into `papers/`.

Current likely candidate:

- `docs/papers/macroir-elife-2025/`

## Working rule

Decisions recorded here should be reflected, when relevant, in:

- `PROJECT.md`
- `PROJECT.yaml`
- `TARGET.md`
- `MIGRATION.md`
