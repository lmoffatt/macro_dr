# Migration Plan

This document defines a phased migration from the current `macro_dr` layout to
the target program-oriented organization.

## Target structure

The target structure currently assumed by the framework is:

```text
program/
theory/
code/
projects/
papers/
```

The migration is incremental. It should not be attempted as a one-step rewrite.

## Phase 1 - Governance scaffolding

Goal: establish the new organizing layer without breaking existing work.

Actions:

1. Keep `projects/program-project-framework/` as the canonical source for the
   framework.
2. Keep `projects/macro-dr-reorganization/` as the canonical source for the
   migration effort.
3. Add classification, migration, and mapping documents before moving other
   material.

Outcome:

- the framework exists
- the migration effort is tracked
- no content has been lost or prematurely moved

## Phase 2 - Program identity and canonical docs

Goal: define where the program-level baseline identity lives.

Actions:

1. Create `program/` and seed it with:
   - program identity
   - roadmap
   - baseline/release concepts
2. Move or copy only the most clearly canonical governance material first.
3. Leave legacy and mixed areas in place until their destination is explicit.

Outcome:

- program identity no longer depends on scattered docs

## Phase 3 - Separate code baseline from mixed docs

Goal: make the implementation baseline legible as code plus code-facing docs.

Actions:

1. Decide whether `code/` will wrap the current source tree or whether a softer
   mapping document is preferred first.
2. Consolidate canonical code-facing docs under the code-facing documentation
   layer.
3. Keep generated workspaces like `build/` outside the conceptual structure.

Outcome:

- implementation and code docs are treated as part of the baseline

## Phase 4 - Separate papers from docs

Goal: move publication-facing material out of the generic docs surface.

Actions:

1. Promote active paper packs into top-level `papers/`.
2. Separate active paper material from archival manuscript versions.
3. Keep paper-local archives inside each paper folder.

Outcome:

- papers become clearly publication-facing deliverables

## Phase 5 - Separate active theory from theory archive

Goal: prevent exploratory and archival theory from being mistaken for baseline
theory.

Actions:

1. Use the accepted two-family theory structure:
   - `theory/macroir/`
   - `theory/scientific-software/`
2. Within each family, distinguish:
   - `docs/` for curated baseline-facing theory
   - `notes/` for active exploratory work
   - `archive/` for superseded theory
3. Promote the first clearly baseline-facing documents into each theory family.
4. Reclassify large mixed theory areas such as `docs/theoretical results/`
   incrementally inside `theory/macroir/`.
5. Keep draft predecessors, `archive_*` files, and generated theory artifacts
   out of the curated `docs/` surface.

Outcome:

- theory gains explicit family structure and an internal distinction between
  baseline, notes, and archive

## Phase 6 - Clean archive boundaries

Goal: keep history without creating a mirrored second repository tree.

Actions:

1. Use local `archive/` folders for project-, paper-, and theory-specific
   history.
2. Use only a small top-level archive for orphaned or imported historical
   material.
3. Avoid reproducing the full active tree under a global archive.

Outcome:

- archive is local and traceable rather than duplicative

## Migration principles

- prefer incremental moves over a big-bang rewrite
- classify before moving
- preserve traceability
- do not treat generated output as release material
- keep generated and operational directories such as `build/` outside the
  curated program structure
- distinguish canonical, active draft, and archive at every step
- do not mix content editing with structural migration unless necessary

## Immediate next moves

1. Write a canonical-vs-archive mapping for the main docs areas.
2. Define the first version of `program/`.
3. Classify the first MacroIR theory families into `docs/`, `notes/`, and
   `archive`.
4. Identify the first active paper pack to promote into `papers/`.
