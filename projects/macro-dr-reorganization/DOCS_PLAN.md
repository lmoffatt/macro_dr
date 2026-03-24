# Docs Cleanup Plan

This document proposes the next cleanup steps for the remaining `docs/` tree
after the first conservative promotions into `program/`, `theory/`, and
`papers/`.

## Goal

Reduce `docs/` to its intended role:

- code-facing documentation
- architecture and ADRs
- specs and guides
- implementation-facing reference material

while moving or archiving the remaining mixed material.

## Current remaining `docs/` categories

### Keep in `docs/`

- `docs/architecture/`
- `docs/adr/`
- `docs/spec/`
- `docs/cli.md`
- `docs/testing.md`
- `docs/testing_cli_recipes.md`
- `docs/dev/`
- `docs/perf/`
- `docs/math/`
  - keep only if these notes are implementation-facing or code-supporting
- `docs/blog/`
  - keep for now as general repository writing
  - if this area grows into a recurring workstream, consider promoting it into
    a dedicated project such as `projects/software-patterns/`

### Move or classify outside `docs/`

- top-level roadmap/vision material still in `docs/`
  - classify into `program/`, `projects/`, or archive as appropriate

Completed migrations:

- former `docs/eLife 2025/` material has been split into:
  - active manuscript sources at
    `papers/macroir-elife-2025/docs/manuscript-drafts/`
  - compiled/history residue at
    `papers/macroir-elife-2025/archive/version-inicial/`
- former `docs/theoretical results/` material has been promoted into
  `theory/macroir/`
- former `docs/audios/` material has been promoted into
  `program/source-notes/audios/`
- `docs/Project Roadmap MacroIR.md` has been promoted into `program/notes/`
- redundant `.docx` version of the MacroDR vision was removed
- stale legacy docs were moved into `program/archive/legacy-docs/`
- `C++ Template Instantiation report.md` was moved into `docs/dev/`

### Archive or remove

- dated `docs/2025_*.txt`
  - already moved into `program/source-notes/audios/`
- cross-cutting legacy docs
  - already moved into `program/archive/legacy-docs/`

### Reference/support material needing explicit decision

- `docs/bibliography/`

These should either remain as support/reference material under `docs/` or move
to a more explicit archival/reference location.

## Proposed next sequence

1. Clean obsolete documentation:
   - continue removing stale duplicated or superseded docs now that promoted
     homes exist
2. Classify the remaining top-level non-code docs:
   - review whether any remaining top-level docs should still leave `docs/`
3. Re-evaluate support folders:
   - `bibliography/`
   - `math/`
   - `perf/`
   - `dev/`
4. Keep `docs/blog/` for now unless it grows into a broader standalone writing
   area or a dedicated project such as `projects/software-patterns/`

## Guiding rule

If a document mainly helps understand, use, test, or implement the current
codebase, it may remain in `docs/`.

If it mainly defines theory, paper narrative, program vision, or archive
history, it should gradually leave `docs/`.
