# MacroIR / MacroDR

This repository is the working repository of the MacroIR / MacroDR program.

It contains several distinct kinds of work at once:

- the evolving code baseline
- the evolving theoretical baseline
- bounded execution and governance projects
- publication-facing paper material
- generated operational workspaces

The repository is currently being reorganized so those layers are easier to
distinguish.

## Current top-level structure

### Program identity

- `program/`
  - program-level identity, governance, notes, source notes, and archive

### Theory

- `theory/`
  - destination for curated baseline theory, theory notes, and theory archive

### Papers

- `papers/`
  - destination for publication-facing material

### Projects

- `projects/`
  - bounded efforts such as reorganization, framework work, and execution
    projects

### Code baseline

The code baseline remains in the existing source tree:

- `include/`
- `src/`
- `tests/`
- `tools/`
- `legacy/`

### Canonical code-facing docs

- `docs/`
  - architecture, ADRs, specs, and guides

### Generated or operational state

- `build/`
- `Testing/`
- `.cache/`
- `stats/`

These are not releases. They are operational outputs and remain outside the
curated program structure.

## Migration status

The repository is in an incremental reorganization.

Current migration-governance projects:

- `projects/program-project-framework/`
- `projects/macro-dr-reorganization/`

The current repo-specific target model is documented in:

- `projects/macro-dr-reorganization/TARGET.md`

## Working rule during migration

Until the migration is further along:

- `docs/` is now mostly code-facing documentation plus a small number of
  support/reference areas
- `program/`, `theory/`, and `papers/` carry most non-code material
- the source tree for code remains in place

The goal is to improve structure without breaking current work.
