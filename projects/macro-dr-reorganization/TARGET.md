# Target Structure for `macro_dr`

This document defines the first practical target structure for the `macro_dr`
repository.

It adapts the general organization model from Project 1 to the actual
constraints of this repository.

## Pragmatic principle

The generic framework proposed this top-level structure:

```text
program/
theory/
code/
projects/
papers/
```

For `macro_dr`, the conceptual distinction still holds, but the implementation
should be less disruptive:

- `program/`, `theory/`, and `papers/` should become explicit top-level
  folders
- `projects/` already exists and should remain the home of bounded efforts
- the code baseline should continue to live in the current source tree
  (`include/`, `src/`, `tests/`, `tools/`, `legacy/`) rather than being
  physically moved under a new `code/` directory
- `docs/` should remain, but its scope should narrow toward canonical
  code-facing and program-facing documentation

This avoids unnecessary disruption to CMake, include paths, and current work.

## Repository target model

The first target model for `macro_dr` is:

```text
README.md
program/
theory/
papers/
projects/
docs/
include/
src/
tests/
tools/
legacy/
```

Generated and operational directories remain outside the conceptual structure:

```text
build/
Testing/
.cache/
stats/
```

These directories are not part of the program architecture. They are generated
or operational state.

## Meaning of each top-level area

## Generated and operational layer

The following directories remain outside the curated program structure:

- `build/`
- `Testing/`
- `.cache/`
- `stats/`

In particular:

- `build/` is a local compilation workspace
- it is regenerable and machine-specific
- it is not a `release`
- it should not be interpreted as part of `program/`, `theory/`, `papers/`, or
  `projects/`

For now it should remain at repo root, since current CMake and workflow
conventions already expect that layout.

### `README.md`

Root orientation for the repository.

It should answer:

- what the repository is
- how the program is currently organized
- where to look for program identity, theory, code, projects, and papers
- which areas are still in migration

### `program/`

Program-level identity and governance.

Expected contents:

- program definition
- program roadmap
- release/readiness concepts
- links to current manifests and active projects

This folder is about the program as a whole, not about one bounded effort.

### `theory/`

Theory that belongs to the evolving program baseline.

Recommended internal structure:

```text
theory/
  README.md
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

Use:

- `theory/macroir/` for MacroIR-specific mathematical and algorithmic theory
- `theory/scientific-software/` for scientific software theory that belongs to
  the program's conceptual foundation
- each family's `docs/` for curated baseline theory
- each family's `notes/` for active exploratory theory
- each family's `archive/` for superseded theory material

### `papers/`

Publication-facing material.

Recommended internal structure:

```text
papers/
  README.md
  <paper-name>/
```

Each paper folder may then contain its own:

- `README.md`
- `outline.md`
- `figures.md`
- `sources.md`
- `drafts/`
- `archive/`

### `projects/`

Bounded efforts that validate, extend, reorganize, or apply the baseline.

Examples already present:

- `projects/p2x2/`
- `projects/eLife_2025/`
- `projects/program-project-framework/`
- `projects/macro-dr-reorganization/`

### `docs/`

Canonical repository documentation that remains code-facing or
architecture-facing.

In the target model, `docs/` should primarily contain:

- architecture
- ADRs
- specs
- canonical guides such as CLI and testing
- documents that explain how the codebase relates to or realizes theory

It should gradually stop being the catch-all home of theory drafts, papers,
and raw archives.

### Code baseline in place

The code baseline remains where it is:

- `include/`
- `src/`
- `tests/`
- `tools/`
- `legacy/`

This is a deliberate deviation from the generic `code/` folder proposal.

The conceptual category is still `code`, but the physical structure should
remain stable for now.

## Current-to-target mapping

### Program

- current source:
  - scattered governance material
  - `projects/program-project-framework/PROGRAM.yaml`
- target:
  - `program/`

### Theory

- current source:
  - selected theory docs under `docs/`
  - mixed exploratory material under `docs/theoretical results/`
  - source notes in `docs/audios/` and dated note files
- target:
  - `theory/macroir/docs/`
  - `theory/macroir/notes/`
  - `theory/macroir/archive/`
  - `theory/scientific-software/docs/`
  - `theory/scientific-software/notes/`
  - `theory/scientific-software/archive/`

Important clarification:

- theory remains theory even when it is central or practically important
- Information Distortion Matrix material belongs in `theory/`
- only implementation-facing explanatory material about how the code realizes
  that theory belongs in `docs/`

### Papers

- current source:
  - legacy paper material from old `docs/` locations
- target:
  - `papers/`
  - with active paper-local docs under `papers/<name>/docs/`
  - and paper-local history under `papers/<name>/archive/`

### Projects

- current source:
  - `projects/`
- target:
  - remains `projects/`

### Code-facing docs

- current source:
  - `docs/architecture/`
  - `docs/adr/`
  - `docs/spec/`
  - `docs/cli.md`
  - `docs/testing.md`
  - `docs/testing_cli_recipes.md`
- target:
  - remains `docs/`

## First applied moves

The first structural moves should be:

1. Create `README.md` at repo root.
2. Create `program/` as the destination for program identity.
3. Create `theory/` scaffolding without moving theory content yet.
4. Create `papers/` scaffolding without moving paper content yet.
5. Keep the code tree in place.

This makes the intended structure visible before any large content move.

## What this target model does not do yet

- It does not move `include/` and `src/` under `code/`.
- It does not immediately rewrite existing `docs/` content.
- It does not collapse all old content into an archive.
- It does not decide every final destination for every historical file.

Those will be handled incrementally.
