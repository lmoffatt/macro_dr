# Organization Model

This document proposes a practical top-level organization for the `macro_dr`
repository.

## Core principle

The repository is a technical container for a program. It is not assumed to be
equivalent to a single project.

The organization should separate:

- the program definition
- the evolving baseline theory
- the evolving baseline code
- bounded projects
- publication-facing papers
- archival material

## Proposed top-level structure

```text
program/
theory/
code/
projects/
papers/
```

### `program/`

Contains the program-level identity and governance documents:

- purpose
- roadmap
- decisions
- baseline and readiness concepts

### `theory/`

Contains theory that belongs to the evolving baseline of the program.

Recommended substructure:

- `README.md`
- `docs/` for curated theory
- `notes/` for active exploratory work
- `archive/` for superseded theory documents

### `code/`

Contains the evolving baseline implementation and its code-facing
documentation.

Recommended substructure:

- `README.md`
- `docs/`
- `include/`
- `src/`
- `tools/`
- `tests/`

### `projects/`

Contains bounded efforts that operate on, validate, extend, or exploit the
baseline.

Recommended substructure for each project:

```text
projects/<name>/
  README.md
  PROJECT.md
  PROJECT.yaml
  theory/
  code/
  results/
  runs/
  archive/
```

Not every project needs all subfolders, but this is the default pattern.

### `papers/`

Contains publication-facing material.

Recommended substructure:

- `README.md`
- `outline.md`
- `figures.md`
- `sources.md`
- `drafts/`
- `archive/`

## Archive rule

Archive material should be organized primarily by ownership and context, not by
replicating the full active repository tree inside a single global archive.

### Preferred rule

- project-specific historical material belongs in `projects/<name>/archive/`
- paper-specific historical material belongs in `papers/<name>/archive/`
- theory-specific historical material belongs in the owning theory family, for
  example `theory/macroir/archive/` or
  `theory/scientific-software/archive/`

This keeps historical material close to the active workstream it came from and
avoids creating a second competing navigation tree.

### Top-level archive

A small top-level archive may still exist, but only for material that is:

- orphaned
- imported from older repositories or vaults
- cross-cutting and not owned by a single active project or workstream

It should not be used to mirror the active tree.

## Build vs release

The repository may contain operational directories such as `build/`, but those
are not conceptually the same as releases.

- `build/` is local, regenerable, machine-specific output
- a `release` is a curated, validated snapshot of the baseline for other people

## Document splitting rule

Start with one `README.md` per concept or folder.

Split it when one or more of these becomes true:

- it exceeds roughly 150 lines
- it serves more than one purpose
- users repeatedly need to link to one section independently
- multiple subtopics evolve separately

When splitting creates three or more companion documents, treat the concept as
a folder with an index `README.md`.

## Immediate implication for `macro_dr`

This framework does not require a one-step big-bang move. It can be applied
incrementally:

1. define the vocabulary and manifests
2. identify existing material by category
3. move active canonical material first
4. isolate archive material
5. reorganize code and theory only as needed to improve coherence
