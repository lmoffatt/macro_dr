# Repository Classification

This document classifies the current `macro_dr` repository according to the
Program/Project Framework defined in Project 1.

## Purpose

The goal is to identify what each major area of the current repository is doing
today and where it should conceptually belong in the target organization.

The classification uses these categories:

- `program`
- `theory`
- `code`
- `project`
- `paper`
- `archive`
- `generated`

## Classification summary

### Program-level material

These files define the current baseline structure or repository governance and
should eventually move under `program/` or under canonical code/theory docs
once the new structure exists.

- `AGENTS.md`
  - role: program/code governance
  - status: canonical
- `projects/program-project-framework/`
  - role: project defining the framework used by the program
  - status: canonical active project
- `projects/macro-dr-reorganization/`
  - role: project applying the framework to this repository
  - status: canonical active project

### Code baseline

These areas are clearly part of the evolving baseline implementation and should
ultimately belong under `code/`.

- `include/`
  - role: code
  - status: canonical
- `src/`
  - role: code
  - status: canonical
- `tests/`
  - role: code validation
  - status: canonical
- `tools/`
  - role: code support
  - status: canonical
- `legacy/`
  - role: code
  - status: canonical but legacy-scoped

### Canonical code-facing documentation

These documents are baseline-facing and should ultimately be attached to
`program/` or `code/docs/`.

- `docs/architecture/`
  - role: canonical code/program docs
  - status: canonical
- `docs/adr/`
  - role: canonical decisions
  - status: canonical
- `docs/spec/`
  - role: canonical behavioral/docs baseline
  - status: canonical, though some files are rough drafts
- `docs/cli.md`
  - role: canonical code-facing guide
  - status: canonical
- `docs/testing.md`
  - role: canonical code-facing guide
  - status: canonical
- `docs/testing_cli_recipes.md`
  - role: canonical code-facing guide
  - status: canonical
- `include/macrodr/cmd/README.md`
  - role: code-facing documentation
  - status: canonical
- `src/core/README.md`
  - role: code-facing documentation
  - status: canonical

In this framework, `docs/` is the software-facing interpretive layer:

- architecture
- ADRs
- specs
- guides
- explanations of how the code realizes the theory

### Theory baseline vs theory archive

The theory material is mixed. Some documents are close to the active conceptual
baseline; many others are exploratory or archival.

Likely baseline-facing theory:

- `docs/macro_dr_vision_and_realization.md`
- `docs/Proof_Oriented_Design_Manifesto.md`
- `docs/general_scientific_software_design_theory.md`
- `docs/provenance_types_hott.md`
- `docs/algebraic_CLI.md`
- `docs/contagiously_safe_design.md`
- `docs/composite_safety_domains.md`
- `docs/hott_contagion_principle.md`
- `docs/runtime_vs_semantic_safety_narrative.md`
- `docs/safety_explosion_limit.md`

Mixed theory area requiring internal classification:

- `docs/theoretical results/`
- dated note files under `docs/2025_*.txt`
- much of `docs/audios/`

These should eventually be split into:

- active baseline theory
- active theory notes
- archive

Important refinement:

- `docs/theoretical results/` should not be treated as a single archive bucket
- it contains a mixture of:
  - core baseline theory
  - active theory notes
  - archival draft variants
- Information Distortion Matrix material there is theory, not code-facing docs

### Papers

These are publication-facing and should ultimately live under `papers/`.

- `papers/macroir-elife-2025/`
  - role: active paper planning pack
  - status: canonical active paper material
- `papers/macroir-elife-2025/docs/manuscript-drafts/`
  - role: active manuscript draft sources
  - status: canonical active paper-local docs
- `papers/macroir-elife-2025/archive/version-inicial/`
  - role: compiled manuscript residue and older local history
  - status: paper-local archive
- `papers/archive/behind-the-paper/`
  - role: published paper-adjacent narrative material
  - status: archive

### Execution projects

These are bounded efforts that apply the code and should remain under
`projects/`, though they need cleanup and clearer internal structure.

- `projects/p2x2/`
  - role: execution project
  - status: canonical active project
- `projects/eLife_2025/`
  - role: execution project / paper support project
  - status: canonical active project

### Archive-heavy material

These areas are largely archival or raw source material rather than baseline
truth.

- `program/source-notes/audios/`
  - role: source notes pending extraction
- `papers/macroir-elife-2025/archive/elife-template-import/`
  - role: manuscript support / imported template
- `docs/CODEBASE_MAP.md`
  - role: obsolete/stale generated overview
- `docs/command_registry.md`
  - role: legacy-oriented documentation
- `docs/command_registry_advanced.md`
  - role: legacy-oriented documentation
- duplicate or former-version files such as:
  - `docs/composite_safety_domains (1).md`
  - multiple `behind the paper` variants
  - duplicated audio transcription files

### Generated / operational material

These areas are not conceptual categories of the framework and should not be
modeled as releases.

- `build/`
  - role: generated operational workspace
- `.cache/`
  - role: generated operational workspace
- `Testing/`
  - role: generated operational workspace
- `stats/`
  - role: generated/measurement workspace

## Immediate classification conclusions

1. The repo already behaves like a program repository, not a single-project
   repository.
2. `docs/` is overloaded: it contains canonical docs, theory, paper material,
   and archive at once.
3. `projects/` already contains real execution projects and is worth keeping.
4. `papers/` exists inside `docs/`, but should conceptually become top-level.
5. Theory needs a stronger distinction between active baseline theory and
   archive/exploratory notes.
6. `build/` is generated operational state, not a release.
