# Theory

This folder is the destination for theory that belongs to the evolving program
baseline or to active theory work.

It is intended to separate:

- MacroIR-specific theory
- scientific software theory that is part of the program's conceptual
  foundation
- curated baseline-facing theory
- active exploratory theory notes
- superseded or archival theory material

## Structure

- `macroir/`
  - MacroIR-specific theory family
- `scientific-software/`
  - scientific software theory family

Each theory family should have:

- `docs/`
  - curated baseline-facing theory
- `notes/`
  - active exploratory theory work, synthesis drafts, comparison notes, and
    provisional formulations that are not yet stable enough for `docs/`
- `archive/`
  - superseded theory material

Theory remains theory even when it is central to implementation or validation.
`docs/` is reserved for software-facing material such as architecture, specs,
guides, and explanations of how the code realizes the theory.

## Migration status

Most theory material still lives under `docs/`:

- selected baseline-facing theory documents
- mixed exploratory work under `docs/theoretical results/`
- archival/source-note material under `docs/audios/` and dated note files

This folder exists so theory can be promoted out of the mixed `docs/` surface
incrementally.

## What belongs in `notes/`

`notes/` is for theory work that is active but not yet part of the curated
baseline.

Typical contents:

- exploratory derivations
- alternative formulations
- comparison notes between competing ideas
- synthesis drafts that are still changing
- questions, open issues, and partial resolutions
- bridge notes from raw source material toward curated theory documents

What does not belong in `notes/`:

- stable baseline theory that should be promoted into `docs/`
- purely historical theory that should go into `archive/`
