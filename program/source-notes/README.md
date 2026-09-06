# Source Notes

This folder is for raw or lightly processed source material that still needs
analysis and extraction before it can be promoted into more structured
documentation.

Current contents include:

- `audios/`
  - transcripts of exploratory spoken thinking about the repository and the
    program

## Status

These materials are not canonical documentation.

They are source notes pending analysis. After extraction, they may lead to:

- updates in `program/`
- updates in `projects/`
- updates in `theory/`
- archival storage

## Audio transcription tooling

Audio transcription tooling exists outside this repository under:

- `/home/lmoffatt/Projects/scripts/`

This repository keeps the resulting transcripts as source material, not as
final documentation.

## Idea distillation pipeline

`scripts/` (extract_ideas.py, consolidate_ideas.py, mentor.py) distills the
audio transcripts into idea ecosystems under `audios/ideas/`, driven by the
prompts in `prompts/`. Adapted 2026-09-01 from `~/Projects/luthier/scripts/`;
see `scripts/README.md`.
