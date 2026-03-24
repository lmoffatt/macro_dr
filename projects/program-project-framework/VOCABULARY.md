# Vocabulary

This document is the canonical source of truth for the core organizational
terms used by the MacroIR / MacroDR framework.

## Program

A long-lived research and engineering mission that evolves over time and is
expected to contain multiple bounded efforts, code changes, theory updates,
execution campaigns, and papers.

MacroIR / MacroDR is a program.

## Baseline

The current accepted state of the program: theory, code, and canonical
documentation taken together at a given time.

The baseline is not every idea in the repository. It is the working reference
state that projects test, extend, validate, or reorganize.

## Project

A bounded effort with a specific goal, explicit deliverables, success criteria,
and exit criteria.

A project may update the baseline, validate it, exploit it, or reorganize it.
A paper is one possible outcome of a project, but not the only one.

## Paper

A curated publication-facing deliverable that communicates claims, results,
methods, or interpretation derived from one or more projects.

Papers summarize selected work; they are not the canonical source of truth for
the baseline.

## Release

A curated, versioned, and validated snapshot of the baseline intended for use
by other people.

A release is analogous to a paper, but for the program baseline rather than for
scientific claims. A release is not a local compiler build and is not defined
by the existence of `build/`.

## Run

A concrete execution of code, analysis, experiment, simulation, or validation
procedure within a project.

Runs generate results and other artifacts.

## Artifact

Any output produced by work in the program. Examples include:

- binaries
- build outputs
- figures
- benchmark tables
- datasets
- reports
- manuscripts
- release notes

Artifacts may be ephemeral or curated.

## Archive

Material intentionally kept for historical or traceability reasons but not part
of the active baseline.

Examples include superseded drafts, raw transcriptions, obsolete notes,
deprecated outputs, and historical versions of documents.

## Build

A local, regenerable, machine-specific compilation/output workspace.

`build/` contains technical outputs of compilation and testing. It is an
operational workspace, not a release.

## Relationship summary

- A program has a baseline.
- Projects act on, validate, or extend the baseline.
- Runs occur inside projects.
- Runs produce artifacts.
- Papers and releases are curated artifact classes.
- Archives preserve history without redefining the baseline.
