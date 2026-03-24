# Lifecycle And Checks

This document defines the minimal stage model and validation logic for both
projects and programs.

## Project stages

Projects move through these stages:

1. `idea`
   - question exists, but scope is still fluid
2. `framed`
   - purpose, scope, deliverables, and success criteria are written
3. `designed`
   - approach, validation strategy, and dependencies are defined
4. `implemented`
   - main artifacts exist in working form
5. `validated`
   - required checks pass
6. `reproducible`
   - another person or agent can rerun the work from documented instructions
7. `delivered`
   - the intended deliverable has been produced
8. `archived`
   - project closed or superseded

Projects do not need to end in a paper. Delivery may be code, a benchmark, a
technical assessment, a report, a demo, or a negative result.

## Project checks

Project checks are evidence that a bounded effort is actually progressing.

Recommended categories:

- `docs`
  - purpose, scope, deliverables, success criteria, current status exist
- `design`
  - validation plan and dependencies are explicit
- `build`
  - code compiles, where code is part of the project
- `tests`
  - relevant automated tests pass
- `science`
  - project-specific scientific or mathematical validations pass
- `repro`
  - runbook exists and outputs can be reproduced
- `delivery`
  - intended deliverable exists and is documented

## Program dimensions

A program is assessed less by a single stage and more by four dimensions:

- `maturity`
  - how far the program has progressed as an organized body of work
- `health`
  - whether the current baseline builds, tests, and runs
- `coherence`
  - whether code, theory, docs, projects, and papers are aligned
- `readiness`
  - whether the current baseline is suitable for release or external use

## Program checks

Program checks evaluate the health of the baseline and the system of work.

Recommended categories:

- `baseline_build`
  - baseline code compiles
- `baseline_tests`
  - core tests pass
- `baseline_docs`
  - canonical docs exist and are current enough to navigate the program
- `baseline_theory_alignment`
  - theory/code mapping has been reviewed or documented
- `project_registry`
  - active projects are declared and tracked
- `project_status`
  - active projects have explicit stage and checks
- `release_readiness`
  - release criteria are defined and assessed
- `repro_infrastructure`
  - reproducible execution path exists for core workflows

## Exit and readiness logic

### Project

A project should not be treated as complete merely because something was
implemented. Completion should normally require:

- deliverable exists
- required checks pass
- current status is documented
- exit criteria are satisfied

### Program

A program baseline should not be treated as ready merely because code exists.
Program readiness should normally require:

- baseline compiles
- baseline core tests pass
- core docs exist
- active projects are tracked
- release criteria are defined

## Rule of thumb

- A project asks: "did this bounded effort succeed?"
- A program asks: "is the evolving whole coherent, healthy, and ready to
  absorb more work?"
