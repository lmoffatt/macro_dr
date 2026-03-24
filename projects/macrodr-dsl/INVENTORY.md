# DSL Inventory

This document records current DSL pain points and candidate improvements.

The goal is to start from real blocking or confusing cases in active projects,
not from abstract language redesign.

## Inventory template

For each case, record:

- case id
- owning project
- script or file
- symptom
- current behavior
- what is missing
- candidate DSL-side improvement
- diagnostics linkage

## Case 1 - `figure_1` failure context

- Case id: `dsl-001`
- Owning project: `eLife_2025`
- Script or file:
  - `projects/eLife_2025/ops/local/figure_1.macroir`
- Symptom:
  - the `figure_1` workflow produces failures that are hard to diagnose from
    current reporting alone
- Current behavior:
  - the active workflow can fail without sufficiently clear script-level context
    about where and why evaluation broke
- What is missing:
  - clearer DSL-side structure for locating the failing expression or command
  - a better notion of evaluation context that diagnostics can surface
  - easier identification of which command output is expected to feed later
    commands
- Candidate DSL-side improvement:
  - inventory the structure of the script and identify which constructs need
    clearer naming, staging, or source-location tracking
- Diagnostics linkage:
  - this case is also the first acceptance case for `projects/macrodr-diagnostics/`

## Next inventory steps

1. inspect the `figure_1.macroir` workflow in detail
2. identify the concrete DSL constructs involved in the failure
3. separate DSL-side issues from diagnostics-only issues
4. add more cases from active scripts once the first case is clear
