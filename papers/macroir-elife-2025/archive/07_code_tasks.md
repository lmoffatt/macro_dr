# Paper 2 (MacroIR / eLife 2025) — Code Tasks (ARCHIVED)

> **Archived 2026-07-14.** Its live items moved to `09_carve_plan.md`, "Freeze preconditions", where they belong: they must land in the frozen commit *before* any rerun, because the binary stamps its git hash into every CSV.
>
> Why it was retired and not just updated: its opening rule was "we do not touch code until the plan docs are reviewed", which made it the owner of exactly the code that has to be touched. A self-blocking node. Kept for history; do not read as current.

This file is intentionally *plan-only*. We do not touch code until the plan docs are reviewed.

## Known gaps (from repo inspection)

- Some FIM-related functions appear declared but not fully implemented in the likelihood command surface.
- There is a known valgrind log in `projects/eLife_2025/` indicating an invalid read (needs investigation).

## Candidate code tasks (to schedule later)

- ~~Implement missing FIM-from-predictions helpers~~ — DONE: the Gaussian-Fisher FIM family is built and run.
- Add a small CLI/regression test that validates:
  - score mean ≈ 0 at θ\*
  - two FIM estimators agree within tolerance in a “good regime”
- Ensure diagnostic CSVs expose the necessary primitives cleanly for the R pipeline.
- Fix or mitigate the valgrind issue and add a minimal reproducer.

## Rule

Any code work for Paper 2 should be done in small, isolated commits with a clear purpose and a matching update to:

- `papers/macroir-elife-2025/00_master_plan_v2.md` (§9)
- `papers/macroir-elife-2025/02_decision_log.md`
