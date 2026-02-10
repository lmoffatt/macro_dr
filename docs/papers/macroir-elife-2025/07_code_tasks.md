# Paper 2 (MacroIR / eLife 2025) — Code Tasks (Future, Not Now)

This file is intentionally *plan-only*. We do not touch code until the plan docs are reviewed.

## Known gaps (from repo inspection)

- Some FIM-related functions appear declared but not fully implemented in the likelihood command surface.
- There is a known valgrind log in `projects/eLife_2025/` indicating an invalid read (needs investigation).

## Candidate code tasks (to schedule later)

- Implement missing FIM-from-predictions helpers (ensure definitions match paper notation).
- Add a small CLI/regression test that validates:
  - score mean ≈ 0 at θ\*
  - two FIM estimators agree within tolerance in a “good regime”
- Ensure diagnostic CSVs expose the necessary primitives cleanly for the R pipeline.
- Fix or mitigate the valgrind issue and add a minimal reproducer.

## Rule

Any code work for Paper 2 should be done in small, isolated commits with a clear purpose and a matching update to:

- `docs/papers/macroir-elife-2025/01_workboard.md`
- `docs/papers/macroir-elife-2025/02_decision_log.md`

