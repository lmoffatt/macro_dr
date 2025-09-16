# MacroDR Core — Implementations

- Purpose: final implementations behind CMD; orchestrate Models, Inference, IO, Math, Utils.
- What goes here: `.cpp` files only; no DSL/CLI code. Keep boundaries thin and typed.
- May wrap `legacy/` components during migration; behavior should remain unchanged.
- Tests: prefer unit tests here and CLI smoke tests that exercise CMD surfaces.
- See `docs/architecture/modules.md` and ADRs 001–002 for design and contracts.

