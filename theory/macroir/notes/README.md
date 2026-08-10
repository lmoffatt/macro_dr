# MacroIR Theory Notes

Active exploratory MacroIR theory, derivations, and comparisons should live
here until they are stable enough for `docs/` or old enough for `archive/`.

Current promoted note families include:

- `Adaptive_MacroR/`
- `Gmean_ij_gvarij/`
- `Macro_Taylor/`
- `Macro_TaylorIR/`
- `MacroLogit/`

These are active theoretical lines that are still exploratory, comparative, or
not yet accepted as stable baseline documentation.

**Status vs the current eLife paper (2026-07):** `Macro_Taylor/` and `Macro_TaylorIR/` are Taylor variance-correction variants, **cut** from the current paper (taylor=false), kept for future. `Gmean_ij_gvarij/` holds the derivations of `gmean_ij` and `gvar_ij`. The gvar_i audit is one level up, `gvar_i_overcount_audit.md`, and its verdict is **historical**: the over-count was fixed in the live kernel (resolution block at the head of that file, 2026-08-06), and what carries forward is that at the same prior `MR` and `IR` predict the same variance, so the difference between them is the gain. Canonical decisions: `papers/1_method/decisions.md` and `papers/_program/decisions.md`; canonical statement of the MR/IR variance question: `papers/1_method/figures_build_plan.md:195-245`. <!-- 2026-08-06: this line pointed the audit at Gmean_ij_gvarij/ (wrong directory) and the decision log at papers/macroir-elife-2025/02_decision_log.md, which no longer exists (the tree was renamed to papers/1_method in the 2026-07-20 split). -->
