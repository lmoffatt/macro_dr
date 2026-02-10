# Paper 2 (MacroIR / eLife 2025) — Workboard

This is the step-by-step checklist. Keep it updated as tasks complete.

## Phase A — Planning docs (now)

- [ ] Read and edit `docs/papers/macroir-elife-2025/00_master_plan.md`
- [ ] Fill/confirm open decisions in `docs/papers/macroir-elife-2025/02_decision_log.md`
- [ ] Edit thresholds + definitions in `docs/papers/macroir-elife-2025/03_metrics_diagnostics.md`
- [ ] Confirm the 4-figure storyboard in `docs/papers/macroir-elife-2025/04_figures_storyboard.md`

## Phase B — Manuscript restructuring (no code changes required)

- [ ] Confirm which LaTeX file is the “source of truth” (`elife-macroir-merged.tex` vs another)
- [ ] Resolve known TODOs in `docs/eLife 2025/version inicial/elife-macroir-merged.tex`
- [ ] Ensure vocabulary consistency: MR vs IR vs “MacroIR”, define all abbreviations early
- [ ] Write a short “When should I use which algorithm?” subsection (ties to Figure 4)

## Phase C — Repro pipeline sanity (use existing scripts first)

- [ ] Run `projects/eLife_2025/ops/local/figure_1.macroir` and regenerate `projects/eLife_2025/figures/data/*`
- [ ] Run `projects/eLife_2025/ops/local/figure_2.macroir` and regenerate `projects/eLife_2025/figures/data/*`
- [ ] Knit:
  - [ ] `projects/eLife_2025/figures/figure_1.Rmd`
  - [ ] `projects/eLife_2025/figures/figure_2.Rmd`
- [ ] Confirm figure intent matches the storyboard (update storyboard or scripts accordingly)

## Phase D — Figure 3 (validation) dataset

- [ ] Decide minimal models + experiments used for validation (likely `scheme_CO` first)
- [ ] For each algorithm variant, generate repeated simulations at true parameters
- [ ] Compute per-algorithm:
  - [ ] Score mean + CI shrink with n
  - [ ] FIM estimator agreement
  - [ ] Residual statistics + autocorrelation
  - [ ] Score-correlation diagnostic
  - [ ] Variance inflation factor

## Phase E — Figure 4 (validity map) dataset

- [ ] Confirm grid definition in `docs/papers/macroir-elife-2025/05_experiment_grid.md`
- [ ] Run grid experiments (Δ/τ_min vs N_ch), store outputs with clear naming
- [ ] Reduce to a compact set of scalar “quality metrics” per grid cell per algorithm
- [ ] Plot validity map + write a short interpretation section

## Phase F — Optional / future work (keep separate; do not block submission)

- [ ] Decide if any MicroIR comparison belongs in Supplementary (default: no)
- [ ] If needed, schedule code tasks from `docs/papers/macroir-elife-2025/07_code_tasks.md` as separate commits

