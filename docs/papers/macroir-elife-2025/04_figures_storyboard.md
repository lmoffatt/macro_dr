# Paper 2 (MacroIR / eLife 2025) — Figures Storyboard

This is the “spec” for the main figures. Update it before updating code.

## Figure 1 — Update timing: MR vs IR

**Question:** What is the conceptual difference between MacroMR and MacroIR?

**Panels (proposed):**

- (A) Schematic timeline: interval, Bayes update(s), Markov propagation
- (B) Example open-probability update plot with arrows (MR vs IR)
- (C) One-sentence punchline: “IR uses information more symmetrically across interval boundaries.”

**Existing prototype:**

- Data generator: `projects/eLife_2025/ops/local/figure_1.macroir`
- Plot: `projects/eLife_2025/figures/figure_1.Rmd`

**Acceptance criteria:**

- A reader can *visually* see “what changes” between MRV and IRV in one example.

## Figure 2 — Score (derivatives) across time

**Question:** Do the derivatives computed by the approximation look stable/consistent?

**Panels (proposed):**

- (A) Per-interval log-likelihood contributions summary vs time (mean ± SE over simulations)
- (B) Per-parameter score contributions vs time (mean ± SE)
- (C) Compare at least two schemes (e.g. `scheme_CO` and `scheme_CCO`) to show robustness

**Existing prototype:**

- Data generator: `projects/eLife_2025/ops/local/figure_2.macroir`
- Plot: `projects/eLife_2025/figures/figure_2.Rmd`

**Acceptance criteria:**

- The figure makes it obvious which algorithm yields biased/noisy gradients (if any).

## Figure 3 — Validation: Score mean + FIM agreement + residuals

**Question:** Does the approximation pass basic statistical identity checks?

**Panels (proposed):**

- (A) Score mean vs number of simulations (CI shrinking, contains 0)
- (B) FIM agreement summary (heatmap of relative errors or eigenvalue ratios)
- (C) Residual diagnostics (mean/var and autocorrelation; optionally PIT histogram)
- (D) Score correlation diagnostic (Cov(sum score) vs sum Cov(score_t))

**Data needed (new, beyond existing prototypes):**

- repeated simulations at true parameters
- extracted score and Hessian (or whatever FIM estimators are available)
- per-interval predictive mean/variance for residuals

**Acceptance criteria:**

- MacroIR passes (or is clearly better) in a regime where others fail.

## Figure 4 — Validity map: Δ/τ_min vs N_ch

**Question:** In which regimes is each approximation reliable?

**Panels (proposed):**

- (A) Validity map per algorithm (or one map with “best algorithm” per cell)
- (B) Same map but using a different metric (robustness check)
- (C) Compute-time cost overlay (optional; secondary)

**Data needed (new):**

- the full grid definition from `docs/papers/macroir-elife-2025/05_experiment_grid.md`
- a metric table (algorithm × grid cell)

**Acceptance criteria:**

- The map gives an actionable rule-of-thumb (“use IR when …”).

