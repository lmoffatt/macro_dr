# Figure 5 family — placement plan and justifications

Audit of the 44 `figure_5_*` candidates, sorted into **main / supplementary / drop**. What matters here is the *justification* for each placement, not the label. Cross-checked by an independent multi-agent read; the two syntheses agreed on everything except the NR correlation panel (resolved below).

**Data-basis caveat (state once in Methods/caption).** Cross-algorithm panels can only come from the numeric-Fisher run (433ed13, the only commit with all five algorithms). The IR-protagonist mechanism panels use the Gaussian-Fisher run (1c2ae6f). Main-text panels are therefore not all on the same F. This split is forced by data availability and does not change any qualitative conclusion, but a reader must be told.

---

## MAIN (5) — the figures that carry the thesis

**1. `distortion_algo_grid`** — five-algorithm information-distortion map over the (interval, N_ch) plane.
*Justification:* the single object that makes the core claim legible: IR sits at distortion ≈ 1 (white) while NR/MNR/R/MR saturate to 10³–10⁴, across the whole design plane. Already curated onto one shared symmetric-log scale with 1.15/1.5 validity iso-lines. Beats `info_distortion_contour_algos`, which tells the identical story but fans out to 10 uncurated PDFs. *Caveat to caption:* numeric-F basis; the extreme naive values are partly finite-difference-Fisher instability, not all signal.

**2. `distortion_decomp_grid`** — sample × correlation decomposition, all five algorithms.
*Justification:* the only figure that shows the *two-mode* failure: NR/MNR blow up in the correlation factor at short Δ and in the sample factor at coarse Δ, while IR stays clean in both. This decomposition is the distinctive mechanistic contribution, not just a restatement of panel 1. *Caveat:* single noise level; caption the naive extremes as above.

**3. `IR_sample_peak_in_Nch`** — the per-sample distortion is non-monotonic in channel number, peaking at N\* = 2σ²/G, with N\* marching 10 → 1000 as noise rises.
*Justification:* the novel, counter-intuitive, verified result. It says the worst case for the Gaussian approximation is an intermediate channel count fixed by the recording noise, not the fewest channels. *Rebuilt (2026-07-08) with the denser low/mid-N_ch runs (noise 0.05/0.2/0.5/1 now filled to N_ch 100–200) and the method fixed:* N\* is now the vertex of a log-log parabola through the three points around the peak, not a raw argmax; peaks pinned at the grid edge (low noise, N\* ≤ 10) are left unmarked rather than plotted as resolved; and a companion panel (`Figure_5_IR_sample_Nstar_scaling.pdf`) tests the scaling law, giving N\* ∝ noise^0.80 (R² 0.59; the low-noise interior peaks 0.2–0.5 sit on the line, the scatter is the vanishing signal at noise ≥ 1). Still worth overlaying the absolute N\* = 2σ²/G curve if the open probability p is pinned down.

**4. `bias_empirical_grid`** — empirical bias (estimate / truth) of all five algorithms.
*Justification:* the correct empirical quantity (θ_pool − θ_sim), not the distortion-induced-bias predictor, which is ≈ 0 at the pooled MLE by design and misses the NR catastrophe. Shows IR unbiased and the naive methods grossly biased in N_ch, in one grid; this is thesis element (b). Use the *pooled* version for main, the *cloud* version as its supplementary twin.

**5. `IR_channel_pooled_loss`** — covariance penalty of concentrating vs spreading channels at a fixed budget.
*Justification:* the crispest single design statement (concentrating channels is ≈ free for k_off, ≈ 4× for k_on, up to ≈ 500× for i), on the honest Gaussian distortion-corrected covariance. Carries thesis element (c). Chosen over `covariance_design`/`covariance_contour`/`channel_loss`, which restate the same content with more panels or as its own per-step derivative.

---

## SUPPLEMENTARY — grouped by the role they play

**Why the naive methods fail (mechanism of the competitor)**
- `NR_correlation_collapse` — the correlation factor of NR collapses cleanly onto Δ/τ (R² ≈ 0.8, up to ≈ 200×, monotonic). *Justification:* the one collapse in the whole set that actually works, because the correlation leg is genuinely monotonic in Δ/τ (temporal over-counting). It explains *why* the naive estimator fails, which a method paper is stronger for showing; kept out of main to keep the main narrative on IR. Soften the "independent of N_ch" claim (holds for N_ch ≥ 100) and clip the sub-1 outliers at Δ/τ = 1.

**IR mechanism, detail**
- `decomp_gaussian` — rigorous Gaussian IR-only sample-vs-correlation decomposition on the design plane. *Fix the stale header comment* (it names the reconstructed correlation but the code uses the direct component; code is right).
- `IR_decomp_lines_G` — the same decomposition as a line view vs interval, with the ±15% self-consistency band front and centre; uses the analytic Gaussian FIM that is reliable exactly at the stressed corner.
- `IR_distortion_ridge` — 2D heatmap of the N\* peak in the (Δ, N_ch) plane. *Remove the argmax "ridge" overlay* (it is an argmax on a near-flat field and jumps 10 → 10⁴); replace with the predicted N\* curve, drop the near-empty i panel.
- `IR_sample_peak_migration` — the interval-axis (Δ\*) twin of the N\* result. Restrict markers to interior peaks; match noise sets across N_ch.
- `IR_sample_noise` — the interval crossover of the per-sample hump, with real bootstrap CIs. The noise-flattening is subtle, so it supports rather than carries.
- `domains_validation` — affine-invariant (Riemannian SPD) distance of the distortion matrix from identity: the reparam-invariant "distance from Gaussian". Validates the domain ordering and the noise-fade; IR-only, so supporting.

**Bias, validation**
- `bias_score_ci_grid` — theoretical (score-based) bias companion for all algorithms, median-FIM stabilised and significance-gated. Disciplines the degenerate i·N_ch cells; supersedes the plain `bias_score_grid`.
- `bias_lines` — IR empirical-vs-DIB vs interval with CIs. *Fix the title clipping to "group_size = 1".*
- `bias_path` — the only bias-*direction* view (the k_on–N_ch / i–N_ch identifiability ridge). Trim to one noise; it is large (16 in).

**Design / elasticity**
- `IR_covariance_design` — per-recording SE landscape vs interval and N_ch. *Fix the stale header* (names `Likelihood_Fisher_Covariance`; actually uses the Gaussian corrected covariance).
- `covariance_contour` — the design-plane map (kinetics flat, i/N_ch degeneracy grows). Regenerate on the Gaussian corrected covariance to kill the noise-10 sandwich artifact.
- `IR_interval_loss` — interval elasticity (oversampling trade-off). *Fix the doubly-stale header* (wrong component and wrong formula) before it feeds a caption.
- `IR_normcov_noise` — channel sub-linearity via N_ch·var. Legible but overlaps the other channel-loss cuts.

**Gaussian-vs-numeric Fisher (methods robustness)**
- `distortion_gauss_vs_numeric_Nch10` — the sharp ≈ 8× variability gap at few channels; the best single justification of the Gaussian-Fisher substitution.
- `distortion_gauss_vs_numeric_lines` — the dense line grid showing numeric ≈ Gaussian everywhere. *Drop/retitle the `relCI_by_Nch` "equal at every N_ch" panel:* it is a median-over-intervals artifact that averages away the short-interval blow-up and contradicts the Nch10 sibling.
- `fisher_variability`, `fisher_ratio_contour_algos` — collapse to a single panel (IR white + NR divergence), not the full 10-PDF fan-out.

**Orientation / schematic**
- `domains_schematic` — cartoon of the three non-Gaussianity sources. *Only usable if the hard-coded N_ch = 30 multinomial boundary is redrawn as the noise-dependent N\* = 2σ²/G line;* as drawn it contradicts the paper's own N\* result.
- `logL_diff_grid` — design-decision orientation (recursion is the big win, the naive M term hurts). Off-thesis metric, so supporting.
- `distortion_contour_masked` — IR-only four-parameter detail with honest reliability masking. Regenerate on the Gaussian basis; keep the information-distortion panel only.

---

## DROP (16) — with the specific reason

| figure | reason |
|---|---|
| `bias_algo_grid` | plots the predicted DIB, which is ≈ 0 at the pooled MLE and misses the NR catastrophe — understates the very contrast it should show; no rendered PDF |
| `bias_contour` | least quantitative IR DIB validation; the same comparison is carried by `bias_lines` (with CIs) and `bias_path` |
| `bias_score_grid` | strictly dominated by `bias_score_ci_grid` (median FIM + significance) |
| `covariance_lines` | exploratory three-view; the contour carries the concept more cleanly |
| `distortion_contour` | unmasked twin of `distortion_contour_masked`; without the mask the noise-10 speckle reads as real distortion |
| `distortion_decomp_grid_by_param` | duplicate layout of `distortion_decomp_grid`, single noise |
| `distortion_diag` | corner/pairwise scatter is a category error for a per-parameter scalar; `lines_bynch` says it better |
| `distortion_lines` | numeric IR-only per-parameter lines, superseded by the Gaussian line view |
| `distortion_lines_bynch` | same, and 9 rows are untenable for print |
| `first` | exploratory DIB calibration, superseded by the later bias figures |
| `IR_channel_loss` | per-step reciprocal of `pooled_loss`; its header's α/L narrative is inverted vs the plotted quantity |
| `IR_decomp_collapse` | four collapse variables, none collapses; division/amplification artifact in panel D |
| `IR_decomp_lines` | numeric twin of `IR_decomp_lines_G`; the Gaussian version is reliable exactly where it matters |
| `IR_sample_collapse_Nnoise` | **wrong premise** — assumes the swept noise scales as N·noise/Δ so N cancels; contradicts the verified y_var = e + N·gSg + N·ms (Current_Noise enters at N⁰). The (D−1)·√N panels are amplified sampling noise |
| `IR_sample_collapse_W` | inconclusive, admittedly mis-parametrized single-variable collapse scan |
| `NR_sample_collapse` | **wrong** — forces a peaked quantity (sample distortion, peak at N\*) onto a monotonic collapse axis; opposite sides of the peak map to the same x. No collapse possible in principle |

---

## Method-WRONG flags (never use as evidence, even in supp, without the fix)

- `IR_sample_peak_in_Nch`: argmax of N\* on a coarse grid → boundary-censored markers. Fixable (interpolate + mark boundaries + overlay predicted law).
- `IR_distortion_ridge`: the argmax ridge overlay is noise, not a ridge. Remove it.
- `distortion_gauss_vs_numeric_lines` `relCI_by_Nch` panel: median-over-intervals artifact; contradicts the Nch10 sibling.
- `domains_schematic`: the static N_ch = 30 multinomial line contradicts N\* ∝ noise.
- Stale header comments in `decomp_gaussian`, `IR_covariance_design`, `IR_interval_loss` name the wrong component/formula. If they feed a caption they will mislead. Correct before use.

---

## Gaps the thesis needs and no figure covers well

1. **Product reconstruction.** Every decomposition figure shows sample and correlation as two panels; none shows sample × correlation = total J/F on the same axes. The "multiplicative" claim is asserted, never demonstrated. One verification panel (the product overlaid on the measured total) closes it.
2. **N\* against its predicted law.** `IR_sample_peak_in_Nch` shows the measured N\* but does not overlay N\* = 2σ²/G. The overlay turns "marches with noise" into "matches the theory".

---

## Suggested arc

- **Figure 5** = (a) `distortion_algo_grid` + (b) `bias_empirical_grid` + (c) `IR_channel_pooled_loss`.
- **Figure 6 (mechanism)** = `distortion_decomp_grid` + `IR_sample_peak_in_Nch`, the novel contribution, given its own figure rather than compressed into Figure 5.
