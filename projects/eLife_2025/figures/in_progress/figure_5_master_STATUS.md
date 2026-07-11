# Figure 5 master — build status & handoff (2026-07-09)

Working files (moved to `figures/paper/` 2026-07-09; this STATUS stays in `in_progress/` as the dev handoff):
- `figure_5_master.Rmd` → `Figure_5_master.pdf` — panel B = **log-det** (signed, diverging).
- `figure_5_master_affine.Rmd` → `Figure_5_master_affine.pdf` — panel B = **affine distance** (magnitude, sequential). DONE 2026-07-09.
- `figure_5_master_frobenius.Rmd` → `Figure_5_master_frobenius.pdf` — panel B = **Frobenius norm** (magnitude, sequential). DONE 2026-07-09.
Panel A (bias) is byte-identical across the three.
Render from `figures/paper/`:
`Rscript -e 'f <- knitr::purl("figure_5_master.Rmd", output=tempfile(fileext=".R"), quiet=TRUE); source(f, echo=FALSE)'`
(paths `../data/...` still resolve because `paper/` is a sibling of `data/` under `figures/`).

## What the figure is
A one-figure **correctness audit** of the 5 algorithms at a **single noise (0.1)**, from the numeric run `../data/433ed13`. Two stacked blocks, each `facet_grid(row ~ algo_f)`, columns = **NR, MNR, R, MR, IR** (established order, IR last), each cell a (Δ·k_off × N_ch) heatmap. Below each block a horizontal stepped colour key.

- **Panel A (bias, 2 rows):** score-based bias = `inv(median Gaussian FIM) · mean-score` at the truth, per parameter **k_off** and **N_ch**. Defined even where the MLE does not converge → far more coverage than empirical θ_pool−θ_sim. Diverging scale (red = over-estimate >1, blue = under <1, white = 1). Rows labelled `bias log10 k_off`, `bias log10 N_ch` (the *parameter* is log10-scaled; the *displayed value* is the linear factor 10^bias). Legend title **"bias factor"**.
- **Panel B (distortion, 3 rows):** program-computed **log-det** of the Information / Sample / Correlation distortion matrices (`log_Det_Likelihood_*`), shown on a **linear det** axis (numbers = exp(log-det)); det=1 = self-consistent. Rows `det information`, `det sample`, `det correlation`. Legend **"det distortion"**.
- x-axis label on BOTH panels: `Δ · k_off` (this is the `interval_in_tau` column; τ_min ≈ 1/k_off). Panel tags **A** and **B** via `labs(tag=)`.

## Data / components (all noise 0.1, 433ed13)
- Bias: `*_battery_sim.csv`, components `Moment_statistics_Sum_dlogL` (score) and `Moment_statistics_Sum_Gaussian_Fisher_Information` (FIM, use probit=quantile q=0.5 = median). Assemble 6×6 FIM, invert (solve→ginv), bias = F⁻¹·score; params index off=1 (→pt[2]), Num_ch_mean=5 (→pt[6]). Valid if finite & non-singular & 95%CI width < CIW_MAX=2 dex.
- Distortion: `*_battery_pool.csv`, `log_Det_Likelihood_{Information,Sample,Correlation}_Distortion`. Scalar per cell (param_index=NA), mean (probit=mean) + 5 bootstrap quantiles (quantile_level 0.025/0.16/0.5/0.84/0.975). `reliable` if 95%CI width ≤ CI_MAX=3 nats.

## VERIFIED facts (don't re-derive)
- **log_Det is NATURAL log (ln)**: linear det = exp(log-det). (ratio logdet/Σln(diag)≈0.97; /log10≈2.3=ln10.)
- **Additivity is EXACT** for IR/R/MR: log-det(info)=log-det(sample)+log-det(correlation), residual 0.000; breaks only in NR/MNR indefinite cells.
- **Singular-FIM cells**: where log-det = exactly 0 the distortion matrix is SINGULAR (condition number Inf, NaN eigenvalues) — a fallback, NOT low distortion. Flag `singular = abs(m) < 1e-9`, mask GREY (`geom_tile` grey85 sized by per-algo `cellw`). Do NOT use condition-number to flag (over-flags: sample/correlation matrices are rank-deficient by design → cond Inf everywhere but log-det still valid).
- Log-det weakness (why we'll compare metrics): near-singular cells cancel (one tiny eigenvalue × one huge → log-det ≈ 0, looks white but is catastrophic). Affine `sqrt(Σ(log λ)²)` doesn't cancel.
- IR score+FIM bias range: k_off factor 0.997–1.056, N_ch 0.984–1.138 (median ~1.002). That's why bias needs fine bands near 1.

## Colour scale logic (settled after long iteration — keep this)
- `RAMP` = the established 13-anchor RdBu (ColorBrewer, colourblind-safe, publishable).
- `bandcols(breaks)`: **perceptually-uniform diverging** — one colour per band = one equal-perceptual step (interpolate `colorRampPalette(rev(RAMP[1:7]), space="Lab")` for blue, `RAMP[7:13]` for red) from the white centre band (`findInterval(0, breaks)`). BOTH sides use the same per-step size (K=max(nb,nr)); band edges (B_LIN/LD_LIN) set the per-algorithm resolution. Blue reaches only as far as data goes; red goes to darkest at ∞.
- Bias edges `B_LIN = c(0.2,0.33,0.5,0.67,0.8,0.9,0.95,0.99,1.01,1.05,1.1,1.25,1.5,2,3,5)` — fine near 1 so IR shows ~3 tones.
- Distortion edges `LD_LIN = c(0.5,0.67,0.8,0.9,1.1,1.25,1.5,2,3,5,10,30,100,300,1e3,1e4,1e5,Inf)`, `LD_BREAKS=log(LD_LIN)`. Inner band 0.9–1.1; up to 1e5 then **Inf**.
- Colour key is drawn MANUALLY (`colorbar()` = a strip of tiles) so **numbers sit AT the boundaries between colours**; labels via `parse(text=fmt_big(...))`; `fmt_big` = superscript `10^n` for ≥1000, `infinity` for Inf (eLife format). Legends at the bottom, one row, full width, thin.
- Contour lines: **bias solid 1.05** (IR's p99 precision), **dashed 1.15**; **distortion solid det 1.15** (VALID from other figs), **dashed 1.5**. (`geom_contour(breaks=c(log10(1/x),log10(x)))` for bias; `log(...)` for distortion.) The threshold values are ALSO marked on each colour key in their line style (solid/dashed vertical marks rising above the strip at those value positions) via `colorbar(..., sol=, dsh=)`; those values must be present in B_LIN/LD_LIN. Colour-key TITLE is on the left side (horizontal, `axis.title.y`), space-saving.
- Points: `dot` (grey15, small) = valid/reliable data cell; `cross` (×) = invalid/unreliable cell.

## Still PENDING on Fig 5
1. Header `<!-- -->` comment still says "empirical bias = theta_pool-theta_sim from pool_runs" — STALE, update to score+FIM.
2. Caption in a .md (not written yet).
3. Colour-key high-end labels can crowd slightly; acceptable now.

## Range check (2026-07-08, noise 0.1, Info distortion) — LD_LIN edges CONFIRMED good
log-det ln median → det: IR 0.06→1.06, R 0.78→2.2, MR 2.20→9, NR 15.0→3e6. Each algo resolves into several
distinct bands (IR at white, clean gradient to ∞). Keep LD_LIN as is.

## AS-BUILT: affine & Frobenius variants (2026-07-09) — the build spec below (kept for history) was WRONG on data availability
Both files built and render clean. **Bias panel A is IDENTICAL** (score+FIM). Both metrics are MAGNITUDES ≥ 0
(0 = self-consistent) → SEQUENTIAL white→red, NO blue/sign; `bandcols()` handles one-sided when the first edge is
0 (findInterval(0,breaks)=1 → nb=0 → all red from white). `*_BREAKS = the linear edges directly` (NOT log),
`m = the value directly`. Additivity does NOT hold (dropped the additivity note).

### CORRECTED FACTS (the original spec assumed Sample/Correlation had pre-computed affine/eigenvalue — they DON'T)
- Pre-computed `Affine_Invariant_Distance_*` and `Eigenvalue_Spectrum_*` exist ONLY for `Information` and
  `Gaussian_Fisher` distortion. **NOT for Sample/Correlation.** So the 3-row info/sample/correlation decomposition
  cannot use pre-computed scalars.
- What IS there: the full 6×6 matrices `Likelihood_{Information,Sample,Correlation}_Distortion` (value_row/value_col
  0..5, probit=mean + 5 bootstrap quantiles). Both variants ASSEMBLE the central (probit=mean) matrix per cell and
  compute the metric from it: affine = `sqrt(Σ(log λ_i)²)`, Frobenius = `sqrt(Σ_ij (D_ij − δ_ij)²)`.
- **KEY SUBTLETY (verified):** metric-of-the-mean-matrix ≠ mean-of-per-replicate-metric (the pre-computed scalars,
  e.g. log_Det, ARE per-replicate-averaged). They MATCH for IR/R/MR (low variance) and DIVERGE for NR/MNR (e.g. NR
  nch100 int1 log-det: from-mean-matrix 7.73 vs pre-computed 4.22; affine 8.47 vs 9.97). This is a nonlinearity
  (Jensen), not an error. We use the mean matrix ("distortion of the typical matrix"). Justification: the contour
  FILL uses one value per cell, the bootstrap CI is NOT used for these plots (confirmed w/ Luciano), so a single
  central summary is fine. The median (q=0.5) matrix was tried and gives MORE indefinite cells (worse), so mean it is.
- **Metric trade-off, now SHOWN in the figures:** affine takes logs → **undefined (NaN) for indefinite central
  matrices** → those NR cells are GREY. Frobenius needs no log → **always finite** → same cells are coloured (large)
  and marked with an **×** (non-PSD). Frobenius also has far better dynamic range for NR/MNR (0→10⁴→∞) so they
  resolve into structure instead of saturating; affine's NR/MNR (~5–15) saturate to dark red.

### Edges as built
- Affine: `AFF_LIN = c(0,0.05,0.1,0.2,0.3,0.5,0.7,1,1.5,2,3,5,8,13,20,Inf)`. Lines: solid 0.2 / dashed 0.5. Grey =
  `!is.finite(m)` (indefinite). Legend "affine distance", rows "affine information/sample/correlation".
- Frobenius: `FROB_LIN = c(0,0.1,0.2,0.3,0.5,0.7,1,1.5,2,3,5,10,30,100,300,1e3,1e4,Inf)`. Lines: solid 0.2 /
  dashed 0.5. × on `indef` cells (finite but non-PSD). Legend "Frobenius norm", rows "Frobenius information/…".
- Ranges (median from-mean-matrix, noise 0.1): affine Info IR 0.26 / R 0.59 / MR 1.19 / NR 8.1 / MNR 8.8;
  Frobenius Info IR 0.27 / R 0.70 / MR 1.74 / NR 115 / MNR 183. NR indefinite cells: 3 (Info), 3 (Sample), 1 (Corr).

### Story decision still open
affine/Frobenius keep no over/under SIGN (all red); log-det keeps the sign (blue=under=safe) but cancels/fails at
near-singular cells (grey). Likely: **log-det as the main panel B (keeps the safe/unsafe sign)** + a supplement with
Frobenius (robust, no grey, shows the NR/MNR magnitude cleanly). Affine is the weakest of the three here (greys the
very cells that matter most). NOTE: for a per-replicate-averaged (not mean-matrix) affine/Frobenius decomposition,
the program would need to emit those scalars for Sample & Correlation — deferred; the mean-matrix version is enough
for the metric-comparison question.

## ORIGINAL BUILD SPEC (SUPERSEDED — assumed pre-computed Sample/Correlation affine/eigenvalue that do not exist)
Copy `figure_5_master.Rmd` → `_affine.Rmd` / `_frobenius.Rmd`. **Bias panel A is IDENTICAL** (score+FIM). Only
panel B changes. Both metrics are MAGNITUDES ≥ 0 (0 = self-consistent) → SEQUENTIAL scale white→red, NO blue,
NO sign. So: `bandcols()` already handles one-sided (findInterval(0,breaks) → ctr=1, nb=0 → all red from white)
IF the first edge is 0. Set the distortion `LD_BREAKS = the linear edges directly` (NOT log) and `m = the value
directly` (no ln/exp). Labels = the value (fmt_big handles Inf). Additivity does NOT hold (distance, not additive)
so the info/sample/correlation rows are independent magnitudes (drop the additivity note).
- **AFFINE** (pre-computed, has bootstrap): components `Affine_Invariant_Distance_Likelihood_{Information,Sample,
  Correlation}_Distortion`. Ranges (Info): IR 0.10–0.65, R 0.42–0.78, MR 0.61–1.49, NR 5–13. Edges e.g.
  `AFF_LIN = c(0, 0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5, 8, 13, Inf)`. Self-consistency lines ~ affine 0.2
  (solid) / 0.5 (dashed) (single values, no reciprocal). Singular flag: `!is.finite(m)` (NaN), not m==0 (affine 0
  = genuine self-consistency). Legend title "affine distance", rows "affine information/sample/correlation".
- **FROBENIUS** (not pre-computed as scalar): compute per cell `sqrt(sum((lambda_i - 1)^2))` from
  `Eigenvalue_Spectrum_Likelihood_{Information,Sample,Correlation}_Distortion` (populated: 6 eigenvalues per cell,
  168–336 rows/algo). Assemble like the FIM in load-bias. Also a magnitude ≥ 0, sequential. Get its ranges, then
  design edges the same way. Row title "Frobenius ...", legend "‖D − I‖_F".
- The story trade-off to state when comparing: affine/Frobenius keep no over/under-confidence SIGN (all red); the
  log-det keeps the sign (blue = under = safe) but cancels/fails at near-singular cells (the grey cells). Decide
  which panel B tells the story best; likely log-det for the sign + a supp with affine to show it's robust.

## OLD NEXT (superseded by build spec above): affine & Frobenius variants of panel B
Make the SAME figure with panel B using **affine** and **Frobenius** instead of log-det, then decide which tells the story best.
- Affine: pre-computed `Affine_Invariant_Distance_Likelihood_{Information,Sample,Correlation}_Distortion` (with bootstrap). It is a **magnitude ≥ 0** (0 = self-consistent) → SEQUENTIAL scale (white→red), NO blue/sign. Doesn't cancel like log-det.
- Frobenius: not a pre-computed scalar; compute `sqrt(Σ(λ_i−1)²)` from `Eigenvalue_Spectrum_Likelihood_*` (eigenvalues). Also a magnitude ≥ 0, sequential.
- Trade-off to state: affine/Frobenius keep no over/under sign (all red); log-det keeps the sign but cancels at near-singular cells. The bias block (score+FIM, diverging) stays identical across all three.

## Broader figure plan (see `figure_5_PLAN.md` in this folder)
Arc converged to (eLife allows many figures; norm ~7 main + figure-supplements; use multipanel):
- **Bloque 1 (¿confiar?):** Fig 5 = bias (A) + distortion/sample/correlation (B) — THIS master figure. Fig 6 optional split.
- **Bloque 2 (precisión):** covariance/design (IR_channel_pooled_loss + IR_interval_loss + covariance_contour), using `Gaussian_Distortion_Corrected_Covariance`. Depends logically on distortion (honest covariance).
- **Bloque 3 (por qué fallan + mapa):** decomposition mechanism; the **N\* peak** result (`figure_5_IR_sample_peak_in_Nch.Rmd` — non-monotonic in N_ch, N\*=2σ²/G ∝ noise, verified; needs new low-N runs already partly landed); regime schematic `domains_schematic` (must redraw the static N_ch=30 line as the noise-dependent N\*).
- All figure_5 candidates were moved to `figures/in_progress/`; the audit/main-supp-drop sort is in `figure_5_PLAN.md`.

## Data-coverage note (Nch × noise, 1c2ae6f gaussian, for the mechanism figures)
noise 0.05: Nch 10,20,50,100 · 0.1: full 10–10000 · 0.2: 10,20,50,100 · 0.5: 10,20,50,100 · 1: 10,20,50,100,200,1000,10000 · 10: 10,100,1000,10000. Still worth: noise 0.5 at 200; noise 1/10 intermediate Nch to resolve the N* bells.
