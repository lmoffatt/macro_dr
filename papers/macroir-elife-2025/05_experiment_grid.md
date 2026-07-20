# Paper 2 (MacroIR / eLife 2025) — Experiment Grid Spec

> Updated: 2026-07-20. Detail layer beneath `00_master_plan_v2.md` §6. Defines the regime sweep for the usage map.

## 1) Control variables
From the two-state `scheme_CO`:
- **N_ch** — effective number of channels. Map axis 1.
- **instrumental noise** — map axis 2, and the axis that carries the two crossovers (§2a). Swept as the dimensionless label of `decisions/D-2_parameter_units.md` (label = 10·ν, ν = `Current_Noise`·k_off/g²; the dispatcher maps `Current_Noise = label/1000`).
- **interval Δ / τ** — Δ = n_samples / fs; τ = 1 / max_k|Re(λ_k)| of the generator Q at the agonist condition. Carried *inside* every cell as 7 values of `interval_in_tau`, not as a separate map axis. Standardize how τ is recorded in run metadata.
- **K_off** — **not swept.** Fixed at 100 on disk (`figure_3_mle_G.macroir:46`), the output path has no K_off component so cells would silently overwrite each other, and building the axis was costed at ~35,000 CPU-hours. Any document still promising an N_ch × K_off map is stale.

## 2) The map
- Primary: **N_ch × noise** heatmaps (bias, distortion), per method, with the interval swept inside each cell.
- The map's job is to place each method's validity boundary against the two **predicted** crossovers below, so the reader can read off the cheapest adequate method for their recording.

## 2a) The two crossovers that structure the noise axis

The instrumental noise has two reference scales, separated by exactly a factor N_ch: the gating variance of one channel, and that of the whole population. In label units, valid across the production interval range (Δ ≤ τ, where the within-interval averaging does not yet suppress the gating term):

| Boundary | Label | Meaning |
|---|---|---|
| A/B | `10 · interval` | instrumental noise reaches the single-channel gating scale |
| B/C | `10 · N_ch · interval` | instrumental noise reaches the total gating noise |

Band A is gating-dominated and finely resolved (expect IR); band B is gating-dominated and coarsely resolved (expect R may suffice); band C is instrumental-dominated (expect LSE and the non-recursive members to be adequate, and the interval filter to buy nothing). Derivation and the coverage audit: `00_master_plan_v2.md` §1a.

## 2b) Coverage

- **Before 2026-07-20:** labels {0.05 … 10}. Band A at the headline cell, reaching band B at short intervals, grazing the B/C boundary only at N_ch 10 / Δ = 0.01 τ / label 10. Band C was systematically absent at canonical N_ch and interval, which is why the map favoured the gating-aware likelihoods by construction.
- **Dispatched 2026-07-20:** `N_NOISE="100 1000 10000 100000"` paired index-wise with `NCHS="10 100 1000 10000"`, i.e. `label = 10·N_ch`, the B/C boundary at Δ = τ, with the interval axis carrying each cell down into band C. `nonlinearsqr` via `dispatch_figure_3_LSE.sh`; {IR, R, MR, NMR} via `dispatch_figure_3_G.sh`. `n_sims = 1000`, `GROUP_SIZE="10 100"`.
- **Not covered:** NR anywhere in bands B/C; the interior of band B at canonical N_ch (the fill jumps from label ≤ 10 to label = 10·N_ch, leaving a gap for N_ch ≥ 100).

**Three divergent label→value maps exist in ops** and any change to the noise axis must touch all three or the divergence deepens: 11 levels in `dispatch_figure_3_G.sh:143-155` and `dispatch_figure_3_LSE.sh:148-160`, 4 levels in `dispatch_figure_2.sh:110-115` and `dispatch_figure_3.sh:167`, 9 levels in `run_figure_3_G_local.sh:125-134`. Figure 2 currently cannot express the band-C labels at all.

## 3) Model
- `scheme_CO` only (two states). `scheme_CCO` (3-state) is out of scope (>2 states = later program component).

## 4) Methods
- **Method zero**, off the lattice: `nonlinearsqr` (display `LSE`), classical nonlinear least squares on the mean current, `family_approximation = 2`, its own dispatcher `dispatch_figure_3_LSE.sh`.
- **The five likelihood approximations:** `NR`, `NMR`, `R`, `MR`, `IR`. Drop the Taylor variance-correction variants (MNRV/MRV/IRV — cut, taylor=false). Naming: MNR → NMR.
- The label `nonlinearsqr` must appear verbatim end to end (`.macroir` axis label → CSV `algorithm` cell → R `ALGOS` entry); any mismatch silently drops the rows, the same failure class as the MNR/NMR bug.
- LSE needs two parameters `Fixed` that the others leave free: `unitary_current` (the amplitude ridge, N and g are unidentifiable from the mean alone) and `Current_Noise` (a pure-variance direction, so its score row is identically zero and the Fisher would be singular). Details in `theory/macroir/notes/nonlinearsqr_lse_plan.md` §7.

## 5) Per-cell protocol
For each cell and algorithm:
1. Fixed non-stationary experiment (single concentration jump / pulse + washout).
2. Simulate n_simulations datasets at θ* (reuse the same samples across intervals/algorithms where possible).
3. Compute the score, the Gaussian Fisher H, and per-interval predictive mean/variance (residuals). Numerical FD-Fisher only to gauge H.
4. MLE / Gauss-Newton local max → empirical parameter covariance.
5. Reduce to the scalar metrics below.

## 6) Scalar metrics per cell (per method)
- score bias norm (max |mean(score_i)/sd(score_i)|) and the projected DIB;
- distortion matrix C = H^(−1/2) J H^(−1/2): det (→ ½ log det C) and its correlation vs sample/geometric split;
- residual mean, var, max ACF (lag>0);
- empirical-vs-sandwich covariance ratio.

Always store the reduced table (so the maps plot cheaply); raw data stays out of git (see `09_carve_plan.md`).

**Two of these do not mean for LSE what they mean for the five.** State it wherever an LSE cell is plotted beside a likelihood cell:
- **r̄²_std ≡ 1 is a tautology for LSE**, since Σ r_std² = n identically once σ̂² = SSE/n. Any residual-variance panel reads "calibrated" for LSE by construction and carries no information.
- **F = Var(score) requires homoscedasticity**, which LSE assumes and the exact simulator violates (the gating variance N·gSg peaks in the transient and vanishes at the plateaus). The identity therefore *breaks* for LSE, and the size of the break is a result, not an artifact: it is the classical method's overconfidence, measured. Related: the shared (n/SSE) factor makes LSE's per-interval scores non-martingale through a rank-1 common mode.

**n_sims must be uniform within a panel.** All these scalars carry a Jensen bias in n_sims; the band-B/C fill is 1000 and the band-A cells are 10000.

## 7) Exact grid values
The definitive values live in the `ops/` dispatch configs (`figure_2*.macroir`, the `figure_3*` Gaussian-Fisher dispatch) — the source of truth once the Gaussian rerun is final. Do not hardcode a second copy here.
