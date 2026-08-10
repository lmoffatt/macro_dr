# D-4 — The calibration verdict

> **BANNER 2026-07-31 — read before quoting any `NMR` number from this file.** Every `NMR` below is
> the **defective build**: the non-recursive interval member with the `N·ms` interval-variance term
> missing, which is what ran between `a3e0a89` and `1f7138b` and produced the `macro_NMR` data. The
> corrected member is now `INR` (`../../_program/nomenclature.md`), it was re-run on 2026-07-31, and
> **no `INR` number exists in this file**. In particular §5's near-identity with `NR` measures the
> missing term, not the method, so it may not be restated as "the interval-mean conductance buys
> nothing without recursion". The `NR`-only numbers are unaffected.

> **Updated 2026-07-28, coverage corrected 2026-07-29.** Luciano's call: `433ed13` is the
> numerical-Fisher demo, so every number that fed the paper's verdict has been recomputed on
> **`1c2ae6f` + `87889e6` + `0ffbda7`** (the freeze, the D-0 fill, and the high-noise columns;
> multi-commit provenance per
> `../../_program/decisions.md` §4) against the **Gaussian Fisher**, which is the program's declared
> anchor. The 2026-07-14 version of this file, computed on `433ed13` against
> `Likelihood_Fisher_Covariance`, is preserved in git history. Its verdict **survives**; the
> magnitudes move by five to ten percent and one headline range narrows.
>
> Scripts are committed this time, not left in a scratchpad: `recompute/d4_headline_and_nch_trend.py`
> and `recompute/d4_distortion_envelopes.py`. Re-run them and this file must reproduce.
>
> **NMR is a SUPPLEMENT member** (Luciano, 2026-07-29), not dropped. The 2026-07-28 record said
> dropped, partly for "no literature attribution", which is false: `NMR` **is** the published
> `MacroINR`, the control of Comm Biol 2025. What survives is that it is measurably redundant with
> `NR` (§5), which is why it is a supplement member and not a body column.
>
> Abbreviations: emp = the covariance of the MLE cloud; Fisher = the covariance the likelihood reports
> from its own Gaussian Fisher information; D = the diagonal `Likelihood_Gaussian_Information_Distortion`
> at the pooled fit θ_pool; **over-confident** = the reported uncertainty is too tight (ratio > 1);
> **conservative** = too loose (ratio < 1).

---

## 0. What changed relative to the 433ed13 version

Two changes, both deliberate, and neither is a bug fix:

1. **The commit.** `433ed13` exists to show that the Gaussian Fisher ≈ the finite-difference Fisher
   (`D-0_freeze_and_rerun_scope.md`). It is a demo, not the production basis. The freeze is `1c2ae6f`,
   and the non-IR noise columns landed on `87889e6`.
2. **The anchor variable.** On the freeze the reported covariance is `Gaussian_Fisher_Covariance` and
   the distortion is `Likelihood_Gaussian_Information_Distortion`. The old file read
   `Likelihood_Fisher_Covariance` and `Likelihood_Information_Distortion`, which are the numerical
   family and exist on the freeze only as a handful of rows.

That the verdict is stable across both is itself a result, and it is the result `433ed13` was run to
produce. Say it once in Methods; do not quote `433ed13` numbers anywhere else.

---

## 1. The verdict: how much of the design space each method gets right

**This is the paper's quantity, and it is a map quantity, not a ranking one.** For every algorithm,
over every measured cell (N_ch × noise × the seven intervals) and every one of the four reported
parameters, ask whether the reported uncertainty is within ±15% of the truth.

| Algo | cells | (param, cell) points | **within ±15%** | D envelope | Reading |
|---|---:|---:|---:|---|---|
| **IR** | 525 | 2100 | **94%** | 0.65 → 1.70 | Calibrated over almost the whole measured plane; both excursions are in the few-channel corner |
| **R** | 224 | 896 | **70%** | 0.33 → 4.29 | Two-sided, and it is **the high-noise columns that rescue it**: over-confident at low noise, calibrated once instrumental noise dominates, which is region 3 of the usage map |
| **MR** | 84 | 336 | **20%** | 0.74 → 3.09 | Worse than R almost everywhere it is measured. The one-endpoint recursion misreports more, not less |
| **NR** | 84 | 336 | **0%** | 1.21 → 3.5 × 10⁴ | Never calibrated at any measured cell |

**Three commits, not two** (corrected 2026-07-29). The scan must read `1c2ae6f` + `87889e6` +
**`0ffbda7`**. The first pass of this table read only the first two and under-reported both IR (93%)
and, badly, R (51%): `0ffbda7` carries the high-noise columns out to a noise label of 10⁷, and R is
calibrated across most of them. Anyone re-deriving these numbers over two directories will reproduce
the wrong pair. The abstract in `../docs/manuscript-drafts/elife_paper.tex` already quotes the
three-commit figure.

Coverage is honestly unequal and must be stated with the number: IR spans N_ch 5 to 10⁴ at eleven
noise levels; R spans seven N_ch values at nine noise levels; MR and NR span four N_ch decades at
three noise levels, **all of them low**. So **MR's 20% and NR's 0% are measured only where the gating
signal is strong**, which is the regime that flatters them least; do not present those two as
plane-wide fractions without saying so. IR's 94% and R's 70%, by contrast, are measured across the
crossover.

The monotone 94 / 70 / 20 / 0 is the sentence the paper wants, and it replaces every "sole survivor"
phrasing (`../../_program/decisions.md` §6). **R's 70% is a better number for the paper than 51% was**,
because the region map's whole point is that there is a region where the classical error bar is honest.
A method that is never calibrated would contradict the map.

## 1a. The headline cell

N_ch 100, noise label 0.1, Δ = 0.1 τ, group_size 10, n_sims 10 000, 1000 MLEs in the cloud. This is
the cell Figure 2 plots.

| Algo | emp/Fisher ellipse-area ratio, geom-mean (per-pair range) | joint 4-param area ratio | emp/corrected |
|---|---|---:|---|
| **NR** | **11.41** (10.04 – 14.83) | 154.65 | 1.03 (1.01 – 1.07) |
| **MR** | **1.63** (1.40 – 1.97) | 3.07 | 0.90 (0.87 – 0.93) |
| **R** | **1.17** (1.08 – 1.32) | 1.56 | 0.95 (0.92 – 0.96) |
| **IR** | **1.01** (0.98 – 1.05) | 1.12 | 1.05 (1.01 – 1.07) |

Diagonal D at the same cell (k_on / k_off / i / N_ch):
NR 13.24 / 17.32 / 10.04 / 9.91 · MR 2.80 / 1.82 / 1.56 / 1.62 · R 1.39 / 1.51 / 1.16 / 1.16 ·
IR 1.05 / 1.13 / 1.00 / 0.90.

**The sandwich correction returns every method to ≈ 1** (0.87 to 1.07 across all four). That is the
machinery working, and it is the same claim the 433ed13 version made.

**Direction convention.** The ratio is `area(emp) / area(Fisher)` with `area = sqrt(det)`. Ratio > 1
means the cloud is wider than the reported ellipse, i.e. the likelihood under-reports the parameter
covariance, i.e. **over-confident**. `../../_program/machinery.md` §4 states this three times with two
different signs and has never been checked against the producer. **That check is still owed**, and
every boundary of the region map inherits it.

## 1b. The N_ch trend, which is new and is a mechanism result

Diagonal D at noise 0.1, Δ = 0.1 τ:

| Algo | N_ch 10 | 100 | 1000 | 10⁴ |
|---|---|---|---|---|
| **NR** (k_on) | 8.77 | 13.24 | 220.77 | 1583.07 |
| **MR** (k_on) | 2.28 | 2.80 | 2.89 | 2.87 |
| **R** (k_on) | 1.36 | 1.39 | 1.46 | 1.45 |
| **IR** (k_on) | 1.20 | 1.05 | 0.99 | 1.00 |

Three separate behaviours, and the contrast is the argument:

- **IR converges to calibration as channels increase.** Its residual error is the Gaussian closure, and
  the closure gets better with N_ch. This is what an approximation failing in the right direction
  looks like.
- **R and MR are flat.** More channels do not rescue them. `../decisions.md` gives the mechanism: the
  distortion splits into a sample part that goes to 1 with N_ch exactly as IR's does, and a
  correlation part that stays elevated and rises. R and IR do not differ in per-sample fidelity, they
  differ in temporal decorrelation.
- **NR diverges.** Two orders of magnitude over four decades, consistent with the `∝ N_ch²`
  conditioning claim in `D-3_novelty_claim.md`.

---

## 2. The two formerly contested cells

### 2.1 MR's direction: over-confident, and the old prose had it backwards

MR's emp/Fisher parameter-covariance ratio is **1.63** at the headline cell and never below 0.74
anywhere measured, i.e. > 1 almost everywhere, i.e. MR **under-reports** the parameter covariance and
is **over-confident**. The three prose copies said MR "overestimates variance", which reads as
conservative. **The copies inverted MR.** Delete that phrasing wherever it survives.

The category error behind the confusion is real and worth one sentence in Methods: "over-confident"
is about the **parameter covariance**; "overestimates variance" was about the **predicted per-interval
observable variance** `y_var`, a different object. The tidy story that MR overestimates the observable
variance because it drops a subtractive boundary term does **not** hold for the production algorithms:
that subtractive term lives only in the cut Taylor branch.

**SETTLED 2026-08-06, and this closes the §4 item.** The `y_var` comparison was run. In the production
path the boundary term IR *adds* to `gSg` is exactly the one it *removes* from `ms`, so the two cancel
and **at the same prior MR and IR predict the same observable variance**: on the figure-1 dumps, same
recording, `MR = IR = 1.048475625791748` at every interval where they still share a prior. Along a
recording they diverge, because the predictive variance divides the gain and the gains differ: MR then
runs +69% to +81% above IR. `VR` is strictly below both from any state, by `μᵀVar_j[gmean_ij|i]`
(0.378085 at that same interval). So the observable-variance direction for MR against IR is: equal at
equal state, above along a recording. Verified three ways (term by term against
`legacy/qmodel.h:4568-4619`; numerically on random fields, difference 0.0 and structural; and on the
dumps). Canonical: `../figures_build_plan.md:195-245`. The earlier "IR *adds* a term, which points the
other way" is withdrawn: it read one half of the cancellation.

### 2.2 IR's corner: two-sided, and the old number was one-sided

IR is two-sided in the few-channel corner. Over the full measured grid its D runs from **0.645**
(N_ch 10, noise 0.05, Δ 0.01 τ, on N_ch: conservative, the safe direction) to **1.701** (N_ch 5,
noise 0.1, Δ 1 τ, on k_off: over-confident). The retired "up to ~1.3 in the corners" both omitted the
conservative side and understated the over-confident one.

On the freeze the low side is **less** extreme than the 433ed13 version reported (0.645 against 0.503)
and the high side is **more** (1.701 against 1.484). Quote the freeze numbers.

---

## 3. The overconfidence factor: 10 to 15

The non-recursive overconfidence quoted in the abstract is the emp/Fisher **ellipse-area** ratio for
NR at the Figure 2 cell: per-pair **10.04 to 14.83**, geom-mean 11.41.

**Quote "10 to 15", not "10 to 16" and not "14 to 21".** The old 10-16 came from pooling NR with NMR
on `433ed13`; with NMR dropped and the anchor moved, the honest range is 10 to 15. The 14-21 figure
was never the ellipse-area factor at all: it is the single-parameter variance distortion on the two
rate directions, a different measure of a different object.

---

## 4. What this file still cannot settle

- ~~**The observable per-interval variance direction for MR against IR.**~~ **CLOSED 2026-08-06**, on
  the figure-1 per-interval dumps rather than the time-resolved ones: equal at equal prior, MR +69% to
  +81% along a recording, VR below both. See §2.1 above and `../figures_build_plan.md:195-245`.
- **The sign convention** (§1a). A code read against the producer, roughly an hour, and every region
  boundary depends on it.
- **The ~100-channel threshold (D-J).** It does not fall out of this table. IR's k_off distortion at
  noise 0.1 is 1.32 / 1.10 / 1.00 / 1.00 across N_ch 10 / 100 / 1000 / 10⁴, so at exactly 100 channels
  IR is still about 10% off, and Figure 6's closure boundary sits at N_ch ≈ 70-150 and moves with
  noise. **Channels or openings cannot be settled by measurement**: P_open is fixed at 0.5 in every
  run and never swept, so the two differ by a constant factor of two and the design cannot distinguish
  them. That is a definitional call, not a recompute.

---

## 5. NMR: measurably redundant with NR, and published — VOID as evidence about the method (2026-07-31)

> **The redundancy is a property of `NMR`, the defective build, and does NOT transfer to `INR`**
> (2026-08-01). Read the strike-through at the foot of this section before citing anything in it.

NMR was in the fill and is on `87889e6` at all three noise levels. Recomputed, it is
**numerically indistinguishable from NR**: envelope 1.226 to 35 738 against NR's 1.209 to 35 358, and
0% of 336 points within ±15% for both. It measures the same failure NR measures, at the same size,
with no literature attribution and no mechanistic role of its own.

**But the original "no attribution" reason was false** (corrected 2026-07-29): `NMR` is `MacroINR`,
the published control of Comm Biol 2025, whose systematic underestimation of the evidence for schemes
with conformational intermediates is that paper's central methodological claim. So what is measured
here is redundancy with `NR`, not absence of a role.

~~Read the other way, the redundancy is a **result**: the interval-mean conductance buys nothing without
recursion, so recursion is the step that matters at the bottom of the ladder and interval conditioning
is the step that matters at the top.~~ **STRUCK 2026-08-01. The redundancy was the bug.** Everything
above measures `NMR`, the build missing the `N·ms` interval-variance term, and the section title's
"measurably redundant with NR" holds for `NMR` and for nothing else. The corrected `INR`
(`figures/data/1f7138b/`, 17 cells at `nsim` 10⁴) separates from `NR` in both moments: median |bias|
in `N_ch` 0.003 against 0.099 log10, and `k_off` information distortion at noise 0.1 / N_ch 10⁴ of
22.2 / 1.00 / 22.1 (total / sample / correlation) against 78.9 / 47.2 / 21.2. **The interval-mean
conductance buys the first moment and the per-sample fidelity**; what it leaves is pure correlation
distortion, and that is recursion's job. So the two steps are orthogonal margins acting on different
moments, not two ends of one ladder. Recompute: `recompute/d3_interval_vs_recursion_2x2.py`; roster
consequence: Q-5 in `../decisions.md`. `../../_program/decisions.md` §2 holds the options.

**Still owed from this section, and not answered by that script:** §3's ×10–15 overconfidence factor
and the `N_ch²` conditioning scaling are standard-error and Fisher-spectrum statistics, both computed
against `NMR`. Rerun D-4's own recipe on `macro_INR` before either enters the paper.

---

## 6. What the paper says now

The window-ignoring approximations over-state what the data know (NR by one to four orders of
magnitude and growing with channel count, MR by ~1.6×, R by ~1.2×, all over-confident), the interval
likelihood is within ±15% over 94% of the measured design space, and where it departs it mostly
departs in the conservative direction. The sandwich correction returns all of them to ≈ 1, which is
what makes the diagnostic a measurement rather than a verdict.

---

## 7. Recompute provenance

Data: `projects/eLife_2025/figures/data/{1c2ae6f,87889e6,0ffbda7}/figure_3_G_nch_{N}_nsim_10000_macro_{ALGO}_noise_{S}_*.csv`,
skip row 1 (the engine git hash). Filters, exactly as the figure notebooks:

- **Cloud**: `_mle_cloud_runs.csv`, `variable == Model_Parameters_Hat`, `statistic == value`,
  **`group_size == 10`** (the file mixes 10 and 100 and pooling them averages garbage),
  `interval_in_tau == 0.1`; one MLE per (simulation_index, sample_index, sub_index, segment_index).
- **Model covariance**: `_battery_pool_G.csv`, `variable == Gaussian_Fisher_Covariance`,
  `probit == mean`, `statistic == value`, `calculus == primitive`, `operation == probit`, indexed by
  (`param_index`, `param_col`), scaled by `cov_scale = 1/group_size = 1/10`.
- **Distortion**: `variable == Likelihood_Gaussian_Information_Distortion`, same filters, diagonal
  only (`param_index == param_col`).
- Parameters reported: indices 0 `k_on`, 1 `k_off`, 2 `unitary_current`, 5 `Num_ch_mean`.
- `area(S) = sqrt(prod(eigvalsh(S)))`.

Traps respected: "noise 0.1" is the label, `Current_Noise = 1e-4`; the interval sweep lives *inside*
each battery file (seven values); n_sims is held at 10 000 throughout, so no Jensen mixing.

Scripts: `recompute/d4_headline_and_nch_trend.py`, `recompute/d4_distortion_envelopes.py`.
Run 2026-07-28.
