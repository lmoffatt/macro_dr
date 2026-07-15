# D-4 — The ranking verdict

> Prepared for Luciano as a 10-minute yes/no. Task **B-4** of `01_writing_plan.md` §3/§5.
> Every number below was either read from a file (with `% src:`) or recomputed from the production
> CSVs under `projects/eLife_2025/figures/data/433ed13/` by re-running the exact aggregation the
> figure notebooks use. The recompute scripts are transcribed at the end (§7) so the arithmetic is
> auditable. Nothing outside this file was edited; the `02_decision_log.md` entry is written only
> after you say yes.
>
> Abbreviations: emp = empirical (the covariance of the MLE cloud); Fisher = the covariance the
> likelihood reports from its own Fisher information; D = the `Likelihood_Information_Distortion`
> scalar (per-parameter, evaluated at the pooled fit θ_pool); "over-confident" = the reported
> uncertainty is **too tight** (emp/Fisher > 1); "under-confident/conservative" = too loose
> (emp/Fisher < 1). The five algorithms and the NMR↔MNR data/display relabel are in
> `docs/figure_provenance.md` §3.

---

## 0. The one thing to decide

Approve the table in §1 as the paper's verdict. It changes exactly two cells relative to the copies
now in the repo (**MR's direction** and **IR's corner bound**), plus one number that was loose
(**R**), and it settles the overconfidence factor at **10 to 16** (§4). If you say yes, §6 lists the
edits an agent then makes to the four drifted copies; you approve the table here, not the edits.

---

## 1. The recomputed ranking (the canonical quantity is the parameter-space distortion)

**One quantity runs the whole table:** the ratio of the parameter uncertainty the estimates *actually
have* (emp) to the parameter uncertainty the likelihood *reports* (Fisher). It is what Figure 2 / S3
plot as an ellipse-area ratio and what Figure 5 plots as the scalar `Likelihood_Information_Distortion`
D. **I verified these two are the same object** (same sign, same magnitude at the shared cell — §2.2),
so "variance distortion", "overconfidence factor", "information distortion" and the S3 area ratio are
four names for one number. Headline cell = N_ch 100, noise label 0.1, Δ = 0.1 τ, group 10, 10 000
recordings (the cell Figure 2/S3 plot; `figure_S3.Rmd:45-50`). Ranges are across the whole grid.

| Algo | Bias | Parameter-covariance distortion (emp/Fisher) — direction and size | Verdict |
|---|---|---|---|
| **IR** | none `% src: Figure_S3_caption.md#L7 (truth cross on cloud); recompute §7` | **≈ 1** at the headline cell (**0.93–1.14** over the 4 params `% src: recompute, 433ed13 …macro_IR…battery_pool.csv#Likelihood_Information_Distortion`); two-sided only in the few-channel corner: **down to ~0.5** (under-confident, conservative) and **up to ~1.5** (over-confident) `% src: recompute, IR D envelope min 0.503 @ N_ch10 Δ0.01τ N_ch-param / max 1.484 @ N_ch10 Δ1τ k_off` | Sole survivor; calibrated across the practical regime; self-consistent (±15 %) over **0.88 → 0.99 → 1.00** of the plane as noise rises `% src: figure_5_distortion_algo_grid.Rmd#L23` |
| **R** | mild bias in the amplitude pair (N_ch, i), predicted by the distortion theory `% src: Figure_S3_caption.md#L7` | **over-confident, ~1.2–1.4** (geom-mean 1.26; per-pair 1.18–1.42) `% src: Figure_S3_caption.md#L7; recompute §7` — **not "factor ~2"** | Usable; correlation distortion at short Δ |
| **MR** | mild bias in the amplitude pair (N_ch, i), predicted by the distortion theory `% src: Figure_S3_caption.md#L7` | **over-confident, ~1.5–2.1** (geom-mean 1.80; per-pair 1.52–2.07) `% src: Figure_S3_caption.md#L7; recompute §7` — **same direction as NMR, not the opposite** | Strawman intermediate; the one-endpoint recursion still misreports the information |
| **NR** | biased `% src: Figure_S3_caption.md#L6 ("Fisher is overconfident for the non-recursive pair"); master_plan §4` | **over-confident, ~10–16**, and **grows with N_ch**: 9–18 at N_ch 100 → 530–1590 at N_ch 10 000 `% src: recompute, 433ed13 …macro_NR…battery_pool.csv#Likelihood_Information_Distortion` | Fails; often does not reach the MLE `% src: master_plan §4 (not re-tested here)` |
| **NMR** | none `% src: Figure_S3_caption.md#L6` | **over-confident, ~10–16**, grows with N_ch (same as NR) `% src: recompute; Figure_S3_caption.md#L6` | Only other unbiased method; possible speed niche (non-recursive, parallelizable) `% src: master_plan §4` |

Mechanism sentence (unchanged, and correct as *mechanism*): every member makes the same two Gaussian
approximations; they differ only in how much of the acquisition interval the hidden state is
conditioned on (no endpoints → one → both — the endpoint ladder, `nomenclature.md`), and the
distortion falls monotonically along that ladder `% src: discussion_plan.md#L38`. That is why NR/NMR
(no conditioning) are ~10–16, R (one-sided recursion) is ~1.3, MR (one endpoint) is ~1.8 — worse than
R, hence "strawman" — and IR (both endpoints) is ≈ 1.

**Direction convention, verified once against the producer** (this is the thing every verdict hangs
on): the S3 ratio is `area(emp cov) / area(Fisher cov)` with `area = sqrt(det)`
(`figure_S3.Rmd:120-122`). Ratio **> 1** means the cloud is **wider** than the Fisher ellipse, i.e.
the likelihood **under-reports** the parameter covariance, i.e. **over-confident**. Every non-IR
member sits at ratio > 1, so **the whole family other than IR is over-confident in parameter space**;
they differ only in how much.

---

## 2. The two contested cells, resolved

### 2.1 MR's variance direction — the copies have the sign backwards

**Claim in all three prose copies:** MR "overestimates variance" (`00_master_plan_v2.md:91`,
`02_decision_log.md:26`, `discussion_plan.md:32`). Read in the table's own "Variance distortion"
column that reads as **under-confident** (reports too-*wide* uncertainty).

**Data:** MR's emp/Fisher parameter-covariance ratio is **1.5–2.1** (geom-mean 1.80), i.e. **> 1**,
i.e. the cloud is wider than MR's Fisher ellipse, i.e. MR **under-reports** the parameter covariance
= **over-confident** `% src: Figure_S3_caption.md#L7; recompute §7`. That is the **same direction as
NMR** ("underestimates variance", which the copies got right), just milder (1.8 vs ~12). So in the
column the table actually reports (parameter covariance), MR's cell should read like NMR's, not its
opposite. **The copies inverted MR.**

**Is the observable-vs-parameter category error real for MR?** Partly, and it is the source of the
confusion, but it does not rescue "overestimates variance":

- The two words name two genuinely different variances. "Over-confident" (1.8) is about the
  **parameter covariance**. "Overestimates variance" was meant to be about the **predicted per-interval
  observable output variance** `y_var` — a different object, `results_plan.md:29-36`.
- The producer confirms MR and IR *do* compute a different observable `y_var`: they split on the
  averaging flag, `av = 1` (MR) vs `av = 2` (IR), and differ in the `gSg` term
  `% src: legacy/qmodel.h#L3856-3869 (MR: gSg = gmean_iᵀ·SmD·gmean_i) vs #L3826-3838 (IR: + P_mean·(gtotal_ij⊙gmean_ij)·u)`.
- **But the tidy story "MR overestimates observable variance because it drops a subtractive boundary
  term N·γᵀΣγ that IR keeps" does not hold for the production algorithms.** That clean *subtractive*
  term `− N·zeta·sSg²` lives only in the `variance_correction = true` **Taylor branch**
  `% src: legacy/qmodel.h#L3789`, which is the IRT/MRT path **cut from the paper** (`taylor = false`;
  see `02_decision_log.md:16`). In the production path IR **adds** `P_mean·(gtotal_ij⊙gmean_ij)·u`,
  which would make IR's observable variance *larger*, i.e. MR would *under*-estimate the observable
  variance — the opposite of the reconciliation `results_plan.md:33` proposes.

**Conclusion for MR.** In the paper's quantity (parameter covariance) MR is **over-confident, ~1.8**;
this is what goes in Table 1. Drop "overestimates variance". Keep the endpoint-ladder mechanism
(MR drops boundary conditioning that IR keeps, so it claims information it does not deliver). **Do not
assert that MR overestimates the observable variance** — the production code points the other way, and
settling the observable-variance direction needs the per-interval `y_var` comparison in §5.

### 2.2 IR's corner bound — ≤ 1.3 vs ~0.5 is NOT a category error; it is a one-sided report

**Claim in copies 1–3:** IR distortion "≈ 1 in the Gaussian regime, up to ~1.3 in the corners"
(`00_master_plan_v2.md:89`, `02_decision_log.md:24`, `discussion_plan.md:30`). **Observation from
Figure 5:** IR distortion "down to ~0.5 in the low-N_ch corner" (`results_plan.md:42`;
`figure_5_distortion_algo_grid.Rmd:24`).

**These are the same quantity, measured at different corners — not two different quantities.** Proof:
Figure 5's D is `Likelihood_Information_Distortion` read diagonally at θ_pool
(`figure_5_distortion_algo_grid.Rmd:80-85`); the S3 ratio is `area(emp)/area(Fisher)` at the same
θ_pool anchor. At the shared headline cell they agree in sign and size (IR: S3 area ratio 0.98; D per
param 0.93–1.14). So the suspected **observable-variance-vs-parameter-covariance category error is NOT
present here** — both "1.3" and "0.5" are the parameter-space information distortion.

What actually happened: the ranking reported **only the over-confident side** and rounded it (IR's
over-confident excursion is D ≈ 1.2 for k_on and reaches **1.48 for k_off** at N_ch 10) while
**omitting the under-confident side** (D **down to 0.503** for N_ch at N_ch 10, Δ 0.01 τ)
`% src: recompute, IR D full envelope 0.503–1.484 across 433ed13 grid, §7`. IR is **two-sided in the
few-channel corner**, and the missing side (~0.5, conservative) is the *safe* direction and the better
sentence for the paper (`results_plan.md:40-42`).

**Recomputed IR cell:** ≈ 1 across the practical regime (0.93–1.14 at the headline cell); in the
few-channel corner it runs from **~0.5** (under-confident, conservative — N_ch, N_ch = 10) to **~1.5**
(over-confident — k_off, N_ch = 10). The old "up to ~1.3" both **omits the 0.5 side** and mildly
**understates the over side** (1.48).

---

## 3. How badly the four copies drifted (diff against §1)

Symbols: **✗** = disagrees with the data; **~** = incomplete/loose but not wrong-signed; **✓** = agrees.

| Cell | `00_master_plan_v2.md §4` (L87-93) | `02_decision_log.md` (L23-28) | `discussion_plan.md` (L28-34) | `elife_paper.tex` fill-hint (L65) | Recomputed (§1) |
|---|---|---|---|---|---|
| **IR** | "up to ~1.3 in the corners" **~** (one-sided) | "≤ ~1.3 in the corners" **~** | "up to ~1.3 in the corners" **~** | "calibrated default" **✓** (no number) | ≈1; corner 0.5–1.5, two-sided |
| **R** | "~factor 2 residual" **✗** | "~factor-2 residual" **✗** | "residual factor ~2" **✗** | "R ~factor 2" **✗** | ~1.2–1.4 (over-conf) |
| **MR** | "overestimates" **✗** (sign) | "overestimates variance" **✗** | "overestimates variance" **✗** | "MR strawman" (no dir); MR absent from the "NR/NMR inflation" list **~** | ~1.5–2.1, **over-conf** |
| **NR** | "inflation ∝ N" **✓** | "inflation ∝ N" **✓** | "inflation ∝ N" **✓** | "NR/NMR inflation" **✓** | ~10–16, grows w/ N_ch |
| **NMR** | "underestimates" **✓** (= over-conf) | "underestimates variance" **✓** | "underestimates" **✓** | "unbiased-but-underdispersed" **✓** | ~10–16, over-conf |
| overconf. factor | — | — | — | — | **10–16** (§4) |

Read-out: **all four copies agree on NR and NMR and are wrong or silent on MR; three of the four
overstate R as "factor 2"; and the IR corner is one-sided in every copy that gives a number.** The
`.tex` fill-hint is the least wrong only because it commits to the fewest numbers. The one genuinely
inverted cell is **MR**; the rest is loose or incomplete, not sign-flipped.

Note also a **bias-column drift** not in the two flagged cells: copies 1–3 mark R and MR bias as
"none", but `Figure_S3_caption.md:7` says R and MR carry a real bias in the amplitude pair (N_ch, i)
that the distortion theory predicts. Minor, but the "none" is an oversimplification; §1 records it.

---

## 4. The overconfidence factor: 10 to 16 (14 to 21 is wrong)

**Right number: 10 to 16.** `% src: Figure_S3_caption.md#L6` and `abstract_draft.md#L158`
(the correction line). It is the emp/Fisher **ellipse-area** ratio (sqrt-det over each parameter pair)
for the non-recursive pair at the Figure 2/S3 cell. **Recomputed directly from the CSVs** and it
reproduces the caption: NR per-pair **9.95–14.18**, MNR per-pair **9.68–15.69** → combined **≈ 10 to
16** `% src: recompute §7, 433ed13 …macro_{NR,NMR}…{mle_cloud_runs,battery_pool}.csv`. The
sandwich-corrected ratio is ≈ 1 for all `% src: recompute §7 (emp/corrected 0.94–1.08)`, confirming
the caption's last clause.

**Why 14 to 21 is wrong.** `abstract_draft.md:162` (the older line, and `discussion_plan.md`/
`abstract_draft.md:124` attribute it to a `figure_2_pub` script) does **not** reproduce as the
Figure 2/S3 overconfidence factor on this cell under any faithful aggregation: the per-pair area ratio
is 9.68–15.69, and the *joint* four-parameter area ratio is ~145–154 `% src: recompute §7`. The only
construction that lands near "14 to 21" is the **single-parameter variance** distortion D for NR/NMR
on the two rate directions (k_on ≈ 11–14, k_off ≈ 18–23) `% src: recompute, 433ed13 battery_pool
Likelihood_Information_Distortion diagonal` — a *different* measure (one parameter's variance, not the
ellipse area). So 14–21 is a real number about a different object; the abstract's sentence is about the
**Figure 2/S3 ellipse-area factor, which is 10 to 16.** Quote 10–16 and cite Figure S3.

This is the one item B-4 also records in `02_decision_log.md` — pending your yes.

---

## 5. What the data cannot settle here, and the run that would

The **observable per-interval variance** direction for MR vs IR (the other half of §2.1's
reconciliation) cannot be settled from the parameter-covariance CSVs, and the production code points
against the "MR overestimates observable variance" story (§2.1). The clean settle is a direct
comparison of the predicted `y_var` per interval for MR vs IR on the same ensemble — the columns are
already in the **Figure 3/4 time-resolved dumps** (`figures/data/figure_3_time_dlik_{MR,IR}.csv`,
`docs/figure_provenance.md` §6), no new run needed; it is an `awk` over the `y_var` column at matching
intervals (dedup the ×2 row duplication first, `figure_provenance.md` §6 trap 1). Until then, the
paper should **not** claim MR overestimates the observable variance; the parameter-space verdict
(over-confident, 1.8) stands on its own and is all Table 1 needs.

Traps respected (`docs/figure_provenance.md`): the S3/Figure-2 data are run 2 = `figure_3_mle.macroir`
→ `433ed13/figure_3_nch_*`, **not** a "figure_2" file (§1 naming trap); "noise 0.1" is the label,
`Current_Noise = 1e-4` (§1 trap 3); the mle_cloud file **mixes group_size 10 and 100** and I filtered
to **group_size = 10** to match Figure 2 (`project_cloud_group_size_column`).

---

## 6. Recommendation

**Approve the §1 table as the verdict.** Concretely, yes means:

1. **MR** cell: over-confident, ~1.5–2.1 (like NMR, milder). Delete "overestimates variance"; keep the
   endpoint-ladder mechanism.
2. **IR** corner cell: two-sided, ~0.5 (conservative) to ~1.5 (over-confident) at few channels, ≈ 1
   elsewhere. Replace "up to ~1.3".
3. **R** cell: ~1.2–1.4 over-confident. Replace "factor ~2".
4. **Overconfidence factor**: 10 to 16 everywhere; strike 14 to 21.
5. Column header: "parameter-covariance distortion (emp/Fisher)", not the ambiguous "Variance
   distortion".

**What changes in the paper if you approve:** the bidirectionality claim — the abstract's and
Discussion's central directional sentence — moves onto a carrier the data support: *the
window-ignoring approximations over-state what the data know (NR/NMR ~10–16×, R ~1.3×, MR ~1.8×, all
over-confident), and the interval likelihood, where it departs, under-states it (IR down to ~0.5,
conservative, at few channels)*. That is one directional story, sign-checked against the producer,
and it retires the "MR overestimates / MR sign discrepancy" flag that currently blocks the abstract
(`abstract_draft.md:156-164`) and the Discussion (`discussion_plan.md:36`).

---

## 7. Recompute provenance (auditable arithmetic)

All from `projects/eLife_2025/figures/data/433ed13/figure_3_nch_100_nsim_10000_macro_{ALGO}_noise_0.1_*.csv`,
skip = 1 (git-hash provenance row), filters exactly as `figure_S3.Rmd`: cloud from
`_mle_cloud_runs.csv` (`Model_Parameters_Hat`, statistic `value`, Num_ch 100, group_size 10, interval
≈ 0.1); model covariances from `_battery_pool.csv` (`Likelihood_Fisher_Covariance` /
`_Distortion_Corrected_Covariance` / `_Information_Distortion`, probit `mean`, param_index == param_col
for the diagonal D), all × cov_scale = 1/group_size = 1/10; `area(S) = sqrt(prod(eigen(S)))`.

**emp/Fisher ellipse-area ratio, geom-mean over the 6 pairs (min–max):**
NR 11.12 (9.95–14.18) · MNR 11.83 (9.68–15.69) · R 1.26 (1.18–1.42) · MR 1.80 (1.52–2.07) ·
IR 0.98 (0.93–1.02). emp/corrected ≈ 1 for all (0.94–1.08).

**IR `Likelihood_Information_Distortion` D across the full 433ed13 grid** (N_ch ≥ 10, noise 0.1/1/10,
7 intervals, 4 params): min **0.503** (N_ch 10, Δ 0.01 τ, N_ch param), max **1.484** (N_ch 10, Δ 1 τ,
k_off); headline cell (N_ch 100, Δ 0.1 τ, noise 0.1): N_ch 0.93 / k_off 1.14 / k_on 1.04 / i 1.02.

**NR/NMR D grows with N_ch** (k_on direction): N_ch 100 → 13.5 / 11.3; N_ch 1000 → 216 / 229;
N_ch 10 000 → 1590 / 1663.

**Joint 4-parameter area ratio** (why 14–21 is not the area factor): NR 145, MNR 154, R 1.79,
MR 3.58, IR 1.07.

Scripts: `s3_ratios.R`, `dist_grid.R`, `ir_env.R`, `joint4.R` (run 2026-07-14; kept in the B-4
scratchpad, not committed).
