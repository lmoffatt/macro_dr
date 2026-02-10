# Paper 2 (MacroIR / eLife 2025) — Master Plan (English)

## 0) Goal in one sentence

Write a clear, reproducible paper that answers: **when should a reader use MacroIR vs MacroMR vs simpler MacroR-like approximations**, and how can they *diagnose* when an approximation is failing, using **score/FIM and residual-based tests**.

## 1) Scope (what this paper is / is not)

### In scope

- **Macroscopic** (many-channel) current traces observed as **interval averages**.
- Likelihood approximations family currently used in this repo (naming used in scripts):
  - `NR`, `R`, `MNR`, `MR`, `MNRV`, `MRV`, `IR`, `IRV`
- **Correctness/validity characterization** using:
  - **Score (gradient) mean test** at true parameters.
  - **Fisher Information Matrix (FIM) consistency** tests (two estimators).
  - **Residual diagnostics** (normalized residuals, autocorrelation; optionally PIT).
  - **Score correlation diagnostics** (Cov(sum score) vs sum Cov(score)).
  - **Variance inflation / sandwich factor** to quantify parameter uncertainty distortion.
- A **validity map** (main figure) over an interpretable 2D regime.

### Out of scope (explicitly excluded to reduce attack surface)

- **MicroIR / microscopic-recursive “gold standard”** comparisons in the main paper.
  - Rationale from the audio notes: it expands complexity (requires additional conditional variance / boundary-state machinery) and can derail the paper into implementation details.
  - We can keep MicroIR as future work or a separate companion note.

## 2) “Story arc” of the paper (4-figure narrative)

We keep the paper anchored to a simple narrative: **problem → mechanism → validation → regime map**.

### Figure 1 — Problem + solution sketch

- Show the **core problem**: inference from *interval-averaged* observations is not the same as point sampling; naïve approximations can bias inference.
- Show a minimal schematic of **MR vs IR update timing**:
  - MR: Bayesian update at the interval boundary, Markov propagation across interval.
  - IR: Bayesian “interval-recursive” update that integrates start/end information more symmetrically.
- Candidate implementation artifact already exists:
  - `projects/eLife_2025/ops/local/figure_1.macroir`
  - `projects/eLife_2025/figures/figure_1.Rmd`

### Figure 2 — Mechanics (what MacroIR is doing)

- This is **not** a full microscopic derivation; it is the “just enough” algebra and algorithm flow:
  - Boundary-state concept (start/end states), marginalizations, and the key quantities used to avoid explicit covariance explosion.
  - Explain where each approximation enters (recursive vs non-recursive, averaging, variance correction, Taylor correction).
- The end state of this figure: a reader can follow the algorithm and see what it approximates.

### Figure 3 — Validation (efficacy before efficiency)

Primary tests (must-pass):

1. **Score mean**: at true parameters, the expected gradient is 0.
2. **FIM equality (two estimators)**: for a well-specified likelihood, `Cov(score)` should match `-E[Hessian]` (up to sampling error).

Secondary but very diagnostic (strongly recommended):

3. **Residual tests**: normalized residuals have mean 0, variance 1, and weak temporal correlation.
4. **Score correlation diagnostic**: difference between `Cov(sum score)` and `sum Cov(score_t)` reveals hidden temporal correlation / missing memory.
5. **Variance inflation factor** (sandwich): quantify how a misspecified likelihood inflates/deflates parameter uncertainty (and implications for evidence).

### Figure 4 — Validity map (where approximations are valid)

Main plot axes (already selected earlier):

- **x-axis:** dimensionless interval length, **Δ / τ_min**
- **y-axis:** ensemble size, **N_ch**

Each cell summarizes “quality” per algorithm using a small set of metrics (from Figure 3), producing a regime map that answers:

- “MacroIR is necessary here.”
- “MacroMR is sufficient here.”
- “Even simpler approximations are ok here.”

## 3) Data + reproducibility anchor (existing artifacts)

We already have a working prototype pipeline in-repo:

- Local MacroIR scripts:
  - `projects/eLife_2025/ops/local/figure_1.macroir`
  - `projects/eLife_2025/ops/local/figure_2.macroir`
- Plotting:
  - `projects/eLife_2025/figures/figure_1.Rmd`
  - `projects/eLife_2025/figures/figure_2.Rmd`
- Example outputs:
  - `projects/eLife_2025/figures/data/*.csv`
- Runs/provenance (where present):
  - `projects/eLife_2025/runs/run-*/meta.json`

This paper plan assumes we evolve from these prototypes into:

1. A **clean minimal reproducible pipeline** for the 4 main figures.
2. A **repeatable grid runner** for the validity map (even if initially manual/local).

## 4) Manuscript editing targets (where the draft lives)

Main LaTeX draft directory (as currently in the repo):

- `docs/eLife 2025/version inicial/`

In particular, we should track and resolve the TODOs in:

- `docs/eLife 2025/version inicial/elife-macroir-merged.tex`

Example TODOs currently in that file include section-label checks, a term in the variance expression, a citation key, and consistency of scaling with `Σ` and `N_ch`.

## 5) What we need to decide (explicit open items)

1. **Case study placement**: keep a real “biological” case study in **Supplementary** (default) or promote to **Main text**.
2. **Which model schemes** anchor the main validity map:
   - Default: `scheme_CO` as primary, plus `scheme_CCO` as a robustness check/supplement.
3. **Thresholds** for “valid/invalid” in the validity map:
   - Define quantitative cutoffs (or rank-based) for each metric.
4. **How to track audio sources** in git:
   - Current preference noted: keep **both** `.ogg` and transcript `.md`.
   - Optional improvement: Git LFS for audio if repo size becomes painful.

## 6) Deliverables checklist

- Planning docs (this folder) kept up to date in git.
- Revised manuscript draft (LaTeX) with the new 4-figure narrative.
- Repro scripts + figure outputs sufficient to regenerate the main figures.
- A single “validity map” dataset + plot for the chosen grid (Δ/τ_min vs N_ch).

## 7) Next action (human workflow)

1. You read and edit:
   - `docs/papers/macroir-elife-2025/00_master_plan.md`
   - `docs/papers/macroir-elife-2025/01_workboard.md`
2. We incorporate your corrections into these docs.
3. Only then do we decide which (if any) code changes are needed and schedule them as separate tasks.

