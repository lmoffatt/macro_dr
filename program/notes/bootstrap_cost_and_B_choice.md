# Diagnostics bootstrap: cost model and how to choose B

Hard-won during the 2026-05 figure_2 runs. **Read this before raising
`n_boostrap_samples` (B), `max_lag`, or `nsim` "because the bootstrap is cheap."
It is not cheap.**

## The cost model (measured + derived from the code)

The diagnostics bootstrap (`calculate_Likelihood_diagnostics_preset_f`, called B
times by `bootstrap_it_two[_paired]`) does **not** re-aggregate cheap per-sim
totals. Each replicate **recomputes the full per-sample evolution/correlation
from the raw `dy`**:
- `calculate_Likelihood_diagnostics_evolution_impl`: `for t in 0..T` × sum over
  the `nsim` resampled sims → O(nsim·T).
- `build_series_reports_integral_impl` → `Series_Moment_statistics(max_lag)` per
  observable → O(nsim·T·max_lag).
- cross-sim covariance (CDM / Fisher_Covariance) per sample → O(nsim·T·p²).

So **per replicate ≈ O(nsim · T · max_lag · p²)**, and total

> **cost ≈ B · nsim · T · max_lag · p²**  — linear in every factor at once.

## Why there is no free lunch (precompute can't save it)

The bootstrap needs the **cross-simulation variance/covariance at each sample**,
over the **resampled set, which changes every replicate**. Two facts kill any
precompute-once scheme:
1. it's a **variance/covariance** (not a plain mean → not a simple additive
   reduction), and
2. the **resampled set changes per replicate**, so nothing carries over.

The per-sim Gram `G_i = Σ_t s_{i,t}s_{i,t}ᵀ` *can* be precomputed and reweighted
(O(nsim·p²)/replicate), and the within-sim autocorrelation (the `max_lag` series)
*can* be precomputed per sim — but the **per-sample mean/variance over the
changing set is irreducibly O(nsim·T) per replicate** (collapsing its T would
need nsim² cross-sim pair Grams, which is worse). Floor: **B·nsim·T**. This is
the same wall as the CDM discussion — confirmed, not a missed optimization.

## The numbers that bit us (2026-05-28/29, dirac)

- nsim=4096: dlikelihood+Fisher = **68 min**, bootstrap = **21 min** (ratio 0.31).
- Then B was raised **~20×** and nsim → 16384.
- Both phases are linear in nsim, so **nsim cancels in the ratio** — only the ×20
  in B moves bootstrap-vs-dlikelihood: new ratio = 0.31 × 20 = 6.18.
- nsim=16384 dlikelihood = ~6 h ⇒ **predicted bootstrap ≈ 6 h × 6.18 ≈ 37 h.**
- Observed: one nsim=16384 cell at 32 cores (3156% CPU), 22 GB / 125 GB RAM, swap
  idle — **not memory, not a parallelism bug, pure compute** — still unfinished
  >10 h into the bootstrap. Matches the 37 h prediction.

## Choosing B (the actual guidance)

MC error of a bootstrap quantile scales as **1/√B**. For a well-conditioned
(≈Gaussian) direction, at **B=100**:

| quantile | MC SE | relative |
|---|---|---|
| 0.5 median | 0.13σ | — |
| 0.16/0.84 (±1σ) | 0.15σ | ~15% |
| 0.025/0.975 (95%) | 0.27σ | ~14% |

- **B=100 is fine** for figure error bars / qualitative comparison (~14–15%
  endpoint noise). `0.975` is NOT specially bad in *relative* terms vs ±1σ.
- **B=200 is the worst value** — 1/√B means 100→200 buys only ×1.4 for ×2 cost.
- To halve the noise you need ×4 B; ~3% needs ~×20 B = the ~37 h run.
- **Default B=100.** Bump B for a *single* cell only if a quantitative claim
  needs a razor tail. For the heavy-tailed / ill-conditioned directions no
  feasible B helps anyway (see `diagnostics_golden_test_plan.md` and the κ(F)/√n
  argument) — report those wide/flagged, don't chase them with B.
- At large nsim the well-conditioned CIs are essentially Gaussian → an analytic
  ±z·σ from the per-sim variance is ~free; the bootstrap earns its keep only on
  the non-Gaussian directions.
- Structural option: **split IDM (per-sim-reducible → high B cheap) from CDM
  (cross-sim → irreducible → low B)** if you want different resolution per class.

## How we diagnosed it (reuse these)
- `sstat -j <jobid>.batch --format=JobID,MaxRSS,MaxVMSize` — running-job peak RAM.
- `ssh <node> free -h` — swap line: GBs used + ~0 available ⇒ thrashing; else not memory.
- `ssh <node> 'top -bn1 | grep -i macrodr'` — %CPU: ~3200% = all 32 cores (compute-bound); ~100% = serial bug.
