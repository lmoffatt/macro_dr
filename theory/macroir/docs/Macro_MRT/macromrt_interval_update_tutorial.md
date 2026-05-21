# MacroMRT Interval Update: Tutorial

This document explains the MacroMRT interval update in three layers,
parallel to `macroir_interval_update_tutorial.md`:

1. **Developer Tutorial (C++/MacroMRT)**
2. **Supplementary Material (Full Derivations)** — see
   `macromrt_macromrt_supplement.tex`
3. **Main Paper Section** — see `macromrt_macromrt_paper_section.md`

---

## Dispatch (CLI flags)

MacroMRT is reached via the standard likelihood-model dispatch with:

| Flag                       | Value           |
|----------------------------|-----------------|
| `adaptive_approximation`   | (model choice)  |
| `recursive_approximation`  | `true`          |
| `averaging_approximation`  | `1`             |
| `variance_approximation`   | `true`          |
| `taylor_qdt_approximation` | (Qdt path; default fine) |
| `micro`                    | `false`         |

The companion macro_MR is the same dispatch with
`variance_approximation = false`. The companion macro_IRT is
`averaging_approximation = 2`. So the two pairwise contrasts that
matter for the IR-vs-MR question are:

- **MR vs MRT**: hold `averaging=1`, toggle `variance` — isolates the
  Taylor σ² correction.
- **MRT vs IRT**: hold `variance=true`, toggle `averaging` — isolates
  the boundary-state lift of MacroIR.

In `figure_2.macroir`, add `macro_MRT` to the labels list and set the
matching `taylor_approximation` (a.k.a. `variance_approximation`) entry
to `true`. The full list becomes
`[macro_NR, macro_R, macro_MNR, macro_MR, macro_MRT, macro_IR, macro_IRT]`
with `taylor_approximation = [false, false, false, false, true, false, true]`.

---

## A.1 Inputs from Qdtm

MacroMRT consumes per-`i_0` boundary moments from the Qdtm (moments)
path:

| Symbol                      | C++ name           | Source                         |
|-----------------------------|--------------------|--------------------------------|
| `gbar(i_0)`                 | `gmean_i[i0]`      | Qdtm marginalised over `i_t`   |
| `sbar(i_0)` (open-channel)  | `gvar_i[i0]`       | Qdtm marginalised over `i_t`   |
| `P(t)`                      | `P_t(i0,it)`       | matrix exponential `exp(Qt)`   |
| `mu0`, `Sigma0`             | `mu_0`, `Sigma_0`  | prior at start of interval     |

The same `gmean_i`/`gvar_i` that MacroMR already uses are sufficient
for MacroMRT; no new boundary-state machinery beyond Qdtm is needed.

> **Verification gate before wiring MRT.** A memory note flags that
> `calc_Qdt_taylor`'s `gmean*`/`gvar*` outputs were placeholder
> formulas that broke micro_MR/IR figure_2 logL. Confirm that the
> macro-side `gmean_i`/`gvar_i` consumed by MacroMR is the validated
> path (not the broken micro path) before enabling MRT in dispatch.
> Cf.\ `project_qdt_taylor_seed_bug.md`.

### Agonist-change handling within an interval

When the agonist concentration changes within an interval, two
distinct `Q` matrices apply on the two sub-intervals. The boundary
kernel `P(t)` and the boundary moments
`Gamma_bar(i_0,i_t)`, `V_bar(i_0,i_t)` must be computed by
composing the two sub-interval contributions before marginalising
over `i_t` in equations \eqref{eq:gbar}--\eqref{eq:sbar} of the
supplement. This is the reason MacroMR cannot bypass the
boundary-conditioned moments even though the measurement is
"instantaneous" — the within-interval transition kernel is non-trivial
on agonist-change steps.

---

## A.2 Predictive Moments

```cpp
// Predictive mean
double y_pred = N_ch * dot(mu_0, gbar);

// Predictive variance (no boundary cross-correlation term — that's IR's)
double V = eps2_0t + N_ch * dot(mu_0, sbar);
```

Compared with MacroIR, the `N_ch * tilde_gamma_Sigma_gamma_scalar`
term is **dropped**. This is the structural M-vs-I difference; MacroMRT
inherits it from MacroMR.

---

## A.3 Effective Direction and Scalar

```cpp
double delta = y_obs - y_pred;
Vector v = gbar + (delta / V) * sbar;       // effective direction
double s = dot(v, Sigma_prop * v);          // effective scalar
```

When `sbar = 0` (no open-channel noise), `v == gbar` and the rest of
the update collapses to the Moffatt 2007 MacroR update.

---

## A.4 Propagation

```cpp
RowVector mu_prior = mu_0 * P_t;
Matrix Sigma_prop  = P_t.transpose() * (Sigma_0 - diag(mu_0)) * P_t
                   + diag(mu_prior);
```

Identical to MacroIR.

---

## A.5 Measurement Update (rank-1 quasi-Laplace, working scheme)

```cpp
// Cross-covariance vector (no tilde — direct multiplication)
Vector Sp_v = Sigma_prop * v;

// Sherman-Morrison rank-1 covariance update
Matrix Sigma_post = Sigma_prop
    - (N_ch / (V + N_ch * s)) * outer_product(Sp_v, Sp_v);

// Mean update (Newton step from prior mean)
RowVector mu_post = mu_prior
    + (delta / (2.0 * V)) * (gbar + v) * Sigma_post;
```

The pattern matches the Macro_IRT update with two substitutions:
`tilde_gamma_Sigma_vec → Sp_v` (i.e. `Sigma_prop * gbar`) and
`tilde_v_Sigma_vec   → Sigma_prop * v`. Reusing the IRT update
template with these substitutions is the cleanest way to share code
between MRT and IRT.

---

## A.6 Trust-Region Simplex Shrink

Identical to MacroIR. Compute α* as the largest scalar in [0,1] keeping
every component of `mu_post` inside `[p_min, p_max]`; scale both the
mean and covariance updates by α*. Record α* as `trust_coefficient`.
The likelihood `log L` uses the unmodified `V` (the variance-inflation
interpretation was discarded; see
`macro_ir_variance_inflation_correction.tex`).

---

## A.7 Log-Likelihood

```cpp
double logL = -0.5 * (std::log(2 * M_PI * V) + delta * delta / V);
```

Identical in form to MacroIR/MacroMR; only `V` differs (per A.2).

---

## B Reduction Cross-Check (sanity test for the implementation)

Set `gvar_i[*] = 0` in the model. With `variance_approximation = true`
and `gvar_i = 0`:

- `v == gbar`,
- `s == gbar^T Sigma_prop gbar`,
- the update reduces to the Moffatt 2007 form.

This must produce identical per-trace `logL` to running with
`variance_approximation = false` on the same data — a one-line check
after wiring MRT.

---

## C Cross-References

- Math: `macromrt_macromrt_supplement.tex`
- Paper section: `macromrt_macromrt_paper_section.md`
- Cheatsheet: `macromrt_cheatsheet.md`
- Tilde operator (used by IRT, dropped here): `MacroIR_tilde.md`
- Trust region: `macro_ir_variance_inflation_correction.tex`
