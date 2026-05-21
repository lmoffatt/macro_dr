# MacroIRT Interval Update: Tutorial

This document explains the MacroIRT interval update in three layers,
parallel to `macroir_interval_update_tutorial.md`:

1. **Developer Tutorial (C++/MacroIRT)**
2. **Supplementary Material (Full Derivations)** — see
   `macroirt_macroirt_supplement.tex`
3. **Main Paper Section** — see `macroirt_macroirt_paper_section.md`

---

## Dispatch (CLI flags)

MacroIRT is reached via the standard likelihood-model dispatch with:

| Flag                       | Value           |
|----------------------------|-----------------|
| `adaptive_approximation`   | (model choice)  |
| `recursive_approximation`  | `true`          |
| `averaging_approximation`  | `2`             |
| `variance_approximation`   | `true`          |
| `taylor_qdt_approximation` | (Qdt path; default fine) |
| `micro`                    | `false`         |

The companion macro_IR is the same dispatch with
`variance_approximation = false`. The companion macro_MRT is
`averaging_approximation = 1`. So the two pairwise contrasts that
matter for the IR-vs-MR question are:

- **IR vs IRT**: hold `averaging=2`, toggle `variance` — isolates the
  Taylor σ² correction in the interval-averaged regime.
- **MRT vs IRT**: hold `variance=true`, toggle `averaging` — isolates
  the boundary-state lift of MacroIR.

---

## A.1 Inputs from Qdt

MacroIRT consumes the full boundary-conditioned moments from the Qdt
(non-moments) path:

| Symbol                         | C++ name              | Source                       |
|--------------------------------|-----------------------|------------------------------|
| `Gamma_bar(i_0,i_t)`           | `Gamma_bar(i0,it)`    | full boundary mean conductance |
| `V_bar(i_0,i_t)` (open-channel)| `V_bar(i0,it)`        | full boundary open-channel noise|
| `P(t)`                         | `P_t(i0,it)`          | matrix exponential `exp(Qt)` |
| `mu0`, `Sigma0`                | `mu_0`, `Sigma_0`     | prior at start of interval   |

Unlike MacroMRT, MacroIRT requires both indices `(i_0,i_t)` to be
preserved so the tilde operator can lift to boundary-state space.

Per-`i_0` collapses are still convenient and shared with MR/MRT:

```cpp
for (int i0 = 0; i0 < K; ++i0) {
    gbar[i0]  = 0; sbar[i0] = 0;
    for (int it = 0; it < K; ++it) {
        gbar[i0] += P_t(i0,it) * Gamma_bar(i0,it);
        sbar[i0] += P_t(i0,it) * V_bar(i0,it);
    }
}
```

---

## A.2 Predictive Moments

```cpp
// Predictive mean (same as MacroIR/MacroMR)
double y_pred = N_ch * dot(mu_0, gbar);

// Boundary-lifted cross-covariance scalar (tilde scalar)
// Same routine MacroIR uses
double tilde_gSg = tilde_gamma_Sigma_gamma_scalar(
    Gamma_bar, P_t, mu_0, Sigma_0);

// Per-i_0 intrinsic noise sum (third variance term — present in IR/IRT, dropped in MR/MRT)
double intrinsic = N_ch * dot(mu_0, sbar);

// Total predictive variance — three terms
double V = eps2_0t + N_ch * tilde_gSg + intrinsic;
```

The `tilde_gSg` and the third term (intrinsic) are the two pieces
MacroMRT drops; they are what justifies IR's existence over MR.

---

## A.3 Effective Direction and Scalar (the IRT extension)

```cpp
double delta = y_obs - y_pred;
Vector v = gbar + (delta / V) * sbar;       // effective direction (same as MRT)

// Effective scalar via the tilde scalar, applied to v (not gamma)
// Uses U_bar_v(i_0,i_t) = Gamma_bar(i_0,i_t) + (delta/V) * V_bar(i_0,i_t)
double s = tilde_v_Sigma_v_scalar(
    Gamma_bar, V_bar, P_t, mu_0, Sigma_0, delta, V);
```

When `V_bar = 0` (no open-channel noise), `v == gbar` and `s == tilde_gSg`;
the update collapses to MacroIR's standard scalar Kalman correction.

---

## A.4 Cross-Covariance Vector (Tilde Vector)

```cpp
// Tilde vector for gamma (same as MacroIR)
Vector tilde_g_Sigma = tilde_gamma_Sigma_vec(
    Gamma_bar, P_t, mu_0, Sigma_0);

// Tilde vector for the effective direction v
Vector tilde_v_Sigma = tilde_v_Sigma_vec(
    Gamma_bar, V_bar, P_t, mu_0, Sigma_0, delta, V);
```

For implementations that already provide `tilde_gamma_Sigma_vec`, the
`tilde_v_Sigma_vec` routine is a one-line variant: replace
`Gamma_bar(i_0,i_t)` with `Gamma_bar(i_0,i_t) + (delta/V) * V_bar(i_0,i_t)`.

---

## A.5 Propagation

```cpp
RowVector mu_prior = mu_0 * P_t;
Matrix Sigma_prop  = P_t.transpose() * (Sigma_0 - diag(mu_0)) * P_t
                   + diag(mu_prior);
```

Identical to MacroIR/MacroMRT.

---

## A.6 Measurement Update (rank-1 quasi-Laplace, working scheme)

```cpp
// Sherman-Morrison rank-1 covariance update with tilde lift
Matrix Sigma_post = Sigma_prop
    - (N_ch / (V + N_ch * s))
      * outer_product(tilde_v_Sigma, tilde_v_Sigma);

// Mean update (Newton step from prior mean)
RowVector mu_post = mu_prior
    + (delta / (2.0 * V)) * (tilde_g_Sigma + tilde_v_Sigma) * Sigma_post;
```

Compared with MacroIR's rank-1 update, the differences are:

1. `gamma` replaced everywhere by `v` (effective direction);
2. `tilde_v_Sigma` used instead of `tilde_g_Sigma` in the covariance
   downdate;
3. mean step gains the `(tilde_g_Sigma + tilde_v_Sigma)` symmetrised
   gradient direction.

The factor of `1/2` in the mean update is the Newton-step factor for
the rank-1 quasi-Laplace energy `δ²/V + P` (not `2δ²/V + 2P` —
i.e.\ the unscaled negative log-posterior); cf.\ MRT supplement.

---

## A.7 Trust-Region Simplex Shrink

Identical to MacroIR/MacroMRT. Compute α* keeping `mu_post` inside the
simplex floor; scale both updates by α*; record `trust_coefficient`;
`log L` uses unmodified `V`. See
`macro_ir_variance_inflation_correction.tex`.

---

## A.8 Log-Likelihood

```cpp
double logL = -0.5 * (std::log(2 * M_PI * V) + delta * delta / V);
```

Identical in form to MacroIR/MacroMRT; only `V` differs (per A.2).

---

## B Reduction Cross-Checks (sanity tests)

1. **IRT → IR**: set `V_bar = 0` in the model. With
   `variance_approximation = true` and `V_bar = 0`, the update must
   produce identical per-trace `logL` to MacroIR.

2. **IRT → MRT (algebraic only)**: not a flag toggle. Run the IRT
   code path on a model where `Sigma_0 - diag(mu_0)` is replaced by
   zero (so the tilde operator's first term vanishes), and check that
   `tilde_v_Sigma` collapses to `(V_bar∘P)^T mu_0` matching the MRT
   per-`i_0` form. This is an algebraic test, not a regression test.

---

## C The Macro Universe (seven-method experimental grid)

To attribute the IRT gain over baseline MR (Moffatt 2007), repeat the
figure_2 sweep over the full set:

```
labels = ["macro_NR", "macro_R", "macro_MNR",
          "macro_MR",  "macro_MRT",
          "macro_IR",  "macro_IRT"]
taylor_approximation = [false, false, false,
                        false, true,
                        false, true]
averaging_approximation = [0, 0, 1,
                           1,  1,
                           2,  2]
recursive_approximation = [false, true, false,
                           true,  true,
                           true,  true]
```

The IRT-vs-MRT contrast (last two) isolates the boundary-state lift;
the MRT-vs-MR contrast (4 vs 5) isolates the Taylor σ² correction.

---

## D Cross-References

- Math: `macroirt_macroirt_supplement.tex`
- Paper section: `macroirt_macroirt_paper_section.md`
- Cheatsheet: `macroirt_cheatsheet.md`
- MR-side counterpart (this tutorial's twin): `Macro_MRT/macromrt_interval_update_tutorial.md`
- Tilde operator (load-bearing for IRT): `MacroIR_tilde.md`
- Trust region: `macro_ir_variance_inflation_correction.tex`
