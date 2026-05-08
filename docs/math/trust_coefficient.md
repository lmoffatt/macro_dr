# Trust coefficient α in the recursive Macro update

The recursive Macro algorithm performs a Kalman-style Bayesian update of the
posterior channel-state distribution at every interval. Two coupled quantities
are updated:

```
μ_new = μ + α · chi · gS                            (mean)
Σ_new = Σ_pre − (α · N / y_var) · gSᵀ gS            (covariance)
```

with

```
chi   = (y − y_mean) / y_var                       (innovation, normalized)
gS    = gᵀ · Σ_pre  (or the averaging=2 variant)   (gain block)
Σ_pre = Aᵀ · SmD · A + diag(μ · A)                 (Markov-propagated prior)
```

A naive (α = 1) Kalman update can drive the posterior off the simplex (negative
or > 1 entries in μ) or off the PSD cone (negative diagonal in Σ). The
algorithm therefore introduces a trust coefficient α ∈ (0, 1] that scales both
updates by the same amount, picked as the largest α satisfying both the
simplex and PSD constraints.

## α_μ — simplex bound on the mean

For each component i, the mean update `μ_i + α · d_i` (with d = chi · gS) must
stay in [0, 1]:

```
d_i > 0 :   α  ≤  (1 − μ_i) / d_i        (stay ≤ 1)
d_i < 0 :   α  ≤  −μ_i      / d_i        (stay ≥ 0)
d_i = 0 :   no constraint
```

α_μ is the smallest such bound across i (clamped at 1). Implemented in
`Macro_DMR::calculate_trust_coefficient` (legacy/qmodel.h).

## α_Σ — PSD bound on the covariance

The Σ down-date subtracts a rank-1 matrix `β · gSᵀ gS` (with β = α · N / y_var).
Its diagonal at index i is `β · gS_i²`. For each diagonal of Σ_new to be
non-negative — a *necessary* PSD condition — we need

```
Σ_pre_(i,i) − β · gS_i²  ≥ 0
α  ≤  Σ_pre_(i,i) / ((N / y_var) · gS_i²)        whenever gS_i ≠ 0
```

α_Σ is the smallest such bound across i (clamped at 1). Implemented in
`Macro_DMR::calculate_psd_trust_coefficient` (legacy/qmodel.h, added 2026-05-08).

If any pre-update diagonal `Σ_pre_(i,i)` is non-positive at an index where
`gS_i ≠ 0`, the function returns 0 (the recursive step would push an
already-degenerate diagonal further negative). The caller can choose to refuse
the step (current behavior — `to_Covariance_Probability` rejects α = 0 outputs)
or fall back to a non-recursive propagation for that interval.

## Combined α

Both updates are scaled by the same coefficient:

```
α = min(α_μ, α_Σ)
```

Choosing the more restrictive bound guarantees both the simplex constraint on
μ_new and the diagonal-PSD constraint on Σ_new. A safety factor (currently
`trust_multiplying_factor = 0.9`) is applied when α < 1, leaving 10 % margin
against the boundary.

## Necessary vs sufficient PSD

`α_Σ` is a necessary condition for PSD, not a sufficient one. A symmetric
matrix can have non-negative diagonal yet still fail PSD due to off-diagonal
correlations (e.g. `[[1, 2], [2, 1]]` has positive diagonal but eigenvalues
3, −1). The full PSD bound for the rank-1 down-date `Σ_pre − β · vᵀv` is

```
β  ≤  1 / (vᵀ · Σ_pre⁺ · v)
```

where `Σ_pre⁺` is the pseudoinverse on the simplex tangent space (Σ_pre has
the all-ones vector in its null space, so it is rank-deficient by construction).
Computing this bound exactly is O(n³) per step (Cholesky / SVD), versus O(n)
for the diagonal bound. In practice the diagonal bound is sufficient for the
regimes encountered in patch-clamp inference; if it ever fails to prevent
non-PSD outputs, `to_Covariance_Probability` will catch the violation
downstream.

## Joseph-form alternative

A stronger remedy is to replace the simple subtraction with the Joseph form,
which is unconditionally PSD-preserving for any α:

```
Σ_new = (I − K H) · Σ_pre · (I − K H)ᵀ + K R Kᵀ
```

with K, H, R the appropriate Kalman gain, observation operator, and observation
noise. This costs O(n³) per step (vs O(n²) for the rank-1 down-date) and
requires deriving (K, H, R) for each averaging branch. We have not switched
because the diagonal-PSD trust coefficient has been sufficient, but the option
remains open.

## When α_Σ is the binding constraint

Empirically, α_Σ < α_μ in regimes with:

- low observation noise (large `N / y_var`), so each interval injects much
  information per unit prior covariance;
- accumulated information from prior intervals that have already squeezed
  Σ_pre below its initial bare-cov values;
- few channels (small N), so the prior covariance scale is small relative to
  the per-interval information gain;
- conductance vectors with most mass on a single state, so `gSᵀgS` is
  rank-1-like and concentrates information narrowly.

A persistent α_Σ < α_μ across many intervals is itself a diagnostic signal:
the algorithm is operating at the edge of its PSD-preserving regime, and a
Joseph-form update or a finer-grained averaging may be appropriate.

## History

Before the SymmetricMatrix half-population fix on 2026-05-08, `gS` was
silently truncated in one direction (one component effectively zero), making
`gSᵀgS` rank-1 in only one tangent direction. The Σ down-date never decremented
all diagonals, and the PSD overshoot now caught by α_Σ was hidden. After the
fix, the correct full `gS` properly cancels around `gS · u = 0`, both diagonals
of Σ are decremented, and the previously latent PSD-overshoot regime became
visible. `calculate_psd_trust_coefficient` was added in the same commit as the
fix to address it without falling back to Joseph form.

## See also

- `legacy/qmodel.h` :: `Macro_DMR::calculate_trust_coefficient`
- `legacy/qmodel.h` :: `Macro_DMR::calculate_psd_trust_coefficient`
- `legacy/qmodel.h` :: `Macro_DMR::safely_calculate_Algo_State_recursive`
  (call site that combines α_μ and α_Σ)
- `legacy/qmodel_types.h` :: `to_Covariance_Probability` (downstream PSD canary
  on the resulting Σ_new)
