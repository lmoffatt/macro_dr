# Trust coefficient α in the recursive Macro update

The recursive Macro algorithm performs a Kalman-style Bayesian update of the
posterior channel-state distribution at every interval. Two coupled quantities
are updated:

```
μ_new = μ + α · chi · gS                            (mean)
Σ_new = Σ_pre − (α · N / y_var) · gSᵀ gS            (covariance)
```

(The two `α`'s are **not** the same coefficient — the mean uses α_μ and the
covariance uses α_Σ, decoupled; and in the standard rank-1 branch the covariance
trust is just α = 1. See *Decoupled α* below.)

with

```
chi   = (y − y_mean) / y_var                       (innovation, normalized)
gS    = Cov(X_end, y)  as row vector               (endpoint-frame gain)
Σ_pre = Aᵀ · SmD · A + diag(μ · A)                 (Markov-propagated prior)
```

`gS` lives at the *endpoint* frame for all averaging values, so that the rank-1
down-date `XTX(gS)` is in the same frame as `Σ_pre`. The concrete form depends
on when the observation enters:

| averaging | observation depends on | `gS` formula                                       |
| --------- | ---------------------- | -------------------------------------------------- |
| 0         | X_mid via instantaneous g | `gᵀ · Σ_mid · P_half`                           |
| 1         | X_start via gmean_i      | `gmean_iᵀ · Σ_start · P`                        |
| 2         | full interval (integrated) | `gmean_iᵀ · SmD · P + p · gtotal_ij`           |

For avg=0/1 the trailing `· P` (or `· P_half`) propagates the start/mid-frame
gain through the remaining Markov dynamics. Without it the down-date is in the
wrong frame — manifested historically as a Distortion-Induced-Bias spike in
macro_R at long intervals / large Num_ch (fixed 2026-05-10).

A naive (α = 1) Kalman update can drive the posterior off the simplex (negative
or > 1 entries in μ). The algorithm therefore introduces a mean trust
coefficient α_μ ∈ (0, 1], picked as the largest α satisfying the simplex
constraint on μ.

The covariance is treated **separately** (the two trusts are *decoupled* —
there is no shared α scaling both updates). In the standard rank-1 branch the
undamped (α = 1) down-date is PSD by construction (see *α_Σ is redundant in the
standard rank-1 branch* below), so Σ takes the full step. A separate PSD trust
α_Σ is computed but is genuinely needed only in the variance_correction (Taylor)
branches; see those sections below.

α is computed in the **smooth (C∞) form**

```
α(alfa_p) = softmin(1, factor · alfa_p; ε)
softmin(a, b; ε) = (a + b − √((a − b)² + ε²)) / 2
```

so that α is C∞ in θ everywhere — no kink at the boundary `factor·alfa_p = 1`,
hence no Heaviside step in `∂α/∂θ`. Previous formulations had:
- A **discontinuous form** (`α = 1` if `alfa_p ≥ 1`, else `α = factor · alfa_p`)
  with a `1 − factor` jump at `alfa_p = 1` — δ-function in `∂α/∂θ`,
  catastrophic score variance at high N.
- A **C0-continuous form** `α = min(1, factor · alfa_p)` — no jump in α, but a
  Heaviside step in `∂α/∂θ` at `factor·alfa_p = 1` — still pumped variance into
  the score whenever the system hovered near the boundary.

The softmin form removes the step entirely. `ε` controls the smoothing band:
- `ε → 0` recovers the C0 form (and its kink).
- `ε ≈ 1e-4` (current default) gives a smooth band of width ε around the
  boundary; residual bias `α ≤ min(1, factor·alfa_p)` of at most `ε / 2` ≈ 5e-5,
  small vs the `1 − factor = 0.1` safety margin.

Two softmins are applied:
1. `α_μ = softmin(1, factor · alfa_p_μ)` — inside `calculate_trust_coefficient`
2. `α_Σ = softmin(1, factor · alfa_p_Σ)` — inside `calculate_psd_trust_coefficient`

The mean and covariance trusts are **decoupled**: the mean update uses α_μ
alone, the covariance update uses α_Σ alone (only in the variance_correction
branches; the standard branch leaves Σ undamped). There is **no** shared
`α = softmin(α_μ, α_Σ)` — that earlier combined form is gone, removing the
α_μ = α_Σ kink and, crucially, keeping the simplex-trust singularity (α_μ's
`(1 − μ_i)/d_i` bound at a residual zero-crossing `d → 0`) out of the
covariance down-date, where it had destabilized the per-sample numerical Fisher.

`factor ∈ (0, 1)` (typically 0.9) is the safety margin: when the constraint
binds the result is ≈ `factor · alfa_p`, leaving a `1 − factor` margin off the
simplex/PSD boundary.

## α_μ — simplex bound on the mean

For each component i, the mean update `μ_i + α · d_i` (with d = chi · gS) must
stay in [0, 1]:

```
d_i > 0 :   α  ≤  (1 − μ_i) / d_i        (stay ≤ 1)
d_i < 0 :   α  ≤  −μ_i      / d_i        (stay ≥ 0)
d_i = 0 :   no constraint
```

α_μ is the smallest such bound `alfa_p` across i, post-processed via
`α = min(1, factor · alfa_p)`. Implemented in
`Macro_DMR::calculate_trust_coefficient` (legacy/qmodel.h).

## α_Σ — PSD bound on the covariance

The Σ down-date subtracts a rank-1 matrix `β · gSᵀ gS` (with β = α · N / y_var).
Its diagonal at index i is `β · gS_i²`. For each diagonal of Σ_new to be
non-negative — a *necessary* PSD condition — we need

```
Σ_pre_(i,i) − β · gS_i²  ≥ 0
α  ≤  Σ_pre_(i,i) / ((N / y_var) · gS_i²)        whenever gS_i ≠ 0
```

α_Σ is the smallest such bound `alfa_p` across i, post-processed via
`α = min(1, factor · alfa_p)`. Implemented in
`Macro_DMR::calculate_psd_trust_coefficient` (legacy/qmodel.h, added 2026-05-08;
made continuous and AD-aware 2026-05-09).

The trust coefficient's role is to *reduce* gS so the down-date preserves
diagonal positivity — never to abort the step. Indices where `Σ_pre_(i,i) ≤ 0`
are already degenerate and no α > 0 brings them positive, so they are
**skipped** in the argmin: their constraint contributes no information to α,
and other indices set the bound. If no index binds (or only degenerate ones
exist), α = 1 and the full Bayesian step is taken; downstream PSD checks
(`to_Covariance_Probability`) catch any residual problem.

## α_Σ is redundant in the standard rank-1 branch

In the standard rank-1 branch the undamped (α = 1) down-date
`Σ_pre − (N / y_var) · gSᵀ gS` is **positive-semidefinite by construction**, so
α_Σ never binds and has been **dropped** — the covariance takes the full α = 1
step, guarded only by the `to_Covariance_Probability` canary.

The reason is the structure of `y_var = e + N · gᵀΣg`: the dominant PSD margin
is the Markov-step / multinomial covariance (the Löwner term
`Σ_pre ⪰ Aᵀ · Cov · A`), which is present even at residual `e = 0`. In closed
form for the 2-state av=1 case, the single nonzero eigenvalue of the α = 1
down-date is

```
λ = (e·A + N·B) / (2 (e + N·gSg))        with A, B ≥ 0   ⇒   λ ≥ 0  identically.
```

PSD first fails only at `α* = 1 + e/(N·gSg) > 1`, which the `cap ≤ 1` code never
reaches.

Note the margin `e / (e + N·gSg)` **grows as N shrinks** (`e` does not scale with
N), so α_Σ is even *less* needed at few channels. At low N the binding
constraint is the simplex/mean trust α_μ, not the PSD trust. (The earlier claim
that small N makes α_Σ bind is backwards for this branch; it holds only for the
variance_correction branches and/or n > 2 states — see below.)

## Decoupled α

The mean and covariance trusts are applied **independently** — they are *not*
combined into a single `min(α_μ, α_Σ)` scaling both updates:

```
μ_new : scaled by α_μ          (simplex trust, calculate_trust_coefficient)
Σ_new : scaled by α_Σ          (PSD trust, only in the variance_correction branches;
                                the standard rank-1 branch leaves Σ undamped, α = 1)
```

Why decouple. At a likelihood-residual zero-crossing `d = chi·gS → 0`, α_μ's
bound `(1 − μ_i)/d_i` is singular, but the *mean* step `α_μ·chi·gS ≡ 0` there,
so the singularity is harmless on the mean. Sharing a single α carried that
singularity into the Σ down-date — which is chi-free and scaled by `N` (≈ 10⁴)
— and destabilized the per-sample numerical Fisher. Decoupling removes that
path. (Keeping the two trusts on the same iterate is justified by the
tempered-likelihood view, not by a Kalman gain identity: the down-date here is
linear in α and carries `N`, not the α² Joseph form.)

A safety factor (currently `trust_multiplying_factor = 0.9`) is applied to each
trust when it binds (< 1), leaving 10 % margin against the respective boundary.

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
for the diagonal bound.

The per-diagonal form is in fact the wrong form regardless: in the safe 2-state
regime it is over-conservative (it can only spuriously bind, since the α = 1
down-date is already PSD by construction — see *α_Σ is redundant in the standard
rank-1 branch*), and in the dangerous n ≥ 3 regime it gives zero protection
against off-diagonal PSD breaches. Where genuine protection is required (the
variance_correction branches, see below), a det/eigen-sign check on the relevant
block — or the Joseph form — is the correct guard. The final PSD safety net is
the `to_Covariance_Probability` canary downstream.

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

**Scope.** This applies only to the variance_correction (Taylor) branches and/or
n > 2 states. In the standard rank-1 branch α_Σ never binds (the α = 1 down-date
is PSD by construction), so it is not applied at all there.

Within those branches, α_Σ < α_μ in regimes with:

- low observation noise (large `N / y_var`), so each interval injects much
  information per unit prior covariance;
- accumulated information from prior intervals that have already squeezed
  Σ_pre below its initial bare-cov values;
- conductance vectors with most mass on a single state, so `gSᵀgS` is
  rank-1-like and concentrates information narrowly.

Note that *few channels* (small N) does **not** make α_Σ bind — the PSD margin
`e/(e + N·gSg)` grows as N shrinks, so the PSD trust is *less* needed at low N.
At low N it is the simplex/mean trust α_μ that binds, not α_Σ.

A persistent α_Σ < α_μ across many intervals is itself a diagnostic signal:
the algorithm is operating at the edge of its PSD-preserving regime, and a
Joseph-form update or a finer-grained averaging may be appropriate.

## Where α_Σ is genuinely needed: the variance_correction branches

α_Σ is retained (and necessary) only in the variance_correction (Taylor)
branches — MacroMRT (av=1) and MacroIRT (av=2). Those use a **rank-2 Woodbury**
down-date whose middle block `m22 = c − 2V²/N` can go negative (indefinite K),
so the self-limiting PSD guarantee of the standard rank-1 branch structurally
breaks.

Even there, the per-diagonal α_Σ form is *necessary but not sufficient*: for
n ≥ 3 it false-passes off-diagonal PSD breaches. The correct guard is a
det/eigen-sign check on the 2×2 middle block `K`, or a Joseph form. (These
Taylor branches are cut from the May-2026 paper.)

## History

Before the SymmetricMatrix half-population fix on 2026-05-08, `gS` was
silently truncated in one direction (one component effectively zero), making
`gSᵀgS` rank-1 in only one tangent direction. The Σ down-date never decremented
all diagonals, and the PSD overshoot now caught by α_Σ was hidden. After the
fix, the correct full `gS` properly cancels around `gS · u = 0`, both diagonals
of Σ are decremented, and the previously latent PSD-overshoot regime became
visible. `calculate_psd_trust_coefficient` was added in the same commit as the
fix to address it without falling back to Joseph form.

In retrospect (2026-06), the "PSD overshoot" this was meant to catch was not a
genuine exit from the PSD cone in the standard rank-1 branch: that down-date is
PSD by construction (see *α_Σ is redundant in the standard rank-1 branch*), so
the per-diagonal α_Σ there could only ever bind spuriously. α_Σ has since been
dropped from the standard branch and is retained only in the variance_correction
(Taylor) branches, where the down-date is a rank-2 Woodbury form that can
genuinely lose PSD.

## See also

- `legacy/qmodel.h` :: `Macro_DMR::calculate_trust_coefficient`
- `legacy/qmodel.h` :: `Macro_DMR::calculate_psd_trust_coefficient`
- `legacy/qmodel.h` :: `Macro_DMR::safely_calculate_Algo_State_recursive`
  (call site; applies α_μ to the mean and leaves Σ undamped in the standard
  branch — α_μ and α_Σ are decoupled, no shared α)
- `legacy/qmodel_types.h` :: `to_Covariance_Probability` (downstream PSD canary
  on the resulting Σ_new)
