# Bayesian-prior regularization for $\bar\gamma_{i\to j}$ and $\mathrm{Var}(\bar\gamma_{i\to j})$

## Companion to `gmean_ij_gvar_by_blocks.tex`

The companion .tex derives, given the integrated conductance
$A(t) = \int_0^t \gamma(X_s)\,ds$ on a Markov chain with generator $Q$ and
state-dependent conductance vector $\boldsymbol\gamma$, the conditional
moments

$$
\bar\gamma_{i\to j}(t)
\;=\;
\frac{U_{ij}(t)}{t\,P_{ij}(t)},
\qquad
\mathrm{Var}\bigl(\bar\gamma_{i\to j}(t)\bigr)
\;=\;
\frac{W_{ij}(t)}{t^2 P_{ij}(t)}
\;-\;
\bigl(\bar\gamma_{i\to j}(t)\bigr)^2,
$$

where $U_{ij}$ and $W_{ij}$ are the first two $\alpha$-derivatives of the
tilted transition matrix at $\alpha=0$. In the codebase the rescaled
quantities $\texttt{gtotal\_ij}$, $\texttt{gtotal\_var\_ij}$ stand in for $U$
and $W$ (modulo factors of $t$, $n$), so the working definitions are

$$
\texttt{gmean\_ij}_{i,j} \;=\; \texttt{gtotal\_ij}_{i,j} \,/\, P_{ij},
\qquad
\texttt{gvar\_ij}_{i,j}  \;=\; \texttt{gtotal\_var\_ij}_{i,j} \,/\, P_{ij}.
$$

These definitions are well-posed in exact arithmetic for every $P_{ij} > 0$,
but **numerically pathological** when $P_{ij} \ll 1$: both numerator and
denominator approach zero, and the FP error of the quotient diverges. This
note describes the methodology that handles those entries without
introducing a discontinuity, without losing differentiability, and without
silently fabricating physically wrong values.

## The pathology

For a continuous-time Markov chain with rate-scale $\lambda$ and integration
window $t$, off-diagonal transition probabilities scale as
$P_{ij}(t) \approx (Q_{ij}\,t)$ for $\lambda t \ll 1$. For entries
representing **multi-hop transitions** ($k$ hops apart, with no direct
$i\to j$ rate), $P_{ij}(t) \approx O((\lambda t)^k)$ — vanishing
geometrically in the time-step ratio.

For such an entry:

1. The conditional mean $\bar\gamma_{i\to j}$ is a property of an event that
   essentially never happens in window $t$.
2. The Monte-Carlo / numerical estimator
   $\texttt{gtotal\_ij}/P_{ij}$ is the ratio of two tiny numbers, each at
   FP-rounding scale. The quotient is dominated by rounding noise; its
   "value" carries no physical content.
3. Worse, propagating that noise downstream (into the IRT/MRT Newton
   updates, which contract $\texttt{gtotal\_ij}$ against $\texttt{gmean\_ij}$
   to form variance kernels) injects spurious derivatives along directions
   that have no physical analogue.

A regularization is needed. The question is *which*.

## What the codebase did, and why it's unsatisfying

The current implementation (multiple sites in `legacy/qmodel.h`) uses

```cpp
gmean_ij = elemDivSoftAbs(gtotal_ij, P, min_P);
         // = gtotal_ij / √(P_ij² + min_P²)
```

This is a smooth replacement of $1/P$ that avoids the divide-by-zero, but
it has no physical interpretation:

- The denominator $\sqrt{P^2 + \varepsilon^2}$ is the L2 norm of the pair
  $(P,\varepsilon)$. It is not a probability, not a count, not anything one
  reasons about a-priori.
- In the limit $P_{ij} \to 0$, the regularized quotient goes to
  $\texttt{gtotal\_ij}/\varepsilon \to 0$ (because $\texttt{gtotal\_ij} \to
  0$ even faster). The implicit "fallback estimate" is $\bar\gamma_{i\to j}
  = 0$ — i.e. "if you can't see the transit, assume zero conductance during
  it." That is not what the physics dictates.
- The derivative
  $\partial\bigl(\texttt{gtotal\_ij}/\sqrt{P^2+\varepsilon^2}\bigr)/\partial\theta$
  is the soft-abs derivative, which for small $P$ deviates from the true
  $\partial(\texttt{gtotal\_ij}/P)/\partial\theta$. The deviation is
  $O((\varepsilon/P)^2)$ — small for well-behaved entries, but for the very
  entries where the regularizer matters, it dominates.
- It opens up a non-physical asymmetric bias whose magnitude tracks
  `min_P` in a way nobody is keeping track of.

In short: it was a stopgap that papers over the singularity with a smooth
kernel, but it does so by smuggling in a fictitious zero-conductance prior
without flagging it as such.

## The physically motivated fallback (v_g endpoint average)

Condition on the rare event $X_t = j \mid X_0 = i$ in the small-$P_{ij}$
limit. The transit happens on the FP-noise scale of the eigendecomposition,
and on that scale the chain spends an essentially uniform fraction of
the window at the two endpoint states $i$ and $j$ — intermediate
residence is short by the same factor. The time-averaged conductance
over the transit is then the endpoint average:

$$
\boxed{\;
\bar\gamma_{i\to j}^{\text{prior}}
\;=\;
\tfrac{1}{2}\bigl(\gamma_i + \gamma_j\bigr),
\;}
$$

and the corresponding second moment (averaging $\gamma^2$ over the same
uniform endpoint dwell) is

$$
\boxed{\;
\overline{\gamma^2}_{i\to j}^{\text{prior}}
\;=\;
\tfrac{1}{2}\bigl(\gamma_i^2 + \gamma_j^2\bigr).
\;}
$$

The variance prior follows algebraically:

$$
\mathrm{Var}^{\text{prior}}\!\bigl(\bar\gamma_{i\to j}\bigr)
\;=\;
\overline{\gamma^2}_{i\to j}^{\text{prior}}
\;-\;
\bigl(\bar\gamma_{i\to j}^{\text{prior}}\bigr)^2
\;=\;
\biggl(\frac{\gamma_i - \gamma_j}{2}\biggr)^{\!2}.
$$

These priors depend only on the per-state conductance vector
$\boldsymbol\gamma$ — they are model parameters, not derived quantities.

### Why not Empirical Bayes from marginals?

An Empirical-Bayes alternative — building the prior from row/column
sums of $\texttt{gtotal\_ij}$ itself — was considered and rejected on a
numerical-stability ground specific to the eigendecomposition path.

When $V$ is ill-conditioned (large $\kappa(V)$), the reconstruction
$P = V \cdot \mathrm{diag}(e^{\lambda\,dt}) \cdot W$ is FP-contaminated
at level $\varepsilon_{\text{mach}}\cdot\kappa(V)$. This is exactly the
regime in which $\texttt{gtotal\_ij}$ entries fall below the noise
floor and the prior is supposed to dominate. The row/column marginals
are then sums of contaminated entries, and the EB prior built from
those marginals feeds the contamination back into the shrinkage target
— the regularization defeats itself in the regime it was designed for.

The v_g prior is conservative against this failure mode: it depends on
no eigendecomposition, only on per-state model parameters. The cost is
losing EB's ability to incorporate multi-hop intermediate-state
structure automatically; but that loss matters most for moderately
rare entries, not for the FP-pathological ones the regularization is
*for*. An iterative scheme (v_g first pass for de-noising, EB second
pass on the cleaner matrix) would formally recover the multi-hop
benefit, but the first pass already discriminates the well-conditioned
regime from the regularized regime, so the second-pass gain is bounded
and not worth the complexity.

## Bayesian shrinkage

Treat the observed quantities as a posterior update under a conjugate
prior. The natural parameterization:

- "Data": $\texttt{gtotal\_ij}$ is the integrated conductance summed over a
  pseudo-count of $P_{ij}$ paths.
- "Prior": a pseudo-count of $\varepsilon_{\text{prior}}$ fictitious paths,
  each contributing
  $\bar\gamma_{i\to j}^{\text{prior}} = (\gamma_i + \gamma_j)/2$.

The posterior mean is the count-weighted average:

$$
\boxed{\;
\texttt{gmean\_ij}_{i,j}
\;=\;
\frac{\texttt{gtotal\_ij}_{i,j} \;+\; \texttt{gmean\_ij\_prior}_{i,j}\cdot
       \texttt{min\_P\_prior}}
     {P_{i,j} \;+\; \texttt{min\_P\_prior}},
\;}
$$

and analogously for the variance using $\texttt{gtotal\_var\_ij}$ and the
endpoint-variance prior:

$$
\boxed{\;
\texttt{gvar\_ij}_{i,j}
\;=\;
\frac{\texttt{gtotal\_var\_ij}_{i,j} \;+\; \texttt{gvar\_ij\_prior}_{i,j}\cdot
       \texttt{min\_P\_prior}}
     {P_{i,j} \;+\; \texttt{min\_P\_prior}}.
\;}
$$

with

$$
\texttt{gmean\_ij\_prior}_{i,j} \;=\; \tfrac{1}{2}(\gamma_i + \gamma_j),
\qquad
\texttt{gvar\_ij\_prior}_{i,j}  \;=\; \tfrac{1}{4}(\gamma_i - \gamma_j)^2.
$$

## On the symmetry of the prior

The prior $(\gamma_i + \gamma_j)/2$ is symmetric in $i \leftrightarrow j$,
but the shrunken estimator is not — and shouldn't be. Rewriting the
posterior mean as a convex combination,

$$
\texttt{gmean\_ij}_{i,j}
\;=\;
\underbrace{\frac{P_{i,j}}{P_{i,j}+\varepsilon}}_{w_{\text{data}}}\cdot
\underbrace{\frac{\texttt{gtotal\_ij}_{i,j}}{P_{i,j}}}_{\text{data estimate}}
\;+\;
\underbrace{\frac{\varepsilon}{P_{i,j}+\varepsilon}}_{w_{\text{prior}}}\cdot
\underbrace{\tfrac{1}{2}(\gamma_i+\gamma_j)}_{\text{prior}},
$$

makes the structure transparent. For an irreversible chain
($P_{ij}\neq P_{ji}$, e.g. any 3+ state model with non-detailed-balance):

- The **data estimate** $\texttt{gtotal\_ij}/P_{ij}$ is generically
  asymmetric in $i \leftrightarrow j$; it carries the directionality of
  the chain.
- The **data weight** $P_{ij}/(P_{ij}+\varepsilon)$ is also asymmetric:
  in the direction with larger $P$, the data dominates more than in the
  reverse direction.
- The **prior weight** $\varepsilon/(P_{ij}+\varepsilon)$ is the
  complementary asymmetric piece.
- Only the **prior value** itself is symmetric, and that is correct: it
  is the default in the absence of data, where there is no directional
  information to encode.

In the limit $P_{ij}, P_{ji} \to 0$ (rare transit in either direction),
the estimator falls back to the same symmetric value
$(\gamma_i + \gamma_j)/2$ — the right behaviour, because we have
observed neither transit and there is no basis for preferring one
direction over the other.

## Properties

| Property | Old `elemDivSoftAbs` | New Bayesian-shrinkage |
|---|---|---|
| Limit $P_{ij} \gg \varepsilon$ | $\to \texttt{gtotal\_ij}/P_{ij}$ | $\to \texttt{gtotal\_ij}/P_{ij}$ |
| Limit $P_{ij} \to 0$ | $\to 0$ (no physical content) | $\to \texttt{gmean\_ij\_prior}$ (endpoint average) |
| Continuity & smoothness in $P_{ij}$ | $C^\infty$ | $C^\infty$ |
| Denominator | $\sqrt{P^2+\varepsilon^2}$ (non-physical) | $P + \varepsilon$ (pseudo-count) |
| Sign behaviour | Always positive | Always positive (both terms $\ge 0$) |
| Derivative formula | Soft-abs (biased near zero) | Quotient rule (exact) |
| Physical interpretation | None | Conjugate prior, posterior mean |
| Parameter $\varepsilon$ | "regularization floor" | "prior pseudo-count" |

The denominator $P + \varepsilon$ is **strictly positive** for any
$\varepsilon > 0$, so the quotient is well-defined and FP-stable without
the soft-abs trick.

The derivative is the textbook quotient rule:

$$
\frac{\partial\,\texttt{gmean\_ij}}{\partial \theta}
\;=\;
\frac{1}{P + \varepsilon}\,
\frac{\partial \texttt{gtotal\_ij}}{\partial\theta}
\;+\;
\frac{\varepsilon}{P + \varepsilon}\,
\frac{\partial \texttt{gmean\_ij\_prior}}{\partial \theta}
\;-\;
\frac{\texttt{gmean\_ij}}{P + \varepsilon}\,
\frac{\partial P}{\partial\theta}.
$$

No soft-abs, no $\sqrt{\cdot}$ derivative. Clean AD propagation.

## Choice of $\texttt{min\_P\_prior}$

The pseudo-count $\varepsilon = \texttt{min\_P\_prior}$ controls the
threshold at which the data and prior contribute equally. Roughly:

- $P_{ij} \gg \varepsilon$: data dominates; behaviour identical to current
  code in the well-behaved regime.
- $P_{ij} \sim \varepsilon$: equal blend.
- $P_{ij} \ll \varepsilon$: prior dominates; the fallback is whatever
  endpoint-physics dictates.

A reasonable choice: $\varepsilon = 10^{-12}$ (the current `min_P` value).
This puts the "prior dominates" regime well below FP-noise on $P$ for any
realistic Markov chain, so the change is observable only on entries that
were already nonsensical under the old code. Tuning could be revisited if
specific schemes have legitimate small-$P_{ij}$ structure.

## Migration

1. Rename `min_P` → `min_P_prior` everywhere it appears in the model files
   and in the Vector_Space-typed Qdt / Qdtm / Qdtg / micro_Qdt* containers
   (`legacy/qmodel_types.h`, `legacy/micro_types.h`, `models_*.h`).
2. Replace every `elemDivSoftAbs(gtotal_ij, P, min_P)` and
   `elemDivSafe(gtotal_var_ij, P, min_P)` site in `legacy/qmodel.h`
   (~10 sites) with the Bayesian-shrinkage expression. Same for the
   corresponding `gvar_ij` line in each block.
3. Construct the prior matrices `gmean_ij_prior` and `gvar_ij_prior` once
   per Qdt evaluation from $\boldsymbol\gamma = $ `get<g>(m)`. They depend
   only on the per-state conductance vector, so a single helper

   ```cpp
   auto compute_g_endpoint_priors(const auto& v_g) {
       // returns pair<gmean_ij_prior_matrix, gvar_ij_prior_matrix>
   }
   ```

   suffices.
4. Drop `elemDivSoftAbs` from `legacy/matrix.h` if no other call sites
   remain.

## Appendix: adaptive choice of $\texttt{min\_P\_prior}$ from $\kappa(V)$

The fixed value $\texttt{min\_P\_prior}=10^{-12}$ is convenient but
arbitrary. A more principled choice is to tie the pseudo-count to the
actual numerical noise floor of $P_{ij}$, which depends on the
conditioning of the eigendecomposition used to construct $P$.

### Where the noise floor comes from

In the Padé / eigendecomposition path,

$$
P \;=\; V \cdot \mathrm{diag}\bigl(e^{\lambda\,dt}\bigr) \cdot W,
\qquad W = V^{-1},
$$

the reconstruction of $P$ from finite-precision $V$, $\lambda$, $W$
introduces componentwise error of order

$$
\bigl|\,\delta P_{ij}\,\bigr|
\;\lesssim\;
\varepsilon_{\text{machine}} \cdot \|P\|_\infty \cdot \kappa(V),
$$

where $\kappa(V) = \|V\|\,\|V^{-1}\|$ is the condition number of the
eigenvector matrix. An entry $P_{ij}$ whose true value falls below this
threshold is, after reconstruction, **indistinguishable from rounding
noise** — it carries no information regardless of what magnitude it
takes. For such an entry the data-side estimator
$\texttt{gtotal\_ij}/P_{ij}$ is meaningless and the prior should
dominate.

### Adaptive pseudo-count

Set

$$
\boxed{\;
\texttt{min\_P\_prior}
\;=\;
\varepsilon_{\text{machine}}\cdot\kappa(V).
\;}
$$

No tuning constants. The right-hand side *is* the FP noise floor of
$P_{ij}$ after the $V\cdot D\cdot W$ reconstruction, so setting the
pseudo-count equal to it places the data-vs-prior crossover at
signal-to-noise ratio $= 1$ — the canonical Bayesian choice. In the
limit $\kappa(V) = 1$ (perfectly conditioned), $\texttt{min\_P\_prior} =
\varepsilon_{\text{mach}}$, which is precisely the residual rounding of
a perfectly conditioned reconstruction.

A textbook-strict error analysis of the matrix product $V\cdot D\cdot W$
adds an $O(N)$ factor for accumulated rounding [Higham2002, §3.5]. In
practice $\kappa_F$ already overestimates $\kappa_2$ by enough to absorb
it. If a small-$N$ low-$\kappa$ case ever shows under-shrinkage,
multiplying by $N$ (the matrix dimension) is the principled tightening,
still parameter-free:
$$\texttt{min\_P\_prior} \;=\; \varepsilon_{\text{machine}}\cdot N\cdot\kappa(V).$$

The Frobenius bound is cheap and adequate:

$$
\kappa_F(V) \;=\; \|V\|_F \, \|W\|_F
\;=\;
\sqrt{\textstyle\sum_{i,j} V_{ij}^2}
\,\cdot\,
\sqrt{\textstyle\sum_{i,j} W_{ij}^2}.
$$

Two $O(N^2)$ passes over matrices already held by the calling Qdt code.
This bounds the spectral condition number from above:
$\kappa_2(V) \le \kappa_F(V)$, so using $\kappa_F$ in the formula is
conservative (gives a slightly larger pseudo-count than strictly needed,
which is the safe direction).

### Behaviour across regimes

| Regime | $\kappa(V)$ | Adaptive $\texttt{min\_P\_prior}$ |
|---|---|---|
| Symmetric / detailed-balanced chain | $\approx 1$ | $\approx \varepsilon_{\text{mach}} \sim 10^{-16}$ — essentially raw data |
| Generic non-reversible chain | $10^2 - 10^4$ | $10^{-12}$ to $10^{-10}$ — mild prior on rare entries |
| Near-defective (closely spaced eigenvalues) | $10^8$+ | $\sim 10^{-7}$ — prior dominates correctly |

The third regime is exactly when the eigendecomposition is numerically
suspect; the adaptive prior shrinks the unreliable entries toward the
physical endpoint average, where a fixed $10^{-12}$ would let the
numerical garbage through.

### Treatment in derivative-aware computation

$\kappa(V)$ depends on the model parameters $\theta$ through $V(\theta)$.
**Compute $\kappa$ from the primitive only and treat the resulting
$\texttt{min\_P\_prior}$ as constant in $\theta$.** This is the simplest
choice and is correct in the sense that:

- $\texttt{min\_P\_prior}$ is a *hyperparameter* of the regularization,
  not a free parameter of the model. Its role is to mark the threshold
  below which data is unreliable; that threshold is a property of the
  numerical reconstruction at the current $\theta$, not a quantity whose
  derivative is meaningful.
- Propagating $\partial \kappa / \partial\theta$ would require tracking
  $\partial \|V\|_F / \partial\theta$ and $\partial \|W\|_F /
  \partial\theta$ through the chain. The added complexity buys nothing
  at the FP-noise scale this constant operates at.
- Within the Bayesian-shrinkage expression, the *form* of the
  denominator $P_{ij} + \varepsilon$ still propagates derivatives
  correctly via the quotient rule; only the constant $\varepsilon$ is
  frozen.

### Implementation entry point

In each `calc_Qdt*_eig` (and the corresponding `_taylor` path that does
not use eigendecomposition — see open questions), inject

```cpp
double kappa_F     = get<kappa_V>(t_Eigs)();           // pre-computed in calc_eigen
double min_P_prior = std::numeric_limits<double>::epsilon() * kappa_F;
```

once per Qdt evaluation, then use `min_P_prior` in the shrinkage formula
for `gmean_ij` and `gvar_ij`. The Taylor / uniformization path does not
have a $\kappa(V)$ to reference; for those paths fall back to a fixed
$\texttt{min\_P\_prior}$ (the old $10^{-12}$) or use a Taylor-specific
error estimate from the truncation order.

## Open questions deferred to implementation

- **Should the prior use marginal `gmean_i` instead of per-state `v_g[i]`?**
  Marginal `gmean_i = Σ_j P_ij · gmean_ij` already integrates over the
  destination distribution, so using it as the prior introduces circular
  feedback (the prior depends on the very quantity being regularized).
  Per-state `v_g[i]` is parameter-driven and circularity-free, so it is the
  natural choice. The numerical difference is small except in pathological
  schemes.
- **Should the same shrinkage apply to `gtotal_ij` itself (the numerator
  *before* division)?** Current proposal: no. The numerator is well-defined
  for all $P_{ij}$, including zero, so no regularization is needed there.
  The shrinkage acts only on the ratio, where the FP pathology lives.
- **Interaction with `force_gmean_in_range`.** The current pipeline applies
  `force_gmean_in_range` after the division. With the shrinkage in place,
  `gmean_ij` is bounded by construction (it lies between
  $\bar\gamma_{i\to j}^{\text{prior}} \in [\gamma_{\min}, \gamma_{\max}]$
  and $\texttt{gtotal\_ij}/P_{ij}$, the latter of which `force_gtotal_in_range`
  already pins to that interval), so `force_gmean_in_range` may become
  redundant. Decide during implementation; conservative move is to keep it
  as a no-op check.
