# Cheap `gvar_ij` from gtotal_ij marginals (Empirical-Bayes shortcut)

Sibling note to
[bayesian_prior_regularization_of_Qdt.md](bayesian_prior_regularization_of_Qdt.md).
That note settles on v_g (per-state, "Poisson endpoint") priors for the
production path because they are robust against ill-conditioned
eigendecompositions. The present note describes an Empirical-Bayes
(EB) alternative that is **cheaper to compute** and may be preferred in
contexts where:

- the full second-moment integral `gtotal_sqr_ij` is not needed,
- `gvar_ij` is consumed only via marginals (`gvar_i`),
- or the eigendecomposition is known to be well-conditioned (small
  `κ(V)`).

The two approaches are not mutually exclusive: the primary
`calc_Qdt_eig` path uses v_g, but `calc_Qdtm_eig` (or any other code
path that needs `gvar_ij` cheaply and tolerates a mildly less informative
prior) can use the EB shortcut described here.

## What gets skipped

The full pipeline (v_g + back-conversion) for `gvar_ij` requires the
second-moment integral `gtotal_sqr_ij`:

$$
\texttt{gtotal\_sqr\_ij}
\;=\;
2 \cdot V \cdot \mathrm{WgV\_E3} \cdot W,
\qquad
\mathrm{WgV\_E3}[k_1,k_3]
\;=\;
\textstyle\sum_{k_2}
W_{k_1, k_2}\,\gamma_{k_2}\,V_{k_2, k_3} \;\cdot\;
E_3(\lambda_{k_1}, \lambda_{k_2}, \lambda_{k_3})_{dt}.
$$

This is a triple-loop **O(N³)** construction whose inner term is the
divided difference $E_3$ (a delicate Padé-cancellation expression), then
two **O(N³)** matrix products on top. Skipping it is the saving.

The EB shortcut avoids `gtotal_sqr_ij` entirely. The only inputs are
quantities already produced by the first-moment path:

- $\texttt{gtotal\_ij}\;=\;V\cdot\mathrm{WgV\_E2}\cdot W$ (already computed),
- $P_{ij}$ (already computed).

## Derivation

In the no-data limit ($P_{ij}\to 0$), an EB prior on $\bar\gamma_{i\to j}$
shrinks the entry toward the average of the two marginal mean
conductances seen by the chain:

$$
\overline{g}_i^{\text{(out)}} \;=\; \textstyle\sum_k \texttt{gtotal\_ij}[i,k]
\;=\; \mathbb{E}\bigl[\bar\gamma \mid X_0=i\bigr],
$$

$$
\overline{g}_j^{\text{(in)}}  \;=\; \textstyle\sum_k \texttt{gtotal\_ij}[k,j]
\;=\; \mathbb{E}\bigl[\bar\gamma \mid X_t=j\bigr],
$$

$$
\bar\gamma_{i\to j}^{\text{prior}} \;=\; \tfrac{1}{2}\bigl(
   \overline{g}_i^{\text{(out)}} + \overline{g}_j^{\text{(in)}}\bigr).
$$

The corresponding variance prior, by the **law of total variance** applied
to the two-marginal mixture, has two contributions: the mean of the
conditional variances, and the variance of the conditional means. The
**second contribution alone** is constructible from these two marginal
means — no second-moment integral required:

$$
\boxed{\;
\mathrm{Var}^{\text{prior, cheap}}\!\bigl(\bar\gamma_{i\to j}\bigr)
\;=\;
\biggl(\frac{\overline{g}_i^{\text{(out)}}
              - \overline{g}_j^{\text{(in)}}}{2}\biggr)^{\!2}.
\;}
$$

This captures "how much the two marginal mean estimates disagree" at
endpoint pair $(i,j)$. In the limit where the marginal means agree (the
chain looks essentially the same from either direction), the prior
variance is zero; the larger the disagreement, the larger the prior
uncertainty.

## What this misses

The full LTV formula includes a "mean of conditional variances" term
proportional to $\sum_k\texttt{gtotal\_var\_ij}[i,k] + \sum_k\texttt{gtotal\_var\_ij}[k,j]$.
That term requires `gtotal_var_ij`, which requires the second-moment
integral. The cheap prior drops it. The omission is benign when:

- the typical within-state conductance variance is small relative to
  the inter-state variance (states are conductance-distinct),
- the consumer of `gvar_ij` only needs `gvar_i` marginals (where the
  full LTV reassembles correctly anyway via $\texttt{gvar\_i} = \texttt{gsqr\_i} - \texttt{gmean\_i}^2$),
- or the prior weight is small ($P_{ij}\gg\varepsilon$) so this term is
  damped out before it influences the posterior.

The omission is significant when:

- states have substantial within-state conductance dispersion,
- the consumer reads `gvar_ij` entries directly in the small-$P_{ij}$
  regime,
- or shrinkage is heavy (large $\varepsilon/P_{ij}$ ratio) and the
  dropped term would dominate the kept one.

For the macro IRT/MRT and micro paths in this codebase, the dominant
consumer of stored `gvar_ij` is the `gvar_i` marginal, where the full
LTV is reconstructed downstream — so the cheap prior is a reasonable
fit.

## Putting it in the Bayesian update

Same conjugate-prior posterior as the main note:

$$
\texttt{gvar\_ij}_{i,j}
\;=\;
\frac{\texttt{gtotal\_var\_ij}_{i,j} \;+\;
      \mathrm{Var}^{\text{prior, cheap}}\!\bigl(\bar\gamma_{i\to j}\bigr)
      \cdot \varepsilon}
     {P_{i,j} + \varepsilon}.
$$

In the small-$P_{ij}$ limit this falls back to
$((\overline{g}_i^{\text{(out)}} - \overline{g}_j^{\text{(in)}})/2)^2$ —
the EB-cheap fallback. In the large-$P_{ij}$ limit it recovers the raw
quotient unchanged. The same pseudo-count $\varepsilon = \texttt{min\_P\_prior}$
applies (eig path: $\varepsilon_{\text{mach}}\cdot\kappa(V)$; non-eig:
$10\sqrt{N}\,\varepsilon_{\text{mach}}$).

## Compute summary

| Path | Construction cost | What's stored |
|---|---|---|
| Full v_g + back-conversion | one E2 build + one E3 build = **2×O(N³)** plus three matrix products | `gtotal_ij`, `gtotal_sqr_ij`, `gmean_ij`, `gvar_ij`, marginals |
| **EB-cheap variance prior** | **one** E2 build (E3 skipped) + one matrix product | `gtotal_ij`, `gmean_ij`, `gvar_ij` *(EB-shrunken)*, marginals |

For typical channel-model sizes ($N \le 10$) the absolute saving is
small; the practical gain is in keeping `calc_Qdtm_eig` lean and
matching its existing interface (no `gtotal_sqr_ij` field in `Qdtm`).

## Caveat: κ(V) sensitivity

The marginals
$\overline{g}_i^{\text{(out)}} = \sum_k \texttt{gtotal\_ij}[i,k]$ are
sums of FP-noise-contaminated entries when $V$ is ill-conditioned
(see the main note for the $\kappa(V)$ argument). The contamination
partially cancels in the sum — each entry contributes only $1/N$ — but
the cancellation is incomplete in pathological regimes. Apply this
shortcut only on call paths where the eigendecomposition's condition
number is known to be modest, or where the prior contribution is small
relative to the data ($\varepsilon \ll P_{ij}$) so contamination is
damped out.

## Suggested helpers (not yet implemented)

Parallel to the existing v_g helpers at
[qmodel.h:1132-1165](../../../legacy/qmodel.h#L1132):

```cpp
// EB-cheap prior on the conditional variance, built from row/col
// marginals of raw gtotal_ij. Skips the gtotal_sqr_ij computation
// entirely — see theory/macroir/notes/Gmean_ij_gvarij/cheap_gvar_ij_via_EB.md.
template <class C_gtotal_ij>
    requires U<C_gtotal_ij, gtotal_ij>
static auto gvar_ij_eb_prior_cheap(const C_gtotal_ij& t_gtotal_ij) {
    const std::size_t N = t_gtotal_ij().nrows();
    Matrix<double> u(N, 1, 1.0);
    Matrix<double> uT(1, N, 1.0);
    auto row_marg = t_gtotal_ij() * u;        // (N×1)  Σ_k gtotal_ij[i,k]
    auto col_marg = uT * t_gtotal_ij();       // (1×N)  Σ_k gtotal_ij[k,j]
    auto diff = row_marg * uT - u * col_marg;  // (N×N) (i,j) = out_i − in_j
    return build<gvar_ij>(elemMult(diff, diff) * 0.25);   // ((out − in)/2)²
}
```

Callers feed this into the existing `calc_g_ij_bayes` overload that
accepts `(gtotal_sqr_ij, gvar_ij)` — except in the EB-cheap path the
"data" argument is `gtotal_var_ij` (which still needs `gmean_ij` to
build, but does not need `gtotal_sqr_ij` if the LTV-mean-of-variances
term is omitted on the data side too). For a fully gsqr-free pipeline,
take the data term to be zero on the entries below the noise floor —
i.e., let the prior dominate without contest in the regularized regime,
and use the raw quotient elsewhere.
