# Two Paths for Summarizing Score Autocorrelation

Status: research note, companion to
[ADR-005](../adr/ADR-005-likelihood-diagnostic-presets.md) and the diagnostic
refactor in [legacy/moment_statistics.h](../../legacy/moment_statistics.h).

## 1. Setup

Let $V_t$ be a leaf observable of the diagnostic pipeline at time step $t$
($V \in$ {`r_std`, `dlogL`, `y_mean`, `y_var`, `trust_coefficient`, `logL`,
`elogL`, `Gaussian_Fisher_Information`}; scalar, vector, or matrix-valued). Let

$$
C_\tau(i) \;=\; \operatorname{Cov}(V_i, V_{i+\tau})
$$

be the lag-$\tau$ per-sample covariance block at time $i$, computed over the
Monte Carlo sample dimension; $C_0(i)$ is the per-sample variance block, the
lag-zero entry of the variance kernel.

Two distinct, commonly-conflated scalar summaries of the autocorrelation
structure are useful for driving the corrections outlined in
[score-autocorrelation-and-distortion-correction.md](score-autocorrelation-and-distortion-correction.md):

- **Path 1 (forward sum of normalized correlations).** A per-sample quantity
  that measures how far forward in time each sample is still correlated with
  itself.
- **Path 2 (Kish–Sokal sandwich / integral correlation lag).** A variance
  inflation factor applied to the total sum, in the form
  $C_0^{-1/2}\,\Omega\,C_0^{-1/2}$.

These are **different quantities**, not alternate computations of the same
summary. Both are stored in the report types; previous code conflated them
into one field labeled "evolution" which obscured the interpretation.

## 2. Path 1: `cross_correlation_lag_forward`

Per sample $i$, normalize the cross-correlation at lag $\tau \geq 0$:

$$
\rho_\tau(i) \;=\; C_0(i)^{-1/2}\,C_\tau(i)\,C_0(i+\tau)^{-1/2}
$$

so that $\rho_0(i) = I$ by construction. Then define the **forward sum**

$$
F(i) \;=\; \sum_{\tau=0}^{\text{max\_lag}} \rho_\tau(i).
$$

For scalar $V$, $F(i)$ is a scalar measuring how many effective independent
samples-per-step at time $i$. For vector / matrix $V$, $F(i)$ is a block of
the same shape as $\rho_\tau(i)$.

**Interpretation.** $F(i)$ tells you, for a single time step $i$, the total
correlation carried into the future. It is naturally per-sample (one $F$ per
$i$); aggregating over $i$ is not the point. Used for detecting **where in the
protocol** the score fails to decorrelate (e.g. agonist transitions, for
`macro_MRV`).

Computed by `series_cross_correlation_lag_forward<Id>(maybe_corr.value())`;
stored as the `cross_correlation_lag_forward<V>` member of `Report_local_var`,
`Report_local_cov`, and `Report_cross`.

## 3. Path 2: `integral_correlation_lag` (Kish–Sokal form)

Sum the cross-correlation blocks over the **symmetric lag window**, sandwiched
by the total variance:

$$
\Omega \;=\; \sum_{i} \sum_{\tau=-\text{max\_lag}}^{\text{max\_lag}} C_\tau(i),
\qquad C_{\text{tot}} \;=\; \sum_i C_0(i),
$$

$$
K \;=\; C_{\text{tot}}^{-1/2}\,\Omega\,C_{\text{tot}}^{-1/2}.
$$

**Interpretation.** $K$ is the variance-inflation matrix for the sum
$\sum_i V_i$: if the $V_i$ were iid, $\operatorname{Var}(\sum V_i) = C_{\text{tot}}$
and $K = I$; correlation between samples gives $K \neq I$. This is the Sokal
autocorrelation window applied to vector/matrix-valued series, and is the form
that enters Laplace-approximated Bayes-factor corrections.

Computed by `series_integral_correlation_lag<Id>(stats, rtol, atol)`; stored
as the `integral_correlation_lag<V>` member of `Report_integral`,
`Report_local_var`, `Report_local_cov`, and `Report_cross`. Present in every
preset since the base reports all need it.

Path 2 is **bidirectional** (symmetric lag window, including $\tau < 0$) while
Path 1 is **forward-only**. For a time-reversal-symmetric process both give
related values; in general they carry different information, and both are
kept.

## 4. The diagonal-variance scatter for matrix-valued $V$

For vector-valued $V$ (e.g. `dlogL`, a p-vector), the per-sample variance
block $C_0(i)$ is a p×p symmetric matrix and `Report_local_var` stores its
diagonal as a p-vector — a factor-of-p saving over `Report_local_cov`, which
stores the full p×p block.

For matrix-valued $V$ (`Gaussian_Fisher_Information`, a p×p SymPosDef), the
outer-product path inside `prod_XY<true>` first flattens via `to_vector_half`
to an m-vector (m = p(p+1)/2, packed upper triangle), then squares via
`Lapack_Product_Self_Transpose_vectorized` to an **m×m** SymPosDef block. The
diagonal of that m×m block is m scalar variances — one per unique entry of
$V$.

To keep `Report_local_var<GFI>`'s stored type matching
`mean_value_type_t<GFI>` = `SymPosDefMatrix<double>` of shape p×p, we
**scatter** those m diagonal values back into a p×p symmetric matrix:

$$
\mathrm{diag\_mat}(i,j) \;=\; \text{(m×m block)}(k,k)
\quad\text{where } k = \text{packed index of } (i,j).
$$

The packed index follows `to_vector_half`: for $0 \leq i \leq j < p$,

$$
k \;=\; \frac{i(2p - i + 1)}{2} + (j - i),
\qquad p \;=\; \frac{\sqrt{8m+1} - 1}{2}.
$$

The result stores p(p+1)/2 unique values in natural p×p shape — no
zero-padding, no information loss, same type as the non-matrix reports expect.

Implemented in `extract_series_diagonal_variance` in
[legacy/moment_statistics.h](../../legacy/moment_statistics.h) as the third
constexpr branch (detected by `value_type(size_t, size_t, double)` and `.set`
being available on the entry type).

## 5. How this maps onto the presets

| Preset | Path 1 stored | Path 2 stored | Per-sample variance storage |
| --- | --- | --- | --- |
| `basic` | — | — (via base `integral_correlation_lag` is not attached; the base carries distortion matrices only) | — |
| `series_var` | yes (per observable) | yes (per observable) | diagonal only (`Report_local_var`) |
| `series_cov` | yes | yes | full p×p block for vector V (`Report_local_cov`); p×p scatter for GFI |
| `series_kernel` | yes | yes | full kernel over lags (`Report_cross`) |
| `series_kernel_full` | yes | yes | kernel + GFI scatter + per-sample derived |

## 6. Why both paths matter

Path 2 drives the **correction of the total log-likelihood** in the
Laplace-Bayes-factor sense: the inflation factor $K$ says how much to widen
the posterior once the approximate likelihood's score has lost martingale
orthogonality. Path 1 drives **diagnosis of where** the correlation is coming
from — which sample, which protocol segment, and thus which approximation
choice is responsible. Reporting one without the other would either (a) apply
a global correction with no insight into its origin, or (b) localize the
problem without being able to translate that into a corrected evidence.
