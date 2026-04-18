# Score Autocorrelation, the Distortion Matrix, and Principled Correction of Approximate Likelihoods

Status: research note, companion to figure 2 of the eLife manuscript. Formulae are of
textbook-standard form; the specialization to the Macro_* family is where this note
contributes.

## 1. What we observed

Running the cross-correlation diagnostic
(`likelihood_derivative_cross_correlation_diagnostics`) on `dlikelihood_predictions`
for `scheme_CO` across the (algorithm, `interval_in_tau`, noise, `Num_ch`)
lattice reveals that the autocorrelation of the **score**
$\tilde{s}_t = \partial_\theta \log \tilde p(y_t \mid y_{1:t-1},\theta)$
is not the Dirac $\delta(\tau)$ that exact likelihood theory predicts.
Representative findings:

- `r_std` cross-correlation decays in ~3 steps at `interval_in_tau = 0.1`.
- `r_std` cross-correlation carries a lag-1 coefficient for `macro_IRV` at
  `interval_in_tau = 1`.
- `macro_MRV` shows a structural residual artefact near the agonist-concentration
  transition (at sample ≈ 500 of figure 2's protocol) that is absent in IRV, IRTV,
  and MNRV.
- The diagonal of the cross-correlation kernel is 1 at $\tau = 0$ by construction;
  the off-diagonal at $\tau = 0$ is non-zero and encodes the (normalized) Fisher
  Information matrix.

This note builds the theoretical framework, relates the kernel to the existing
`Information_Distortion_Matrix`, and proposes correction strategies.

## 2. Theoretical framework

### 2.1 The true score is a martingale difference sequence

For a correctly specified model with likelihood
$L(\theta) = \prod_{t=1}^{T} p(y_t \mid y_{1:t-1}, \theta)$,
the sample-wise score

$$
s_t(\theta) \;=\; \partial_\theta \log p(y_t \mid y_{1:t-1}, \theta)
$$

satisfies, under the data-generating law at $\theta$:

$$
\mathbb{E}[s_t \mid \mathcal{F}_{t-1}] \;=\; 0,
\qquad
\operatorname{Cov}[s_t, s_{t+\tau}] \;=\; \mathcal{I}_t(\theta)\,\delta(\tau),
$$

where $\mathcal{F}_{t-1} = \sigma(y_{1:t-1})$ and $\mathcal{I}_t$ is the per-sample
Fisher Information. Summing over $t$ gives the total Fisher Information
$\mathcal{I}(\theta) = \sum_t \mathcal{I}_t(\theta)$.

The martingale property is structural: it is forced by the fact that $p$ is a
proper conditional density. **Any measured deviation from
$\operatorname{Cov}[s_t, s_{t+\tau}] = 0$ at $\tau \neq 0$ implies the likelihood
is mis-specified or approximated.**

### 2.2 Approximate score decomposition

Let $\tilde p$ be the Macro_* family's approximation to the true conditional density.
Write

$$
\tilde s_t(\theta) \;=\; s_t(\theta) + \epsilon_t(\theta),
$$

where $\epsilon_t$ is the approximation-induced error in the score. $\epsilon_t$ is a
deterministic functional of the true hidden state trajectory $x_{1:t}$ and the
observed history $y_{1:t}$, because the algorithm's internal quantities
$\hat y_t, \hat V_t, \hat P_t$ depend on the history and approximate their exact
counterparts with an error determined by the algorithm's truncation scheme.

Crucially, $\epsilon_t$ inherits the memory of the hidden state: the state $x_t$ is
Markov with autocorrelation time $\tau_{\text{kin}}$ (the slowest kinetic time
constant), so $\epsilon_t$ and $\epsilon_{t+\tau}$ are correlated for
$|\tau| \lesssim \tau_{\text{kin}}/\Delta t$.

### 2.3 The cross-correlation kernel

Define the **score cross-correlation kernel** at lag $\tau$:

$$
K_{ij}(\tau) \;=\; \operatorname{Cov}[\tilde s_{i,t},\, \tilde s_{j,t+\tau}].
$$

Expanding the decomposition:

$$
K_{ij}(\tau) \;=\;
\underbrace{\delta(\tau)\, \mathcal{I}_{ij}(\theta)}_{\text{true Fisher}} \;+\;
\underbrace{\operatorname{Cov}[\epsilon_{i,t}, \epsilon_{j,t+\tau}]}_{\text{distortion kernel}} \;+\;
\underbrace{\operatorname{Cov}[s_{i,t}, \epsilon_{j,t+\tau}] + \operatorname{Cov}[\epsilon_{i,t}, s_{j,t+\tau}]}_{\text{cross-terms}}.
$$

Under weak assumptions (the algorithm's error is mean-zero conditional on
$\mathcal{F}_{t-1}$, i.e. $\mathbb{E}[\epsilon_t \mid \mathcal{F}_{t-1}] = 0$, which
holds for Macro_* at first order) the cross-terms vanish and

$$
K_{ij}(\tau) \;=\; \delta(\tau)\, \mathcal{I}_{ij}(\theta) \;+\;
\operatorname{Cov}[\epsilon_{i,t}, \epsilon_{j,t+\tau}].
$$

The measured normalized kernel (what `series_cross_correlation_kernel` emits, scaled
to $[-1, 1]$ by the unit-variance diagonal) is

$$
\rho_{ij}(\tau) \;=\;
\frac{K_{ij}(\tau)}{\sqrt{K_{ii}(0)\, K_{jj}(0)}}.
$$

The diagonal $\rho_{ii}(0) = 1$ by construction; the off-diagonal $\rho_{ij}(0)$ is
the Fisher-correlation coefficient between gradients w.r.t. parameters $i$ and $j$.

## 3. Relation to the existing distortion matrix

The `Information_Distortion_Matrix` in the codebase is the **zero-frequency (DC)
component** of the kernel's Fourier transform, i.e. the sum over all lags:

$$
\mathcal{J}_{ij}(\theta) \;=\; \sum_{\tau=-T+1}^{T-1} K_{ij}(\tau)
\;=\; \mathcal{I}_{ij}(\theta) \;+\;
\sum_{\tau \neq 0} \operatorname{Cov}[\epsilon_{i,t}, \epsilon_{j,t+\tau}].
$$

$\mathcal{J}$ is the quantity that actually drives the asymptotic variance of
estimators built from $\tilde L$. The existing diagnostic output reports
$\mathcal{J}$ but not $K(\tau)$ directly.

The **lag structure of $K$ carries strictly more information** than
$\mathcal{J}$: two algorithms can have identical $\mathcal{J}$ but very different
$K(\tau)$ profiles (one concentrating error energy at short lags, another at long
lags), and this matters for:

- Effective sample size (a short-memory kernel inflates variance by a small factor;
  a long-memory kernel by a large factor, even with the same integrated strength).
- Bootstrap validity (block-bootstrap block length should track $K$'s support).
- Model-selection: two algorithms with equal $\mathcal{J}$ but different $K$
  lag structure differ in *where* their approximation fails along the trace.

## 4. Principled corrections

Three distinct corrections, each consuming different functionals of $K$.

### 4.1 Sandwich / HAC variance of the MLE

When the score is autocorrelated, the asymptotic variance of the approximate MLE
$\hat\theta = \arg\max \tilde L(\theta)$ is **not** the inverse of any single
Fisher-like quantity. It is the sandwich estimator:

$$
\operatorname{Var}[\hat\theta] \;=\; \mathcal{J}(\theta)^{-1} \;\Omega(\theta)\; \mathcal{J}(\theta)^{-1},
\qquad
\Omega(\theta) \;=\; \sum_{\tau = -(T-1)}^{T-1} K(\tau).
$$

In time-series econometrics this is the **Newey–West / HAC estimator**. In
well-specified MLE, $\Omega = \mathcal{J} = \mathcal{I}$ and the sandwich collapses
to $\mathcal{I}^{-1}$. In the approximate case the three differ and only the
sandwich is correct.

An unbiased kernel-sum estimator with Bartlett weights is

$$
\hat\Omega \;=\; \hat K(0) + \sum_{\tau=1}^{m} w(\tau)\,\bigl(\hat K(\tau) + \hat K(\tau)^{\top}\bigr),
\qquad w(\tau) = 1 - \tau/(m+1),
$$

with bandwidth $m \sim T^{1/3}$ for consistency. The existing bootstrap output
supplies $\hat K(\tau)$ directly.

**Practical consequence.** Naive confidence intervals using
$\tilde{\mathcal I}^{-1}$ are narrower than the truth by a factor
$\sqrt{\lambda_{\min}(\Omega \mathcal{I}^{-1})}$ along the most-affected parameter
direction. For conditions where $\text{interval\_in\_tau} \to 1$ with algorithms
other than IRTV, this factor is meaningful (≥ 10% or more).

### 4.2 Effective sample size

For a scalar statistic the naive vs. effective sample sizes relate by

$$
N_{\text{eff}} \;=\; \frac{N}{1 + 2\sum_{\tau=1}^{\infty}\bar\rho(\tau)},
$$

where $\bar\rho(\tau)$ is a scalar summary (e.g. average of normalized diagonal
autocorrelations) of $K$. Useful for block-bootstrap design and cross-validation
fold sizing.

### 4.3 Score prewhitening

If $K(\tau)$ is known, construct the covariance operator
$\mathbf{K} \in \mathbb{R}^{Tp \times Tp}$ (block-Toeplitz, with block
$K(\tau)$ on the $\tau$-th block diagonal) and its Cholesky factor $L$:

$$
\mathbf{K} \;=\; L L^{\top},
\qquad
w_t \;=\; (L^{-1}\, \tilde s)_t.
$$

The whitened score $w_t$ is iid by construction. Inference on $w$ recovers the
correct Fisher Information and, crucially, can **improve point estimates** (not
just CIs) when the fit objective uses iterative score-matching or natural-gradient
steps, because the preconditioner is now the true covariance rather than the
inflated $\tilde{\mathcal I}$.

For short-memory kernels (support $\ll T$), $\mathbf{K}$ is banded and inversion
is $O(Tp^2 m)$ with bandwidth $m$ — tractable for the scales in figure 2.

## 5. Theoretical kernel for the Macro_* family

The Macro_* family differs in how each algorithm approximates the within-bin
state evolution. Schematically, let $\Delta t$ be the sample interval and
$\tau_{\text{kin}}$ the slowest kinetic time constant of the CTMC. Let
$h \equiv \Delta t / \tau_{\text{kin}}$.

| Algorithm | Mean treatment | Variance treatment | Leading error |
|---|---|---|---|
| `MNRV` | steady-state, no recursion | no recursion | $O(1)$ at all $h$ |
| `MRV`  | recursive, no within-bin integration | recursive | $O(h)$ |
| `IRV`  | recursive with 2nd-order within-bin integration | recursive | $O(h^3)$ |
| `IRTV` | IRV + Taylor variance correction | 4th-order | $O(h^4)$ or beyond |

For each algorithm, $\epsilon_t$ has a leading-order expression as a functional of
the exact filter residuals and the truncated Taylor remainder. For IRV the
canonical form is

$$
\epsilon_t^{\text{IRV}} \;=\; -\frac{\Delta t^3}{6}\, \partial_\theta
\bigl[\,\mathrm{tr}(Q^3 P_t) + \text{data-dependent terms}\bigr] \;+\; O(h^4),
$$

where $Q$ is the CTMC generator and $P_t$ is the conditional state distribution
at the start of bin $t$. The kernel

$$
\operatorname{Cov}[\epsilon_t^{\text{IRV}}, \epsilon_{t+\tau}^{\text{IRV}}]
$$

is then a function of the two-time joint distribution of
$(P_t, P_{t+\tau})$, which for a CTMC decays as
$\exp(-|\tau|\Delta t / \tau_{\text{kin}})$. This immediately predicts

$$
K_{ij}^{\text{IRV}}(\tau) \approx K_{ij}^{\text{IRV}}(0) \cdot
e^{-|\tau|\,h},
\qquad
K_{ij}^{\text{IRV}}(0) \propto h^{6}.
$$

**Testable prediction.** The empirical kernel at fixed algorithm should collapse
onto a single curve when plotted against $\tau \cdot h$ (time-lag in units of
kinetic time constant), and its amplitude should scale as $h^{6}$ for IRV,
$h^{2}$ for MRV, and $h^{8}$ for IRTV. Verifying this collapse across your
lattice of `interval_in_tau` values is a direct test of the theoretical derivation.

Analogous expressions hold for MRV and IRTV with different leading orders in $h$.
Writing them out explicitly (with the data-dependent terms) is the work that
converts the current "measure and report" diagnostic into a predictive theory.

## 6. Paper-level contribution

What you have empirically:

- Validated bootstrap kernel $\hat K(\tau)$ across a dense lattice of
  $(\text{algorithm}, \text{interval\_in\_tau}, \text{noise}, \text{N}_{\text{ch}})$.
- The integrated distortion $\hat{\mathcal J}$ (as `Information_Distortion_Matrix`).
- Clear algorithm-wise separation: IRTV ≳ IRV ≳ MNRV ≳ MRV at large $h$, with MRV
  failing specifically at agonist-concentration transitions.

What the theoretical framework adds:

- **A priori prediction** of $K^{\text{theory}}(\tau)$ per algorithm, so users
  choose the right Macro_* variant *before* running the expensive bootstrap.
- **Unified explanation** of the four algorithms as successive orders in a single
  expansion parameter $h$. The paper gets to say "Macro_IRV is the $O(h^2)$
  truncation; here is the next-order correction, and here is the empirical regime
  where each is valid".
- **HAC-corrected confidence intervals** as the recommended reporting standard,
  with the existing `Information_Distortion_Matrix` as a special case
  ($\Omega = \mathcal J$ would be the naive read).

## 7. Immediate practical steps

Ranked by cost / payoff:

1. **HAC variance from existing bootstrap.** Compute
   $\hat\Omega = \sum_\tau \hat K(\tau)$ with Bartlett weights from the kernels you
   already write to CSV. Compare to $\hat{\mathcal J}$ (the diagonal-only /
   zero-lag version) per (algorithm, interval, noise, N_ch). Publish the ratio.
   No new C++ code beyond an R or Python post-processor on the existing CSV.

2. **Dimensional collapse plot.** Plot $\hat K(\tau) / \hat K(0)$ against
   $\tau \cdot h$ for each algorithm, overlaying all `interval_in_tau` values.
   If the curves collapse, the leading-order expansion is valid and you have a
   one-parameter family per algorithm.

3. **Amplitude scaling plot.** Plot $\hat K(0)$ or $\operatorname{tr}\hat\Omega$
   against $h$ on log-log axes per algorithm. Slopes should match the theoretical
   orders $\{2, 6, 8\}$ for $\{\text{MRV}, \text{IRV}, \text{IRTV}\}$.

4. **Analytic derivation.** Write out $\epsilon_t^{\text{MRV}}$,
   $\epsilon_t^{\text{IRV}}$, $\epsilon_t^{\text{IRTV}}$ to leading order, and
   compute their two-time covariances explicitly for `scheme_CO`. Verify
   agreement with the empirical kernels.

5. **Score prewhitening experiment.** Re-fit `scheme_CO` with the prewhitened
   score on a small subset of conditions. Compare MLE bias and variance to the
   naive fit. This is the strongest form of the correction — if the prewhitened
   fits are more accurate than the naive ones, the kernel is not just a
   diagnostic but a correction lever for the estimator itself.

Steps 1–3 are post-processing of existing data and should land within a day. Step
4 is the theoretical core of the paper extension. Step 5 is more ambitious and
may deserve its own short companion paper.

## Appendix A. Glossary

| Symbol | Meaning |
|---|---|
| $y_t$ | observed current at sample $t$ |
| $x_t$ | hidden channel state at time $t$ |
| $\theta$ | parameter vector (on, off, unitary_current, noise, baseline, N_ch) |
| $s_t, \tilde s_t$ | exact and approximate per-sample score |
| $\epsilon_t$ | approximation error in the score |
| $\mathcal I$ | Fisher Information matrix |
| $\mathcal J$ | integrated distortion matrix ($= \sum_\tau K(\tau)$) |
| $K(\tau)$ | score cross-correlation kernel at lag $\tau$ |
| $\rho(\tau)$ | normalized kernel ($K$ rescaled to diagonal 1 at $\tau=0$) |
| $\Omega$ | HAC long-run variance estimator |
| $h$ | $\Delta t / \tau_{\text{kin}}$, the dimensionless sampling ratio |
| $\tau_{\text{kin}}$ | slowest kinetic time constant of the CTMC |

## Appendix B. Assumptions used above

- The data-generating model is correctly specified (in particular, the CTMC
  structure of `scheme_CO` is assumed exact); all error is algorithmic, not
  structural. If there is additional structural misspecification, $\epsilon_t$
  picks up an additional contribution with its own kernel; the corrections above
  still apply but the decomposition into "Fisher + distortion" breaks.
- $\mathbb{E}[\epsilon_t \mid \mathcal F_{t-1}] = 0$ to leading order. This holds
  for the Macro_* family by construction (each algorithm preserves the conditional
  mean exactly; the error enters through variance and higher moments). It breaks
  for crude approximations (e.g. mean-field), where an additional bias correction
  is needed.
- The experimental design is stationary within each agonist pulse regime. Near
  concentration transitions, the kernel has a non-stationary component
  (visible as the MRV artefact at sample ≈ 500). Stationary HAC theory needs a
  localized variant there — see Robinson (1998) or Dahlhaus (2012) for
  locally-stationary extensions.
