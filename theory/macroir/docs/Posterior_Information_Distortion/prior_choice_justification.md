# Prior Choice Justification — Posterior Information Distortion Framework

## Scope

This document justifies the **structural role** of a proper Gaussian prior in
the Posterior Information Distortion (PID) framework, the practical default
of a **one-decade prior in log10 transformed coordinates**, and the operational
consequences for parameter inference and algorithm validation in the eLife
2025 manuscript.

It does NOT propose a new prior class, does NOT relitigate the choice of
log10 transformation (already standard in `legacy/parameters.h` via
`Log10_Tr`), and does NOT claim the prior is informative in a
subjective-Bayes sense. The prior is **always present** as a structural
element of the framework, not as an ad-hoc regularizer applied selectively.

## 1. The framework is Posterior, not Likelihood

The PID framework treats the **posterior** as the primary object of
inference, not the likelihood. Every diagnostic
([[supplement-posterior-main]], [[supplement-posterior-evidence]],
[[supplement-information-gain]]) is defined on
$\mathbf H_{\mathrm{post}} = \mathbf H_{\mathrm{lik}} + \mathbf H_{\mathrm{prior}}$
and $\mathbf{JT}_{\mathrm{post}} = \mathbf{JT}_{\mathrm{lik}} + \mathbf H_{\mathrm{prior}}$,
not on $\mathbf H_{\mathrm{lik}}$ and $\mathbf{JT}_{\mathrm{lik}}$ alone.

This is a structural commitment, not a numerical convenience. The
consequence: **a proper prior is always required**. The "Likelihood-only"
framework documented in
`theory/macroir/docs/Likelihood_Information_Distortion/` exists as the
historical precursor and as a degenerate limit
$\mathbf H_{\mathrm{prior}} \to \mathbf 0$ of the present framework, but is
not the operational target.

The reason is empirical and structural at once. The likelihood Hessian
$\mathbf H_{\mathrm{lik}}$ computed via finite differences of the AD score
(`calculate_mnumerical_fisher_information` in
[`src/core/likelihood.cpp`](../../../../src/core/likelihood.cpp)) can be
indefinite at some parameter regions — not because of numerical noise
(see `calculate_mnumerical_fisher_information` audit in this session),
but because the likelihood surface is genuinely non-concave there for
the algorithm in question. The empirical case of macro_R at
$\Delta t = 1\tau$ where $\lambda_{\min}(\mathbf H_{\mathrm{lik}}) \approx -42$
across 100 bootstrap replicates illustrates this: it is the algorithm's
property, not the estimator's. See
[[supplement-algorithm-validation]] for the empirical evidence.

In that regime, the "what value of curvature should I use" question is
not answerable without auxiliary information. The prior IS that auxiliary
information. Making it always present, rather than only invoked when
$\mathbf H_{\mathrm{lik}}$ fails, makes the framework operationally
uniform: every diagnostic has the same meaning everywhere, with no
case-distinction between identifiable and non-identifiable directions.

## 2. Three roles the prior plays

The prior plays three distinct roles in the PID framework, each
mathematically necessary:

### 2.1 Regularizer of curvature

$\mathbf H_{\mathrm{post}}$ is positive definite for any proper prior
satisfying the dominance condition
[[supplement-posterior-main#sec:dominance-condition]]:
$\mathbf H_{\mathrm{prior}} \succ \mathbf H_{\mathrm{lik}}^{-}$ (where
$\mathbf H_{\mathrm{lik}}^{-}$ is the absolute magnitude of the negative
eigenpart of $\mathbf H_{\mathrm{lik}}$). When this condition holds,
$\mathbf C_{\mathrm{post}}$ is a well-defined SPD-normalized matrix
globally, with no retained-eigenspace machinery required.

The dominance condition has a checkable sufficient form:
$\lambda_{\min}(\mathbf H_{\mathrm{prior}}) > \max(0, -\lambda_{\min}(\mathbf H_{\mathrm{lik}}))$.
This is a scalar test, computable in two eigenvalue evaluations.

### 2.2 Deterministic augmented observation

The prior contributes to the curvature
($\mathbf H_{\mathrm{post}} = \mathbf H_{\mathrm{lik}} + \mathbf H_{\mathrm{prior}}$)
and to the total score covariance
($\mathbf{JT}_{\mathrm{post}} = \mathbf{JT}_{\mathrm{lik}} + \mathbf H_{\mathrm{prior}}$)
additively, propagating the information identity to the posterior
without contributing sampling variance. The prior is, mathematically, a
"single deterministic observation" with infinite weight; see
[[supplement-posterior-main]] Section 3 for the bookkeeping.

### 2.3 Replacement for arbitrary subspace cutoff

The Likelihood-only path requires a numerical tolerance
(`rtol=1e-10`, `atol=0`) to define the active subspace in which the
distortion matrix is computed; see
[[implementation-subspace-information-distortion]]. This tolerance is
arbitrary: it has no statistical interpretation and produces
discontinuous diagnostics across its boundary.

The prior replaces this with an objective scale: $\lambda(\mathbf H_{\mathrm{prior}})$
defines what counts as "the data does not inform this direction".
When the prior dominates and the likelihood does not, the diagnostic
reports the prior — a well-defined statistical statement rather than a
numerical truncation.

## 3. The log10 transformation and "diagonal in transformed coordinates"

All eLife models use the parameter transformation infrastructure of
`legacy/parameters.h`:

- **`Log10_Tr`**: applied to rate constants, conductance amplitudes, and
  other quantities that are positive and span many orders of magnitude
- **`Identity_Tr`** (= Linear): applied to current baseline, leakage,
  and other quantities expected to live on a finite linear scale
- **`Fixed_Tr`**: for parameters held at a fixed value

The prior covariance in the PID framework is **diagonal in the transformed
(log10) coordinates**. This has two consequences:

- $\lambda_{\min}(\mathbf H_{\mathrm{prior}}) = 1 / \sigma_{\max}^2$ where
  $\sigma_{\max}$ is the largest transformed standard deviation — trivially
  computable from the diagonal of the prior covariance, no eigendecomposition
  required.
- $\lambda_{\min}(\mathbf H_{\mathrm{prior}})$ is direction-independent in
  the transformed basis. For the Stage-A dominance test of
  [[supplement-posterior-main#sec:dominance-condition]], this gives a
  scalar bound directly.

The transformation choice itself is justified by physical considerations:
ion-channel rate constants span 4–7 orders of magnitude across kinetic
schemes (the eLife paper covers Markov rates from $10^{-2}$ to $10^{5}\,
\mathrm{s}^{-1}$). A linear prior on a quantity that spans seven decades
would be either uninformative everywhere (large $\sigma$, no
regularization) or pathologically restrictive locally (small $\sigma$,
forces particular values). Log10 transformation lets a single
$\sigma_{\mathrm{transformed}}$ provide uniform, scale-aware
regularization across all rate constants.

## 4. Default: one decade of prior, $\sigma_{\mathrm{transformed}}^2 = 2$

The default for the eLife 2025 manuscript is

$$\sigma_{\mathrm{transformed}}^2 = 2 \quad \Longleftrightarrow \quad
\sigma_{\mathrm{transformed}} \approx 1.41 \quad \Longleftrightarrow \quad
\text{95\% prior interval} \approx \pm 2.77 \text{ decades}.$$

This is what the existing
[`projects/eLife_2025/scheme_1_prior.csv`](../../../../projects/eLife_2025/scheme_1_prior.csv)
and sibling files specify uniformly. Sometimes summarized loosely as
"one decade of prior", though strictly speaking the 95% interval is
closer to three decades.

The choice is justified on three grounds:

### 4.1 Physical magnitude argument

A 95% prior interval of $\pm 2.77$ decades covers the published spread
of kinetic parameters in ion-channel literature for the channels under
study. It is wide enough that the prior does not contradict any value
reported in independent biophysics, but narrow enough that it excludes
unphysical regimes (e.g., rates $> 10^9\,\mathrm{s}^{-1}$ for typical
allosteric transitions).

### 4.2 Dominance-condition argument

With $\sigma_{\mathrm{transformed}}^2 = 2$,
$\lambda_{\min}(\mathbf H_{\mathrm{prior}}) = 1/2$. Empirically, in the
figure_2 cells where $\mathbf H_{\mathrm{lik}}$ is positive definite for
macro_IR, $\lambda_{\min}(\mathbf H_{\mathrm{lik}})$ ranges from
$\approx 5$ (cell 1, $\Delta t = 1\tau$) to $\approx 480$ (cell 7,
$\Delta t = 0.01\tau$). In these cells the prior is dominated by the
likelihood: it contributes regularization but does not shape inference.
The posterior is data-dominated.

In the pathological cell (macro_R at $\Delta t = 1\tau$), where
$\lambda_{\min}(\mathbf H_{\mathrm{lik}}) \approx -42$, the scalar
dominance bound fails: $1/2 < 42$. The MATRIX-level dominance must be
checked direction-by-direction — and the empirical finding that
$\mathbf G_{\mathrm{lik}}$ is singular in the same direction means that
the prior IS what defines the posterior in that direction. This is
exactly what the prior is supposed to do: take over where the data is
silent.

### 4.3 Sensitivity argument

The correct-specification identity guarantees
$\mathbf C_{\mathrm{post}} = \mathbf I$ at $\theta_0$ for ANY proper
prior strength
[[supplement-posterior-main#sec:posterior-c-identity]]. The choice of
$\sigma^2 = 2$ only affects $\mathbf C_{\mathrm{post}}$ in regions of
misspecification, where the prior interacts with the likelihood
mismatch. Section 5 below specifies how to verify that the choice does
not unduly drive inference.

## 5. Sensitivity analysis recipe

Required as a robustness companion to any inference reported with
$\sigma^2 = 2$. The recipe produces a single supplementary figure.

### 5.1 Protocol

For one representative cell (recommended: macro_IR at $\Delta t = 0.1\tau$,
median-noise regime), rerun with the prior variance taking values in
$\{0.5, 1.0, 2.0, 4.0, 8.0\}$. Report, per variance:

1. $\lambda_{\min}(\mathbf H_{\mathrm{post}})$ range across replicates
2. Fraction of replicates in each safety category (see
   [[supplement-safety-categorization]])
3. Median $[\mathbf C_{\mathrm{post}}]_{ii}$ for the headline parameter
4. Median DCC scalar diagnostic (log-det)

### 5.2 Expected behavior

By the limiting regimes in [[supplement-posterior-main]] Section 2.3
("flat-prior limit" and "informative-prior limit"):

- As $\sigma^2 \to \infty$: posterior diagnostics converge to the
  Likelihood IDM diagnostics. If the Likelihood IDM is finite (no
  singularity), the convergence is smooth; if the Likelihood IDM
  diverges (singular $\mathbf H_{\mathrm{lik}}$), the posterior
  diagnostics blow up too.
- As $\sigma^2 \to 0$: $\mathbf C_{\mathrm{post}} \to \mathbf I$ from
  the inside (Section 5.3 of main supplement). Inference contracts to
  the prior mean.

For a well-identified model, the diagnostics should be **approximately
flat** across $\sigma^2 \in [0.5, 8]$ — a two-orders-of-magnitude
sensitivity sweep. Sharp dependence on $\sigma^2$ in this range
indicates the prior is materially shaping inference, which is a flag
for prior dominance (= weak data) in the affected direction.

### 5.3 Empirical companion to identity

This is the empirical version of the population-level guarantee
$\mathbf C_{\mathrm{post}}(\theta_0) = \mathbf I$ for any proper prior.
The sensitivity figure demonstrates how closely the finite-sample
estimator satisfies the same invariance.

## 6. Posterior IDM versus Likelihood IDM at the chosen prior

A practical consequence of $\sigma^2 = 2$: in all figure_2 cells where
$\mathbf H_{\mathrm{lik}}$ is well-conditioned (the vast majority), the
Posterior IDM ($\mathbf C_{\mathrm{post}}$) and the Likelihood IDM
($\mathbf C_{\mathrm{lik}}$) differ by at most a few percent. The two
frameworks agree numerically where both are defined.

The crucial difference is at the boundary: in cells or replicates where
$\mathbf H_{\mathrm{lik}}$ is singular or indefinite, the Posterior IDM
remains defined as an SPD-normalized matrix, while the Likelihood IDM
requires the retained-eigenspace machinery (and even then, may fail).
The Posterior framework is therefore globally applicable — it does NOT
require case distinction on the conditioning of $\mathbf H_{\mathrm{lik}}$.

This is the reason the eLife paper adopts the Posterior framework as
the diagnostic of record. Likelihood IDM remains useful as a
supplementary spectral readout in tables (where defined) but is not the
headline diagnostic.

## 7. Connection to the evidence pipeline

The same prior used for the PID diagnostics is the same prior used by
the parallel-tempering evidence-computation pipeline of the macro-dr
codebase. The prior files in `projects/eLife_2025/scheme_*_prior.csv`
are loaded once per run via
[`load_Prior`](../../../../legacy/parameters_distribution.h) or
constructed inline via the `create_prior` DSL function and used
identically by:

- The PID diagnostics (Likelihood and Posterior IDM, GFD, DCC, DIB)
- The Posterior evidence computation
  ([[supplement-posterior-evidence]])
- The information gain computation
  ([[supplement-information-gain]])
- The MAP optimization (algorithm validation supplement, see Section 8
  below)

This is a coherence statement, not a methodological one: the manuscript
uses one prior across all derived quantities. There is no possibility
of the failure mode where different priors are used for evidence and
for parameter estimation.

## 8. Implication for MAP optimization and per-replicate validation

The shift to Posterior-everywhere has a direct consequence for
algorithm validation: the canonical evaluation point is now the
**MAP** ($\theta_{\mathrm{MAP}}$), not the MLE ($\theta_{\mathrm{MLE}}$)
nor the simulator parameter ($\theta_{\mathrm{sim}}$).

- $\theta_{\mathrm{MAP}}$ is the argmax of $\log L(\theta) + \log\pi(\theta)$
- $\theta_{\mathrm{MLE}}$ is the argmax of $\log L(\theta)$ alone
- $\theta_{\mathrm{sim}}$ is the parameter value used by the simulator

Under the PID framework, all diagnostics are evaluated at
$\theta_{\mathrm{MAP}}$. This has three practical consequences:

1. **Optimization is robust globally**. Even in directions where the
   likelihood is singular or non-concave, the posterior has a unique
   maximum because the prior contributes positive curvature
   ($+\sigma^{-2}$) in every direction. There is no "case distinction"
   between identifiable and non-identifiable parameters at the level
   of the optimization itself.

2. **Non-identifiable parameters become diagnostically visible**. When
   the posterior maximum in some direction is dominated by the prior
   (the likelihood is silent), the posterior covariance in that
   direction equals the prior covariance up to small data corrections.
   This is the operational signature of non-identifiability: the prior
   is "not consumed" by the data in that direction. See
   [[supplement-safety-categorization]] for the per-direction reporting.

3. **The MAP is the operational reference**. In production use, the
   user does not know $\theta_{\mathrm{sim}}$. The user has data,
   computes MAP, and reads diagnostics at MAP. The empirical
   verification we run in simulation (where $\theta_{\mathrm{sim}}$ is
   known) compares the operational answer at $\theta_{\mathrm{MAP}}$
   with the ground truth at $\theta_{\mathrm{sim}}$. This is the
   algorithm-bias validation of
   [[supplement-algorithm-validation]].

## 9. Limitations and explicit non-claims

The prior choice in this framework is explicitly NOT:

- **Informative in the elicitation sense**. The $\sigma^2 = 2$ is a
  computational default, not a statement of subjective belief about
  parameter values. A subjective prior, elicited from
  expert biophysicists, would generally be tighter and
  parameter-specific.

- **Parameter-specific**. A more defensible prior would assign
  different variances to different parameters based on their physical
  meaning (rate constants for fast vs slow processes, conductance vs
  noise parameters). The uniform $\sigma^2 = 2$ across all parameters
  is operationally convenient and matches the empirical practice of
  the lab, but is conceptually conservative — it could be tightened.

- **Joint** (off-diagonal in transformed coordinates). The diagonal
  prior ignores correlations between rate constants that a
  literature-elicited prior would capture (e.g., the detailed-balance
  constraints among kinetic rates of a Markov scheme). For inference
  on detailed balance, a non-diagonal prior would be more appropriate.

- **Robust to all parameter scales**. For very large rate constants
  (e.g., $> 10^{5}\,\mathrm{s}^{-1}$), the log10 transformation may
  amplify numerical error in the AD-propagated derivatives. In such
  regimes the $\sigma^2 = 2$ choice may be operationally optimistic.

In the eLife 2025 manuscript scope (2-state and small-state allosteric
schemes at conventional rate ranges), none of these limitations affect
the headline conclusions. Replicates affected by edge-case parameter
scales will show up as MARGINAL or UNRELIABLE in the per-category
safety report of [[supplement-safety-categorization]], which is the
operational safeguard.

## Cross-references

- [[supplement-posterior-main]] — the mathematical framework, including
  dominance condition (Section 2.1), limiting regimes (Section 2.3),
  and the consistency identity (Section 4)
- [[supplement-posterior-evidence]] — how the prior enters the
  evidence correction and the peak-term correction
- [[supplement-information-gain]] — how the prior interacts with the
  information gain budget
- [[supplement-algorithm-validation]] — empirical validation at
  $\theta_{\mathrm{MAP}}$ using the per-replicate covariance test
- [[supplement-safety-categorization]] — per-replicate and per-cell
  reporting policy
- [[implementation-posterior-distortion]] — how the framework is
  realized in `src/core/likelihood.cpp` and
  `legacy/parameters_distribution.h`
