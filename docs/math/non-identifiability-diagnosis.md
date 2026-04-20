# Diagnosing Parameter Non-Identifiability from FC and DCC

Status: research note, companion to the Fisher-diagnostics pipeline described
in [correlation-lag-paths.md](correlation-lag-paths.md) and
[score-autocorrelation-and-distortion-correction.md](score-autocorrelation-and-distortion-correction.md).
Introduces the identifiability diagnostic bundle added to
`Analisis_derivative_diagnostic_base`:

- `Correlation_Of<Fisher_Covariance>`, `Correlation_Of<Distortion_Corrected_Covariance>`
- `Eigenvalue_Spectrum<·>`, `Effective_Rank<·>`, `Null_Space_Projector<·>`

Correlations are computed in **spectral form** — directly from H's
eigendecomposition — so `|ρ| ≤ 1` is guaranteed by Cauchy–Schwarz on
unit vectors; the numerically fragile `basis · reduced · basisᵀ` round-trip
that could leak `|ρ| > 1` out of the old implementation is bypassed.

## 1. Why correlation matrices?

The Fisher information matrix $\mathcal{I}(\theta)$ and its inverse — the
Fisher covariance $C_F = \mathcal{I}^{-1}$ — tell you the asymptotic covariance
of the MLE. The **Pearson correlation matrix** of $C_F$,

$$
R_{ij} \;=\; \frac{(C_F)_{ij}}{\sqrt{(C_F)_{ii}\,(C_F)_{jj}}},
$$

normalizes away scale and makes the off-diagonals directly interpretable:
$R_{ij}\approx 1$ means a linear combination $\theta_i + c\,\theta_j$ is
almost invariant in the likelihood — the two parameters trade off in a
near-singular direction. This is a **structural identifiability** diagnosis:
it tells you the *shape* of the ridge, not just that one exists.

The same construction on the Distortion-Corrected Covariance $C_{DCC}$
measures correlation **after the model's own distortion has been unwound** —
so near-unity entries in $R_{DCC}$ are identifiability problems that persist
even when the approximate likelihood's score-correlation bias has been
accounted for.

## 2. What each diagnostic reports

| Field | Type | Meaning |
| --- | --- | --- |
| `Correlation_Of<C>` | p×p SPD (ParameterIndexed) | Full $R$, parameter-labeled rows/cols; `|ρ|≤1` by construction |
| `Eigenvalue_Spectrum<C>` | p-vector (sorted ↓) | H's eigenvalues; κ = λ_max/λ_min + gap location |
| `Effective_Rank<C>` | `size_t` | #{λ_k > rtol·λ_max}; null defect = p − effective_rank |
| `Spectrum_Condition_Number<C>` | `double` | κ = λ_max/λ_min scalar; lattice-scale sloppiness screen |
| `Null_Space_Projector<C>` | p×p SPD (ParameterIndexed) | `Π = Σ_{k : λ_k ≤ tol} v_k v_kᵀ`; formal nulls only; empty for sloppy-not-null |
| `Worst_Subspace_Projector<C>` | p×p SPD (ParameterIndexed) | Projector onto the smallest-eigenvalue cluster (gap-based); populated whenever κ ≥ 10⁴ |

**Null vs. Worst projector** — they answer different questions:

- `Null_Space_Projector` uses a strict threshold (`rtol · λ_max`, default `1e-10`) so
  it's non-empty **only if some direction is numerically null**. Use this
  for model reduction ("can I drop a parameter?") and for formal rank analysis.
- `Worst_Subspace_Projector` is always populated whenever `κ ≥ 10⁴` (the
  "sloppy" threshold from Section 4). It projects onto the smallest
  eigenvalue's eigenvector, extending upward through near-degenerate
  neighbors (within 10× of `λ_min`). Use this for everyday sloppiness
  diagnosis ("where is my worst direction?"). In your figure_3 stationary
  case with κ=10⁹ and no formal null, `Null_Space_Projector` is empty
  while `Worst_Subspace_Projector` correctly localizes the sloppy
  kinetic-rate combination.

Both have the same ParameterIndexed payload shape, so the CSV column
layout is identical; only the `operation` tag differs
(`null_space_projector` vs `worst_subspace_projector`).

Why the projector instead of a "smallest eigenvector":
individual eigenvectors are ambiguous under rotation within a degenerate
subspace, so reporting them is misleading when null defect ≥ 2. The
projector is basis-invariant within its span and carries the same
information — plus it composes cleanly when you have more than one null
direction. `Π_{ii} ∈ [0,1]` is the fraction of parameter *i*'s variance
living in the non-identifiable subspace: `Π_{ii}=0` means fully
identifiable individually; `Π_{ii}=1` means fully unidentifiable on its
own. Off-diagonals reveal entanglement.

**Note on the null projector for FC vs DCC.** Because `DCC = H⁻¹ J H⁻¹`,
any `v` with `Hv = 0` satisfies `vᵀ DCC v = 0` — so DCC inherits H's null
subspace. We emit the same projector (derived from H's eigendecomposition)
under both `variable=Fisher_Covariance` and `variable=Distortion_Corrected_Covariance`
for schema consistency; the two rows will be identical for a given cell.
Likewise `Eigenvalue_Spectrum` is H's spectrum (sorted descending); FC's
spectrum is `1/λ_H` (reciprocation preserves the null mask and gap location).

## 3. Why the three-quantity bundle (rather than correlation alone)

Correlation carries only pairwise-ridge information and **misses genuine
non-identifiability in two important cases**:

1. **Individual-coordinate nulls.** If a parameter θ_i has infinite marginal
   variance on its own — i.e. the null direction is the coordinate axis e_i
   — all `ρ(i, j)` can be 0 even though θ_i is unidentifiable. The signal
   is in `FC_{ii} → ∞`, equivalently `Π_{ii} ≈ 1`. Correlation alone says
   "no coupling"; the null projector correctly flags θ_i.

2. **Higher-order collinearities.** A 3-way combo like `θ_1 + θ_2 − 2θ_3 ≈
   const` distributes loading across 3 parameters. Each pairwise `|ρ|` can
   be < 1, yet the direction is a rigorous null. The `Eigenvalue_Spectrum`
   picks up the near-zero eigenvalue and the null projector's off-diagonals
   reveal the 3-way coupling.

**Theorem (pairwise).** `|ρ(i, j)| = 1` iff the 2×2 block of the null
projector at rows/cols `{i, j}` has rank 1. Rank-2 block (= "the whole
(i,j) plane is null, but on-off are not coupled along a ridge") gives
`ρ(i, j) = 0`. So correlation saturates at ±1 **only when the null
structure happens to be a simple parameter-pair difference**; the
projector catches the non-saturating cases.

## 4. Thresholds that are actionable in practice

| Diagnostic | Threshold | Interpretation |
| --- | --- | --- |
| `Effective_Rank / p` | = 1 | All directions identified; no rank deficit. |
| | < 1 | Null defect = p − Effective_Rank; there are unidentifiable directions. |
| `Eigenvalue_Spectrum` | κ = λ_max/λ_min < 10⁶ | Well-conditioned; MCMC converges without reparameterization. |
| | 10⁶ – 10¹⁰ | Sloppy; reparameterize along the smallest eigenvectors, or scale proposals by √λ. |
| | > 10¹⁰ | Effectively rank-deficient; an eigenvalue is at numerical-noise floor. |
| `Π_{ii}` | < 0.01 | Parameter *i* is essentially fully identified individually. |
| | 0.01 – 0.9 | Partial entanglement; look at `Π_{ij}` and `ρ(i, j)` to identify partners. |
| | > 0.9 | Parameter *i* is essentially unidentifiable; fix, re-parameterize, or profile. |
| `\|ρ(i, j)\|` | > 0.99 | Strong pairwise ridge; candidate for reparameterization to (θ_i + θ_j) / (θ_i − θ_j). |

The thresholds are protocol-dependent. `figure_3.macroir` was set up to
show this: stationary conditions collapse the kinetic contribution, so
non-recursive algorithms should exhibit `Effective_Rank < p` for rate-pair
parameters whose effects are only distinguishable through kinetic
relaxation. Recursive variants that integrate the kinetic history recover
the missing rank.

## 5. Workflow for diagnosing a suspected non-identifiability

Concretely, when you suspect a ridge (e.g. `macro_MNRV` on `figure_3`'s
stationary protocol):

1. **Run the `basic` preset** over the (algorithm, interval_in_tau,
   noise, Num_ch) lattice. The bundle you need lives in `base`.
2. **Rank-check at a glance.** Plot `Effective_Rank<Fisher_Covariance>` as
   a heatmap across the lattice. Cells with `Effective_Rank < p` are
   flagged for follow-up. This is your lattice-scale yes/no screen.
3. **Compare algorithms.** For a given protocol cell, `Effective_Rank` near
   `p` on `macro_MRV` / `macro_IRV` / `macro_IRTV` but `< p` on `macro_MNRV`
   means the non-identifiability is an **artifact of the non-recursive
   algorithm** (loses transient information). If all algorithms show
   `Effective_Rank < p`, the protocol itself is the bottleneck.
4. **Spectrum eyeball.** For flagged cells, plot
   `Eigenvalue_Spectrum<Fisher_Covariance>` on a log scale. A large gap
   confirms rank deficiency; a slow decay says "sloppy, not degenerate."
5. **Worst-direction localization.** Read off `Worst_Subspace_Projector`
   diagonal `Π_{ii}`: parameters with `Π_{ii}` near 1 dominate the worst
   direction; off-diagonals `Π_{ij}` tell you *which* parameters are
   co-involved. This is populated whenever κ ≥ 10⁴, so it's always
   available for flagged cells. For cells where `Effective_Rank < p`
   (formal null defect), cross-reference with `Null_Space_Projector` —
   it will give the same answer when the worst direction *is* the null
   direction, and will differ only when the null defect involves multiple
   directions of different magnitudes.
6. **Pairwise check via correlation.** `Correlation_Of<Fisher_Covariance>`
   closes the loop. Entries near ±1 that align with `Π_{ii}` loadings
   identify the specific ridge combination. (Correlation at 0 despite
   `Π_{ii}` loading signals a 3+-way higher-order collinearity — see
   Section 3.)
7. **Compare FC vs DCC.** Both share H's null subspace (so `Π` and
   `Effective_Rank` rows match) — but the correlations can differ if the
   distortion alters the informative-subspace posterior shape. If
   DCC correlations drop significantly versus FC's, the apparent ridge is
   partly distortion-driven and removable by correction.
8. **Cross-reference with IDM and CDM.** `log_Det<Information_Distortion_Matrix>`
   far from 0 signals information loss; `Correlation_Distortion_Matrix`
   eigenvalues far from 1 flag per-sample vs total covariance disagreement.
   A non-identifiable direction should show either a near-zero IDM
   eigenvalue or a CDM eigenvalue far from 1.

## 6. What to do with a confirmed ridge

- **Reparameterize.** If $\theta_i + c\,\theta_j$ is well-determined but
  $\theta_i - c\,\theta_j$ is not, fit to the sum and impose a prior on the
  difference.
- **Extend the protocol.** Add a segment that excites the direction the
  ridge hides in — for stationary protocols, add an agonist step; for
  step protocols, try a different amplitude.
- **Add a restraint.** If two rate constants appear only as a ratio, fix
  the ratio from auxiliary data.
- **Accept the ridge and report the joint.** Marginal CIs on the two
  ridge-coupled parameters are misleading; emit the full 2D posterior
  projection or the reparameterized summary.

## 7. Known limitations of the diagnostic

- **Local, not global.** The bundle is built from H at the current $\theta$
  — the Hessian of the log-likelihood at the MLE or a reference point.
  Bimodal posteriors, non-identifiable subsets separated by saddle points,
  or near-flat directions that only emerge far from $\theta$ are invisible
  to it. Use a profile-likelihood scan for those.
- **Threshold-sensitive.** `Effective_Rank` depends on `rtol`; we use
  `rtol = 1e-10` by default (very strict). Tuning it looser flags more
  directions as null. The `Eigenvalue_Spectrum` is threshold-free and is
  the right scalar to eyeball when the cutoff matters.
- **Correlation alone is insufficient** (theorem in Section 3). The
  bundle's other members — `Eigenvalue_Spectrum`, `Effective_Rank`,
  `Null_Space_Projector` — exist precisely to catch cases correlation misses.
- **Pairwise correlation.** $|\rho_{ij}|$ saturates at 1 only for
  parameter-pair-difference null directions. Use `Π_{ii}` (projector
  diagonal) and `Eigenvalue_Spectrum` to catch higher-order collinearities.
  `log_Det<Fisher_Covariance>`, which is $-\log\det\mathcal{I}$ up to sign,
  is the complementary scalar that
  catches this: it diverges to $-\infty$ when any eigenvalue of
  $\mathcal{I}$ vanishes, regardless of whether the null direction aligns
  with parameter pairs.
- **Bootstrap distribution.** Every diagnostic in the bundle is wrapped in
  `Probit_statistics` via the preset bootstrap; the CSV emits mean + CIs
  per cell, not just point estimates. Use the upper CI for classifying
  `Π_{ii}` and `|ρ|`, and the lower CI for `Effective_Rank` — i.e. the
  conservative bound for the "is this direction identified?" question.
- **`Effective_Rank` is integer-valued per bootstrap sample.** Its
  bootstrap mean is continuous (fraction of samples where each mode is
  retained); interpret it as the expected rank, not a count.

## 8. Appendix: CSV column layout

The bundle lands in the following CSV columns (long format):

```
variable,                         operation,                   row, col, statistic, value
Fisher_Covariance,                probit.spectrum,             0,      , mean,      <λ_1>
Fisher_Covariance,                probit.spectrum,             ...,    , ...,       ...
Fisher_Covariance,                probit.effective_rank,          ,    , mean,      <r>
Fisher_Covariance,                probit.correlation_of,       on, off, mean,       <ρ>
Fisher_Covariance,                probit.null_space_projector, on, off, mean,       <Π_{on,off}>
Distortion_Corrected_Covariance,  probit.correlation_of,       ...,    , ...,       ...
```

`variable` tags the source covariance (FC or DCC). `operation` discriminates
the four member types (`correlation_of`, `spectrum`, `effective_rank`,
`null_space_projector`). Rows/cols are parameter names for ParameterIndexed
payloads (correlation, null projector) and integer mode indices for the
spectrum; `effective_rank` is a scalar with no row/col. `statistic` holds
the bootstrap moment (`mean`) or probit quantiles (`probit_0.025`,
`probit_0.975`).
