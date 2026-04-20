# ADR-005 — Likelihood Diagnostic Presets

Status: Accepted. Implemented in
[include/macrodr/cmd/likelihood.h](../../include/macrodr/cmd/likelihood.h),
[src/core/likelihood.cpp](../../src/core/likelihood.cpp), and
[legacy/moment_statistics.h](../../legacy/moment_statistics.h).

## Context

`Analisis_derivative_diagnostic` used to be a class template keyed by two bools:

```cpp
template <bool include_evolution, bool include_series_cross_correlation>
struct Analisis_derivative_diagnostic { … };
```

Three DSL commands wrapped the three meaningful instantiations
(`…<false,false>`, `…<true,false>`, `…<false,true>`). The compile surface grew
without a natural way to add intermediate trade-offs between **compute/storage
cost** and **diagnostic depth**, and the two-bool flags conflated orthogonal
choices (time-series moments vs. full cross-correlation kernel) while forcing
matrix-valued quantities like `Gaussian_Fisher_Information` into a single tier.

## Decision

Replace the template with **five named presets**, each a `var::Vector_Space`
alias, plus a shared base that every preset extends:

| Preset | Extra content beyond the base |
| --- | --- |
| `basic` | — (base scalars + distortion matrices + Σ moments only) |
| `series_var` | 6 observables as `Report_local_var<V>` + `Per_sample_derived_diagnostics` |
| `series_cov` | 7 observables as `Report_local_cov<V>`, GFI as `Report_local_var` + `Per_sample_derived_diagnostics` |
| `series_kernel` | 5 observables as `Report_cross<V>` (no auxiliaries, no per-sample) |
| `series_kernel_full` | 7 observables as `Report_cross<V>`, GFI as `Report_local_var` + `Per_sample_derived_diagnostics` |

The shared `Analisis_derivative_diagnostic_base` always carries:
`Information_Distortion_Matrix`, `Sample_Distortion_Matrix`,
`Correlation_Distortion_Matrix`, `Fisher_Covariance`,
`Distortion_Corrected_Covariance`, their `log_Det<…>` (for Laplace-approximated
Bayes-factor corrections), `Distortion_Induced_Bias`, Σ moments for the seven
observables, plus an **identifiability diagnostic bundle** on `Fisher_Covariance`
and `Distortion_Corrected_Covariance`:

- `Correlation_Of<V>` — p×p Pearson correlation matrix, computed via the
  spectral form so `|ρ| ≤ 1` is guaranteed by construction.
- `Eigenvalue_Spectrum<V>` — H's eigenvalues sorted descending (κ and gap
  readable from the plot).
- `Effective_Rank<V>` — count of eigenvalues above `rtol·λ_max`; null defect
  = p − effective rank.
- `Null_Space_Projector<V>` — rotation-invariant p×p projector onto the
  non-identifiable subspace; `Π_{ii}` gives per-parameter non-identifiability
  fraction.

See [docs/math/non-identifiability-diagnosis.md](../math/non-identifiability-diagnosis.md)
for interpretation and workflow.

### Report families (defined in `legacy/moment_statistics.h`)

Three dedicated report types for the 6–7 leaf observables, picked per preset:

- `Report_integral<V>` — `count + integral_correlation_lag` only (Path 2 summary).
- `Report_local_var<V>` — `count + series_mean + series_diagonal_variance +
  cross_correlation_lag_forward + integral_correlation_lag`. Compact per-sample
  storage (p values at each time step for vector V, p(p+1)/2 for matrix V).
- `Report_local_cov<V>` — full lag-0 block per sample (p×p for vector V).
- `Report_cross<V>` — full kernel over lags for V (most expensive).

### `Per_sample_derived_diagnostics`

Unlike the leaf-observable reports, `Sample_Distortion_Matrix` and
`Distortion_Induced_Bias` are **Evolution-wrapped per-sample derived quantities**
(one per time step, carrying themselves — not moment statistics over samples).
They stay outside the `Series_Moment_report` family and are concatenated
separately:

```cpp
using Per_sample_derived_diagnostics =
    Probit_statistics<Evolution_of<Vector_Space<Sample_Distortion_Matrix,
                                                Distortion_Induced_Bias>>>;
```

The concatenation pattern in C++ requires nesting because `var::concatenate_t`
is binary:

```cpp
using Analisis_derivative_diagnostic_series_cov =
    var::concatenate_t<
        var::concatenate_t<Analisis_derivative_diagnostic_base, Per_V_Reports>,
        var::Vector_Space<Per_sample_derived_diagnostics>>;
```

Runtime assembly mirrors this: `concatenate(base, series)` returns a
`Vector_Space`, then `push_back_var(combined, per_sample)` appends the
`Evolution_of<…>` — `push_back_var` requires a `Vector_Space` as its first arg,
not the `Evolution_of`.

### DSL surface

Five commands registered in
[src/cli/command_manager.cpp](../../src/cli/command_manager.cpp):

- `likelihood_derivative_basic_diagnostics`
- `likelihood_derivative_series_var_diagnostics`
- `likelihood_derivative_series_cov_diagnostics`
- `likelihood_derivative_series_kernel_diagnostics`
- `likelihood_derivative_series_kernel_full_diagnostics`

All take the same signature: `(dlikelihood_predictions, n_boostrap_samples,
probits, seed, max_lag)`. Every preset now depends on `max_lag` because
`integral_correlation_lag` (Path 2) is part of the base reports.

The generic `write_csv(var::Indexed<T> const&, std::string)` template at
[include/macrodr/cmd/likelihood.h:314](../../include/macrodr/cmd/likelihood.h#L314)
covers indexed variants automatically — no per-preset `Indexed<…>` write_csv
overloads needed.

### Shared PSD decomposition

Every distortion computation (`IDM`, `SDM`, `DCC`, `FC`, `DIB`) and every
kernel normalization pass share one eigendecomposition of the anchor matrix
(`H`, `J_sample`, `SDM`, or lag-0 covariance blocks) via `lapack::PSDDecomposition`
in [legacy/lapack_headers.h](../../legacy/lapack_headers.h). A single
`compute_psd_decomp` per anchor produces basis + retained/full spectrum +
diagonal inverse/inv-sqrt/sqrt scalings; downstream helpers (`apply_normalized_congruence`,
`apply_inverse_congruence`, `apply_inverse_vector`, `apply_inverse_as_matrix`,
`apply_psd_whitening`, `fc_correlation_from_decomp`, `dcc_correlation_from_decomp`,
`null_space_projector`) consume it without re-decomposing. This replaces what
was previously 5 eigendecompositions of H per bootstrap sample (and T × max_lag
eigendecompositions in the kernel loop) with one per anchor, and gives all
identifiability outputs for free from the same data.

## Why not alternatives

- **Keep one template, add more bools.** Each new axis (kernel vs. moments,
  cov vs. var, include per-sample) doubles the instantiations. The 2-bool form
  was already hard to pick between at the call site; three bools would be
  worse.
- **Union of everything in one type.** Defeats the point: the five tiers exist
  *because* full cross-correlation kernel × all observables × all samples is
  too expensive for routine runs. The preset names let the caller pick their
  budget.
- **Let DSL pick dynamically.** Type-erasure would force heap allocation and
  lose the CSV-writing generic machinery that currently relies on the member
  list being known statically.

## GFI as a matrix-valued special case

`Gaussian_Fisher_Information` is a p×p `SymPosDefMatrix`. Including it as
`Report_local_cov` would cost p²×p² per sample — too much. It only appears as
`Report_local_var` in `series_cov` and `series_kernel_full`, where its
per-sample "diagonal variance" is scattered from the m×m packed Sokal–Kish
block (m = p(p+1)/2) back into a p×p symmetric matrix holding Var(GFI[i,j]).
See [docs/math/correlation-lag-paths.md](../math/correlation-lag-paths.md) for
the derivation.

## Consequences

- Removing the 2-bool template broke three call sites
  ([command_manager.cpp](../../src/cli/command_manager.cpp) Indexed aliases
  and three write_csv registrations); these were replaced with the new preset
  write_csv overloads.
- Scripts that referenced the old DSL commands
  (`likelihood_derivative_diagnostics`, `…_evolution`, `…_cross_correlation`)
  must be updated to pick a preset. No backwards-compat aliases.
- CSV schemas gained two orthogonal columns: `variable` (the leaf observable V)
  and `operation` (the report-type discriminator, e.g.
  `integral_correlation_lag`, `cross_correlation_lag_forward`).
- Adding a new leaf observable now means adding it to the 1–2 presets that
  should carry it, rather than touching every `<bool,bool>` instantiation.
