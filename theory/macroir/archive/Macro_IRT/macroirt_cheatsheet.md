# MacroIRT Developer Cheat Sheet

MacroIRT is **MacroIR + state-dependent (open-channel) noise** handled
via a Taylor/Laplace correction in boundary-state space. Compared with
MacroMRT, the σ² treatment is identical; the difference is the
**boundary-state lift** (full tilde operator) replacing per-\(i_0\)
scalars.

Dispatch: `averaging=2`, `recursive=true`, `variance=true`.

## State Variables (C++ notation)
- `mu0`        → \( \boldsymbol{\mu}_0 \)
- `Sigma0`     → \( \boldsymbol{\Sigma}_0 \)
- `P_t`        → \( \mathbf{P}(t) \)
- `Gamma_bar`  → \( \overline{\Gamma}_{i_0\to i_t} \) — boundary mean conductance
- `V_bar`      → \( \overline{V}_{i_0\to i_t} \)      — boundary open-channel noise

## Boundary Contractions
\[
\overline{\gamma}_0(i_0) \, = \sum_{i_t} P_{i_0\to i_t}(t)\,\overline{\Gamma}_{i_0\to i_t},
\]
\[
\overline{\sigma^2}_0(i_0) \, = \sum_{i_t} P_{i_0\to i_t}(t)\,\overline{V}_{i_0\to i_t}.
\]

## Predictive Mean
\[
\overline{y}^{\mathrm{pred}} = N_{\mathrm{ch}}\,(\boldsymbol{\mu}_0\cdot\overline{\boldsymbol{\gamma}}_0).
\]

## Unified Tilde Operator (lift–modulate–collapse)

Vector form (boundary cross-covariance):
\[
\widetilde{u^T\Sigma} \, = \overline{u}_0^T(\boldsymbol{\Sigma}_0 - \mathrm{diag}\,\boldsymbol{\mu}_0) + (\overline{U}\circ \mathbf{P})^T\boldsymbol{\mu}_0.
\]

Scalar form:
\[
\widetilde{u^T\Sigma w} \, = \overline{u}_0^T(\boldsymbol{\Sigma}_0 - \mathrm{diag}\,\boldsymbol{\mu}_0)\overline{w}_0 + \boldsymbol{\mu}_0^T[(\overline{U}\circ \mathbf{P})\circ(\overline{W}\circ \mathbf{P})]\mathbf{1}.
\]

Specialised:
- `tilde_gamma_Sigma_vec`         ↔ \(\widetilde{\gamma^T\Sigma}\)
- `tilde_gamma_Sigma_gamma_scalar`↔ \(\widetilde{\gamma^T\Sigma\gamma}\)
- `tilde_v_Sigma_vec`             ↔ \(\widetilde{v^T\Sigma}\)
- `tilde_v_Sigma_v_scalar`        ↔ \(\widetilde{v^T\Sigma v}\)

## Predictive Variance
\[
V \, = \epsilon^2_{0\to t} + N_{\mathrm{ch}}\,\widetilde{\gamma^T\Sigma\gamma} + N_{\mathrm{ch}}\sum_{i_0}\mu_0(i_0)\,\overline{\sigma^2}_0(i_0).
\]

The **third term is what MRT drops**: the boundary cross-correlation
between channel state and current. Holding the Taylor σ² treatment
constant, this term plus the boundary lift of the cross-covariance
vector is the entire IR-vs-MR delta.

## Innovation and Effective Direction
\[
\delta = \overline{y}^{\mathrm{obs}} - \overline{y}^{\mathrm{pred}},
\]
\[
\mathbf{v} = \overline{\boldsymbol{\gamma}}_0 + \frac{\delta}{V}\,\overline{\boldsymbol{\sigma}^2}_0.
\]

## Effective Variance Scalar
\[
s = \widetilde{v^T\Sigma v}.
\]

## Propagation

Mean:
\[
\boldsymbol{\mu}^{\mathrm{prior}} = \boldsymbol{\mu}_0\,\mathbf{P}(t).
\]

Covariance:
\[
\boldsymbol{\Sigma}^{\mathrm{prop}} = \mathbf{P}(t)^T(\boldsymbol{\Sigma}_0 - \mathrm{diag}\,\boldsymbol{\mu}_0)\mathbf{P}(t) + \mathrm{diag}(\boldsymbol{\mu}^{\mathrm{prior}}).
\]

## Measurement Update (rank-1 quasi-Laplace)

Covariance:
\[
\boldsymbol{\Sigma}^{\mathrm{post}} = \boldsymbol{\Sigma}^{\mathrm{prop}} - \frac{N_{\mathrm{ch}}}{V + N_{\mathrm{ch}}\,s}\,(\widetilde{v^T\Sigma})^T(\widetilde{v^T\Sigma}).
\]

Mean:
\[
\boldsymbol{\mu}^{\mathrm{post}} = \boldsymbol{\mu}^{\mathrm{prior}} + \frac{\delta}{2V}\,\boldsymbol{\Sigma}^{\mathrm{post}}\,(\widetilde{\gamma^T\Sigma} + \widetilde{v^T\Sigma}).
\]

## Reduction to MacroIR (variance_taylor_correction = false)

Setting \(\mathbf{v} = \overline{\boldsymbol{\gamma}}_0\),
\(\widetilde{v^T\Sigma} = \widetilde{\gamma^T\Sigma}\),
\(s = \widetilde{\gamma^T\Sigma\gamma}\):
\[
\boldsymbol{\Sigma}^{\mathrm{post}}_{\mathrm{IR}} = \boldsymbol{\Sigma}^{\mathrm{prop}} - \frac{1}{V}\,(\widetilde{\gamma^T\Sigma})^T(\widetilde{\gamma^T\Sigma}),
\]
\[
\boldsymbol{\mu}^{\mathrm{post}}_{\mathrm{IR}} = \boldsymbol{\mu}^{\mathrm{prior}} + \frac{1}{V}\,\boldsymbol{\Sigma}^{\mathrm{prop}}\,\widetilde{\gamma^T\Sigma}\,\delta.
\]
This is the standard MacroIR scalar Kalman update.

## Reduction to MacroMRT (averaging_approximation = 1)

Replace each tilde with its per-\(i_0\) collapse:
- \(\widetilde{\gamma^T\Sigma}\to \boldsymbol{\Sigma}^{\mathrm{prop}}\overline{\gamma}_0\),
- \(\widetilde{v^T\Sigma}\to \boldsymbol{\Sigma}^{\mathrm{prop}}\mathbf{v}\),
- \(\widetilde{v^T\Sigma v}\to \mathbf{v}^T\boldsymbol{\Sigma}^{\mathrm{prop}}\mathbf{v}\),
- drop the third variance term \(N_{\mathrm{ch}}\widetilde{\gamma^T\Sigma\gamma}\).

This recovers the MacroMRT cheatsheet formulas exactly.

## Log-Likelihood

\[
\log L = -\tfrac{1}{2}\bigl[\log(2\pi V) + \delta^2/V\bigr].
\]

The trust-region simplex shrink \(\alpha^\star\) (cf. macroir variance
inflation correction) applies identically and only modifies the state
update, not \(\log L\).

## Mapping to MacroMRT (the comparison axis)

| Quantity                   | MacroMRT                                                  | MacroIRT                                                     |
|----------------------------|-----------------------------------------------------------|--------------------------------------------------------------|
| Cross-covariance vector    | \(\boldsymbol{\Sigma}^{\mathrm{prop}}\overline{\gamma}_0\) | \(\widetilde{\gamma^{T}\Sigma}\)                              |
| Effective direction        | \(\mathbf{v}\) (per-\(i_0\))                              | \(\widetilde{v^{T}\Sigma}\) (full boundary lift)              |
| Effective scalar           | \(\mathbf{v}^{T}\Sigma^{\mathrm{prop}}\mathbf{v}\)        | \(\widetilde{v^{T}\Sigma v}\)                                 |
| Variance third term        | dropped                                                   | \(N_{\mathrm{ch}}\widetilde{\gamma^{T}\Sigma\gamma}\)         |

Holding the Taylor σ² treatment fixed by comparing MRT vs IRT
isolates the boundary-state contribution of MacroIR.
