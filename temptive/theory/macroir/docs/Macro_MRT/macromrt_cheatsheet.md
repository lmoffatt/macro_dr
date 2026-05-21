# MacroMRT Developer Cheat Sheet

MacroMRT is the **instantaneous-current** Bayesian update with
**state-dependent (open-channel) noise** handled via a Taylor/Laplace
correction. It is the controlled baseline against which MacroIRT is
compared: same Taylor σ² treatment, but no boundary-state lift.

Dispatch: `averaging=1`, `recursive=true`, `variance=true`.

## State Variables (C++ notation)
- `mu0`        → \( \boldsymbol{\mu}_0 \)         — prior mean state distribution
- `Sigma0`     → \( \boldsymbol{\Sigma}_0 \)      — prior per-channel covariance
- `P_t`        → \( \mathbf{P}(t) \)              — transition matrix \(e^{Qt}\)
- `gbar`       → \( \overline{\boldsymbol{\gamma}}_0 \) — per-\(i_0\) marginalised mean conductance (from Qdtm)
- `gvar`       → \( \overline{\boldsymbol{\sigma}^2}_0 \) — per-\(i_0\) marginalised open-channel noise (from Qdtm)

## Boundary Contractions (from Qdtm)

Per-start-state mean conductance (marginalised over end state):
\[
\overline{\gamma}_0(i_0) \, = \sum_{i_t} P_{i_0\to i_t}(t)\,\overline{\Gamma}_{i_0\to i_t}.
\]

Per-start-state open-channel noise (marginalised over end state):
\[
\overline{\sigma^2}_0(i_0) \, = \sum_{i_t} P_{i_0\to i_t}(t)\,\overline{V}_{i_0\to i_t}.
\]

These are the **same boundary-conditioned moments** that MacroIR uses,
collapsed over \(i_t\) to per-\(i_0\) scalars. No tilde operator. They are
reduced to pure per-state \(\gamma_i,\sigma^2_i\) only when \(Q\) is
constant over \([0,t]\) and the recording window is negligible.

## Predictive Mean
\[
\overline{y}^{\mathrm{pred}} = N_{\mathrm{ch}}\,(\boldsymbol{\mu}_0\cdot\overline{\boldsymbol{\gamma}}_0).
\]

## Predictive Variance

State-dependent variance:
\[
V \, = \epsilon^2_{0\to t} + N_{\mathrm{ch}}\,(\boldsymbol{\mu}_0\cdot\overline{\boldsymbol{\sigma}^2}_0).
\]

(MR drops the boundary cross-covariance term that IR carries via
\(\widetilde{\gamma^{T}\Sigma\gamma}\); MRT inherits this drop and only
adds the Taylor σ² correction.)

## Innovation and Effective Direction
\[
\delta = \overline{y}^{\mathrm{obs}} - \overline{y}^{\mathrm{pred}},
\]
\[
\mathbf{v} = \overline{\boldsymbol{\gamma}}_0 + \frac{\delta}{V}\,\overline{\boldsymbol{\sigma}^2}_0.
\]

Reduces to \(\mathbf{v} = \overline{\boldsymbol{\gamma}}_0\) when
\(\overline{\boldsymbol{\sigma}^2}_0 = 0\) — i.e. MRT collapses to MR.

## Effective Variance Scalar
\[
s = \mathbf{v}^{T}\boldsymbol{\Sigma}^{\mathrm{prop}}\mathbf{v}.
\]

## Propagation

Mean:
\[
\boldsymbol{\mu}^{\mathrm{prior}} = \boldsymbol{\mu}_0\,\mathbf{P}(t).
\]

Covariance:
\[
\boldsymbol{\Sigma}^{\mathrm{prop}} = \mathbf{P}(t)^T(\boldsymbol{\Sigma}_0-\mathrm{diag}\,\boldsymbol{\mu}_0)\mathbf{P}(t) + \mathrm{diag}(\boldsymbol{\mu}^{\mathrm{prior}}).
\]

## Measurement Update (rank-1 quasi-Laplace; main implementation form)

**Convention for \(V\):** all formulas below use \(V = V_{\mathrm{obs}} = \varepsilon^2 + N_{\mathrm{ch}}\,\boldsymbol{\mu}\!\cdot\!\overline{\boldsymbol{\sigma}^2}_0\) (measurement-noise piece at \(r=\boldsymbol{\mu}\), what the energy/Hessian derivation produces), **not** \(V_{\mathrm{pred}} = V_{\mathrm{obs}} + N_{\mathrm{ch}}\,s_{\boldsymbol{\gamma}}\). In code: `V_obs = r_y_var - N*gSg`.

### Start-frame form (canonical, before propagation through P)

Covariance:
\[
\boldsymbol{\Sigma}^{\mathrm{post}}_{\mathrm{start}} = \boldsymbol{\Sigma}_p - \frac{N_{\mathrm{ch}}}{V_{\mathrm{obs}} + N_{\mathrm{ch}}\,s}\,\boldsymbol{\Sigma}_p\mathbf{v}\mathbf{v}^{T}\boldsymbol{\Sigma}_p, \qquad s := \mathbf{v}^{T}\boldsymbol{\Sigma}_p\mathbf{v}.
\]

Mean (single Newton step from \(p_0 = \boldsymbol{\mu}\)):
\[
\boldsymbol{\mu}^{\mathrm{post}}_{\mathrm{start}} = \boldsymbol{\mu} + \frac{\delta}{2V_{\mathrm{obs}}}\,\boldsymbol{\Sigma}^{\mathrm{post}}_{\mathrm{start}}\,(\overline{\boldsymbol{\gamma}}_0 + \mathbf{v}).
\]

### Endpoint-frame form (what the code computes; av=1)

Code stores everything in endpoint frame: `mu_prior_end = μ·P`, `sigma_pre = PᵀSmDP + diag(μ_end)`. Propagating the start-frame Newton step through `P` and expanding the SM down-date gives the implementable row-form expressions. Define
- \(\mathtt{gS} := \overline{\boldsymbol{\gamma}}_0^{T}\,\boldsymbol{\Sigma}_p\,\mathbf{P}\),
- \(\mathtt{vS} := \mathbf{v}^{T}\,\boldsymbol{\Sigma}_p\,\mathbf{P}\),
- \(\mathtt{vSv} := \mathbf{v}^{T}\,\boldsymbol{\Sigma}_p\,\mathbf{v} = s\),
- \(\mathtt{b'} := \overline{\boldsymbol{\gamma}}_0^{T}\,\boldsymbol{\Sigma}_p\,\mathbf{v} = s_{\boldsymbol{\gamma}} + (\delta/V_{\mathrm{obs}})\,b_{\mathrm{tilde}}\), where \(b_{\mathrm{tilde}} := \overline{\boldsymbol{\gamma}}_0^{T}\,\boldsymbol{\Sigma}_p\,\overline{\boldsymbol{\sigma}^2}_0\),
- \(\mathtt{sm} := N_{\mathrm{ch}}/(V_{\mathrm{obs}} + N_{\mathrm{ch}}\,\mathtt{vSv})\).

Covariance:
\[
\boldsymbol{\Sigma}^{\mathrm{post}}_{\mathrm{end}} = \boldsymbol{\Sigma}^{\mathrm{prop}}_{\mathrm{end}} - \mathtt{sm}\,\mathtt{vS}^{T}\mathtt{vS}.
\]

Mean (row form, the SM expansion is explicit; **not** \((\mathtt{gS}+\mathtt{vS})\cdot\boldsymbol{\Sigma}^{\mathrm{post}}\)):
\[
\boldsymbol{\mu}^{\mathrm{post}}_{\mathrm{end}} = \boldsymbol{\mu}^{\mathrm{prior}}_{\mathrm{end}} + \frac{\delta}{2V_{\mathrm{obs}}}\,\Bigl[(\mathtt{gS} + \mathtt{vS}) - \mathtt{sm}\,(\mathtt{b'} + \mathtt{vSv})\,\mathtt{vS}\Bigr].
\]

At \(\mathbf{v}=\overline{\boldsymbol{\gamma}}_0\) the bracket reduces to \(2\,\mathtt{gS}\,(V_{\mathrm{obs}}/V_{\mathrm{pred}})\), recovering \((\delta/V_{\mathrm{pred}})\,\mathtt{gS}\) — the canonical Moffatt-2007 endpoint Kalman update. This collapse is structural, not coincidental.

(Rank-2 form including the \(\tfrac12\log V\) likelihood term is given in the supplement; for ensembles with large \(N_{\mathrm{ch}}\) the rank-1 form is the working scheme.)

## Reduction to MR (variance_taylor_correction = false)

Setting \(\mathbf{v}=\overline{\boldsymbol{\gamma}}_0\) and dropping the
σ² Taylor correction:
\[
\boldsymbol{\Sigma}^{\mathrm{post}}_{\mathrm{MR}} = \boldsymbol{\Sigma}^{\mathrm{prop}} - \frac{N_{\mathrm{ch}}}{V + N_{\mathrm{ch}}\,\overline{\gamma}_0^T\boldsymbol{\Sigma}^{\mathrm{prop}}\overline{\gamma}_0}\,\boldsymbol{\Sigma}^{\mathrm{prop}}\overline{\gamma}_0\overline{\gamma}_0^T\boldsymbol{\Sigma}^{\mathrm{prop}},
\]
\[
\boldsymbol{\mu}^{\mathrm{post}}_{\mathrm{MR}} = \boldsymbol{\mu}^{\mathrm{prior}} + \frac{\delta}{V}\,\boldsymbol{\Sigma}^{\mathrm{post}}_{\mathrm{MR}}\,\overline{\gamma}_0.
\]
This is the Moffatt 2007 update.

## Log-Likelihood

\[
\log L = -\tfrac{1}{2}\bigl[\log(2\pi V) + \delta^2/V\bigr].
\]

The trust-region simplex shrink \(\alpha^\star\) (cf. MacroIR variance
inflation correction) applies identically to MRT and only modifies the
state update, not \(\log L\).

## Mapping to MacroIRT

| Quantity                   | MacroMRT                                                  | MacroIRT                                                     |
|----------------------------|-----------------------------------------------------------|--------------------------------------------------------------|
| Cross-covariance vector    | \(\boldsymbol{\Sigma}^{\mathrm{prop}}\overline{\gamma}_0\) | \(\widetilde{\gamma^{T}\Sigma}\)                              |
| Effective direction        | \(\mathbf{v}\) (per-\(i_0\))                              | \(\widetilde{v^{T}\Sigma}\) (full boundary lift)              |
| Effective scalar           | \(s = \mathbf{v}^{T}\Sigma^{\mathrm{prop}}\mathbf{v}\)    | \(s = \widetilde{v^{T}\Sigma v}\)                             |
| Variance third term        | dropped                                                   | \(N_{\mathrm{ch}}\widetilde{\gamma^{T}\Sigma\gamma}\)         |

The IR "win" — if real — must come from the third variance term and the
boundary-state lift of the cross-covariance, holding the Taylor σ²
correction constant.
