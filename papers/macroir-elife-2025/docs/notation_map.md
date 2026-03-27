# MacroIR Code-to-Manuscript Notation Map

This note maps the diagnostic/intermediate state variables in `legacy/qmodel.h` to the notation used in the MacroIR manuscript drafts.

Scope:

- interval-recursive (`IR`) state moments
- boundary-state quantities over one interval `[0,t]`
- pre-observation (`y0`) and post-observation (`y1`) objects

Primary code reference:

- `legacy/qmodel.h`

Primary manuscript references:

- `papers/macroir-elife-2025/docs/manuscript-drafts/elife-macroir-revised3.tex`
- `papers/macroir-elife-2025/docs/manuscript-drafts/elife-template.tex`
- `papers/macroir-elife-2025/docs/manuscript-drafts/elife-macroir-merged.tex`

## 1. Naming conventions used in code

- `0t`: object involving the start of the interval (`0`) and the end of the interval (`t`)
- `y0`: before conditioning on the observed interval average
- `y1`: after conditioning on the observed interval average
- `t11`, `t10`, `t20`, `t2`: historical labels for propagated or updated one-time moments

The manuscript notation is cleaner than the code notation. The code keeps extra intermediate objects for debugging and diagnostics.

## 2. Main mapping

### Start-of-interval one-time prior

- Code: `P_mean`
- Manuscript: `\boldsymbol{\mu}^{\mathrm{prior}}_0`
- Meaning: prior mean occupancy vector at the start of the interval

- Code: `P_Cov`
- Manuscript: `\Sigma^{\mathrm{prior}}_0`
- Meaning: prior covariance of the occupancy vector at the start of the interval

### End-of-interval one-time prior before observation update

- Code: `P_mean_t11_y0` for averaging mode 2, `P_mean_t2_y0` for averaging mode 1
- Manuscript: `\boldsymbol{\mu}^{\mathrm{prior}}_0 P(t)` or equivalently `\boldsymbol{\mu}^{\mathrm{prior}}(t)`
- Meaning: propagated prior mean occupancy at the end of the interval before conditioning on the observation

- Code: `P_Cov_t11_y0` for averaging mode 2, `P_Cov_t2_y0` for averaging mode 1
- Manuscript: `\Sigma^{\mathrm{prior}}_t`
- Meaning: propagated prior covariance at the end of the interval before conditioning on the observation

The propagation formula in the manuscript is the standard one:

`\Sigma^{\mathrm{prior}}_t = P(t)^{\mathrm T} (\Sigma^{\mathrm{prior}}_0 - \mathrm{diag}(\mu^{\mathrm{prior}}_0)) P(t) + \mathrm{diag}(\mu^{\mathrm{prior}}_0 P(t))`

This appears in the manuscript drafts and is implemented directly in `qmodel.h`.

### End-of-interval posterior after observation update

- Code: `P_mean_t10_y1` for averaging mode 2, `P_mean_t1_y1` for averaging mode 1
- Manuscript: intermediate posterior one-time mean after assimilating `y_{0\to t}^{\mathrm{obs}}`
- Meaning: posterior mean occupancy at the end of the interval

- Code: `P_Cov_t10_y1` for averaging mode 2, `P_Cov_t1_y1` for averaging mode 1
- Manuscript: intermediate posterior one-time covariance after assimilating `y_{0\to t}^{\mathrm{obs}}`
- Meaning: posterior covariance at the end of the interval

### Boundary-state prior mean over the interval

- Code: `P_mean_0t_y0`
- Manuscript: `\mu^{\mathrm{prior}}_{0\to t}` from Eq. `boundary_mean_prior`
- Meaning: prior mean of the boundary-state distribution over ordered pairs `(i_0,i_t)`

Elementwise:

`(P_mean_0t_y0)_{i_0,i_t} = P(X_0=i_0, X_t=i_t) = (\mu^{\mathrm{prior}}_0)_{i_0} P_{i_0\to i_t}(t)`

Code formula:

`P_mean_0t_y0 = diag(\mu^{\mathrm{prior}}_0) P(t)`

This is the exact matrix form of the manuscript’s boundary-state prior mean.

### Boundary-state posterior mean over the interval

- Code: `P_mean_0t_y1`
- Manuscript: `\mu^{\mathrm{post}}_{0\to t}`
- Meaning: posterior mean of the boundary-state distribution after conditioning on the observed interval average

This is the matrix-form reduced representation of the posterior over start/end pairs.

### Reduced start-to-end cross-covariance

- Code: `P_cross_cov_0t_y0`
- Manuscript: no dedicated symbol in the current revised draft; derived reduced object associated with the full boundary covariance `\Sigma^{\mathrm{prior}}_{0\to t}`
- Meaning: cross-covariance between occupancy indicators at the start and end of the interval before conditioning on the observation

Elementwise:

`(P_cross_cov_0t_y0)_{i,j} = \mathrm{Cov}(x_0^{(i)}, x_t^{(j)})`

where `x_0` and `x_t` are the occupancy-indicator vectors at the start and end of the interval.

Code formula:

`P_cross_cov_0t_y0 = (\Sigma^{\mathrm{prior}}_0 - diag(\mu^{\mathrm{prior}}_0)) P(t) + diag(\mu^{\mathrm{prior}}_0) P(t)`

which simplifies to:

`P_cross_cov_0t_y0 = \Sigma^{\mathrm{prior}}_0 P(t)`

This is a reduced `k x k` object. It is not the full boundary covariance over pair states, which would be `k^2 x k^2`.

### Reduced start-to-end cross-covariance after observation update

- Code: `P_cross_cov_0t_y1`
- Manuscript: reduced posterior analogue of the previous object; no dedicated symbol in the current revised draft
- Meaning: cross-covariance between start and end occupancies after conditioning on the observed interval average

For averaging mode 2, the code updates it by a Kalman-like rank-1 correction.

## 3. Relation to the full boundary covariance in the manuscript

The old manuscript draft writes the full boundary-state covariance explicitly:

`(\Sigma^{\mathrm{prior}}_{0\to t})_{(i_0\to i_t)(j_0\to j_t)}`

This is a covariance over ordered pairs of states, so it lives in a formal `k^2`-dimensional space.

The code does not construct this object explicitly in the main MacroIR recursion. Instead, it keeps reduced contractions of it:

- `P_mean_0t_y0`: the boundary-state mean as a `k x k` matrix
- `P_cross_cov_0t_y0`: the reduced start/end cross-covariance as a `k x k` matrix
- `GS`, `gS0`, `d_GS`: contractions needed for the Kalman-like update

This is exactly the computational shortcut described in the manuscript when it says that MacroIR avoids explicit construction of the full `k^2` boundary covariance.

## 4. What is easy to confuse

### `P_mean_0t_y0` versus `P_cross_cov_0t_y0`

- `P_mean_0t_y0` is a joint probability table over start and end states.
- `P_cross_cov_0t_y0` is a covariance matrix between start occupancy indicators and end occupancy indicators.

They have the same shape, but they are not the same object.

For a single channel:

`P_cross_cov_0t_y0 = P_mean_0t_y0 - \mu_0^{\mathrm T} \mu_t`

For the Gaussian macroscopic occupancy recursion used in MacroIR, the code stores the reduced cross-covariance directly.

### `P_Cov_t11_y0` versus `P_cross_cov_0t_y0`

- `P_Cov_t11_y0` is a same-time covariance at the end of the interval.
- `P_cross_cov_0t_y0` is a cross-time covariance from interval start to interval end.

## 5. Suggested manuscript-symbol aliases for internal use

If we want a stable internal notation in comments or figures, the following aliases are the most consistent:

- `P_mean` -> `\mu_0`
- `P_Cov` -> `\Sigma_0`
- `P_mean_t11_y0` or `P_mean_t2_y0` -> `\mu_t^-`
- `P_Cov_t11_y0` or `P_Cov_t2_y0` -> `\Sigma_t^-`
- `P_mean_0t_y0` -> `M_{0,t}^-`
- `P_mean_0t_y1` -> `M_{0,t}^+`
- `P_cross_cov_0t_y0` -> `C_{0,t}^-`
- `P_cross_cov_0t_y1` -> `C_{0,t}^+`

Here:

- `M_{0,t}` denotes the boundary-state mean matrix over ordered pairs
- `C_{0,t}` denotes the reduced start/end cross-covariance
- `-` means pre-observation
- `+` means post-observation

These aliases are not yet used in the manuscript, but they are mathematically unambiguous.

