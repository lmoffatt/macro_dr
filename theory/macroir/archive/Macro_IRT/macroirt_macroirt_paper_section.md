\section{Interval-Averaged Bayesian Inference with State-Dependent Noise (MacroIRT)}

MacroIRT extends MacroIR by accounting for state-dependent open-channel
noise via a Taylor/Laplace correction in boundary-state space. The
boundary-conditioned moments $\overline{\Gamma}_{i_0\to i_t}$ and
$\overline{V}_{i_0\to i_t}$ enter the same lift--modulate--collapse
tilde operator that MacroIR uses, but the effective sensitivity vector
incorporates the per-state intrinsic variance.

The predicted interval-averaged measurement is
\[
\overline{y}^{\mathrm{pred}}_{0\to t}
= N_{\mathrm{ch}}\,\mu_0\cdot\overline{\gamma}_0,
\]
unchanged from MacroIR. The predictive variance is
\[
V = \epsilon^2_{0\to t}
+ N_{\mathrm{ch}}\,\widetilde{\gamma^{T}\Sigma\gamma}
+ N_{\mathrm{ch}}\sum_{i_0}\mu_{0,i_0}\,(\overline{\sigma^2}_0)_{i_0},
\]
where the third term is the boundary-conditioned open-channel variance
already present in MacroIR. The IRT correction does not enter the
predictive variance directly; it enters the measurement update via an
effective direction
\[
v = \overline{\gamma}_0 + \frac{\delta}{V}\,\overline{\sigma^2}_0,
\qquad
s = \widetilde{v^{T}\Sigma v},
\]
where $\delta=\overline{y}^{\mathrm{obs}}_{0\to t}-\overline{y}^{\mathrm{pred}}_{0\to t}$
and the effective scalar $s$ is computed via the same tilde operator
applied to $v$ rather than $\gamma$.

The IRT measurement update is a rank-1 quasi-Laplace correction in
the boundary-state representation:
\[
\Sigma^{\mathrm{post}}
= \Sigma^{\mathrm{prop}}
- \frac{N_{\mathrm{ch}}}{V + N_{\mathrm{ch}}s}\,
(\widetilde{v^{T}\Sigma})^{T}(\widetilde{v^{T}\Sigma}),
\]
\[
\mu^{\mathrm{post}}
= \mu^{\mathrm{prior}}
+ \frac{\delta}{2V}\,\Sigma^{\mathrm{post}}\,
(\widetilde{\gamma^{T}\Sigma} + \widetilde{v^{T}\Sigma}),
\]
with propagation $\mu^{\mathrm{prior}}=\mu_0 P(t)$ and
$\Sigma^{\mathrm{prop}}=P(t)^{T}(\Sigma_0-\mathrm{diag}\,\mu_0)P(t)+\mathrm{diag}(\mu^{\mathrm{prior}})$.
Setting $\overline{\sigma^2}_0=0$ gives $v=\overline{\gamma}_0$,
$\widetilde{v^{T}\Sigma}=\widetilde{\gamma^{T}\Sigma}$,
$s=\widetilde{\gamma^{T}\Sigma\gamma}$, and the update collapses to the
standard MacroIR scalar Kalman correction.

The structural difference between MacroIRT and MacroMRT is the use of
the tilde operator in the cross-covariance vector and the effective
scalar: MRT uses $\Sigma^{\mathrm{prop}}\overline{\gamma}_0$ and
$v^{T}\Sigma^{\mathrm{prop}}v$ in place of $\widetilde{\gamma^{T}\Sigma}$
and $\widetilde{v^{T}\Sigma v}$, and drops the third predictive
variance term. The IRT-vs-MRT contrast therefore isolates the value of
the boundary-state lift, holding the σ² Taylor correction fixed; the
MRT-vs-MR contrast isolates the value of the σ² correction, holding the
measurement-model class fixed. Together they decompose the IRT vs.\
MR gain into its two contributions.

Full derivations of $\widetilde{v^{T}\Sigma}$, $\widetilde{v^{T}\Sigma v}$,
the rank-2 Hessian including the $\tfrac{1}{2}\log V$ term, the
Woodbury inversion, and the Newton-MAP iteration are provided in
Supplementary Material Section Sx.
