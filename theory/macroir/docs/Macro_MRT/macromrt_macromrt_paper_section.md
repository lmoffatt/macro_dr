\section{Instantaneous-Current Bayesian Inference with State-Dependent Noise (MacroMRT)}

MacroMRT performs Bayesian updating using \emph{instantaneous} current
measurements (sampled at the end of each kinetic interval), accounting
for state-dependent open-channel noise via a Taylor/Laplace correction.
For an interval of duration $t$, the prior state distribution is
$(\mu_0,\Sigma_0)$ and the transition matrix is $P(t)=e^{Qt}$. Each
state contributes a per-state mean conductance and an open-channel
noise variance; over the recording window MacroMRT uses the
boundary-conditioned moments $\overline{\Gamma}_{i_0\to i_t}$ and
$\overline{V}_{i_0\to i_t}$ marginalised over the end state to yield
per-$i_0$ scalars
\[
(\overline{\gamma}_0)_{i_0}
= \sum_{i_t} P_{i_0\to i_t}(t)\,\overline{\Gamma}_{i_0\to i_t},
\qquad
(\overline{\sigma^2}_0)_{i_0}
= \sum_{i_t} P_{i_0\to i_t}(t)\,\overline{V}_{i_0\to i_t}.
\]
These are the same boundary-conditioned moments used by MacroIR,
collapsed to per-$i_0$ scalars by marginalising over $i_t$ rather than
carrying both indices through a tilde operator. They reduce to pure
per-state $\gamma_i,\sigma^2_i$ when $Q$ is constant over $[0,t]$ and
the recording window is negligible relative to the kinetic timescales.

The predicted instantaneous measurement is
\[
\overline{y}^{\mathrm{pred}}
= N_{\mathrm{ch}}\,\mu_0\cdot\overline{\gamma}_0,
\]
and its predictive variance is
\[
V = \epsilon^2_{0\to t}
+ N_{\mathrm{ch}}\sum_{i_0}\mu_{0,i_0}\,(\overline{\sigma^2}_0)_{i_0},
\]
i.e.\ instrument noise plus the ensemble-averaged open-channel
contribution. Compared with MacroIR, the boundary cross-correlation
term $N_{\mathrm{ch}}\widetilde{\gamma^{T}\Sigma\gamma}$ is absent;
this is the central structural difference between the M and I families.

The MRT measurement update is a rank-1 quasi-Laplace correction. With
innovation $\delta=\overline{y}^{\mathrm{obs}}-\overline{y}^{\mathrm{pred}}$
and effective direction
\[
v = \overline{\gamma}_0 + \frac{\delta}{V}\,\overline{\sigma^2}_0,
\qquad
s = v^{T}\Sigma^{\mathrm{prop}}v,
\]
the posterior covariance and mean at the end of the interval are
\[
\Sigma^{\mathrm{post}}
= \Sigma^{\mathrm{prop}}
- \frac{N_{\mathrm{ch}}}{V + N_{\mathrm{ch}}s}\,
\Sigma^{\mathrm{prop}}vv^{T}\Sigma^{\mathrm{prop}},
\]
\[
\mu^{\mathrm{post}}
= \mu^{\mathrm{prior}}
+ \frac{\delta}{2V}\,\Sigma^{\mathrm{post}}\,(\overline{\gamma}_0 + v),
\]
with the standard propagation
$\mu^{\mathrm{prior}}=\mu_0 P(t)$ and
$\Sigma^{\mathrm{prop}}=P(t)^{T}(\Sigma_0-\mathrm{diag}\,\mu_0)P(t)+\mathrm{diag}(\mu^{\mathrm{prior}})$.
Setting $\overline{\sigma^2}_0=0$ recovers the Moffatt (2007) MacroR
update; a rank-2 form including the $\tfrac{1}{2}\log V$ term of the
Gaussian likelihood is given in the supplement and is required only
when $N_{\mathrm{ch}}$ is small enough that the heteroscedastic
log-determinant influences the MAP.

MacroMRT is the controlled baseline against which MacroIRT is
evaluated: identical Taylor σ² treatment, identical boundary-conditioned
moment inputs, but no boundary-state lift. Any improvement of IRT
over MRT is therefore attributable to the boundary cross-correlation
captured by the tilde operator and not to the σ² correction itself.

Full derivations of the gradient, rank-2 Hessian, Woodbury inversion,
and Newton-MAP iteration are provided in Supplementary Material
Section Sx.
