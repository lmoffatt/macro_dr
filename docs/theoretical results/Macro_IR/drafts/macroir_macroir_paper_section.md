\section{Interval-Averaged Bayesian Inference with MacroIR}

MacroIR performs Bayesian updating using interval-averaged current
measurements. For an interval of duration $t$, the state distribution at the
beginning is $(\mu_0,\Sigma_0)$ and the transition matrix is $P(t)=e^{Qt}$.
The interval-averaged current depends on the ``boundary states''
$(i_0\to i_t)$ consisting of start and end states. The mean count of channels
in each boundary state is
\[
\mu_{0\to t,(i_0\to i_t)} = \mu_{0,i_0}P_{i_0\to i_t}(t),
\]
and the corresponding covariance, derived by multinomial splitting, is
\[ \Sigma_{0\to t,(i_0\to i_t)(j_0\to j_t)}= P_{i_0\to i_t}(t)\left[\Sigma_{0,i_0j_0} - \delta_{i_0j_0}\mu_{0,i_0} \right]P_{j_0\to j_t}(t)+ \delta_{i_0j_0}\delta_{i_tj_t}\,\mu_{0,i_0}P_{i_0\to i_t}(t).\]

To avoid forming the full boundary covariance, MacroIR analytically
marginalizes over initial and final states. The mean interval current given
start state $i_0$ is
\[
(\overline{\gamma}_0)_{i_0}
= \sum_{i_t} P_{i_0\to i_t}(t)\,\overline{\Gamma}_{i_0\to i_t},
\]
yielding the predicted measurement
\[
\overline{y}^{\mathrm{pred}}_{0\to t}
= N_{\mathrm{ch}}\;\mu_0\cdot\overline{\gamma}_0.
\]

The predicted variance is
\[
\sigma^2_{\overline{y}^{\mathrm{pred}}_{0\to t}} =
\epsilon^2_{0\to t} + N_{\mathrm{ch}}\widetilde{\gamma^{T}\Sigma\gamma}+ N_{\mathrm{ch}}\sum_{i_0}\mu_{0,i_0}(\sigma^2_{\overline{\gamma}_0})_{i_0},
\]
where $\widetilde{\gamma^{T}\Sigma\gamma}$ is a collapsed quadratic form
obtained without constructing the boundary covariance. The cross-covariance
vector between the measurement and the end state is denoted
$\widetilde{\gamma^{T}\Sigma}\in\mathbb{R}^K$.

The end-of-interval covariance is then
\[
\Sigma(t) =
P(t)^T(\Sigma_0 - \mathrm{diag}(\mu_0))P(t) + \mathrm{diag}(\mu(t))-
\frac{1}{\sigma^2_{\overline{y}^{\mathrm{pred}}_{0\to t}}}
\widetilde{\gamma^{T}\Sigma}^T
\widetilde{\gamma^{T}\Sigma}.
\]

Full derivations of
$\widetilde{\gamma^{T}\Sigma\gamma}$ and $\widetilde{\gamma^{T}\Sigma}$,
including index-wise expansions and proofs of the marginalization identities,
are provided in Supplementary Material Section Sx.
