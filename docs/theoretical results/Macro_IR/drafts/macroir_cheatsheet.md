# MacroIR Developer Cheat Sheet

## State Variables (C++ notation)
- `mu0`        → \( \boldsymbol{\mu}_0 \)
- `Sigma0`     → \( \boldsymbol{\Sigma}_0 \)
- `P_t`        → \( \mathbf{P}(t) \)
- `Gamma_bar`  → \( \overline{\Gamma}_{i_0\to i_t} \)
- `V_bar`      → \( \overline{V}_{i_0\to i_t} \)

## Boundary Contractions
Start-index contraction:
\[
w0\_bar(i_0)=\sum_{i_t}P_{i_0\to i_t}\,W_{i_0\to i_t}.
\]

Intrinsic variance:
\[
sigma2\_gamma0(i_0)
\, =
\sum_{i_t}P_{i_0\to i_t}\,V_{i_0\to i_t}.
\]

## Predictive Mean
\[
ybar\_pred = N_{\text{ch}}\,(mu0\cdot gamma0\_bar).
\]

## Unified Tilde Operator
Tilde = **one** operator.  
Different expressions yield vector or scalar automatically.

### Vector expression
\[
\widetilde{u^T\Sigma}
\, =
\overline{u}_0^T(\Sigma_0-\mathrm{diag}\,\mu_0)
+(\overline{U}\circ P)^T\mu_0.
\]

### Scalar expression
\[
\widetilde{u^T\Sigma w}
\, =
\overline{u}_0^T(\Sigma_0-\mathrm{diag}\,\mu_0)\overline{w}_0
\, +
\mu_0^T
[(\overline{U}\circ P)\circ(\overline{W}\circ P)]\mathbf{1}.
\]

Key special cases:
- `tilde_gamma_Sigma_vec` ↔ \(\widetilde{\gamma^T\Sigma}\)
- `tilde_gamma_Sigma_gamma_scalar` ↔ \(\widetilde{\gamma^T\Sigma\gamma}\)

## Predictive Variance
\[
sigma2\_ybar =
epsilon2\_0t
\, + N_{\text{ch}}\widetilde{\gamma^T\Sigma\gamma}
\, + N_{\text{ch}}\sum_{i_0}mu0(i_0)\,sigma2\_gamma0(i_0).
\]

## Measurement Update
Innovation:
\[
\delta = ybar\_obs - ybar\_pred.
\]

Mean:
\[
mu\_post=\mu\_{\text{prior}}
+\frac{1}{\sigma^2}\Sigma\_{\text{prop}}\;\widetilde{\gamma^T\Sigma}\;\delta.
\]

Covariance:
\[
Sigma\_post
=Sigma\_{\text{prop}}
-\frac{1}{\sigma^2}
\widetilde{\gamma^T\Sigma}^T\widetilde{\gamma^T\Sigma}.
\]

## Propagation
\[
\mu\_{\text{prior}}=\mu_0 P,
\qquad
Sigma\_{\text{prop}}
\, =
P^T(\Sigma_0-\mathrm{diag}\,\mu_0)P-\mathrm{diag}(\mu\_{\text{prior}}).
\]

**Always remember:**  
Tilde is one operator; shapes differ because of expression structure, not operator type.
`