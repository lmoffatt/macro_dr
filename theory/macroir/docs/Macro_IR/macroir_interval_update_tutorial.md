# MacroIR Interval Update: Tutorial, Supplement, and Paper Section

This document explains the MacroIR interval update in three layers:

1. **Developer Tutorial (C++/MacroIR)**  
2. **Supplementary Material (Full Derivations)**  
3. **Main Paper Section (Condensed; refers to Supplement)**

---

## A.1 Developer Tutorial (MacroIR C++ Implementation Guide)

### Overview

MacroIR performs Bayesian updates using *interval-averaged currents*.  
Instead of tracking all possible boundary states \((i_0 → i_t)\), MacroIR uses analytic marginalization to avoid constructing the full \(K^2 \times K^2\) boundary covariance.

Key quantities:

- `mu_0` — prior mean state distribution at interval start  
- `Sigma_0` — prior covariance  
- `P_t` — transition matrix \(e^{Qt}\)  
- `Gamma_bar(i0, it)` — mean interval current conditional on start/end states  
- `Var_Gamma_bar(i0, it)` — conditional variance  
- `N_ch` — number of channels in the patch

---

## A.1.1 Predicted Mean Interval Current

### Explicit Formula

\[
(\overline{\gamma}_0)_{i_0} = \sum_{i_t} P_{i_0\to i_t}(t)\,\overline{\Gamma}_{i_0\to i_t}.
\]

\[
\overline{y}^{\mathrm{pred}}_{0\rightarrow t}
= N_{\mathrm{ch}}\,\mu_0 \cdot \overline{\gamma}_0.
\]

### C++ Implementation Sketch

```cpp
for (int i0 = 0; i0 < K; ++i0) {
    gamma_bar0[i0] = 0.0;
    for (int it = 0; it < K; ++it)
        gamma_bar0[i0] += P_t(i0,it) * Gamma_bar(i0,it);
}

double y_pred = 0.0;
for (int i0 = 0; i0 < K; ++i0)
    y_pred += mu_0[i0] * gamma_bar0[i0];

y_pred *= N_ch;
```

---

## A.1.2 Predicted Interval Variance

### Components

- Measurement noise:

$$
\epsilon^2_{0\to t} = \epsilon^2/t + \nu^2.
$$

- Conditional variance:

$$
(\sigma^2_{\overline{\gamma}_0})_{i_0}
= \sum_{i_t} P_{i_0\to i_t}(t)\operatorname{Var}(\overline{\Gamma}_{i_0\to i_t}).
$$

-- Main structural term:

$$
\widetilde{\gamma^{T}\Sigma\gamma}
= \overline{\gamma}_0^T(\Sigma_0 - \mathrm{diag}(\mu_0))\overline{\gamma}_0 + \sum_{i_0}\mu_0[i_0]\sum_{i_t}(\Gamma\_bar(i_0,i_t)P\_t(i_0,i_t))^2.
$$


### Combined Variance

$$
\sigma^2_{\overline{y}^{\mathrm{pred}}_{0\to t}}
= \epsilon^2_{0\to t}+ N_{\mathrm{ch}} \widetilde{\gamma^{T}\Sigma\gamma} + N_{\mathrm{ch}} \sum_{i_0} \mu_0[i_0] (\sigma^2_{\overline{\gamma}_0})_{i_0}.
$$


---

## A.1.3 Cross-Covariance Vector \( \widetilde{\gamma^{T}\Sigma} \)

### Explicit Formula
\[
(\widetilde{\gamma^{T}\Sigma})_{i_t}
= \overline{\gamma}_0^T(\Sigma_0 - \mathrm{diag}(\mu_0))P_t[:,i_t]+ \sum_{i_0}\mu_0[i_0]P_{i_0\to i_t}(t)\overline{\Gamma}_{i_0\to i_t}.
\]

### C++ Sketch

```cpp
for (int it = 0; it < K; ++it) {
    double sum1 = 0.0;
    for (int j0 = 0; j0 < K; ++j0) {
        double inner = 0.0;
        for (int i0 = 0; i0 < K; ++i0) {
            double S = Sigma_0(j0,i0);
            if (j0 == i0) S -= mu_0[j0];
            inner += S * P_t(i0,it);
        }
        sum1 += gamma_bar0[j0] * inner;
    }

    double sum2 = 0.0;
    for (int i0 = 0; i0 < K; ++i0)
        sum2 += mu_0[i0] * P_t(i0,it) * Gamma_bar(i0,it);

    tilde_vec[it] = sum1 + sum2;
}
```

---

## A.1.4 Final Covariance Update

\[
\Sigma^{\text{prop}}(t)
= P_t^T(\Sigma_0 - \mathrm{diag}(\mu_0))P_t + \mathrm{diag}(\mu(t)).
\]

\[
\Sigma(t)
= \Sigma^{\mathrm{prop}}(t)- \frac{1}{\sigma^2_{\overline{y}^{\mathrm{pred}}_{0\to t}}}
\widetilde{\gamma^{T}\Sigma}^T\widetilde{\gamma^{T}\Sigma}.
\]

---

## B Supplementary Material (Full Derivations)

### Includes:

- Boundary-state mean derivation  
- Boundary-state covariance derivation (multinomial splitting)  
- Full expansion of  
  \(\widetilde{\gamma^{T}\Sigma\gamma}\)  
  and  
  \(\widetilde{\gamma^{T}\Sigma}\)  
- Marginalization identities  
- Final Bayesian update identity

*(Full LaTeX version provided separately in Document B.)*

---

## C Main Paper Section

### Summary

MacroIR evaluates interval-averaged currents using boundary-state statistics but avoids explicit \(K^2\)-dimensional boundary covariances via analytic marginalization. The key collapsed quantities are:

\[
\overline{\gamma}_0 = P(t)\odot \overline{\Gamma},
\quad
\widetilde{\gamma^{T}\Sigma\gamma},
\quad
\widetilde{\gamma^{T}\Sigma}.
\]

These determine the interval likelihood:
\[
\mathcal{L} = \mathcal{N}(\overline{y}^{\mathrm{obs}};\overline{y}^{\mathrm{pred}},\sigma^2_{\overline{y}^{\mathrm{pred}}}),
\]
and the updated covariance:
\[
\Sigma(t) = P_t^T(\Sigma_0 - D_{\mu_0})P_t + D_{\mu(t)}- \frac{1}{\sigma^2_{\overline{y}^{\mathrm{pred}}}} 
\widetilde{\gamma^{T}\Sigma}^T\widetilde{\gamma^{T}\Sigma}.
\]

All derivations appear in Supplementary Section Sx.
