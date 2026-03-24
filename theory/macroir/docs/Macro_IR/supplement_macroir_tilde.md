
# Supplement: Full Derivations of the MacroIR Boundary-State Operators

This supplement expands all intermediate steps that were omitted in the main text for clarity.

---

# S.1 Boundary-State Geometry

A microscopic channel trajectory over \([0,t]\) is summarized by its start/end states \((i_0,i_t)\).  
Given generator \(Q\in\mathbb{R}^{K\times K}\), the transition matrix is:
\[
P(t)=\exp(Qt).
\]

For every boundary pair:
\[
\overline{\Gamma}_{i_0\to i_t}
= \mathbb{E}\left[\overline{I}\_{0\to t} \mid X(0)=i_0,X(t)=i_t\right],
\]
\[
\overline{V}_{i_0\to i_t}
= \mathrm{Var}\left[\overline{I}\_{0\to t} \mid i_0,i_t\right].
\]

These encode the entire trajectory distribution via conditional path measures (Feynman–Kac representation), but MacroIR uses only \((i_0,i_t)\).

---

# S.2 Collapse From Boundary Space to State Space

We often need expectations conditioned only on \(i_0\).  
Marginalization over \(i_t\) gives:
\[
(\overline{\gamma}_0)_{i_0}
\, =
\sum_{i_t} P_{i_0\to i_t}(t)\,\overline{\Gamma}_{i_0\to i_t}.
\]

This object arises naturally when computing the predictive mean.

Similarly, the intrinsic variance term:
\[
(\sigma^2_{\overline{\gamma}_0})_{i_0}
\, =
\sum_{i_t}P_{i_0\to i_t}(t)\,\overline{V}_{i_0\to i_t}.
\]

---

# S.3 Deriving the Predictive Mean

For \(N_{\text{ch}}\) independent channels:

\[
\overline{y}^{\text{pred}}
\, =
N_{\text{ch}}\,\mathbb{E}[\overline{I}_{0\to t}]
\, =
N_{\text{ch}}\,\sum_{i_0}\mu_{0,i_0}\,(\overline{\gamma}_0)_{i_0}.
\]

Here \(\mu_{0,i_0}\) is the expected number of channels initially in state \(i_0\).

---

# S.4 Deconstructing the Predictive Variance

We seek:
\[
\sigma^2_{\overline{y}}
= \mathrm{Var}\left[\frac{1}{N_{\text{ch}}}\sum_k \overline{I}_{0\to t}^{(k)}\right].
\]

Direct computation yields three terms:

1. Measurement noise  
2. Contribution from the uncertainty in state distribution \((\mu_0,\Sigma_0)\)  
3. Contribution from conditional microscopic variance  

The second term is:
\[
\mathrm{Var}\left[\sum_i N_i(0)(\overline{\gamma}_0)_i\right]
= \widetilde{\gamma^T\Sigma\gamma}.
\]

This is a contraction of a quadratic form against the lifted boundary structure.  
Its explicit form is given in S.6 below.

---

# S.5 Unified View of the Tilde Operator

We want an operator $\widetilde{\cdot}$ satisfying:

- linearity in arguments,
- compatibility with the Markov kernel \(P(t)\),
- correct reduction to known cases (mean propagation, covariance propagation),
- ability to incorporate boundary-conditioned second-order effects.

The operator arises naturally from writing the macroscopic state-space expressions in terms of expectations over the joint distribution of \((i_0,i_t)\).

---

# S.6 Explicit Derivation of \(\widetilde{u^T\Sigma}\)

Consider the random variable:
\[
Z = u^T X, \qquad X\sim\text{state counts at }0.
\]

Propagate correlations through the interval using boundary-conditioned weights:

\[
\mathbb{E}[Z(i_t)]
\, =
\sum_{i_0}
u_{i_0}\,(\Sigma_0-\mathrm{diag}\,\mu_0)_{i_0,*}\,P_{i_0\to i_t}
\, +
\sum_{i_0}
\mu_{i_0}\overline{U}_{i_0\to i_t} P_{i_0\to i_t}.
\]

This yields exactly:
\[
\widetilde{u^T\Sigma}
\, =
\overline{u}_0^T(\Sigma_0-\mathrm{diag}\,\mu_0)
\, + (\overline{U}\circ P)^{T}\mu_0.
\]

The first term captures covariance propagation; the second captures boundary-conditioned fluctuation.

---

# S.7 Deriving the Scalar Expression \(\widetilde{u^T\Sigma w}\)

We expand the quadratic form:

\[
u^T \Sigma_0 w
= \mathbb{E}[(u\cdot X)(w\cdot X)] - u\cdot\mu_0 \, w\cdot\mu_0.
\]

MacroIR modifies this by replacing each endpoint-\(i_t\) contribution with the conditional expectation and variance implied by \(\overline{U}, \overline{W},P\):

\[
\widetilde{u^T\Sigma w}
\, =
(\overline{u}_0)^T(\Sigma_0-\mathrm{diag}\,\mu_0)\,\overline{w}_0
\, +
\mu_0^T
\left[
(\overline{U}\circ P)\circ(\overline{W}\circ P)
\right]\mathbf{1}.
\]

This completes the second-order contraction.

---

# S.8 The MacroIR Update

Given innovation \(\delta\) and scalar predictive variance \(\sigma^2\):

\[
\mu^{\text{post}}
\, =
\mu^{\text{prior}}
\, +
\frac{1}{\sigma^2}
\Sigma^{\text{prop}}
\widetilde{\gamma^T\Sigma}\,
\delta,
\]
\[
\Sigma^{\text{post}}
\, =
\Sigma^{\text{prop}}
\, -
\frac{1}{\sigma^2}
\widetilde{\gamma^T\Sigma}^{T}\widetilde{\gamma^T\Sigma}.
\]

These derive from standard Kalman algebra with the “design” given by a tilde-transformed object.

---

# S.9 Final Remarks

Although the tilde operator produces different ranks depending on the expression in which it appears, it is a single operator defined by a consistent lift–modulate–collapse rule.

All MacroIR second-order statistics reduce to tilde-expressions involving:
- \(\overline{\Gamma}\),
- \(\overline{V}\),
- \(P(t)\),
- \(\mu_0\),
- \(\Sigma_0\).

