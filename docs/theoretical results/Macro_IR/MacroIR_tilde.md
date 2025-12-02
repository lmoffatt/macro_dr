# MacroIR: Interval-Based Bayesian Updating Using Boundary-State Lifting

MacroIR performs Bayesian updating of a macroscopic state-occupancy distribution for an ensemble of Markov channels based on **interval-averaged** observations.  
The key difficulty is that an interval-averaged current depends on the entire hidden *trajectory* of each channel over \([0,t]\). MacroIR resolves this through:

\, - **Boundary states** \((i_0,i_t)\), encoding start–end microscopic states.
\, - A **unified lift–modulate–collapse transformation**, called the **tilde operator**, that embeds detailed boundary-state information into tractable macroscopic updates.

This document defines:

1. The interval problem and the boundary-state representation  
2. Predictive mean and variance  
3. The generalized tilde operator and its specializations  
4. The MacroIR interval update  
5. Conceptual summary


# 1. Interval Problem and Boundary-State Geometry

## 1.1 Macroscopic state variables

For an ensemble of \(N_{\text{ch}}\) independent channels with \(K\) microscopic states, let
\[
\boldsymbol{\mu}_0 \in \mathbb{R}^K,
\qquad
\boldsymbol{\Sigma}_0 \in \mathbb{R}^{K\times K}
\]
be the prior mean and covariance at time \(t=0\).

Microscopic state evolution is Markovian with generator \(\mathbf{Q}\) and transition matrix:
\[
\mathbf{P}(t) = \exp(\mathbf{Q}t),
\qquad
P_{i_0\to i_t}(t) = [\mathbf{P}(t)]_{i_0,i_t}.
\]

## 1.2 Boundary states \((i_0,i_t)\)

A **boundary state** specifies a channel’s state at both ends of the interval \([0,t]\):
\[
(i_0,i_t) .
\]

For each boundary pair, precompute:
\[
\overline{\Gamma}_{i_0\to i_t}
\qquad
\text{(mean interval-averaged current)},
\]
\[
\operatorname{Var}(\overline{\Gamma}_{i_0\to i_t})
\qquad
\text{(interval-averaged current variance)}.
\]

These encode the microscopic model’s behavior over the interval.

To obtain a start-state–indexed conditional mean:
\[
(\overline{\gamma}_0)_{i_0} =
\sum_{i_t} P_{i_0\to i_t}(t)\,\overline{\Gamma}_{i_0\to i_t}.
\]


# 2. Predictive Mean and Variance of the Interval-Averaged Observation

## 2.1 Predictive mean

The predicted macroscopic interval-averaged current is:
\[
\overline{y}^{\mathrm{pred}}_{0\to t}
\, =
N_{\text{ch}}\,
\boldsymbol{\mu}_0\cdot \overline{\boldsymbol{\gamma}}_0 .
\]

## 2.2 Predictive variance

Measurement noise:
\[
\epsilon^2_{0\to t} \qquad = \frac{\epsilon^2}{t} + \nu^2.
\]

Intrinsic channel variability consists of:

1. Uncertainty in the initial ensemble state,  
2. Variance of the interval current conditioned on boundary states.

Define:
\[
(\sigma^2_{\overline{\gamma}_0})_{i_0}
\, =
\sum_{i_t}
P_{i_0\to i_t}(t)\,
\operatorname{Var}(\overline{\Gamma}_{i_0\to i_t}).
\]

The full predictive variance is:
\[
\sigma^2_{\overline{y}^{\mathrm{pred}}_{0\to t}}
\, =
\epsilon^2_{0\to t}
\, +
N_{\text{ch}}
\,\widetilde{\gamma^{T}\Sigma\gamma}
\, +
N_{\text{ch}}
\sum_{i_0}
(\boldsymbol{\mu}_0)_{i_0}\,
(\sigma^2_{\overline{\gamma}_0})_{i_0}.
\]

The central term \(\widetilde{\gamma^{T}\Sigma\gamma}\) is a **triple-tilde** scalar defined below.


# 3. The Tilde Operator: Generalized Lift–Modulate–Collapse

The **tilde operator** is a three-stage transformation applied to any state-space object \(X\):

1. **Lift**: reinterpret \(X\) in boundary-state coordinates.  
2. **Modulate**: apply interval-conditioned microscopic weights (e.g.  
   \(\overline{\Gamma}_{i_0\to i_t},\overline{\mathbf{V}},P(t)\)).  
3. **Collapse**: sum over the initial boundary index \(i_0\), yielding a macroscopic object.

Formally, a tilde always has the form:
\[
\widetilde{X}_{i_t}
\, =
\sum_{i_0}
\text{Lift}(X)_{i_0,i_t}
\cdot
\text{Modulate}(i_0,i_t).
\]

The operator depends on the interval:
\[
\widetilde{(\cdot)} \equiv \widetilde{(\cdot)}_{0\to t},
\]
though the interval is usually clear from context.

## 3.1 First-order tilde: mean and boundary expectations (special cases)

Mean propagation is a first-order tilde:
\[
\boldsymbol{\mu}^{\mathrm{prior}}(t)
\, =
\boldsymbol{\mu}_0\,\mathbf{P}(t)
\, =
\widetilde{\boldsymbol{\mu}_0}.
\]

Similarly, the start-state boundary-average
\[
(\overline{\gamma}_0)_{i_0}
\, =
\sum_{i_t}
P_{i_0\to i_t}(t)\,\overline{\Gamma}_{i_0\to i_t}
\]
is a *backward* first-order tilde (or forward tilde with time indices transposed).

These unify neatly under the lift–modulate–collapse schema.

## 3.2 Second-order tilde: covariance propagation (special case)

Covariance propagation follows exactly the same schema applied to a second-order object:
\[
\boldsymbol{\Sigma}^{\mathrm{prop}}(t)
\, =
\mathbf{P}(t)^T\,
(\boldsymbol{\Sigma}_0 - \mathrm{diag}(\boldsymbol{\mu}_0))\,
\mathbf{P}(t)
\, -
\mathrm{diag}(\boldsymbol{\mu}^{\mathrm{prior}}(t)),
\]
which is a first-second–order tilde transformation.

## 3.3 True MacroIR tilde operators (second and third order)

For inference, we require tilde applied to second-order contractions involving the covariance.

### Double-tilde (vector)
For \(u\in\mathbb{R}^K\) with boundary weights \(\overline{\mathbf{W}}\):
\[
\widetilde{u^{T}\Sigma}
\, =
\overline{\mathbf{w}}_0^{T}\,
(\boldsymbol{\Sigma}_0 - \operatorname{diag}(\boldsymbol{\mu}_0))
\, +
(\overline{\mathbf{W}}\circ \mathbf{P}(t))^{T}\boldsymbol{\mu}_0,
\]
where
\[
(\overline{\mathbf{w}}_0)_{i_0}
\, =
\sum_{i_t}
P_{i_0\to i_t}(t)\,\overline{W}_{i_0\to i_t}.
\]

### Triple-tilde (scalar)
\[
\widetilde{u^{T}\Sigma w}
\, =
\overline{\mathbf{u}}_0^{T}
(\boldsymbol{\Sigma}_0 - \operatorname{diag}(\boldsymbol{\mu}_0))
\overline{\mathbf{w}}_0
+
\boldsymbol{\mu}_0^{T}
\left[
(\overline{\mathbf{U}}\circ \mathbf{P}(t))
\circ
(\overline{\mathbf{W}}\circ \mathbf{P}(t))
\right]\mathbf{1}.
\]

For current-based likelihoods:
\, - Triple tilde: \(\widetilde{\gamma^{T}\Sigma\gamma}\).  
\, - Under TaylorIR: mixed forms \(\widetilde{\gamma^{T}\Sigma v}\).


## 3.4 Directionality

Tilde is defined as contraction over the **initial** boundary index \(i_0\).  
If contraction over \(i_t\) is needed, we apply tilde to:
\[
\overline{\mathbf{W}}^{T},
\qquad
\mathbf{P}(t)^{T}.
\]

No independent backward-tilde operator is necessary.


# 4. MacroIR Interval Update

## 4.1 Propagation step

Mean:
\[
\boldsymbol{\mu}^{\mathrm{prior}}(t)
\, =
\boldsymbol{\mu}_0 \mathbf{P}(t).
\]

Covariance:
\[
\boldsymbol{\Sigma}^{\mathrm{prop}}(t)
\ =
\mathbf{P}(t)^T
(\boldsymbol{\Sigma}_0 - \operatorname{diag}\boldsymbol{\mu}_0)
\mathbf{P}(t)
\, -
\operatorname{diag}(\boldsymbol{\mu}^{\mathrm{prior}}(t)).
\]

## 4.2 Measurement update

Let
\[
\widetilde{\gamma^{T}\Sigma} \in \mathbb{R}^K,
\qquad
\sigma^2 \qquad = \sigma^2_{\overline{y}^{\mathrm{pred}}_{0\to t}}.
\]

Define the innovation:
\[
\delta
\ =
\overline{y}^{\mathrm{obs}}_{0\to t}
\, -
\overline{y}^{\mathrm{pred}}_{0\to t}.
\]

Mean update:
\[
\boldsymbol{\mu}^{\mathrm{post}}(t)
\, =
\boldsymbol{\mu}^{\mathrm{prior}}(t)
\, +
\frac{1}{\sigma^2}
\boldsymbol{\Sigma}^{\mathrm{prop}}(t)
\,
\widetilde{\gamma^{T}\Sigma}
\,\delta.
\]

Covariance update:
\[
\boldsymbol{\Sigma}^{\mathrm{post}}(t)
\, =
\boldsymbol{\Sigma}^{\mathrm{prop}}(t)
\, -
\frac{1}{\sigma^2}
(\widetilde{\gamma^{T}\Sigma})^{T}
(\widetilde{\gamma^{T}\Sigma}).
\]

This is a rank-1 scalar-Kalman measurement update in the macroscopic state space.


# 5. Conceptual Summary

MacroIR reduces a trajectory-dependent microscopic observation to a macroscopic Bayesian update through the following structure:

1. **Boundary states** \((i_0,i_t)\) encode interval-conditioned behavior without explicit trajectory enumeration.  
2. **Boundary-conditioned matrices** \(\overline{\Gamma}\), \(\overline{V}\), and \(\mathbf{P}(t)\)  
   capture microscopic physics over \([0,t]\).  
3. The **tilde operator** is a unified lift–modulate–collapse transform:  
   - First-order tilde covers mean and boundary expectations,  
   - Second-order tilde covers covariance propagation,  
   - Second/third-order tilde contractions yield \(\widetilde{\gamma^{T}\Sigma}\) and  
     \(\widetilde{\gamma^{T}\Sigma\gamma}\) for inference.  
4. The update equation reduces to a **scalar Kalman correction** in \(K\)-dimensional macroscopic state space.  
5. This avoids constructing the full \(K^2\times K^2\) boundary covariance while preserving all necessary microscopic detail.

The tilde operator is the core mechanism that makes MacroIR efficient, expressive, and compatible with microscopic kinetics without explicit path integration.
