# MacroIR Interval Update: Unified Theory and Implementation

This document provides a comprehensive exposition of the MacroIR interval-update formalism, combining conceptual clarity with rigorous mathematical foundations.

---

## 1. Problem Formulation and Core Concepts

### 1.1 The Interval Update Problem

We consider an ensemble of $N_{ch}$ identical Markov channels with K microscopic states. Over a time interval $[0,t]$ we observe the interval-averaged macroscopic current $\overline{y}_{0\to t}$. Our goal is to update the Bayesian belief over the channel state-occupancy vector at time $t$ given this observation.

**Core Challenge**: The interval-averaged current depends on the entire trajectory of each channel, not just endpoint states. Direct integration over all trajectories is intractable.

### 1.2 Macroscopic State Space

Let the population of $N_{\text{ch}}$ independent channels be described macroscopically at time $t$ by:

\[
\boldsymbol{\mu}(t) \in \mathbb{R}^K,
\qquad
\boldsymbol{\Sigma}(t) \in \mathbb{R}^{K\times K},
\]

where $K$ is the number of microscopic states. At the start of the interval $[0,t]$, define:

\[
\boldsymbol{\mu}_0 \equiv \boldsymbol{\mu}(0),
\qquad
\boldsymbol{\Sigma}_0 \equiv \boldsymbol{\Sigma}(0).
\]

These represent **expected channel counts**, not probabilities.

### 1.3 Microscopic Dynamics

Microscopic evolution is Markovian with generator $\mathbf{Q}$. Over $[0,t]$:

\[
\mathbf{P}(t) = e^{\mathbf{Q}t},
\qquad
P_{i_0 \to i_t}(t) = [\mathbf{P}(t)]_{i_0,i_t}.
\]

The indices are crucial:
- $i_0$: initial state at time $0$  
- $i_t$: final state at time $t$

---

## 2. Boundary-State Representation: The Core Insight

### 2.1 The Boundary-State Trick

The fundamental insight is to introduce **boundary states** $(i_0,i_t)$ representing start-end pairs. For each boundary pair, we precompute:

- Mean interval-averaged current: $\overline{\Gamma}_{i_0 \to i_t}$
- Variance contribution: $\overline{V}_{i_0 \to i_t}$

These are arranged in $K \times K$ matrices:

\[
\overline{\mathbf{\Gamma}},\quad
\overline{\mathbf{V}}
\qquad \text{indexed by } (i_0,i_t).
\]

### 2.2 Boundary-State Random Variables

For the ensemble, define the boundary count:
\[
N_{i_0\to i_t} = \text{number of channels that start in } i_0 \text{ and end in } i_t.
\]

The **boundary-state mean** is:
\[
\mu^{\text{prior}}_{0\to t,(i_0\to i_t)}
= \mathbb{E}[N_{i_0\to i_t}]
= (\boldsymbol{\mu}_0)_{i_0}\,P_{i_0\to i_t}(t).
\]

The **boundary-state covariance** is:
\[
\Sigma^{\text{prior}}_{0\to t,(i_0\to i_t)(j_0\to j_t)} = 
P_{i_0\to i_t}(t)
\big[(\boldsymbol{\Sigma}_0)_{i_0 j_0} - \delta_{i_0 j_0}(\boldsymbol{\mu}_0)_{i_0}\big]
P_{j_0\to j_t}(t) - \delta_{i_0 j_0}\delta_{i_t j_t}
(\boldsymbol{\mu}_0)_{i_0}P_{i_0\to i_t}(t).
\]

### 2.3 Boundary-Conditioned Expectations

The **start-indexed conditional mean** of the interval current is:

\[
(\overline{\gamma}_0)_{i_0}
= \sum_{i_t} P_{i_0 \to i_t}(t)\,\overline{\Gamma}_{i_0 \to i_t}.
\]

Similarly, intrinsic variance conditioned on $i_0$:

\[
(\sigma^2_{\overline{\gamma}_0})_{i_0}
= \sum_{i_t} P_{i_0 \to i_t}(t)\,\overline{V}_{i_0 \to i_t}.
\]

---

## 3. Predictive Mean and Variance

### 3.1 Predictive Mean

The macroscopic predicted interval-averaged current is:

\[
\overline{y}^{\text{pred}}_{0\to t}
= N_{\text{ch}}\,\boldsymbol{\mu}_0 \cdot \overline{\boldsymbol{\gamma}}_0.
\]

This is the mean across the ensemble given macroscopic uncertainty in the starting state.

### 3.2 Predictive Variance Structure

Total predictive variance:

\[
\sigma^2_{\overline{y}^{\text{pred}}}
= \epsilon^2_{0\to t}
\,+ N_{\text{ch}}\,\widetilde{\gamma^{\!T}\Sigma\gamma}
\,+ N_{\text{ch}} \sum_{i_0} (\mu_0)_{i_0}\,(\sigma^2_{\overline{\gamma}_0})_{i_0}.
\]

Where:

- **Instrumental noise**: $\epsilon^2_{0\to t} = \frac{\epsilon^2}{t} + \nu^2$
- **Uncertainty in initial distribution**: $\widetilde{\gamma^{\!T}\Sigma\gamma}$
- **Intrinsic microscopic variance**: $\sum_{i_0} \mu_{0,i_0}\,(\sigma^2_{\overline{\gamma}_0})_{i_0}$

---

## 4. The Unified Tilde Operator

### 4.1 Conceptual Definition

The tilde operator is a **single mathematical operator**, exactly like matrix multiplication: different expressions yield vectors or scalars based on operand structure, not because there are "different" tildes.

Given:
- Boundary-conditioned matrix $\overline{\mathbf{W}}_{0\to t}$
- Transition matrix $\mathbf{P}(t)$  
- Macroscopic prior $\boldsymbol{\mu}_0, \boldsymbol{\Sigma}_0$

The tilde operator performs:

1. **Lift**: Map state-space vectors into boundary space via $\overline{\mathbf{W}}$
2. **Modulate**: Combine with transition probabilities $P_{i_0\to i_t}(t)$  
3. **Collapse**: Sum out initial index $i_0$, producing objects indexed by final indices (or nothing, for scalars)

### 4.2 Directionality Principle

By definition, tilde contracts over the **initial** boundary index $i_0$:

\[
\widetilde{(\cdot)} = \sum_{i_0} (\cdot).
\]

If backward contraction is needed, apply tilde to **transposed matrices**: $\overline{\mathbf{W}}^{\!T}, \mathbf{P}(t)^{\!T}$.

**Tilde is direction-agnostic**; directionality lives in the matrices, not the operator.

---

## 5. Canonical Tilde Expressions

### 5.1 Start-Index Contraction

For any boundary-conditioned matrix $\overline{\mathbf{W}}$:

\[
(\overline{\mathbf{w}}_0)_{i_0}
= \sum_{i_t} P_{i_0\to i_t}(t)\, \overline{W}_{i_0\to i_t}.
\]

### 5.2 Vector-Valued Tilde: $\widetilde{u^{\!T}\Sigma}$

Yields a **vector** because $u^T\Sigma$ is a row vector in state-space.

\[
\boxed{
\widetilde{u^{\!T}\Sigma}
\,= \overline{\mathbf{u}}_0^{\!T}\,(\boldsymbol{\Sigma}_0 - \operatorname{diag}\boldsymbol{\mu}_0)
\,+ (\overline{\mathbf{U}}\circ \mathbf{P}(t))^{\!T} \boldsymbol{\mu}_0
}
\]

Where:
\[
(\overline{\mathbf{u}}_0)_{i_0} = \sum_{i_t} P_{i_0\to i_t}(t)\,\overline{U}_{i_0\to i_t}.
\]

When $u = \gamma$ (interval current vector): $\widetilde{\gamma^{\!T}\Sigma} \in \mathbb{R}^K$

This is the "design vector" in the MacroIR update.

### 5.3 Scalar-Valued Tilde: $\widetilde{u^{\!T}\Sigma w}$

Yields a **scalar** because $u^T\Sigma w$ is an inner product.

\[
\boxed{
\widetilde{u^{\!T}\Sigma w}
\,= \overline{\mathbf{u}}_0^{\!T} \left(\boldsymbol{\Sigma}_0 - \operatorname{diag}\boldsymbol{\mu}_0\right) \overline{\mathbf{w}}_0
\,+ \boldsymbol{\mu}_0^{\!T} \left[ (\overline{\mathbf{U}}\circ \mathbf{P}(t)) \circ (\overline{\mathbf{W}}\circ \mathbf{P}(t)) \right] \mathbf{1}
}
\]

Special cases:
- $\widetilde{\gamma^{\!T}\Sigma\gamma}$ → central part of predictive variance
- $\widetilde{\gamma^{\!T}\Sigma v}$ → relevant for TaylorIR expansions

---

## 6. Derivation Insights

### 6.1 How the Solution Works

The core difficulty is that interval-averaged current depends on entire trajectories. The solution:

1. **Introduce boundary states** $(i_0\to i_t)$ and compute prior mean/covariance of boundary counts by combining:
   - Prior over initial counts
   - Multinomial splitting induced by $\mathbf{P}(t)$

2. **Precompute conditional statistics** $\overline{\Gamma}_{i_0\to i_t}$ and $\operatorname{Var}(\overline{\Gamma}_{i_0\to i_t})$ using microscopic model (done once per kinetic scheme)

3. **Express macroscopic current** as linear functional of boundary counts plus noise, yielding closed-form expressions

4. **Map back to K-dimensional state space** using matrix identities, avoiding explicit $K^2 \times K^2$ storage

### 6.2 From Initial to Boundary Counts

Given $N_{i_0}(0)$, the number ending in $i_t$ is multinomial:

\[
\mathbb{E}[N_{i_0\to i_t}\mid N_{i_0}(0)] = N_{i_0}(0)\,P_{i_0\to i_t}(t)
\]
\[
\operatorname{Var}[N_{i_0\to i_t}\mid N_{i_0}(0)] = N_{i_0}(0)\,P_{i_0\to i_t}(t)(1-P_{i_0\to i_t}(t))
\]
\[
\operatorname{Cov}[N_{i_0\to i_t}, N_{i_0\to j_t}\mid N_{i_0}(0)] = -N_{i_0}(0)\,P_{i_0\to i_t}(t)\,P_{i_0\to j_t}(t)
\]

Applying law of total expectation/covariance yields the boundary mean/covariance formulas.

### 6.3 Collapsing the Boundary Quadratic Form

Naively, $\boldsymbol{\Sigma}^{\text{prior}}_{0\to t}$ is $K^2\times K^2$. Instead:

1. Insert explicit boundary covariance into $\widetilde{\gamma^{\top}\Sigma\gamma}$
2. Split into terms corresponding to covariance parts  
3. Recognize inner sums produce coarse-grained quantities $(\overline{\gamma}_0)_{i_0}$
4. After algebra, obtain the compact expressions in Section 5

---

## 7. MacroIR Interval Update Algorithm

### 7.1 Propagation (Markov Evolution)

**Mean propagation**:
\[
\boldsymbol{\mu}^{\text{prop}}(t) = \boldsymbol{\mu}_0\,\mathbf{P}(t).
\]

**Covariance propagation**:
\[
\boldsymbol{\Sigma}^{\text{prop}}(t) = \mathbf{P}(t)^{\!T}\,(\boldsymbol{\Sigma}_0 - \operatorname{diag}\boldsymbol{\mu}_0)\,\mathbf{P}(t) - \operatorname{diag}(\boldsymbol{\mu}^{\text{prop}}(t)).
\]

### 7.2 Measurement Update (Kalman-like Correction)

Define innovation:
\[
\delta = \overline{y}^{\text{obs}}_{0\to t} - \overline{y}^{\text{pred}}_{0\to t}.
\]

Define design vector and predictive variance:
\[
\mathrm{g} \equiv \widetilde{\gamma^{\!T}\Sigma} \in \mathbb{R}^K,
\qquad
\sigma^2 \equiv \sigma^2_{\overline{y}^{\text{pred}}}.
\]

**Mean update**:
\[
\boxed{
\boldsymbol{\mu}^{\text{post}}(t)
\,= \boldsymbol{\mu}^{\text{prop}}(t)
\,+ \frac{\delta}{\sigma^2}\,\mathrm{g}
}
\]

**Covariance update**:
\[
\boxed{
\boldsymbol{\Sigma}^{\text{post}}(t)
\,= \boldsymbol{\Sigma}^{\text{prop}}(t)
\,- \frac{1}{\sigma^2}\,\mathrm{g}\,\mathrm{g}^{\!T}
}
\]

This is a **rank-one correction**, analogous to scalar Kalman update.

---

## 8. Conceptual Summary

MacroIR resolves trajectory-dependent interval observations by expressing microscopic physics in **boundary-state matrices** indexed by $(i_0,i_t)$. The Bayesian update combines these with the macroscopic prior via a **single unified operator** — the tilde operator.

**Key principles**:
- Tilde is not "double" or "triple" — it's one operator whose expressions have different algebraic ranks
- $\widetilde{u^{\!T}\Sigma}$: vector | $\widetilde{u^{\!T}\Sigma w}$: scalar
- Just like matrix multiplication yields different tensor ranks depending on input shape

Embedding the tilde operator inside a scalar Kalman update gives MacroIR its power: full microscopic fidelity is preserved while computations remain entirely macroscopic.

The algorithm never constructs the full $K^2 \times K^2$ boundary covariance explicitly, making it computationally feasible while maintaining rigorous Bayesian updating.

