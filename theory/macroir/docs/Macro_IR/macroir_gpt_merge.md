
# MacroIR Interval Update

## Unified Boundary-State & Tilde Operator Spec


## 0. Scope and conventions

We consider an ensemble of $N_{\text{ch}}$ independent Markov channels with $K$ microscopic states. Over an interval $[0,t]$ we observe a scalar interval-averaged macroscopic current $\overline{y}_{0\to t}^{\text{obs}}$. We want the posterior mean and covariance of the macroscopic state at time $t$, and the predictive mean/variance of the interval current.

**Convention.** In what follows:

* $(\boldsymbol{\mu}_0, \boldsymbol{\Sigma}_0)$ describe the **per-channel** occupancy statistics at the start of the interval:

  * $\mu_{0,i} = \mathbb{E}[\text{fraction of channels in state } i]$,
  * $\sum_i \mu_{0,i} = 1$.
* Macroscopic means and variances per ion channel, so they both increase with $N_{\text{ch}}$, the number of channels.

This matches the usage of the explicit $N_{\text{ch}}$ factors in the predictive mean/variance formulas in Moffatt 2007.

You can convert to raw-counts as ${\boldsymbol{n}}_0 = N_{\text{ch}}\boldsymbol{\mu}_0$ if needed; all formulas simply scale accordingly.

---

## 1. Indices, shapes, and operators

* $K$: number of microscopic states.

* Indices:

  * $i_0, j_0$: start (time $0$) states.
  * $i_t, j_t$: end (time $t$) states.
  * $a,b$: generic state indices at time $t$.

* state Vectors are **columns**:

  * $\boldsymbol{\mu}_0 \in \mathbb{R}^K$, $\boldsymbol{\mu}^{\text{prior}}(t)\in\mathbb{R}^K$, etc.


* state properties Vectors are **columns**:

  * $\boldsymbol{\gamma}_0 \in \mathbb{R}^K$, $\boldsymbol{\sigma^2}_0(t)\in\mathbb{R}^K$, etc.

  
* state covariance Matrices are indexed by two states (i_a, j_a) at the same time a:

  * $\boldsymbol{\Sigma}_0, \boldsymbol{\Sigma}^{\text{prop}}_{t}\in\mathbb{R}^{K\times K}$ (symmetric).

* state transition and properties of Matrices are indexed by the state at start time $0$, $i_0$ and at the end time $t$ $i_t$    
  * $\mathbf{P}(t)\in\mathbb{R}^{K\times K}$ with entries $P_{i_0\to i_t}$ (asymmetric;  $\sum_{i_t} P_{i_0\to i_t}=1$  sum of rows is 1)

  * mean and variance of mean conductance  matrices $\overline{\mathbf{\Gamma}},\overline{\mathbf{V}}\in\mathbb{R}^{K\times K}$ indexed by $(i_0 \to i_t)$.

* Elementwise (Hadamard) product: $A\circ B$.

* $\operatorname{diag}(x)$: diagonal matrix from vector $x$.

* $\operatorname{diag}(A)$: vector of diagonal entries of matrix $A$.

* $\mathbf{1}$: all-ones column vector.

---

## 2. Inputs and precomputed microscopic objects

For each interval $[0,t]$ we need:

### 2.1 Macroscopic prior

* $\boldsymbol{\mu}_0$
* $\boldsymbol{\Sigma}_0$

### 2.2 Markov dynamics

* Generator $\mathbf{Q}$ (offline)
* Transition matrix:
  $$
  \mathbf{P}(t) = e^{\mathbf{Q}t},\qquad
  P_{i_0\to i_t}(t) = [\mathbf{P}(t)]_{i_0 i_t}.
  $$

### 2.3 Boundary-conditioned current statistics

* Mean interval current:
  $$
  \overline{\Gamma}_{i_0\to i_t}
  \, =
  \mathbb{E}[\frac{1}{t}{\int_0 ^t y(\tau) d\tau}\mid i_0,i_t].
  $$
* Variance contribution:
  $$
  \overline{V}_{i_0\to i_t}
  \, =
  \operatorname{Var}(\frac{1}{t}{\int_0 ^t y(\tau) d\tau}\mid i_0,i_t).
  $$

### 2.4 Measurement setup

* Number of channels $N_{\text{ch}}$
* Instrument/binning noise:
  $$
  \epsilon^2_{0\to t} = \frac{\epsilon^2}{t} + \nu^2
  $$
* Observation $\overline{y}^{\text{obs}}_{0\to t}$

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
## 3. Boundary-lifted current statistics

Define
$$
\mathbf{G} = \overline{\mathbf{\Gamma}} \circ \mathbf{P}(t),\qquad
G_{i_0 i_t} = \overline{\Gamma}_{i_0\to i_t} P_{i_0\to i_t}(t).
$$

### 3.1 Start-conditioned mean interval current

$$
(\overline{\gamma}_0)_{i_0}
\,=
\sum_{i_t} G_{i_0 i_t}.
$$

Matrix form:
$$
\overline{\boldsymbol{\gamma}}_0 = \mathbf{G}\mathbf{1}.
$$

### 3.2 Start-conditioned intrinsic variance

Define $\mathbf{V}=\overline{\mathbf{V}}\circ\mathbf{P}(t)$, then
$$
(\sigma^2_{\overline{\gamma}_0})_{i_0}
\,=
\sum_{i_t} V_{i_0 i_t}.
$$

Matrix form:
$$
\sigma^2_{\overline{\gamma}_0} = \mathbf{V}\mathbf{1}.
$$

---

## 4. Predictive mean of the interval current

$$
\boxed{
\overline{y}^{\text{pred}}_{0\to t}
\,=
N_{\text{ch}}\,\boldsymbol{\mu}_0 \overline{\boldsymbol{\gamma}}_0
}
$$

---

## 5. The tilde operator

A unified linear operator coupling prior moments, $\mathbf{P}(t)$, and boundary matrices.

### 5.1 Tilde over mean state : $\widetilde{\boldsymbol{\mu}_{0}}$

$$
\boxed{
\widetilde{\boldsymbol{\mu}_{0}}
\,=
\boldsymbol{\mu}_0\, \mathbf{P}(t)
}
$$

### 5.2 Tilde over state Covariance: $\widetilde{\boldsymbol{\Sigma}_{0}}$

$$
\boxed{
\widetilde{\boldsymbol{\Sigma}}
\,=
\mathbf{P}(t)^\top
(\boldsymbol{\Sigma}_0-\operatorname{diag}(\boldsymbol{\mu}_0))\,
\mathbf{P}(t)
\, -
\operatorname{diag}(\boldsymbol{\mu}_0\mathbf\,{P}(t))
}
$$


### 5.3 tilde over bilinear product: $\widetilde{u^\top \Sigma w}$

$$
\boxed{
\widetilde{u^\top\Sigma w}
\,=
\overline{\boldsymbol{u}}_0^\top
(\boldsymbol{\Sigma}_0-\operatorname{diag}\boldsymbol{\mu}_0)
\overline{\boldsymbol{w}}_0
\,+
\boldsymbol{\mu}_0^\top
[(\overline{\mathbf{U}}\circ\mathbf{P})\circ(\overline{\mathbf{W}}\circ\mathbf{P})]\mathbf{1}
}
$$

Special case:
$$
\boxed{
\widetilde{\gamma^\top\Sigma\gamma}
\,=
\overline{\boldsymbol{\gamma}}_0^\top
(\boldsymbol{\Sigma}_0-\operatorname{diag}\boldsymbol{\mu}_0)
\overline{\boldsymbol{\gamma}}_0
+
\boldsymbol{\mu}_0^\top[\mathbf{G}\circ\mathbf{G}]\mathbf{1}
}
$$

### 5.4 Vector tilde: $\widetilde{u^\top\Sigma}$

$$
\boxed{
\widetilde{u^\top\Sigma}
\,=
\mathbf{P}(t)^\top
(\boldsymbol{\Sigma}_0-\operatorname{diag}\boldsymbol{\mu}_0)
\overline{\boldsymbol{u}}_0 +
(\overline{\mathbf{U}}\circ\mathbf{P}(t))^\top \boldsymbol{\mu}_0
}
$$

For the interval current:
$$
\boxed{
\mathbf{g}
\,=
 \widetilde{\gamma^\top\Sigma}=
\mathbf{P}(t)^\top
(\boldsymbol{\Sigma}_0-\operatorname{diag}\boldsymbol{\mu}_0)
\overline{\boldsymbol{\gamma}}_0
\, +
\mathbf{G}^\top\boldsymbol{\mu}_0
}
$$

---

## 6. Predictive variance of the interval current

$$
\boxed{
\sigma^2_{\overline{y}^{\text{pred}}}
\,=
\epsilon^2_{0\to t}
+
N_{\text{ch}}\,\widetilde{\gamma^\top\Sigma\gamma}
+
N_{\text{ch}}\sum_{i_0}\mu_{0,i_0}
(\sigma^2_{\overline{\gamma}_0})_{i_0}
}
$$

---

## 7. Propagation to time $t$

### 7.1 Mean

$$
\boxed{
\boldsymbol{\mu}^{\text{prop}}(t)
\,=
\widetilde{\boldsymbol{\mu}}=
 \boldsymbol{\mu}_0\, \mathbf{P}(t)
}
$$

### 7.2 Covariance

$$
\boxed{
\boldsymbol{\Sigma}^{\text{prop}}(t)
\,=
\mathbf{P}(t)^\top
(\boldsymbol{\Sigma}_0-\operatorname{diag}\boldsymbol{\mu}_0)
\mathbf{P}(t)
\, -
\operatorname{diag}(\boldsymbol{\mu}^{\text{prior}}(t))
}
$$

---

## 8. Measurement update

Let
$$
\delta = \overline{y}^{\text{obs}}_{0\to t} -\overline{y}^{\text{pred}}_{0\to t}.
$$

### 8.1 Mean update

$$
\boxed{
\boldsymbol{\mu}^{\text{post}}(t)
\,=
\boldsymbol{\mu}^{\text{prop}}(t)
\, +
\frac{\mathbf{g}\,\delta}{\sigma^2_{\overline{y}^{\text{pred}}}}
}
$$

### 8.2 Covariance update

$$
\boxed{
\boldsymbol{\Sigma}^{\text{post}}(t)
\,=
\boldsymbol{\Sigma}^{\text{prop}}(t)-
\frac{\mathbf{g}\mathbf{g}^\top}{\sigma^2_{\overline{y}^{\text{pred}}}}
}
$$

---

## 9. Summary of workflow

1. Inputs: $\boldsymbol{\mu}_0,\boldsymbol{\Sigma}_0, \mathbf{P}(t), \overline{\mathbf{\Gamma}}\, \overline{\mathbf{V}}, N_{\text{ch}}, \epsilon^2, \nu^2, \overline{y}^{\text{obs}}_{0\to t}$
2. Compute:

   * $\mathbf{G}=\overline{\mathbf{\Gamma}}\circ\mathbf{P}(t)$
   * $\overline{\boldsymbol{\gamma}}_0=\mathbf{G}\mathbf{1}$
   * $\mathbf{V}=\overline{\mathbf{V}}\circ\mathbf{P}(t)$
   * $\sigma^2_{\overline{\gamma}_0}=\mathbf{V}\mathbf{1}$
3. Predictive current:
   $$
   \overline{y}^{\text{pred}}_{0\to t}
   \,=
   N_{\text{ch}}\,\boldsymbol{\mu}_0\,\overline{\boldsymbol{\gamma}}_0
   $$
   $$
   \sigma^2_{\overline{y}^{\text{pred}}}
   \,=
   \epsilon^2_{0\to t}
   +
   N_{\text{ch}}\widetilde{\gamma^\top\Sigma\gamma}
   +
   N_{\text{ch}}\boldsymbol{\mu}_0\sigma^2_{\overline{\gamma}_0}
   $$
4. Cross-covariance: $\mathbf{g}=\widetilde{\gamma^\top\Sigma}$
5. Propagate:

   * $\boldsymbol{\mu}^{\text{prior}}(t)=\mathbf{P}(t)^\top\boldsymbol{\mu}_0$
   * $\boldsymbol{\Sigma}^{\text{prop}}(t)=\mathbf{P}(t)^\top(\boldsymbol{\Sigma}_0-\operatorname{diag}\boldsymbol{\mu}_0)\mathbf{P}(t)-\operatorname{diag}(\boldsymbol{\mu}^{\text{prior}}(t))$
6. Update:

   * $\boldsymbol{\mu}^{\text{post}}(t)=\boldsymbol{\mu}^{\text{prior}}(t)+\mathbf{g}\delta/\sigma^2_{\overline{y}^{\text{pred}}}$
   * $\boldsymbol{\Sigma}^{\text{post}}(t)=\boldsymbol{\Sigma}^{\text{prop}}(t)-\mathbf{g}\mathbf{g}^\top/\sigma^2_{\overline{y}^{\text{pred}}}$
7. Use $(\boldsymbol{\mu}^{\text{post}}(t),\boldsymbol{\Sigma}^{\text{post}}(t))$ as next-interval prior.



All steps use only $K$- and $K\times K$-objects; no $K^2\times K^2$ boundary covariance is ever instantiated.

## 10. Numerical notes and invariants

* $\boldsymbol{\Sigma}_0$ and $\boldsymbol{\Sigma}^{\text{prop}}(t)$ should be symmetric; small asymmetries from floating point should be symmetrized explicitly.
* $\sigma^2_{\overline{y}^{\text{pred}}}$ must be strictly positive; in practice, impose a lower floor to avoid division by zero in the update.
* $\boldsymbol{\Sigma}^{\text{post}}(t)$ is obtained by a rank-1 subtraction; numerical PSD violations can be mitigated by:

  * enforcing symmetry,
  * clipping tiny negative eigenvalues, if necessary, in post-processing.
* Complexity:

  * Matrix-matrix multiplies: (O(K^3)) (dominated by covariance propagation).
  * Hadamard products and tilde expressions: (O(K^2)).

# Appendix A — Pedagogical narrative (optional)

This appendix compresses the “story” from both documents into a conceptual picture; it is **not** needed for implementation.

### A.1 Why boundary states?

The interval-averaged current of a channel depends on its entire trajectory over ([0,t]), not just its start or end state. Directly integrating over all trajectories is intractable.

MacroIR’s trick:

1. Introduce **boundary states** $(i_0,i_t)$.
2. For each boundary pair, precompute:

   * mean interval current $\overline{\Gamma}_{i_0\to i_t}$,
   * variance $\overline{V}_{i_0\to i_t}$.
3. Note that, given initial counts, channels from each $i_0$ split multinomially into end states with probabilities $P_{i_0\to i_t}(t)$. The random boundary counts (N_{i_0\to i_t}) therefore contain *all* the trajectory information that matters for the interval current.

The macroscopic interval current can then be written as a linear functional of the boundary counts plus additive noise:
$$
\overline{y}_{0\to t}
\approx
\sum_{i_0,i_t}
\overline{\Gamma}_{i_0\to i_t}N_{i_0\to i_t}
\, +
\text{noise}.
$$

### A.2 Why the tilde operator?

Naively, the boundary count covariance is a (K^2\times K^2) object. We never want to build it. Instead, observe:

* The prior information lives in the **state space** at $t=0$: $(\boldsymbol{\mu}_0,\boldsymbol{\Sigma}_0)$.
* The microscopic physics over $[0,t]$ is fully captured by:

  * the transition matrix $\mathbf{P}(t)$,
  * the boundary tables $\overline{\mathbf{\Gamma}},\overline{\mathbf{V}}$.

The tilde operator is the algebraic mechanism that:

1. **Lifts** a state-space vector (like $\overline{\boldsymbol{\gamma}}_0)$ into the boundary space via the boundary matrices.
2. **Modulates** it with the transition probabilities $\mathbf{P}(t)$ and the initial covariance $\boldsymbol{\Sigma}_0$.
3. **Collapses** back to either:

   * a scalar $\widetilde{u^{\top}\Sigma w}$, or
   * a vector $\widetilde{u^{\top}\Sigma}$,

without ever explicitly constructing the full boundary covariance.

Conceptually:

* $\widetilde{\gamma^{\top}\Sigma\gamma}$ is the scalar you’d get by forming the boundary covariance and applying the quadratic form defined by the boundary current weights.
* $\widetilde{\gamma^{\top}\Sigma}$ is the cross-covariance between the macroscopic state at time $t$ and the interval current.

The lift–modulate–collapse pattern is the same in both cases; only the **rank** of the underlying algebraic expression differs, just like $v^{\top} A$ (vector) vs. $v^{\top}Aw$ (scalar) in ordinary matrix algebra.

### A.3 Why the update looks like a scalar Kalman filter

Once you have:

* a prior Gaussian state at time $t$: $(\boldsymbol{\mu}^{\text{prior}}(t),\boldsymbol{\Sigma}^{\text{prop}}(t))$,
* a scalar observation $\overline{y}^{\text{obs}}_{0\to t}$,
* its predictive mean/variance $\overline{y}^{\text{pred}}_{0\to t},\sigma^2_{\overline{y}^{\text{pred}}}$,
* a cross-covariance vector $\mathbf{g}$,

you are exactly in the setting of a 1D Kalman update in a $K$-dimensional state space:

* $\mathbf{g}$ plays the role of the “measurement vector times prior covariance”.
* $\sigma^2_{\overline{y}^{\text{pred}}}$ is the total effective measurement variance.
* The update
  $$
  \boldsymbol{\mu}^{\text{post}}
  \, =
  \boldsymbol{\mu}^{\text{prior}}
  \, +
  \frac{\mathbf{g}}{\sigma^2_{\overline{y}^{\text{pred}}}}\delta,
  \qquad
  \boldsymbol{\Sigma}^{\text{post}}
  \,=\boldsymbol{\Sigma}^{\text{prop}}
  \frac{\mathbf{g}\mathbf{g}^{\top}}{\sigma^2_{\overline{y}^{\text{pred}}}}
  $$
  is precisely the Gaussian conditioning formula for a joint variable $(\boldsymbol{N_{ch}}(t),\overline{y}_{0\to t})$.

MacroIR’s distinctive feature is not the Kalman structure itself, but that $\mathbf{g}$ and $\sigma^2_{\overline{y}^{\text{pred}}}$ are computed from **interval physics** via boundary matrices and the unified tilde operator, rather than from a naive instant-state observation model.

