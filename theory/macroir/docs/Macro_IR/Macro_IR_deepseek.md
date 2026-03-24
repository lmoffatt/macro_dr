# MacroIR Interval Update: Boundary-State Formulation and Implementation Guide

This document explains how to implement the MacroIR interval update using the boundary-state representation. It has three parts: (1) explicit problem formulation and equations, (2) a short “how it was solved” sketch, and (3) a tutorial-style derivation.

## 1. Problem formulation and core equations

We consider an ensemble of identical Markov channels with K microscopic states. Over a time interval \([0,t]\) we observe the interval-averaged macroscopic current \(\overline{y}_{0\to t}\). We want to update the Bayesian belief over the channel state-occupancy vector at time \(t\) given this observation.

### 1.1 State variables and prior

Let  
\[
\boldsymbol{\mu}^{\text{prior}}_0 \in \mathbb{R}^K,\qquad
\boldsymbol{\Sigma}^{\text{prior}}_0 \in \mathbb{R}^{K\times K}
\]
be the prior mean and covariance of the channel-count vector at the start of the interval. Component-wise,
\[
(\boldsymbol{\mu}^{\text{prior}}_0)_i = \mathbb{E}[N_i(0)],\qquad
(\boldsymbol{\Sigma}^{\text{prior}}_0)_{ij}
= \operatorname{Cov}(N_i(0), N_j(0)).
\]

The channel state follows a continuous-time Markov process with generator \(\mathbf{Q}\), so the transition matrix over \([0,t]\) is  
\[
\mathbf{P}(t) = \exp(\mathbf{Q}t),\qquad
P_{i_0\to i_t}(t) = \bigl[\mathbf{P}(t)\bigr]_{i_0 i_t}.
\]

We denote the number of channels by \(N_{\text{ch}}\).

### 1.2 Boundary-state random variables

For a single channel, a “boundary state” is the pair of start and end states \((i_0,i_t)\). For the whole ensemble, define the count
\[
N_{i_0\to i_t} = \text{number of channels that start in } i_0 \text{ and end in } i_t.
\]

The **boundary-state mean** (per ensemble, not normalized) is
\[
\mu^{\text{prior}}_{0\to t,(i_0\to i_t)}
= \mathbb{E}[N_{i_0\to i_t}]
= \mathbb{E}[N_{i_0}(0)]\,P_{i_0\to i_t}(t)
= (\boldsymbol{\mu}^{\text{prior}}_0)_{i_0}\,P_{i_0\to i_t}(t).
\]

The **boundary-state covariance** is
\[
\Sigma^{\text{prior}}_{0\to t,(i_0\to i_t)(j_0\to j_t)} = \operatorname{Cov}\bigl(N_{i_0\to i_t}, N_{j_0\to j_t}\bigr).
\]

Using the law of total covariance, conditioning on the initial counts \(N_{i_0}(0)\) and accounting for the multinomial splitting from each start state, one obtains
\[
\Sigma^{\text{prior}}_{0\to t,(i_0\to i_t)(j_0\to j_t)} = \\ P_{i_0\to i_t}(t)
\bigl[(\boldsymbol{\Sigma}^{\text{prior}}_0)_{i_0 j_0} - \delta_{i_0 j_0}(\boldsymbol{\mu}^{\text{prior}}_0)_{i_0}\bigr]
P_{j_0\to j_t}(t) - \delta_{i_0 j_0}\delta_{i_t j_t}
(\boldsymbol{\mu}^{\text{prior}}_0)_{i_0}P_{i_0\to i_t}(t).
\]

The term \(\boldsymbol{\Sigma}^{\text{prior}}_0 - \mathrm{diag}(\boldsymbol{\mu}^{\text{prior}}_0)\) comes from covariances between initial counts and the negative correlations of the multinomial split within each start state; the term with \(\delta_{i_0 j_0}\delta_{i_t j_t}\) comes from the variance of the multinomial counts when all four indices coincide.

### 1.3 Interval-averaged current and measurement model

For each boundary pair \((i_0\to i_t)\), let
\[
\overline{\Gamma}_{i_0\to i_t}
\]
denote the **conditional mean interval-averaged current** of a *single* channel over \([0,t]\), given that it starts in state \(i_0\) at \(0\) and ends in state \(i_t\) at \(t\).

Define the mean interval current from a given start state \(i\) as
\[
(\overline{\gamma}_0)_i
= \sum_{j} P_{i\to j}(t)\,\overline{\Gamma}_{i\to j}.
\]

Then the mean macroscopic interval-averaged current predicted by the prior is
\[
\overline{y}^{\text{pred}}_{0\to t}
= N_{\text{ch}}\;\boldsymbol{\mu}^{\text{prior}}_0\cdot \overline{\boldsymbol{\gamma}}_0
= N_{\text{ch}}\sum_i (\boldsymbol{\mu}^{\text{prior}}_0)_i\,(\overline{\gamma}_0)_i
= N_{\text{ch}}\sum_i (\boldsymbol{\mu}^{\text{prior}}_0)_i
\sum_j P_{i\to j}(t)\,\overline{\Gamma}_{i\to j}.
\]

The measurement model assumes an effective scalar noise variance
\[
\epsilon^2_{0\to t} = \frac{\epsilon^2}{t} + \nu^2,
\]
combining instrumental noise and additional interval-averaging noise.

The full predictive variance of \(\overline{y}_{0\to t}\) has three pieces:
1. Instrumental and binning noise \(\epsilon^2_{0\to t}\).
2. Uncertainty in the initial ensemble state, propagated through the deterministic function \(\overline{\boldsymbol{\gamma}}_0\).
3. Intrinsic stochastic variability of the interval current given the start and end states.

The second contribution can be written as
\[
N_{\text{ch}}\;\widetilde{\gamma^{\top}\Sigma\gamma},
\]
where
\[
\widetilde{\gamma^{\top}\Sigma\gamma}
= \boldsymbol{\gamma}_{0\to t}^{\top}\,
\boldsymbol{\Sigma}^{\text{prior}}_{0\to t}\,
\boldsymbol{\gamma}_{0\to t}
\]
is a scalar quadratic form over boundary states, with
\[
(\boldsymbol{\gamma}_{0\to t})_{(i_0\to i_t)}
= \overline{\Gamma}_{i_0\to i_t}.
\]

The third contribution is
\[
N_{\text{ch}}\sum_i (\boldsymbol{\mu}^{\text{prior}}_0)_i\,
(\sigma^2_{\overline{\gamma}_0})_i,
\quad
(\sigma^2_{\overline{\gamma}_0})_i
= \sum_j P_{i\to j}(t)\,\operatorname{Var}(\overline{\Gamma}_{i\to j}).
\]

Putting everything together, the predictive variance is
\[
\sigma^2_{\overline{y}^{\text{pred}}_{0\to t}}
= \epsilon^2_{0\to t} + N_{\text{ch}}\;\widetilde{\gamma^{\top}\Sigma\gamma} + N_{\text{ch}}\sum_i (\boldsymbol{\mu}^{\text{prior}}_0)_i
\sum_j P_{i\to j}(t)\,\operatorname{Var}(\overline{\Gamma}_{i\to j}).
\]

### 1.4 Key tilde quantities

The central “tilde” objects for the update are:

Scalar:
\[
\widetilde{\gamma^{\top}\Sigma\gamma}
= \boldsymbol{\gamma}_{0\to t}^{\top}\,
\boldsymbol{\Sigma}^{\text{prior}}_{0\to t}\,
\boldsymbol{\gamma}_{0\to t}.
\]

Vector indexed by end state \(i_t\):
\[
\bigl(\widetilde{\gamma^{\top}\Sigma}\bigr)_{i_t} = \sum_{i_0}
\bigl(\boldsymbol{\gamma}_{0\to t}^{\top}
\boldsymbol{\Sigma}^{\text{prior}}_{0\to t}\bigr)_{(i_0\to i_t)}.
\]

These can be expressed directly in terms of \(\boldsymbol{\mu}^{\text{prior}}_0\), \(\boldsymbol{\Sigma}^{\text{prior}}_0\), \(\mathbf{P}(t)\) and \(\overline{\Gamma}_{i_0\to i_t}\) without explicitly storing the full \(K^2\times K^2\) boundary covariance.

1. Explicit scalar expression:
\[
\widetilde{\gamma^{\top}\Sigma\gamma} = \overline{\boldsymbol{\gamma}}_0^{\top}
\bigl(\boldsymbol{\Sigma}^{\text{prior}}_0 - \mathrm{diag}(\boldsymbol{\mu}^{\text{prior}}_0)\bigr)
\overline{\boldsymbol{\gamma}}_0 + \boldsymbol{\mu}^{\text{prior}}_0\cdot
\Bigl[
\bigl(\overline{\mathbf{\Gamma}}_{0\to t}\circ\mathbf{P}(t)\bigr)
\circ
\bigl(\overline{\mathbf{\Gamma}}_{0\to t}\circ\mathbf{P}(t)\bigr)
\Bigr]\mathbf{1},
\]
where \(\overline{\mathbf{\Gamma}}_{0\to t}\) is the matrix with entries
\[
(\overline{\mathbf{\Gamma}}_{0\to t})_{i_0 i_t}
= \overline{\Gamma}_{i_0\to i_t},
\]
\(\circ\) denotes Hadamard (elementwise) product, and \(\mathbf{1}\) is the all-ones vector over the end-state index.

2. Explicit vector expression:
\[
\bigl(\widetilde{\gamma^{\top}\Sigma}\bigr)_{i_t} = \sum_{j_0}(\overline{\gamma}_0)_{j_0}
\sum_{i_0}\bigl[(\boldsymbol{\Sigma}^{\text{prior}}_0)_{j_0 i_0} - \delta_{j_0 i_0}(\boldsymbol{\mu}^{\text{prior}}_0)_{j_0}\bigr]
P_{i_0\to i_t}(t) + \sum_{i_0}(\boldsymbol{\mu}^{\text{prior}}_0)_{i_0}
\overline{\Gamma}_{i_0\to i_t}P_{i_0\to i_t}(t).
\]

In matrix form,

\[
\widetilde{\gamma^{\top}\Sigma}
= \overline{\boldsymbol{\gamma}}_0^{\top}
\bigl(\boldsymbol{\Sigma}^{\text{prior}}_0
      - \mathrm{diag}(\boldsymbol{\mu}^{\text{prior}}_0)\bigr)
\mathbf{P}(t)
\;+\;
\Bigl(\boldsymbol{\mu}^{\text{prior}}_0 \circ 
      \operatorname{diag}(\overline{\mathbf{\Gamma}}_{0\to t} \circ \mathbf{P}(t))
\Bigr)^{\top}.
\]

where the second term is just another way to write the last sum component-wise.

### 1.5 Covariance propagation and update

The prior mean at time \(t\) is
\[
\boldsymbol{\mu}^{\text{prior}}(t) = \boldsymbol{\mu}^{\text{prior}}_0\,\mathbf{P}(t),
\qquad
\mu^{\text{prior}}_a(t)
= \sum_{i_0}(\boldsymbol{\mu}^{\text{prior}}_0)_{i_0}
P_{i_0\to a}(t).
\]

The **propagated covariance** (without conditioning on the measurement) is
\[
\boldsymbol{\Sigma}^{\text{prior}}_{\text{prop}}(t) = \mathbf{P}(t)^{\top}
\bigl(\boldsymbol{\Sigma}^{\text{prior}}_0 - \mathrm{diag}(\boldsymbol{\mu}^{\text{prior}}_0)\bigr)\mathbf{P}(t) + \mathrm{diag}(\boldsymbol{\mu}^{\text{prior}}(t)) - \mathrm{diag}(\boldsymbol{\mu}^{\text{prior}}(t)),
\]
which simplifies to
\[
\bigl[\mathbf{P}(t)^{\top}
\bigl(\boldsymbol{\Sigma}^{\text{prior}}_0 - \mathrm{diag}(\boldsymbol{\mu}^{\text{prior}}_0)\bigr)\mathbf{P}(t)\bigr]_{ab} - \delta_{ab}\,\mu^{\text{prior}}_a(t).
\]

In index notation,
\[
\Sigma^{\text{prior}}_{\text{prop},ab}(t) = \sum_{i_0,j_0}
P_{i_0\to a}(t)
\bigl[(\boldsymbol{\Sigma}^{\text{prior}}_0)_{i_0 j_0} - \delta_{i_0 j_0}(\boldsymbol{\mu}^{\text{prior}}_0)_{i_0}\bigr]
P_{j_0\to b}(t) - \delta_{ab}\sum_{i_0}(\boldsymbol{\mu}^{\text{prior}}_0)_{i_0}
P_{i_0\to a}(t).
\]

The **Bayesian correction** uses the scalar predictive variance
\(\sigma^2_{\overline{y}^{\text{pred}}_{0\to t}}\) and the vector
\(\widetilde{\gamma^{\top}\Sigma}\). The covariance update is
\[
\boldsymbol{\Sigma}^{\text{prior}}(t) = \boldsymbol{\Sigma}^{\text{prior}}_{\text{prop}}(t) - \frac{1}{\sigma^2_{\overline{y}^{\text{pred}}_{0\to t}}}
\bigl(\widetilde{\gamma^{\top}\Sigma}\bigr)^{\top}
\bigl(\widetilde{\gamma^{\top}\Sigma}\bigr),
\]
or component-wise
\[
\Sigma^{\text{prior}}_{ab}(t) = \Sigma^{\text{prior}}_{\text{prop},ab}(t) - \frac{1}{\sigma^2_{\overline{y}^{\text{pred}}_{0\to t}}}
\bigl(\widetilde{\gamma^{\top}\Sigma}\bigr)_a
\bigl(\widetilde{\gamma^{\top}\Sigma}\bigr)_b.
\]

(The corresponding mean update is the usual scalar-Kalman update in the state space using \(\widetilde{\gamma^{\top}\Sigma}\) and the innovation
\(\overline{y}^{\text{obs}}_{0\to t}-\overline{y}^{\text{pred}}_{0\to t}\).)

## 2. Conceptual “hint” of the solution

The core difficulty is that the interval-averaged current depends on the entire trajectory of each channel over \([0,t]\), not just the states at the endpoints. Directly integrating over all trajectories is intractable.

The trick is:

1. Introduce **boundary states** \((i_0\to i_t)\) and compute the prior mean and covariance of the boundary counts \(N_{i_0\to i_t}\) by combining  
   a) the prior over initial counts, and  
   b) the multinomial splitting induced by \(\mathbf{P}(t)\).
2. For each boundary state, precompute the **conditional interval current statistics** \(\overline{\Gamma}_{i_0\to i_t}\) and \(\operatorname{Var}(\overline{\Gamma}_{i_0\to i_t})\) using the microscopic model (this part is done once for each kinetic scheme and time step).
3. Express the macroscopic interval-averaged current as a **linear functional of the boundary counts**, plus additive noise. This produces closed-form expressions for  
   \(\overline{y}^{\text{pred}}_{0\to t}\),  
   \(\sigma^2_{\overline{y}^{\text{pred}}_{0\to t}}\),  
   and the tilde quantities \(\widetilde{\gamma^{\top}\Sigma\gamma}\) and \(\widetilde{\gamma^{\top}\Sigma}\).
4. Finally, map everything back to the **K-dimensional state space** using matrix identities, so the algorithm never has to store K²-dimensional objects. The update becomes a scalar-Kalman correction in the macroscopic state space with a carefully constructed “measurement vector” given by \(\widetilde{\gamma^{\top}\Sigma}\).

This is why the final implementation only needs \(\boldsymbol{\mu}_0^{\text{prior}}\), \(\boldsymbol{\Sigma}_0^{\text{prior}}\), \(\mathbf{P}(t)\), and the tables \(\overline{\Gamma}_{i_0\to i_t}\), \(\operatorname{Var}(\overline{\Gamma}_{i_0\to i_t})\).

## 3. Tutorial-style derivation

### 3.1 From initial counts to boundary counts

Start from the random vector of initial counts \(N_i(0)\), \(i=1,\dots,K\), with known mean \(\boldsymbol{\mu}_0^{\text{prior}}\) and covariance \(\boldsymbol{\Sigma}_0^{\text{prior}}\).

Given \(N_{i_0}(0)\), the number of channels that end in any state \(i_t\) after time \(t\) is multinomial with probabilities \(P_{i_0\to i_t}(t)\). Therefore,
\[
\mathbb{E}[N_{i_0\to i_t}\mid N_{i_0}(0)]
= N_{i_0}(0)\,P_{i_0\to i_t}(t),
\]
\[
\operatorname{Var}[N_{i_0\to i_t}\mid N_{i_0}(0)]
= N_{i_0}(0)\,P_{i_0\to i_t}(t)\bigl(1-P_{i_0\to i_t}(t)\bigr),
\]
\[
\operatorname{Cov}[N_{i_0\to i_t}, N_{i_0\to j_t}\mid N_{i_0}(0)]
= -N_{i_0}(0)\,P_{i_0\to i_t}(t)\,P_{i_0\to j_t}(t)
\quad\text{for }i_t\neq j_t.
\]

Channels that start in different states \(i_0\neq j_0\) are conditionally independent given the initial counts. Applying the law of total expectation and covariance yields the boundary mean and covariance formulas in section 1.2.

### 3.2 From boundary counts to interval current

The total interval-averaged current from all channels can be written as
\[
\overline{y}_{0\to t} = \frac{1}{t}\int_0^t I_{\text{macro}}(s)\,ds
\approx \sum_{i_0,i_t} \overline{\Gamma}_{i_0\to i_t} N_{i_0\to i_t} + \text{noise}.
\]

By construction,
\[
\mathbb{E}\bigl[\overline{y}_{0\to t}\mid \{N_{i_0\to i_t}\}\bigr] = \sum_{i_0,i_t} \overline{\Gamma}_{i_0\to i_t} N_{i_0\to i_t},
\]
so the unconditional mean is
\[
\mathbb{E}[\overline{y}_{0\to t}] = \sum_{i_0,i_t} \overline{\Gamma}_{i_0\to i_t}
\,\mathbb{E}[N_{i_0\to i_t}],
\]
which simplifies to the expression for \(\overline{y}^{\text{pred}}_{0\to t}\) using the boundary means and the definition of \(\overline{\boldsymbol{\gamma}}_0\).

For the variance, write
\[
\operatorname{Var}(\overline{y}_{0\to t}) = \underbrace{\operatorname{Var}\bigl(\mathbb{E}[\overline{y}_{0\to t}
\mid \{N_{i_0\to i_t}\}]\bigr)}_{\text{variance due to boundary counts}} + \underbrace{\mathbb{E}\bigl[\operatorname{Var}(\overline{y}_{0\to t}
\mid \{N_{i_0\to i_t}\})\bigr]}_{\text{intrinsic interval variance}} + \epsilon^2_{0\to t}.
\]

The first term is precisely the quadratic form \(\widetilde{\gamma^{\top}\Sigma\gamma}\) times \(N_{\text{ch}}\) (because the mean current contributed by a given boundary configuration is linear in \(N_{i_0\to i_t}\)). The second term is obtained by summing the intrinsic variances \(\operatorname{Var}(\overline{\Gamma}_{i_0\to i_t})\) weighted by the expected boundary counts, which yields the third term in the predictive variance formula.

### 3.3 Collapsing the boundary quadratic form

Naively, \(\boldsymbol{\Sigma}^{\text{prior}}_{0\to t}\) is a \(K^2\times K^2\) matrix. Directly constructing it would be infeasible for large K. Instead, we exploit its structure:

1. Insert the explicit boundary covariance into
   \[
   \widetilde{\gamma^{\top}\Sigma\gamma}
   = \sum_{i_0,i_t}\sum_{j_0,j_t}
   \overline{\Gamma}_{i_0\to i_t}\,
   \Sigma^{\text{prior}}_{0\to t,(i_0\to i_t)(j_0\to j_t)}\,
   \overline{\Gamma}_{j_0\to j_t}.
   \]
2. Split the result into two terms \(T_1\) and \(T_2\) corresponding to the two parts of the boundary covariance.
3. Recognize that the inner sums over \(i_t\) and \(j_t\) produce the coarse-grained quantities
   \[
   (\overline{\gamma}_0)_{i_0}
   = \sum_{i_t}P_{i_0\to i_t}(t)\,\overline{\Gamma}_{i_0\to i_t},
   \quad
   (\overline{\gamma}_0)_{j_0}
   = \sum_{j_t}P_{j_0\to j_t}(t)\,\overline{\Gamma}_{j_0\to j_t}.
   \]
4. After some algebra, the first term becomes
   \[
   T_1    = \overline{\boldsymbol{\gamma}}_0^{\top}
   \bigl(\boldsymbol{\Sigma}^{\text{prior}}_0 - \mathrm{diag}(\boldsymbol{\mu}^{\text{prior}}_0)\bigr)
   \overline{\boldsymbol{\gamma}}_0,
   \]
   and the second term becomes
   \[
   T_2    = \boldsymbol{\mu}^{\text{prior}}_0\cdot
   \Bigl[
   \bigl(\overline{\mathbf{\Gamma}}_{0\to t}\circ\mathbf{P}(t)\bigr)
   \circ
   \bigl(\overline{\mathbf{\Gamma}}_{0\to t}\circ\mathbf{P}(t)\bigr)
   \Bigr]\mathbf{1}.
   \]

Exactly the same strategy applied to
\(\widetilde{\gamma^{\top}\Sigma}\) yields its vector form in section 1.4.

### 3.4 Covariance update as a scalar-Kalman correction

At this point, the interval measurement \(\overline{y}_{0\to t}\) is a scalar linear observation of the random vector \(\boldsymbol{N}(t)\) (or equivalently, of the macroscopic occupancy vector), with effective “measurement vector” proportional to \(\widetilde{\gamma^{\top}\Sigma}\) and scalar variance \(\sigma^2_{\overline{y}^{\text{pred}}_{0\to t}}\).

The propagation step \(\boldsymbol{\mu}^{\text{prior}}(t)\), \(\boldsymbol{\Sigma}^{\text{prior}}_{\text{prop}}(t)\) is purely Markovian (no observation). The conditioning on the scalar \(\overline{y}^{\text{obs}}_{0\to t}\) is then identical in structure to a one-dimensional Kalman update in the K-dimensional state space. This directly gives the rank-1 covariance correction
\[
\boldsymbol{\Sigma}^{\text{prior}}(t) = \boldsymbol{\Sigma}^{\text{prior}}_{\text{prop}}(t) - \frac{1}{\sigma^2_{\overline{y}^{\text{pred}}_{0\to t}}}
\bigl(\widetilde{\gamma^{\top}\Sigma}\bigr)^{\top}
\bigl(\widetilde{\gamma^{\top}\Sigma}\bigr),
\]
and the corresponding mean update (not expanded here) in the standard Kalman form.

### 3.5 Implementation recipe (per interval)

For implementation, one time step \([0,t]\) can be coded as:

1. Given \(\boldsymbol{\mu}^{\text{prior}}_0\), \(\boldsymbol{\Sigma}^{\text{prior}}_0\), precomputed \(\mathbf{P}(t)\), \(\overline{\Gamma}_{i_0\to i_t}\), and \(\operatorname{Var}(\overline{\Gamma}_{i_0\to i_t})\), construct  
   a) \(\overline{\boldsymbol{\gamma}}_0\) by summing over end states,  
   b) the Hadamard products \(\overline{\mathbf{\Gamma}}_{0\to t}\circ\mathbf{P}(t)\), etc.
2. Compute \(\overline{y}^{\text{pred}}_{0\to t}\) and \(\sigma^2_{\overline{y}^{\text{pred}}_{0\to t}}\) using the formulas in sections 1.3 and 1.4.
3. Compute the vector \(\widetilde{\gamma^{\top}\Sigma}\) using the matrix or index form in section 1.4.
4. Propagate \(\boldsymbol{\mu}^{\text{prior}}_0\) and \(\boldsymbol{\Sigma}^{\text{prior}}_0\) to time \(t\) via \(\mathbf{P}(t)\) to obtain \(\boldsymbol{\mu}^{\text{prior}}(t)\) and \(\boldsymbol{\Sigma}^{\text{prior}}_{\text{prop}}(t)\).
5. Apply the scalar-Kalman correction to \(\boldsymbol{\Sigma}^{\text{prior}}_{\text{prop}}(t)\) using \(\widetilde{\gamma^{\top}\Sigma}\) and \(\sigma^2_{\overline{y}^{\text{pred}}_{0\to t}}\); similarly update the mean \(\boldsymbol{\mu}^{\text{post}}(t)\) using the innovation term.

This completes the closed-form interval update needed for MacroIR without ever constructing the full boundary-state covariance explicitly.
