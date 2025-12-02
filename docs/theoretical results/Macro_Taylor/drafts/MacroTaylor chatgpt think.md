# Laplace / Taylor Approximation for the Multi-Channel Bandpass Model

We summarize the approximation, the resulting Gaussian posterior (covariance and mean), and the derivations, explicitly keeping the dependence on the number of channels \(N_{ch}\).

---

## 1. Approximation Setup

We observe a single current measurement \(y\) from a population of \(N_{ch}\) identical channels. The parameter vector is \(p \in \mathbb{R}^d\) (channel properties).

- **Mean current (per channel):** \(p \cdot \gamma\)
- **Noise variance contribution (per channel):** \(p \cdot \sigma^2\)
- **Measurement noise floor:** \(\epsilon^2\)

The **total mean** and **total variance** of the current scale with the number of channels:
\[
\begin{aligned}
\mathbb{E}[y \mid p] &= N_{ch} (p \cdot \gamma), \\
\mathrm{Var}(y \mid p) &= \epsilon^2 + N_{ch} (p \cdot \sigma^2).
\end{aligned}
\]

We use a **Gaussian prior** with precision proportional to \(N_{ch}\):
\[
p \sim \mathcal{N}\!\left(\mu,\; \left(N_{ch}\Sigma_p^{-1}\right)^{-1}\right) 
\;=\; \mathcal{N}\!\left(\mu,\; \frac{1}{N_{ch}}\Sigma_p\right).
\]

### Energy function

We define the (unnormalized) negative log-posterior **without the factor \(1/2\)**:
\[
E(p) 
\, = \underbrace{\frac{(y - N_{ch} p \cdot \gamma)^2}{\epsilon^2 + N_{ch} p \cdot \sigma^2}}_{M(p)}
\, + \underbrace{(p - \mu)^T (N_{ch}\Sigma_p^{-1}) (p - \mu)}_{P(p)} .
\]

Introduce:
\[
\begin{aligned}
\delta &:= y - N_{ch} (p \cdot \gamma), \\
V &:= \epsilon^2 + N_{ch} (p \cdot \sigma^2), \\
\mathbf{K} &:= N_{ch}\Sigma_p^{-1}.
\end{aligned}
\]

Then:
\[
E(p) = \frac{\delta^2}{V} + (p - \mu)^T \mathbf{K} (p - \mu).
\]

### Laplace / Taylor approximation

We approximate the posterior around an expansion point \(p_*\) (typically \(p_*=\mu\) or the previous iterate) by a second-order Taylor expansion:
\[
E(p) \approx E(p_*) + \mathbf{g}^T (p - p_*) + \tfrac{1}{2} (p - p_*)^T \mathbf{H} (p - p_*),
\]
where
\[
\mathbf{g} = \nabla E(p_*) ,\qquad \mathbf{H} = \nabla^2 E(p_*).
\]

Under this approximation, the posterior is Gaussian:
\[
p \mid y \;\approx\; \mathcal{N}(p_{\text{post}}, \Sigma_{\text{post}}),
\]
with
\[
\Sigma_{\text{post}} \approx 2\mathbf{H}^{-1}, 
\qquad
p_{\text{post}} \approx p_* - \mathbf{H}^{-1}\mathbf{g}.
\]

---

## 2. Result: Posterior Covariance and Posterior Mean

All quantities below are evaluated at the expansion point \(p_*\) (most naturally \(p_* = \mu\)).

### 2.1 Gradient and Hessian (recap)

- **Prior gradient and Hessian**
\[
\nabla P(p) = 2 N_{ch} \Sigma_p^{-1} (p - \mu), 
\qquad 
\nabla^2 P(p) = 2 N_{ch}\Sigma_p^{-1}.
\]

- **Likelihood gradient**

Using \(\delta = y - N_{ch} (p\cdot\gamma)\), \(V = \epsilon^2 + N_{ch} (p\cdot\sigma^2)\), and
\(\nabla \delta = -N_{ch}\gamma\), \(\nabla V = N_{ch}\sigma^2\):
\[
\nabla M(p) 
\, = -N_{ch} \left( \frac{2\delta}{V}\,\gamma + \frac{\delta^2}{V^2}\,\sigma^2 \right).
\]

- **Likelihood Hessian**
\[
\nabla^2 M(p) =
\frac{2 N_{ch}^2}{V} \gamma \gamma^T 
\, + \frac{2 N_{ch}^2 \delta}{V^2} \left( \gamma (\sigma^2)^T + \sigma^2 \gamma^T \right) 
\, + \frac{2 N_{ch}^2 \delta^2}{V^3} \sigma^2 (\sigma^2)^T.
\]

- **Total gradient and Hessian**
\[
\begin{aligned}
\mathbf{g} &= 2 N_{ch}\Sigma_p^{-1}(p - \mu)
           - N_{ch} \left( \frac{2\delta}{V}\gamma + \frac{\delta^2}{V^2}\sigma^2 \right), \\
\mathbf{H} &= 2 N_{ch}\Sigma_p^{-1}
\, + \frac{2 N_{ch}^2}{V} \gamma \gamma^T
\, + \frac{2 N_{ch}^2 \delta}{V^2} ( \gamma (\sigma^2)^T + \sigma^2 \gamma^T )
\, + \frac{2 N_{ch}^2 \delta^2}{V^3} \sigma^2 (\sigma^2)^T.
\end{aligned}
\]

### 2.2 Factorization of the Hessian

Factor out \(2 N_{ch}\):
\[
\mathbf{H} = 2 N_{ch}\mathbf{A},
\]
with
\[
\mathbf{A} = \Sigma_p^{-1}
\, + \frac{N_{ch}}{V}\gamma\gamma^T
\, + \frac{N_{ch}\delta}{V^2}(\gamma(\sigma^2)^T + \sigma^2\gamma^T)
\, + \frac{N_{ch}\delta^2}{V^3}\sigma^2(\sigma^2)^T.
\]

A key algebraic observation is that the last three terms form a **perfect square**:

Let
\[
\alpha := \frac{\delta}{V}, \qquad
\mathbf{v} := \gamma + \alpha \sigma^2.
\]
Then
\[
\gamma\gamma^T 
\, + \alpha(\gamma(\sigma^2)^T + \sigma^2\gamma^T)
\, + \alpha^2 \sigma^2(\sigma^2)^T
\, = (\gamma + \alpha\sigma^2)(\gamma + \alpha\sigma^2)^T
\, = \mathbf{v}\mathbf{v}^T.
\]

Hence
\[
\mathbf{A} = \Sigma_p^{-1} + \frac{N_{ch}}{V}\,\mathbf{v}\mathbf{v}^T.
\]

Define also
\[
\mathbf{u} := \sqrt{\frac{N_{ch}}{V}}\,\mathbf{v}, 
\quad\text{so that}\quad
\mathbf{A} = \Sigma_p^{-1} + \mathbf{u}\mathbf{u}^T.
\]

### 2.3 Posterior covariance via Sherman–Morrison

The *Sherman–Morrison* formula for a rank-1 update is
\[
(\mathbf{B} + \mathbf{u}\mathbf{u}^T)^{-1}
\, = \mathbf{B}^{-1} - \frac{\mathbf{B}^{-1}\mathbf{u}\mathbf{u}^T\mathbf{B}^{-1}}
                         {1 + \mathbf{u}^T \mathbf{B}^{-1}\mathbf{u}}.
\]

Taking \(\mathbf{B} = \Sigma_p^{-1}\) (so \(\mathbf{B}^{-1} = \Sigma_p\)), we obtain
\[
\mathbf{A}^{-1}
\, = \Sigma_p - \frac{\Sigma_p\mathbf{u}\mathbf{u}^T\Sigma_p}
                   {1 + \mathbf{u}^T\Sigma_p\mathbf{u}}.
\]

Now
\[
\mathbf{u} = \sqrt{\tfrac{N_{ch}}{V}}\,\mathbf{v}
\quad\Rightarrow\quad
\mathbf{u}\mathbf{u}^T = \frac{N_{ch}}{V}\mathbf{v}\mathbf{v}^T,
\]
and
\[
\mathbf{u}^T\Sigma_p\mathbf{u} 
\, = \frac{N_{ch}}{V}\,\mathbf{v}^T\Sigma_p\mathbf{v}
\;=\; \frac{N_{ch}}{V}\,s,
\]
where
\[
s := \mathbf{v}^T\Sigma_p\mathbf{v}.
\]

Thus
\[
\mathbf{A}^{-1} 
\, = \Sigma_p 
\, - \frac{\frac{N_{ch}}{V}\Sigma_p\mathbf{v}\mathbf{v}^T\Sigma_p}
         {1 + \frac{N_{ch}}{V}s}
\, = \Sigma_p - \frac{N_{ch}}{V + N_{ch}s}\,\Sigma_p\mathbf{v}\mathbf{v}^T\Sigma_p.
\]

Because \(\mathbf{H} = 2N_{ch}\mathbf{A}\), we have
\[
\Sigma_{\text{post}} \approx 2\mathbf{H}^{-1}
\, = 2(2N_{ch}\mathbf{A})^{-1}
\, = \frac{1}{N_{ch}} \mathbf{A}^{-1}.
\]

Therefore the **approximate posterior covariance** is:
\[
\boxed{
\Sigma_{\text{post}} 
\;\approx\; \frac{1}{N_{ch}}\Sigma_p 
\, - \frac{1}{V + N_{ch}s}\,\Sigma_p\mathbf{v}\mathbf{v}^T\Sigma_p
}
\]
with
\[
\mathbf{v} = \gamma + \frac{\delta}{V}\sigma^2,
\qquad
s = \mathbf{v}^T\Sigma_p\mathbf{v},
\qquad
\delta = y - N_{ch}(p_* \cdot \gamma),
\qquad
V = \epsilon^2 + N_{ch}(p_* \cdot \sigma^2).
\]

### 2.4 Posterior mean

In general, the Laplace mean is given by a **Newton step**:
\[
p_{\text{post}} \approx p_* - \mathbf{H}^{-1}\mathbf{g}.
\]

We already have \(\mathbf{H}^{-1} = \frac{1}{2N_{ch}}\mathbf{A}^{-1} = \frac{1}{2}\Sigma_{\text{post}}\), because \(\Sigma_{\text{post}} \approx 2\mathbf{H}^{-1}\).

For arbitrary \(p_*\),
\[
\mathbf{g} = 2 N_{ch}\Sigma_p^{-1}(p_* - \mu)
           - N_{ch} \frac{\delta}{V}(\gamma + \mathbf{v})
\]
(using \(\gamma + \mathbf{v} = 2\gamma + \frac{\delta}{V}\sigma^2\) to compress the likelihood term).

This yields a somewhat bulky general formula for \(p_{\text{post}}\) (which you already have in expanded form).

#### Special (and most useful) case: expansion around the prior mean

If we expand at the prior mean, \(p_* = \mu\), then
\[
p_* - \mu = 0 \quad\Rightarrow\quad
\mathbf{g} = -N_{ch} \frac{\delta}{V}(\gamma + \mathbf{v}).
\]

Using \(\mathbf{H}^{-1} = \tfrac12 \Sigma_{\text{post}}\):
\[
\begin{aligned}
p_{\text{post}}
&\approx \mu - \mathbf{H}^{-1}\mathbf{g} \\
&= \mu - \frac{1}{2}\Sigma_{\text{post}}\left( -N_{ch}\frac{\delta}{V}(\gamma + \mathbf{v})\right) \\
&= \mu + \frac{N_{ch}\delta}{2V}\,\Sigma_{\text{post}}(\gamma + \mathbf{v}).
\end{aligned}
\]

So the **posterior mean** for a Laplace expansion at the prior mean is:
\[
\boxed{
p_{\text{post}} \;\approx\;
\mu + \frac{N_{ch}\,\delta}{2V}\,\Sigma_{\text{post}}(\gamma + \mathbf{v}),
}
\]
with
\[
\delta = y - N_{ch}(\mu \cdot \gamma),
\qquad
V = \epsilon^2 + N_{ch}(\mu \cdot \sigma^2),
\qquad
\mathbf{v} = \gamma + \frac{\delta}{V}\sigma^2,
\]
and \(\Sigma_{\text{post}}\) as in the covariance formula above.

---

## 3. Hint of the Derivation (High-Level Path)

1. **Write the energy**  
   Start from
   \[
   E(p) = \frac{\delta^2}{V} + (p-\mu)^T N_{ch}\Sigma_p^{-1} (p-\mu),
   \]
   with \(\delta = y - N_{ch}(p\cdot\gamma)\), \(V = \epsilon^2 + N_{ch}(p\cdot\sigma^2)\).

2. **Compute \(\nabla E\) and \(\nabla^2 E\)**  
   Use \(\nabla\delta = -N_{ch}\gamma\), \(\nabla V = N_{ch}\sigma^2\), and the quotient rule on \(\delta^2/V\). This yields the gradient and Hessian given in §2.1.

3. **Factor the Hessian**  
   Separate a global factor \(2N_{ch}\) to define \(\mathbf{H} = 2N_{ch}\mathbf{A}\). Show that the likelihood contribution in \(\mathbf{A}\) is a structured quadratic form in \(\gamma\) and \(\sigma^2\).

4. **Recognize a perfect square / rank-1 update**  
   Observe that the three likelihood Hessian terms can be written as
   \[
   \frac{N_{ch}}{V}
   \bigl(\gamma + \alpha\sigma^2\bigr)\bigl(\gamma + \alpha\sigma^2\bigr)^T
   \]
   with \(\alpha = \delta/V\). This leads to
   \[
   \mathbf{A} = \Sigma_p^{-1} + \frac{N_{ch}}{V}\mathbf{v}\mathbf{v}^T,
   \quad \mathbf{v} = \gamma + \frac{\delta}{V}\sigma^2,
   \]
   i.e. \(\mathbf{A}\) is a rank-1 update of \(\Sigma_p^{-1}\).

5. **Apply Sherman–Morrison**  
   With \(\mathbf{A} = \Sigma_p^{-1} + \mathbf{u}\mathbf{u}^T\), \(\mathbf{u} = \sqrt{N_{ch}/V}\,\mathbf{v}\), use Sherman–Morrison to obtain a closed form for \(\mathbf{A}^{-1}\), hence for \(\Sigma_{\text{post}} \approx \frac{1}{N_{ch}}\mathbf{A}^{-1}\).

6. **Posterior mean via Newton step**  
   Use the Laplace relation
   \[
   p_{\text{post}} \approx p_* - \mathbf{H}^{-1}\mathbf{g},
   \quad
   \Sigma_{\text{post}} \approx 2\mathbf{H}^{-1}
   \Rightarrow \mathbf{H}^{-1} \approx \tfrac12\Sigma_{\text{post}}.
   \]
   Plug in the gradient. For \(p_*=\mu\), the prior term in \(\mathbf{g}\) drops out, and the update simplifies to the compact formula in §2.4.

---

## 4. Complete Derivation

Below is the derivation in more detail, step by step.

### 4.1 Gradients and Hessians

Recall
\[
E(p) = M(p) + P(p) = \frac{\delta^2}{V} + (p-\mu)^T N_{ch}\Sigma_p^{-1}(p-\mu)
\]
with
\[
\delta = y - N_{ch}(p\cdot\gamma), \quad
V = \epsilon^2 + N_{ch}(p\cdot\sigma^2).
\]

#### Prior term

\[
P(p) = (p-\mu)^T N_{ch}\Sigma_p^{-1}(p-\mu).
\]

- Gradient:
\[
\nabla P(p) 
\, = 2N_{ch}\Sigma_p^{-1}(p-\mu).
\]

- Hessian:
\[
\nabla^2 P(p) = 2N_{ch}\Sigma_p^{-1}.
\]

#### Likelihood term

\[
M(p) = \frac{\delta^2}{V}.
\]

We have
\[
\nabla\delta = -N_{ch}\gamma, \quad
\nabla V    = N_{ch}\sigma^2.
\]

Using the quotient rule in vector form:
\[
\nabla M(p) 
\, = \nabla\left(\frac{\delta^2}{V}\right)
\, = \frac{2\delta (\nabla\delta) V - \delta^2 (\nabla V)}{V^2}.
\]

Substitute:
\[
\nabla M(p) 
\, = \frac{2\delta(-N_{ch}\gamma)V - \delta^2 N_{ch}\sigma^2}{V^2}
\, = -N_{ch}\left(
     \frac{2\delta}{V}\gamma 
\, + \frac{\delta^2}{V^2}\sigma^2
  \right).
\]

Taking another gradient (coordinate-wise) gives the Hessian. It is convenient to write
\[
\nabla M(p) = -N_{ch}\left( a(p)\gamma + b(p)\sigma^2 \right),
\]
with
\[
a(p) = \frac{2\delta}{V}, 
\quad
b(p) = \frac{\delta^2}{V^2}.
\]

Then
\[
\nabla^2 M(p) 
\, = -N_{ch}\left( (\nabla a)\gamma^T + (\nabla b)(\sigma^2)^T\right).
\]

Compute
\[
\nabla a = \nabla\left(\frac{2\delta}{V}\right)
\, = 2\frac{\nabla\delta}{V} - 2\delta\frac{\nabla V}{V^2}
\, = 2\left( \frac{-N_{ch}\gamma}{V} - \frac{\delta N_{ch}\sigma^2}{V^2} \right),
\]
and
\[
\nabla b = \nabla\left(\frac{\delta^2}{V^2}\right)
\, = 2\frac{\delta\nabla\delta}{V^2} - 2\delta^2\frac{\nabla V}{V^3}
\, = 2\left( \frac{\delta(-N_{ch}\gamma)}{V^2} - \frac{\delta^2 N_{ch}\sigma^2}{V^3} \right).
\]

Insert into \(\nabla^2 M(p)\) and collect terms; after simplification one obtains
\[
\nabla^2 M(p)
\, = \frac{2 N_{ch}^2}{V} \gamma \gamma^T
\, + \frac{2 N_{ch}^2 \delta}{V^2} \left( \gamma (\sigma^2)^T + \sigma^2 \gamma^T \right)
\, + \frac{2 N_{ch}^2 \delta^2}{V^3} \sigma^2 (\sigma^2)^T,
\]
as used earlier.

#### Total gradient and Hessian

Summing prior and likelihood:
\[
\begin{aligned}
\mathbf{g} &= \nabla E(p) = \nabla M(p) + \nabla P(p), \\
\mathbf{H} &= \nabla^2 E(p) = \nabla^2 M(p) + \nabla^2 P(p),
\end{aligned}
\]
which matches the formulas in §2.1.

### 4.2 Factorization into \(\mathbf{H} = 2N_{ch}\mathbf{A}\)

Collect the Hessian terms:

\[
\begin{aligned}
\mathbf{H} &= 2 N_{ch}\Sigma_p^{-1}
\, + \frac{2 N_{ch}^2}{V} \gamma \gamma^T
\, + \frac{2 N_{ch}^2 \delta}{V^2} ( \gamma (\sigma^2)^T + \sigma^2 \gamma^T )
\, + \frac{2 N_{ch}^2 \delta^2}{V^3} \sigma^2 (\sigma^2)^T \\
&= 2N_{ch}\Big[
  \Sigma_p^{-1}
\, + \frac{N_{ch}}{V} \gamma\gamma^T
\, + \frac{N_{ch}\delta}{V^2} ( \gamma (\sigma^2)^T + \sigma^2 \gamma^T )
\, + \frac{N_{ch}\delta^2}{V^3} \sigma^2 (\sigma^2)^T
\Big].
\end{aligned}
\]

Define
\[
\mathbf{A} 
\, = \Sigma_p^{-1}
\, + \frac{N_{ch}}{V} \gamma\gamma^T
\, + \frac{N_{ch}\delta}{V^2} ( \gamma (\sigma^2)^T + \sigma^2 \gamma^T )
\, + \frac{N_{ch}\delta^2}{V^3} \sigma^2 (\sigma^2)^T,
\]
so that \(\mathbf{H} = 2N_{ch}\mathbf{A}\).

### 4.3 Rank-1 structure of \(\mathbf{A}\)

Let \(\alpha = \delta/V\), \(\mathbf{v} = \gamma + \alpha\sigma^2\). Then
\[
\mathbf{v}\mathbf{v}^T
\, = \gamma\gamma^T
\, + \alpha(\gamma(\sigma^2)^T + \sigma^2\gamma^T)
\, + \alpha^2\sigma^2(\sigma^2)^T.
\]

Therefore
\[
\begin{aligned}
\mathbf{A}
&= \Sigma_p^{-1} 
\, + \frac{N_{ch}}{V} \left[
    \gamma\gamma^T
\, + \alpha(\gamma(\sigma^2)^T + \sigma^2\gamma^T)
\, + \alpha^2\sigma^2(\sigma^2)^T
  \right] \\
&= \Sigma_p^{-1} + \frac{N_{ch}}{V}\,\mathbf{v}\mathbf{v}^T.
\end{aligned}
\]

Let
\[
\mathbf{u} := \sqrt{\frac{N_{ch}}{V}}\,\mathbf{v},
\]
so that
\[
\mathbf{A} = \Sigma_p^{-1} + \mathbf{u}\mathbf{u}^T,
\]
a **rank-1 update** of \(\Sigma_p^{-1}\).

### 4.4 Sherman–Morrison and \(\Sigma_{\text{post}}\)

Sherman–Morrison gives:
\[
\mathbf{A}^{-1} 
\, = (\Sigma_p^{-1} + \mathbf{u}\mathbf{u}^T)^{-1}
\, = \Sigma_p - \frac{\Sigma_p\mathbf{u}\mathbf{u}^T\Sigma_p}{1 + \mathbf{u}^T\Sigma_p\mathbf{u}}.
\]

Compute
\[
\mathbf{u}\mathbf{u}^T = \frac{N_{ch}}{V}\mathbf{v}\mathbf{v}^T,
\qquad
\mathbf{u}^T\Sigma_p\mathbf{u} = \frac{N_{ch}}{V} \mathbf{v}^T\Sigma_p\mathbf{v} = \frac{N_{ch}}{V}s.
\]

Thus
\[
\mathbf{A}^{-1}
\, = \Sigma_p 
\, - \frac{\frac{N_{ch}}{V}\Sigma_p\mathbf{v}\mathbf{v}^T\Sigma_p}{1 + \frac{N_{ch}}{V}s}
\, = \Sigma_p - \frac{N_{ch}}{V + N_{ch}s}\,\Sigma_p\mathbf{v}\mathbf{v}^T\Sigma_p.
\]

Since
\[
\Sigma_{\text{post}} \approx 2\mathbf{H}^{-1} = \frac{1}{N_{ch}}\mathbf{A}^{-1},
\]
we obtain
\[
\Sigma_{\text{post}} 
\approx \frac{1}{N_{ch}}\Sigma_p 
\, - \frac{1}{V + N_{ch}s}\,\Sigma_p\mathbf{v}\mathbf{v}^T\Sigma_p,
\]
exactly as stated in §2.3.

### 4.5 Posterior mean via Newton step

For the Laplace mean:
\[
p_{\text{post}} \approx p_* - \mathbf{H}^{-1}\mathbf{g}.
\]

We use \(\mathbf{H}^{-1} = \frac{1}{2N_{ch}}\mathbf{A}^{-1} = \frac{1}{2}\Sigma_{\text{post}}\).

The gradient at general \(p_*\) is
\[
\mathbf{g} = 2 N_{ch}\Sigma_p^{-1}(p_* - \mu)
           - N_{ch} \frac{\delta}{V}(\gamma + \mathbf{v}),
\]
where \(\delta, V, \mathbf{v}\) are evaluated at \(p_*\).

Then
\[
\begin{aligned}
p_{\text{post}}
&\approx p_* - \mathbf{H}^{-1}\mathbf{g} \\
&= p_* - \frac{1}{2}\Sigma_{\text{post}}
    \left[
       2 N_{ch}\Sigma_p^{-1}(p_* - \mu)
     - N_{ch} \frac{\delta}{V}(\gamma + \mathbf{v})
    \right].
\end{aligned}
\]

This is the general expression.  

For the important special case \(p_* = \mu\):

- The prior term vanishes: \(2 N_{ch}\Sigma_p^{-1}(p_* - \mu) = 0\).
- The gradient simplifies to 
  \(\mathbf{g} = -N_{ch} \frac{\delta}{V}(\gamma + \mathbf{v})\).

Thus
\[
\begin{aligned}
p_{\text{post}}
&\approx \mu - \frac{1}{2}\Sigma_{\text{post}}\left(-N_{ch} \frac{\delta}{V}(\gamma + \mathbf{v})\right) \\
&= \mu + \frac{N_{ch}\delta}{2V}\,\Sigma_{\text{post}}(\gamma + \mathbf{v}),
\end{aligned}
\]
which matches the compact form in §2.4.
