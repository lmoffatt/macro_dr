# Bayesian Approximation for Channel Noise Models

## 1. The Approximation
To derive closed-form expressions for the posterior statistics, we utilize a **Laplace Approximation** (Gaussian approximation) of the posterior distribution.

We approximate the Energy function $E(p)$ (the negative log-posterior) via a second-order Taylor expansion around the prior mean $\mu$. This implies:
1.  The posterior distribution is treated as a Multivariate Gaussian: $P(p|y) \approx \mathcal{N}(p_{post}, \Sigma_{post})$.
2.  The Posterior Mean $p_{post}$ is estimated using a single **Newton-Raphson** step starting from the prior mean $\mu$.
3.  The Posterior Covariance $\Sigma_{post}$ is approximated by the inverse of the Hessian of the Energy function (scaled by 2, given the definition of $E(p)$).

## 2. The Result
Given the data $y$, the number of channels $N_{ch}$, and prior parameters ($\mu, \Sigma_p$), the approximate posterior statistics are:

### Definitions
* **Residual:** $\delta = y - N_{ch} (\mu \cdot \gamma)$
* **Total Variance:** $V = \epsilon^2 + N_{ch} (\mu \cdot \sigma^2)$
* **Effective Observation Vector:** $\mathbf{v} = \gamma + \frac{\delta}{V} \sigma^2$
* **Projected Variance:** $s = \mathbf{v}^T \Sigma_p \mathbf{v}$

### Posterior Covariance ($\Sigma_{post}$)
$$
\Sigma_{post} \approx \frac{1}{N_{ch}} \Sigma_p - \frac{1}{V + N_{ch} s} \Sigma_p \mathbf{v} \mathbf{v}^T \Sigma_p
$$

### Posterior Mean ($p_{post}$)
$$
p_{post} \approx \mu + \frac{N_{ch}\delta}{2V} \Sigma_{post} (\gamma + \mathbf{v})
$$

---

## 3. Hint of the Derivation
The derivation follows these high-level steps:

1.  **Hessian Calculation:** We compute the Hessian $\mathbf{H}$ of the energy function $E(p)$. We observe that the likelihood contribution to the Hessian involves complex outer products of $\gamma$ and $\sigma^2$.
2.  **Factorization:** We recognize that the likelihood terms in the Hessian form a **perfect square**. This allows us to collapse the complex terms into a single rank-1 update using the vector $\mathbf{v} = \gamma + \frac{\delta}{V}\sigma^2$.
3.  **Sherman-Morrison Formula:** We define the posterior covariance as proportional to $\mathbf{H}^{-1}$. Since $\mathbf{H}$ is now expressed as an invertible matrix plus a rank-1 update, we apply the Sherman-Morrison formula to invert it analytically.
4.  **Newton-Raphson Update:** We calculate the posterior mean using $p_{post} = \mu - \mathbf{H}^{-1}\mathbf{g}$. By substituting the analytical form of $\mathbf{H}^{-1}$ and the gradient $\mathbf{g}$, the expression simplifies significantly into a form dependent on $\Sigma_{post}$.

---

## 4. Complete Derivation

### 4.1 The Energy Function, Gradient, and Hessian

We define the energy function $E(p)$:
$$E(p) = \frac{(y - N_{ch} p \cdot \gamma)^2}{V(p)} + (p - \mu)^T (N_{ch} \Sigma_p^{-1}) (p - \mu)$$
where $V(p) = \epsilon^2 + N_{ch} (p \cdot \sigma^2)$.

**The Gradient ($\mathbf{g}$):**
Evaluated at the prior mean $\mu$:
$$
\mathbf{g} = - N_{ch} \left( \frac{2\delta}{V} \gamma + \frac{\delta^2}{V^2} \sigma^2 \right) = -N_{ch} \frac{\delta}{V} \left( 2\gamma + \frac{\delta}{V}\sigma^2 \right)
$$
*(Note: The prior term vanishes because we evaluate at $p=\mu$.)*

**The Hessian ($\mathbf{H}$):**
$$
\mathbf{H} = 2 N_{ch} \Sigma_p^{-1} + \frac{2 N_{ch}^2}{V} \gamma \gamma^T + \frac{2 N_{ch}^2 \delta}{V^2} \left( \gamma (\sigma^2)^T + \sigma^2 \gamma^T \right) + \frac{2 N_{ch}^2 \delta^2}{V^3} \sigma^2 (\sigma^2)^T
$$

### 4.2 Deriving the Posterior Covariance

We seek $\Sigma_{post} \approx 2 \mathbf{H}^{-1}$. Let us first simplify $\mathbf{H}$.
Factor out $2 N_{ch}$:
$$
\mathbf{H} = 2 N_{ch} \left[ \Sigma_p^{-1} + \frac{N_{ch}}{V} \left( \gamma \gamma^T + \frac{\delta}{V} (\gamma (\sigma^2)^T + \sigma^2 \gamma^T) + \frac{\delta^2}{V^2} \sigma^2 (\sigma^2)^T \right) \right]
$$

**The Perfect Square Factorization:**
The term in the inner parenthesis is a perfect square expansion of the vector $\mathbf{v} = \gamma + \frac{\delta}{V}\sigma^2$:
$$
\mathbf{v}\mathbf{v}^T = \left(\gamma + \frac{\delta}{V}\sigma^2\right)\left(\gamma + \frac{\delta}{V}\sigma^2\right)^T = \gamma \gamma^T + \frac{\delta}{V} (\gamma (\sigma^2)^T + \sigma^2 \gamma^T) + \frac{\delta^2}{V^2} \sigma^2 (\sigma^2)^T
$$

Thus, the Hessian becomes a rank-1 update to the prior precision:
$$
\mathbf{H} = 2 N_{ch} \left[ \Sigma_p^{-1} + \frac{N_{ch}}{V} \mathbf{v}\mathbf{v}^T \right]
$$

**Inversion via Sherman-Morrison:**
We need $\Sigma_{post} = 2\mathbf{H}^{-1} = \frac{1}{N_{ch}} \left[ \Sigma_p^{-1} + \frac{N_{ch}}{V} \mathbf{v}\mathbf{v}^T \right]^{-1}$.
Applying the Sherman-Morrison formula $(A + uv^T)^{-1} = A^{-1} - \frac{A^{-1}uv^TA^{-1}}{1 + v^TA^{-1}u}$:

$$
\left[ \Sigma_p^{-1} + \frac{N_{ch}}{V} \mathbf{v}\mathbf{v}^T \right]^{-1} = \Sigma_p - \frac{\Sigma_p \left(\frac{N_{ch}}{V} \mathbf{v}\mathbf{v}^T\right) \Sigma_p}{1 + \frac{N_{ch}}{V} \mathbf{v}^T \Sigma_p \mathbf{v}}
$$

Let $s = \mathbf{v}^T \Sigma_p \mathbf{v}$. Simplify the fraction:
$$
= \Sigma_p - \frac{\frac{N_{ch}}{V} \Sigma_p \mathbf{v}\mathbf{v}^T \Sigma_p}{\frac{V + N_{ch}s}{V}} = \Sigma_p - \frac{N_{ch}}{V + N_{ch}s} \Sigma_p \mathbf{v}\mathbf{v}^T \Sigma_p
$$

Finally, substitute this back into the expression for $\Sigma_{post}$:
$$
\Sigma_{post} = \frac{1}{N_{ch}} \left( \Sigma_p - \frac{N_{ch}}{V + N_{ch}s} \Sigma_p \mathbf{v}\mathbf{v}^T \Sigma_p \right)
$$
$$
\boxed{\Sigma_{post} = \frac{1}{N_{ch}} \Sigma_p - \frac{1}{V + N_{ch} s} \Sigma_p \mathbf{v} \mathbf{v}^T \Sigma_p}
$$

### 4.3 Deriving the Posterior Mean

We use the Newton-Raphson update step starting from $\mu$:
$$p_{post} = \mu - \mathbf{H}^{-1} \mathbf{g}$$

From the gradient calculation, we rewrite $\mathbf{g}$ using $\mathbf{v}$:
$$
\mathbf{g} = -N_{ch} \frac{\delta}{V} \left( \gamma + \left( \gamma + \frac{\delta}{V}\sigma^2 \right) \right) = -N_{ch} \frac{\delta}{V} (\gamma + \mathbf{v})
$$

Substitute $\mathbf{H}^{-1} = \frac{1}{2} \Sigma_{post}$ and $\mathbf{g}$:
$$
p_{post} = \mu - \left( \frac{1}{2} \Sigma_{post} \right) \left( -N_{ch} \frac{\delta}{V} (\gamma + \mathbf{v}) \right)
$$

Rearranging the scalars yields the final compact result:
$$
\boxed{p_{post} = \mu + \frac{N_{ch}\delta}{2V} \Sigma_{post} (\gamma + \mathbf{v})}
$$