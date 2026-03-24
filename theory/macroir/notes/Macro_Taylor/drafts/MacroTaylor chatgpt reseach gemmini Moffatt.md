
# Bayesian Inference for Ion Channel Kinetics: The Second-Order Approximation

**Nomenclature Convention:**
Following the convention of *Moffatt (2007)*, we define the parameter distribution statistics using **scaled covariance matrices** that are independent of the channel count $N$.
* **Prior:** $P(\mathbf{r}) = \mathcal{N}\left(\mathbf{r}; \boldsymbol{\mu}, \frac{1}{N}\boldsymbol{\Sigma}\right)$
* **Posterior:** $P(\mathbf{r}|y) \approx \mathcal{N}\left(\mathbf{r}; \mathbf{r}_{post}, \frac{1}{N}\boldsymbol{\Sigma}_{post}\right)$

---

## 1. The Result: Recursive Update Equations

Given a macroscopic current observation $y$ and a prior state $(\boldsymbol{\mu}, \boldsymbol{\Sigma})$, the Gaussian posterior is obtained via the following updates.

### 1.1 Definitions at the Expansion Point
We evaluate the following at the prior mean $\mathbf{r} = \boldsymbol{\mu}$:

* **Residual:** $\delta = y - N (\boldsymbol{\mu} \cdot \boldsymbol{\gamma})$
* **Total Variance:** $V = \varepsilon^2 + N (\boldsymbol{\mu} \cdot \boldsymbol{\sigma}^2)$
* **Effective Direction:** $\mathbf{v} = \boldsymbol{\gamma} + \frac{\delta}{V}\boldsymbol{\sigma}^2$
* **Projected Prior Variance:** $s = \mathbf{v}^T \boldsymbol{\Sigma} \mathbf{v}$

### 1.2 Posterior Scaled Covariance ($\boldsymbol{\Sigma}_{post}$)
The covariance update is a rank-1 reduction of the prior covariance:
$$
\boxed{\boldsymbol{\Sigma}_{post} = \boldsymbol{\Sigma} - \frac{N}{V + Ns} \boldsymbol{\Sigma}\mathbf{v}\mathbf{v}^T\boldsymbol{\Sigma}}
$$

### 1.3 Posterior Mean ($\mathbf{r}_{post}$)
The mean update moves along the transformed effective direction:
$$
\boxed{\mathbf{r}_{post} = \boldsymbol{\mu} + \frac{\delta}{2V} \boldsymbol{\Sigma}_{post}(\boldsymbol{\gamma} + \mathbf{v})}
$$

---

## 2. Detailed Derivation

### 2.1 The Energy Function
We define the energy function $E(\mathbf{r})$ as the negative log-posterior.
$$
P(\mathbf{r}|y) \propto \exp\left( -E(\mathbf{r}) \right)
$$
$$
E(\mathbf{r}) = \underbrace{\frac{(y - N\mathbf{r} \cdot \boldsymbol{\gamma})^2}{2(\varepsilon^2 + N\mathbf{r} \cdot \boldsymbol{\sigma}^2)}}_{Likelihood \, M(\mathbf{r})} + \underbrace{\frac{N}{2}(\mathbf{r} - \boldsymbol{\mu})^T \boldsymbol{\Sigma}^{-1} (\mathbf{r} - \boldsymbol{\mu})}_{Prior \, P(\mathbf{r})}
$$

### 2.2 Gradient and Hessian
We approximate $E(\mathbf{r})$ via a second-order Taylor expansion around $\boldsymbol{\mu}$.

**The Gradient ($\mathbf{g}$):**
Using the quotient rule on the likelihood term $M(\mathbf{r})$, and noting that $\nabla \text{Prior}$ vanishes at $\mathbf{r}=\boldsymbol{\mu}$:
$$
\mathbf{g} = \nabla E(\boldsymbol{\mu}) = - \frac{N\delta}{2V} \left( 2\boldsymbol{\gamma} + \frac{\delta}{V}\boldsymbol{\sigma}^2 \right)
$$
Using the definition $\mathbf{v} = \boldsymbol{\gamma} + \frac{\delta}{V}\boldsymbol{\sigma}^2$, we simplify this to:
$$
\mathbf{g} = -\frac{N\delta}{2V}(\boldsymbol{\gamma} + \mathbf{v})
$$

**The Hessian ($\mathbf{H}$):**
Differentiation of the gradient terms yields the Hessian. We observe that the prior contributes $N\boldsymbol{\Sigma}^{-1}$, and the likelihood contributes a structured matrix involving $\boldsymbol{\gamma}$ and $\boldsymbol{\sigma}^2$:
$$
\mathbf{H} = N\boldsymbol{\Sigma}^{-1} + \frac{N^2}{V}\boldsymbol{\gamma}\boldsymbol{\gamma}^T + \frac{N^2 \delta}{V^2} \left( \boldsymbol{\gamma}(\boldsymbol{\sigma}^2)^T + \boldsymbol{\sigma}^2\boldsymbol{\gamma}^T \right) + \frac{N^2 \delta^2}{V^3} \boldsymbol{\sigma}^2(\boldsymbol{\sigma}^2)^T
$$

### 2.3 Factorization (The "Perfect Square" Insight)
The complex likelihood terms in the Hessian form a perfect square expansion of the vector $\mathbf{v}$.
Recall $\mathbf{v}\mathbf{v}^T = (\boldsymbol{\gamma} + \frac{\delta}{V}\boldsymbol{\sigma}^2)(\boldsymbol{\gamma} + \frac{\delta}{V}\boldsymbol{\sigma}^2)^T$. Expanding this product reproduces exactly the likelihood terms in $\mathbf{H}$.

Thus, the Hessian simplifies to a rank-1 update of the prior precision:
$$
\mathbf{H} = N \left[ \boldsymbol{\Sigma}^{-1} + \frac{N}{V} \mathbf{v}\mathbf{v}^T \right]
$$

### 2.4 Deriving the Posterior Scaled Covariance
By definition of the Gaussian approximation, the actual posterior covariance is $\mathbf{H}^{-1}$.
In our scaling convention, the actual covariance is $\frac{1}{N}\boldsymbol{\Sigma}_{post}$. Therefore:
$$
\frac{1}{N}\boldsymbol{\Sigma}_{post} = \mathbf{H}^{-1} = \frac{1}{N} \left[ \boldsymbol{\Sigma}^{-1} + \frac{N}{V} \mathbf{v}\mathbf{v}^T \right]^{-1}
$$
Multiplying by $N$, we seek to invert the bracketed term:
$$
\boldsymbol{\Sigma}_{post} = \left[ \boldsymbol{\Sigma}^{-1} + \mathbf{u}\mathbf{u}^T \right]^{-1} \quad \text{where } \mathbf{u} = \sqrt{\frac{N}{V}}\mathbf{v}
$$

**Applying Sherman-Morrison:**
Using $(A + uv^T)^{-1} = A^{-1} - \frac{A^{-1}uv^TA^{-1}}{1 + v^TA^{-1}u}$:
$$
\boldsymbol{\Sigma}_{post} = \boldsymbol{\Sigma} - \frac{\boldsymbol{\Sigma} \left( \frac{N}{V}\mathbf{v}\mathbf{v}^T \right) \boldsymbol{\Sigma}}{1 + \frac{N}{V}\mathbf{v}^T \boldsymbol{\Sigma} \mathbf{v}}
$$
Substituting $s = \mathbf{v}^T \boldsymbol{\Sigma} \mathbf{v}$ and multiplying numerator/denominator by $V$:
$$
\boldsymbol{\Sigma}_{post} = \boldsymbol{\Sigma} - \frac{N}{V + Ns} \boldsymbol{\Sigma}\mathbf{v}\mathbf{v}^T\boldsymbol{\Sigma}
$$

### 2.5 Deriving the Posterior Mean
We perform a single Newton-Raphson step from $\boldsymbol{\mu}$:
$$
\mathbf{r}_{post} = \boldsymbol{\mu} - \mathbf{H}^{-1}\mathbf{g}
$$
Substitute $\mathbf{H}^{-1} = \frac{1}{N}\boldsymbol{\Sigma}_{post}$ and $\mathbf{g} = -\frac{N\delta}{2V}(\boldsymbol{\gamma} + \mathbf{v})$:
$$
\mathbf{r}_{post} = \boldsymbol{\mu} - \left( \frac{1}{N}\boldsymbol{\Sigma}_{post} \right) \left( -\frac{N\delta}{2V}(\boldsymbol{\gamma} + \mathbf{v}) \right)
$$
The $N$ and $1/N$ cancel, leaving:
$$
\mathbf{r}_{post} = \boldsymbol{\mu} + \frac{\delta}{2V} \boldsymbol{\Sigma}_{post}(\boldsymbol{\gamma} + \mathbf{v})
$$

---

## 3. Comparison with Moffatt (2007)
This derivation is a generalization of the equations found in *Moffatt (2007)*.
If we set state-dependent noise to zero ($\boldsymbol{\sigma}^2 = 0$), then:
1.  $\mathbf{v} \to \boldsymbol{\gamma}$
2.  $V \to \varepsilon^2$
3.  The term $(\boldsymbol{\gamma} + \mathbf{v}) \to 2\boldsymbol{\gamma}$

The mean update becomes:
$$
\mathbf{r}_{post} = \boldsymbol{\mu} + \frac{\delta}{2\varepsilon^2}\boldsymbol{\Sigma}_{post}(2\boldsymbol{\gamma}) = \boldsymbol{\mu} + \frac{\delta}{\varepsilon^2}\boldsymbol{\Sigma}_{post}\boldsymbol{\gamma}
$$
This matches Eq. 50 of the 2007 paper exactly.