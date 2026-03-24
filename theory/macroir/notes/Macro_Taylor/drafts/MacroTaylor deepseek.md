# Bayesian Inference for Ion Channel Kinetics with State-Dependent Noise

## 1. Approximation and Scaling Convention

We employ the **Laplace approximation** (second-order Taylor expansion) of the negative log-posterior around the prior mean. Following Moffatt (2007), we adopt the scaling convention where:

- **Prior distribution**: $P(\mathbf{r}) = \mathcal{N}\left(\mathbf{r}; \boldsymbol{\mu}, \frac{1}{N}\boldsymbol{\Sigma}\right)$
- **Posterior distribution**: $P(\mathbf{r}|y) = \mathcal{N}\left(\mathbf{r}; \mathbf{r}_{post}, \frac{1}{N}\boldsymbol{\Sigma}_{post}\right)$

Here $\boldsymbol{\Sigma}$ and $\boldsymbol{\Sigma}_{post}$ represent the **scaled covariance matrices** that are independent of $N$, with the actual covariance being $\frac{1}{N}\boldsymbol{\Sigma}$.

## 2. Results

### Posterior Scaled Covariance
$$\boxed{\boldsymbol{\Sigma}_{post} = \boldsymbol{\Sigma} - \frac{N}{V + Ns} \boldsymbol{\Sigma}\mathbf{v}\mathbf{v}^T\boldsymbol{\Sigma}}$$

### Posterior Mean
$$\boxed{\mathbf{r}_{post} = \boldsymbol{\mu} + \frac{N\delta}{2V} \left(\frac{1}{N}\boldsymbol{\Sigma}_{post}\right)(\boldsymbol{\gamma} + \mathbf{v}) = \boldsymbol{\mu} + \frac{\delta}{2V} \boldsymbol{\Sigma}_{post}(\boldsymbol{\gamma} + \mathbf{v})}$$

### Definitions
- **Residual**: $\delta = y - N\boldsymbol{\mu} \cdot \boldsymbol{\gamma}$
- **Total variance**: $V = \varepsilon^2 + N\boldsymbol{\mu} \cdot \boldsymbol{\sigma}^2$  
- **Effective direction**: $\mathbf{v} = \boldsymbol{\gamma} + \frac{\delta}{V}\boldsymbol{\sigma}^2$
- **Prior variance of effective observable**: $s = \mathbf{v}^T\boldsymbol{\Sigma}\mathbf{v}$

Where:
- $y$: observed macroscopic current
- $\boldsymbol{\gamma}$: conductance vector per state
- $\boldsymbol{\sigma}^2$: state-dependent noise contributions
- $\varepsilon^2$: instrumental noise variance
- $\boldsymbol{\mu}, \boldsymbol{\Sigma}$: prior mean and scaled covariance of state proportions
- $N$: number of channels

## 3. Derivation

### 3.1 Problem Setup

We begin with the macroscopic current measurement model:

$$P(y|\mathbf{r}) = \mathcal{N}\left(y; N\mathbf{r} \cdot \boldsymbol{\gamma}, \varepsilon^2 + N\mathbf{r} \cdot \boldsymbol{\sigma}^2\right)$$
$$P(\mathbf{r}) = \mathcal{N}\left(\mathbf{r}; \boldsymbol{\mu}, \frac{1}{N}\boldsymbol{\Sigma}\right)$$

The **energy function** (negative log-posterior) is:

$$E(\mathbf{r}) = \frac{(y - N\mathbf{r} \cdot \boldsymbol{\gamma})^2}{2(\varepsilon^2 + N\mathbf{r} \cdot \boldsymbol{\sigma}^2)} + \frac{N}{2}(\mathbf{r} - \boldsymbol{\mu})^T \boldsymbol{\Sigma}^{-1} (\mathbf{r} - \boldsymbol{\mu})$$

### 3.2 Gradient and Hessian at Prior Mean

At $\mathbf{r} = \boldsymbol{\mu}$, define:
- $\delta = y - N\boldsymbol{\mu}\cdot\boldsymbol{\gamma}$
- $V = \varepsilon^2 + N\boldsymbol{\mu}\cdot\boldsymbol{\sigma}^2$
- $\mathbf{v} = \boldsymbol{\gamma} + \frac{\delta}{V}\boldsymbol{\sigma}^2$

**Gradient**:
$$\mathbf{g} = -\frac{N\delta}{V}\boldsymbol{\gamma} + \frac{N\delta^2}{2V^2}\boldsymbol{\sigma}^2 = -\frac{N\delta}{2V}(\boldsymbol{\gamma} + \mathbf{v})$$

**Hessian**:
$$\mathbf{H} = N\boldsymbol{\Sigma}^{-1} + \frac{N^2}{V}\mathbf{v}\mathbf{v}^T$$

### 3.3 Posterior Covariance

The actual posterior covariance is $\mathbf{H}^{-1}$, and by our scaling convention:

$$\frac{1}{N}\boldsymbol{\Sigma}_{post} = \mathbf{H}^{-1} = \left(N\boldsymbol{\Sigma}^{-1} + \frac{N^2}{V}\mathbf{v}\mathbf{v}^T\right)^{-1}$$

Using Sherman-Morrison:

$$\frac{1}{N}\boldsymbol{\Sigma}_{post} = \frac{1}{N}\left(\boldsymbol{\Sigma} - \frac{N}{V + Ns} \boldsymbol{\Sigma}\mathbf{v}\mathbf{v}^T\boldsymbol{\Sigma}\right)$$

Multiplying by $N$:

$$\boldsymbol{\Sigma}_{post} = \boldsymbol{\Sigma} - \frac{N}{V + Ns} \boldsymbol{\Sigma}\mathbf{v}\mathbf{v}^T\boldsymbol{\Sigma}$$

### 3.4 Posterior Mean

Using Newton-Raphson:

$$\mathbf{r}_{post} = \boldsymbol{\mu} - \mathbf{H}^{-1}\mathbf{g} = \boldsymbol{\mu} + \frac{1}{N}\boldsymbol{\Sigma}_{post} \cdot \frac{N\delta}{2V}(\boldsymbol{\gamma} + \mathbf{v})$$

Simplifying:

$$\mathbf{r}_{post} = \boldsymbol{\mu} + \frac{\delta}{2V} \boldsymbol{\Sigma}_{post}(\boldsymbol{\gamma} + \mathbf{v})$$

## 4. Special Case: $\boldsymbol{\sigma}^2 = 0$ (Moffatt 2007)

When state-dependent noise is zero:
- $\mathbf{v} = \boldsymbol{\gamma}$
- $V = \varepsilon^2$
- $s = \boldsymbol{\gamma}^T\boldsymbol{\Sigma}\boldsymbol{\gamma}$

**Posterior scaled covariance**:
$$\boldsymbol{\Sigma}_{post} = \boldsymbol{\Sigma} - \frac{N}{\varepsilon^2 + N\boldsymbol{\gamma}^T\boldsymbol{\Sigma}\boldsymbol{\gamma}} \boldsymbol{\Sigma}\boldsymbol{\gamma}\boldsymbol{\gamma}^T\boldsymbol{\Sigma}$$

**Posterior mean**:
$$\mathbf{r}_{post} = \boldsymbol{\mu} + \frac{\delta}{2\varepsilon^2} \boldsymbol{\Sigma}_{post}(2\boldsymbol{\gamma}) = \boldsymbol{\mu} + \frac{\delta}{\varepsilon^2} \boldsymbol{\Sigma}_{post}\boldsymbol{\gamma}$$

This exactly matches Moffatt (2007) Equations 49-50 when written in the same scaling convention.

## 5. Connection to Moffatt (2007) Recursive Algorithm

The results correspond to one iteration of the macroscopic recursive algorithm:

**Measurement update**:
$$\mathbf{r}_{post} = \boldsymbol{\mu} + \frac{\delta}{V} \boldsymbol{\Sigma}\boldsymbol{\gamma} - \frac{N\delta}{V(V + Ns)} \boldsymbol{\Sigma}\mathbf{v}\mathbf{v}^T\boldsymbol{\Sigma}\boldsymbol{\gamma}$$

**Covariance update**:
$$\boldsymbol{\Sigma}_{post} = \boldsymbol{\Sigma} - \frac{N}{V + Ns} \boldsymbol{\Sigma}\mathbf{v}\mathbf{v}^T\boldsymbol{\Sigma}$$

These maintain the same computational structure as Moffatt (2007) while generalizing to include state-dependent noise.

## 6. Implementation Notes

1. **Initialization**: Start with equilibrium prior $\boldsymbol{\mu}_0$, $\boldsymbol{\Sigma}_0$
2. **Recursive processing**: Update mean and covariance for each measurement
3. **Computational efficiency**: Maintains $O(k^2)$ complexity per time step
4. **Numerical stability**: Sherman-Morrison provides stable covariance updates

The algorithm can estimate kinetic parameters, conductance, and channel number from macroscopic currents while properly accounting for both instrumental noise and state-dependent channel fluctuations.

---

*Using the Moffatt (2007) scaling convention where $\boldsymbol{\Sigma}$ and $\boldsymbol{\Sigma}_{post}$ are $N$-independent scaled covariance matrices, with actual covariances being $\frac{1}{N}\boldsymbol{\Sigma}$ and $\frac{1}{N}\boldsymbol{\Sigma}_{post}$ respectively.*