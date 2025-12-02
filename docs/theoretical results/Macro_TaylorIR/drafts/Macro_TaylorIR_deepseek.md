# MacroTaylorIR: Bayesian Inference for Time-Averaged Currents with State-Dependent Noise

## 1. Problem Formulation and Implementation Equations

### System Definition
We analyze ensembles of \(N_{\text{ch}}\) ion channels with discrete states, observing **time-averaged currents** over successive intervals \([t_n, t_{n+1}]\) while accounting for **state-dependent noise** contributions.

### Core Mathematical Objects
- \(N_{\text{ch}}\): number of ion channels
- \(\boldsymbol{\mu}_0\): initial state probability vector  
- \(\boldsymbol{\Sigma}_0\): initial state covariance matrix
- \(\mathbf{Q}\): Markov transition rate matrix
- \(\mathbf{P}(t) = \exp(\mathbf{Q}t)\): state transition probabilities
- \(\overline{\mathbf{\Gamma}}\): expected currents for state transitions
- \(\overline{\mathbf{V}}\): state-dependent noise contributions

### Complete Implementation Equations

**Predicted average current**:
\[
\overline{y}^{\text{pred}}_{0\to t} = N_{\text{ch}}\,\boldsymbol{\mu}_0 \cdot \overline{\boldsymbol{\gamma}}_0
\]

**Total prediction variance**:
\[
\sigma^2_{\overline{y}^{\text{pred}}} = \epsilon^2_{0\to t} + N_{\text{ch}}\,\widetilde{\gamma^{\!T}\Sigma\gamma} + N_{\text{ch}}\,\boldsymbol{\mu}_0 \cdot \boldsymbol{\sigma}^2_{\overline{\gamma}_0}
\]

**State-dependent noise terms**:
- Innovation: \(\delta = \overline{y}^{\text{obs}}_{0\to t} - \overline{y}^{\text{pred}}_{0\to t}\)
- Total variance: \(V = \sigma^2_{\overline{y}^{\text{pred}}}\)
- Effective direction: \(\mathbf{v} = \boldsymbol{\gamma} + \frac{\delta}{V}\boldsymbol{\sigma}^2\)
- Effective variance: \(s = \widetilde{v^{\!T}\Sigma v}\)

**Bayesian update equations**:

**Propagation**:
\[
\boldsymbol{\mu}^{\text{prior}}(t) = \boldsymbol{\mu}_0\mathbf{P}(t)
\]
\[
\boldsymbol{\Sigma}^{\text{prop}}(t) = \mathbf{P}(t)^{\!T}(\boldsymbol{\Sigma}_0 - \operatorname{diag}\boldsymbol{\mu}_0)\mathbf{P}(t) - \operatorname{diag}(\boldsymbol{\mu}^{\text{prior}}(t))
\]

**Measurement update**:
\[
\boxed{
\boldsymbol{\Sigma}^{\text{post}} = \boldsymbol{\Sigma}^{\text{prop}} - \frac{N_{\text{ch}}}{V + N_{\text{ch}} s}\,(\widetilde{v^{\!T}\Sigma})^{\!T}(\widetilde{v^{\!T}\Sigma})
}
\]
\[
\boxed{
\boldsymbol{\mu}^{\text{post}} =
 \boldsymbol{\mu}^{\text{prop}} 
 \,+ 
  \frac{\delta}{2V}
  \boldsymbol{\Sigma}^{\text{post}} 
  (\widetilde{\gamma^{\!T}\Sigma} + \widetilde{v^{\!T}\Sigma})
  }
\]
\[
\boxed{
\boldsymbol{\mu}^{\text{post}} =
 \boldsymbol{\mu}^{\text{prop}} 
\,+ \frac{N_{\text{ch}}}{V + N_{\text{ch}} s}\,(\widetilde{v^{\!T}\Sigma})^{\!T}  v^\top (p_0 - \mu)
 \,+ 
  \frac{\delta}{2V}
  \boldsymbol{\Sigma}^{\text{post}} 
  (\widetilde{\gamma^{\!T}\Sigma} + \widetilde{v^{\!T}\Sigma})
  }
\]

## 2. Solution Approach Summary

### Key Innovation
This framework unifies two methodological advances:

1. **From MacroIR**: Algebraic marginalization of boundary states via the unified tilde operator
2. **From MacroTaylor**: Laplace approximation with state-dependent noise modeling

### Mathematical Strategy
The core insight applies the Taylor/Laplace approximation framework to the **boundary state representation** while maintaining computational efficiency through the tilde operator.

**Computational Efficiency**: By deriving closed-form expressions for the tilde terms, we maintain \(O(K^3)\) complexity instead of the prohibitive \(O(K^6)\) required for explicit boundary state manipulation.

**Implementation Hint**: The algorithm processes each interval through:
1. Current and variance prediction using precomputed transition statistics
2. Bayesian updates via extended marginalized covariance vectors
3. Uncertainty propagation to subsequent time points

## 3. Detailed Solution Derivation

### 3.1 Theoretical Foundation

We begin with the measurement model incorporating both time-averaging and state-dependent noise:

\[
P(y|\mathbf{r}) = \mathcal{N}\left(y; N\mathbf{r} \cdot \boldsymbol{\gamma}, \varepsilon^2 + N\mathbf{r} \cdot \boldsymbol{\sigma}^2\right)
\]
\[
P(\mathbf{r}) = \mathcal{N}\left(\mathbf{r}; \boldsymbol{\mu}, \frac{1}{N}\boldsymbol{\Sigma}\right)
\]

The energy function (negative log-posterior) is:

\[
E(\mathbf{r}) = \frac{(y - N\mathbf{r} \cdot \boldsymbol{\gamma})^2}{2(\varepsilon^2 + N\mathbf{r} \cdot \boldsymbol{\sigma}^2)} + \frac{N}{2}(\mathbf{r} - \boldsymbol{\mu})^T \boldsymbol{\Sigma}^{-1} (\mathbf{r} - \boldsymbol{\mu})
\]

### 3.2 Laplace Approximation in Boundary State Space

The key innovation is applying the Laplace approximation to the **boundary state representation**. At \(\mathbf{r} = \boldsymbol{\mu}\):

**Gradient**:
\[
\mathbf{g} = -\frac{N\delta}{V}\boldsymbol{\gamma} + \frac{N\delta^2}{2V^2}\boldsymbol{\sigma}^2 = -\frac{N\delta}{2V}(\boldsymbol{\gamma} + \mathbf{v})
\]

**Hessian**:
\[
\mathbf{H} = N\boldsymbol{\Sigma}^{-1} + \frac{N^2}{V}\mathbf{v}\mathbf{v}^T
\]

### 3.3 Boundary State Projection via Tilde Operator

The critical step is projecting these instantaneous quantities through the boundary state using the tilde operator:

**Effective gradient in boundary space**:
\[
\widetilde{\gamma^{\!T}\Sigma} + \widetilde{v^{\!T}\Sigma}
\]

**Effective Hessian in boundary space**:
The second term \(\frac{N^2}{V}\mathbf{v}\mathbf{v}^T\) projects to \(\frac{N^2}{V}(\widetilde{v^{\!T}\Sigma})^{\!T}(\widetilde{v^{\!T}\Sigma})\) in the boundary state measurement update.

### 3.4 Update Equations Derivation

**Covariance Update**:
Using the matrix inversion lemma on the boundary-state-projected Hessian:

\[
\boldsymbol{\Sigma}^{\text{post}} = \boldsymbol{\Sigma}^{\text{prop}} - \frac{N_{\text{ch}}}{V + N_{\text{ch}} s}\,(\widetilde{v^{\!T}\Sigma})^{\!T}(\widetilde{v^{\!T}\Sigma})
\]

**Mean Update**:
The Newton-Raphson step in boundary space:

\[
\boldsymbol{\mu}^{\text{post}} = \boldsymbol{\mu}^{\text{prior}} - \mathbf{H}^{-1}\mathbf{g} = \boldsymbol{\mu}^{\text{prior}} + \frac{1}{N_{\text{ch}}}\boldsymbol{\Sigma}^{\text{post}} \cdot \frac{N_{\text{ch}}\delta}{2V}(\widetilde{\gamma^{\!T}\Sigma} + \widetilde{v^{\!T}\Sigma})
\]


### 3.5 Mathematical Consistency Verification

**Dimensional Analysis**:
- All terms maintain proper units and scaling throughout
- The effective direction \(\mathbf{v}\) combines current sensitivity and noise contributions
- The scaling factor \(\frac{\delta(2V + N_{\text{ch}} s)}{2V(V + N_{\text{ch}} s)}\) properly balances measurement and state uncertainty

**Special Case Recovery**:
- **When \(\boldsymbol{\sigma}^2 = 0\)**: \(\mathbf{v} = \boldsymbol{\gamma}\), recovering exact MacroIR equations
- **When transitions are ignored**: Reduces to simple macroscopic filtering
- **With explicit boundary states**: Mathematically equivalent to \(O(K^6)\) computation

**Computational Complexity**:
- Each tilde expression requires \(O(K^3)\) operations
- The tilde operator reduces boundary state operations to \(K \times K\) matrix operations
- Dominant cost: matrix exponentiation and multiplication (\(O(K^3)\))

### 3.6 Implementation Considerations

**Numerical Stability**:
- Use logarithmic computations for small probabilities
- Ensure covariance matrices remain positive semi-definite
- Employ specialized algorithms for matrix exponentials

**Precomputation Strategy**:
- Cache transition matrices for common interval lengths
- Precompute boundary-conditioned current statistics
- Optimize tilde operator implementations

This unified framework enables efficient Bayesian inference for ion channel kinetics using time-averaged current measurements while properly accounting for both instrumental noise and state-dependent channel fluctuations, achieving mathematical rigor with computational practicality through the elegant tilde operator formalism.

---

## Implementation-Oriented Notes

### Suggested C++ Naming:
- `tilde_gamma_Sigma_vec` ↔ \(\widetilde{\gamma^{\!T}\Sigma}\)  
- `tilde_v_Sigma_vec` ↔ \(\widetilde{v^{\!T}\Sigma}\)
- `tilde_gamma_Sigma_gamma_scalar` ↔ \(\widetilde{\gamma^{\!T}\Sigma\gamma}\)
- `tilde_v_Sigma_v_scalar` ↔ \(\widetilde{v^{\!T}\Sigma v}\)

### State-space objects:
- `mu0`, `Sigma0`
- `mu_prior`, `Sigma_prop`
- `mu_post`, `Sigma_post`
- `P_t`
- `Gamma_bar`, `V_bar`
- `gamma0_bar`, `sigma2_gamma0`

These names remain close to the math while keeping implementation straightforward.

---

*The unified tilde operator convention treats \(\widetilde{\cdot}\) as a single operator whose expressions have different algebraic ranks depending on input structure, maintaining mathematical elegance while enabling computational efficiency.*