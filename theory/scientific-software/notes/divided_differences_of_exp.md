# Numerically stable evaluation of $E_1$, $E_2$, $E_3$, and $E_e$

## Where these functions come from

For a continuous-time Markov chain with generator $Q\in\mathbb{R}^{K\times K}$
and a state-dependent conductance vector
$\boldsymbol\gamma = (\gamma_1,\dots,\gamma_K)^\top$, the first two moments
of the integrated conductance
$A(t) = \int_0^t \gamma(X_s)\,ds$ — both unconditional and conditional on
endpoint states — are linear combinations of so-called
**Van Loan integrals** [VanLoan1978]:

$$
F_1(t) = e^{Q t},
\qquad
F_2(t) = \int_0^t e^{Q s}\,ds,
\qquad
F_3(t) = \int_0^t\!\!\int_0^s e^{Q r}\,\Gamma\,e^{Q(s-r)}\,dr\,ds,
$$

with $\Gamma = \mathrm{diag}(\boldsymbol\gamma)$, and a similar triple
integral $F_4(t)$ for the second moment. After eigendecomposition
$Q = V\Lambda V^{-1} = V\Lambda W$ with $\Lambda = \mathrm{diag}(\lambda_1,\dots,\lambda_K)$,
each Van Loan integral becomes an elementwise combination of the
eigenvectors with **divided differences of the exponential** applied to the
spectrum scaled by $t$:

$$
F_2(t)_{ij}
\;=\;
\sum_{m} V_{im}\,
        \underbrace{\frac{e^{\lambda_m t}-1}{\lambda_m}}_{= t\,E_1(\lambda_m t)}
        \,W_{mj},
$$

$$
[F_3(t)]_{ij}
\;=\;
\sum_{m,n} V_{im}\,
           \underbrace{
             \bigl(\widehat{WgV}\bigr)_{mn}\cdot
             t^2\, E_e(\lambda_m t,\lambda_n t)
           }_{\text{interaction term}}
           \,W_{nj},
$$

with analogous (triple-index) sums for $F_4$ involving the third
divided difference of $\exp$. Here $E_1$, $E_e$, $E_2$, $E_3$ are the
divided-difference quotients defined below.

These divided differences are smooth functions of their arguments (entire,
in fact), but their *naive evaluation* loses precision catastrophically
when two or more arguments approach each other. This note describes the
correct way to evaluate them, with names and references for the underlying
principles.

## Divided differences of $\exp$

For nodes $x_0,x_1,\dots,x_n\in\mathbb{C}$ and a smooth function $f$, the
$n$-th divided difference $f[x_0,\dots,x_n]$ is defined recursively
[deBoor2005]:

$$
f[x_0] = f(x_0),
\qquad
f[x_0,\dots,x_n]
\;=\;
\frac{f[x_1,\dots,x_n] - f[x_0,\dots,x_{n-1}]}{x_n - x_0}.
$$

When some nodes coincide, the formula is interpreted via the smooth
limit: at $x_0 = x_1 = \dots = x_n = x$, $f[\underbrace{x,\dots,x}_{n+1}]
= f^{(n)}(x)/n!$. The function $f[x_0,\dots,x_n]$ is symmetric in its
arguments and entire whenever $f$ is.

Applied to $f(x) = e^x$ (and the related $E_1$), the codebase's notation:

$$
E_1(x)
= \exp[0, x]
= \frac{e^x - 1}{x},
\qquad
E_e(x,y)
= \exp[x, y]
= \frac{e^x - e^y}{x - y},
$$

$$
E_2(x,y)
= E_1[x, y]
= \frac{E_1(y) - E_1(x)}{y - x}
= \exp[0, x, y],
$$

$$
E_3(x,y,z)
= \exp[x, y, z]
= \frac{E_e(y,z) - E_e(x,y)}{z - x}.
$$

The identification $E_2(x,y) = \exp[0,x,y]$ (a second divided difference
of $\exp$) is convenient because it places all four quantities under one
roof: $E_1, E_e, E_2, E_3$ are divided differences of $\exp$ of orders 1
through 3.

## The numerical pathology

Both numerator and denominator of $(e^x - e^y)/(x-y)$ vanish as $y\to x$.
The function is smooth; the *formula* is not. Specifically:

- $e^x$ has relative error $O(\varepsilon)$ where $\varepsilon \approx
  2.2\times 10^{-16}$ is double-precision machine epsilon.
- The subtraction $e^x - e^y$ loses one decimal digit of precision for
  every order of magnitude that $|x-y|$ is smaller than $|x|$. By
  $|x-y| \lesssim \varepsilon\,|x|$ the numerator is pure rounding noise.
- The denominator $x-y$ is small but accurate. Their ratio is therefore
  numerical garbage.

This is the classical phenomenon of **catastrophic cancellation** [Higham2002],
and the recipe is equally classical: rewrite the expression so that the
small-difference quantity is computed directly rather than as a
subtraction of two close numbers.

## The right recipe

### Order 1: $E_1$, $E_e$ via `expm1`

The C standard library function `expm1(z) = exp(z) - 1` is required to
return a result accurate to machine epsilon for *all* $z$, including
$|z|\ll 1$. This is the canonical antidote to the $\exp(z)-1$
cancellation and is available in `<cmath>` (and in every serious numerics
library).

The transformation:

$$
E_e(x,y)
= \frac{e^x - e^y}{x - y}
= e^y\cdot\frac{e^{x-y}-1}{x-y}
= e^y\cdot\frac{\mathrm{expm1}(x-y)}{x-y},
$$

is well-conditioned for *all* $x, y$. The remaining ratio
$\mathrm{expm1}(z)/z$ is itself well-conditioned (its Taylor series
$1 + z/2 + z^2/6 + \dots$ converges everywhere and starts at $1$); for
$|z|$ below a small threshold a direct Taylor expansion can replace it,
but in practice `expm1(z)/z` evaluated as written is already accurate
because `expm1` does the cancellation correctly internally.

The coincidence limit $x = y$ is just $E_e(x,x) = e^x$, obtained from the
formula via $\mathrm{expm1}(0)/0 \to 1$. The implementation handles this
as an exact $x = y$ branch (no FP comparison threshold needed if `expm1`
is used).

Same recipe for $E_1$:

$$
E_1(x) = \frac{e^x - 1}{x} = \frac{\mathrm{expm1}(x)}{x},
\quad
E_1(0) = 1.
$$

### Order $\ge 2$: Opitz form or Stewart's recursion

For higher-order divided differences, the cleanest robust method is the
**Opitz representation** [Opitz1964], also presented as Theorem 1 in
[Najfeld-Havel1995]:

$$
\exp[x_0, x_1, \dots, x_n]
\;=\;
\bigl(e^{J(x_0,\dots,x_n)}\bigr)_{0,n},
$$

where $J$ is the $(n+1)\times(n+1)$ bidiagonal matrix

$$
J(x_0,\dots,x_n)
=
\begin{pmatrix}
x_0 & 1 & & \\
    & x_1 & \ddots & \\
    & & \ddots & 1 \\
    & & & x_n
\end{pmatrix}.
$$

The Opitz form recovers the divided difference (the upper-right corner of
$\exp(J)$) by a single matrix exponential. Crucially, **the algorithm
that computes $\exp(J)$ for a small bidiagonal $J$ does not subtract
nearly-equal exponentials**: any standard Padé / scaling-and-squaring
exponentiator gives the divided difference accurately, regardless of
whether the nodes are well-separated or coincident.

For our $E_3$ in particular, $J$ is a $4\times 4$ matrix; one Padé
matrix-exp per evaluation. For workloads where this is too expensive, the
Taylor / divided-difference recursion of [Stewart1996] or [McCurdy-Ng-Parlett1984]
gives a comparable result at lower cost: it switches between

- the *outer formula* (recursive divided difference) when all node
  distances exceed a threshold;
- a *truncated Taylor series* about the cluster centroid when some
  nodes coincide.

The Stewart algorithm is the recipe used in MATLAB's `funm` for the
"matrix function on a triangular matrix" path, and it is what
[Higham2008] (§4.5, §10.6) recommends for general $f[x_0,\dots,x_n]$
with possibly clustered nodes.

### Threshold for switching

When the formula-vs-Taylor switch is needed (Stewart's path, or any
ad-hoc switch), the **correct switching threshold is relative, not
absolute**. The naive evaluation of $\exp[x,y]$ has the precision

$$
\frac{|\widehat{\exp[x,y]} - \exp[x,y]|}{|\exp[x,y]|}
\;\sim\;
\frac{\varepsilon\,|e^x|}{|e^x - e^y|}
\;\sim\;
\frac{\varepsilon\,\max(|x|,|y|)}{|x - y|}.
$$

So evaluating the formula directly is fine as long as

$$
|x - y| \;\gtrsim\; \sqrt{\varepsilon}\cdot\max(|x|,|y|,1).
$$

For double precision this puts the crossover at relative gap $\sim
10^{-8}$, **not** at any absolute threshold tied to a model constant like
`min_P`. Tying the crossover to `min_P` is dimensionally incorrect:
`min_P` carries a meaning ("probability mass below which entries are
treated as data-noise") that has nothing to do with the FP precision of
divided differences.

For $n$-th order divided differences ($n \ge 2$), the naive recursion
loses an additional factor of $|x-y|$ per recursive step, so the
threshold for using the Taylor branch is roughly
$|x-y| \lesssim \varepsilon^{1/(n+1)}\cdot\max(|x|,|y|)$ [Higham2008, §10.5].
At order 3 (our $E_3$) this is $\varepsilon^{1/4}\approx 1.2\times 10^{-4}$
relative, meaning the Taylor branch is needed surprisingly often. The
Opitz form sidesteps this concern entirely.

## What the current code does, and what to change

The current `Ee`, `E1`, `E2`, `E3` in [`legacy/qmodel.h`](../../../legacy/qmodel.h)
(around lines 761–855) use the **naive divided-difference formula with an
absolute-threshold switch**:

```cpp
if (sqr(primitive(x) - primitive(y)) < eps)
    return exp_x;                       // limit branch
else
    return (exp_x - exp_y) / (x - y);    // naive formula
```

The threshold `eps` is passed as `t_min_P = get<min_P>(m)()` at the call
sites in `calc_Qdt*_eig`. This has three issues:

1. **Overloads `min_P`.** The constant `min_P` is the Bayesian-shrinkage
   pseudo-count from
   [`bayesian_prior_regularization_of_Qdt.md`](../../macroir/notes/Gmean_ij_gvarij/bayesian_prior_regularization_of_Qdt.md);
   it has no relationship to the FP-precision threshold needed here.
   Changing one to fix the other is the wrong adjustment surface.

2. **Absolute threshold is dimensionally wrong.** The naive formula's
   precision loss is *relative* to $\max(|x|,|y|)$, so the right test
   is `|x - y| < eps_rel · max(|x|, |y|, floor)`, not
   `|x-y|² < eps_abs`.

3. **The naive formula is unnecessary.** $E_1$ and $E_e$ should use
   `expm1` and have no switch. $E_2$ and $E_3$ should use Opitz or
   Stewart's recursion, also eliminating the FP-comparison switch (or
   moving it inside a well-justified Taylor cutoff, with the threshold
   $\varepsilon^{1/(n+1)}\cdot\max(|x|,|y|)$ rather than `min_P`).

### Recommended replacements

$E_1$, $E_e$ (orders 0 and 1):

```cpp
double E1(double x) {
    return std::expm1(x) / x;          // accurate for all x including 0
}                                       // (limit x→0 is 1)

double Ee(double x, double y, double ex, double ey) {
    if (x == y) return ex;
    return ey * std::expm1(x - y) / (x - y);
}
```

`expm1(z)/z` is well-conditioned for $|z|\ll 1$. The `x == y` branch
handles the exact-coincidence case; no FP-threshold required.

$E_2$, $E_3$ (orders 2 and 3): either

(a) **Opitz form via a small `expm`**:

```cpp
double E_n(span<const double> nodes) {
    // Build bidiagonal J of size (n+1)x(n+1)
    // Compute expm(J) (Padé scaling-and-squaring, or a tiny
    // hand-rolled Pade-3,3 for n+1 ≤ 4).
    // Return (expm(J))(0, n).
}
```

This is robust for any node configuration but requires a Padé routine.

(b) **Stewart's recursion** with a relative threshold:

```cpp
constexpr double rel_eps_threshold_E2 = std::cbrt(2.2e-16);   // ε^{1/3} for order 2
constexpr double rel_eps_threshold_E3 = std::pow(2.2e-16, 0.25); // ε^{1/4} for order 3

double E2(double x, double y, double ex, double ey) {
    double scale = std::max({std::abs(x), std::abs(y), 1.0});
    if (std::abs(x - y) < rel_eps_threshold_E2 * scale) {
        // Taylor at (x+y)/2 — derivation in note appendix
        double m = 0.5 * (x + y);
        double d = 0.5 * (x - y);
        double em = std::exp(m);
        return em * (1.0 + d*d/6.0 + d*d*d*d/120.0 + ...);
    }
    return (E1(y, ey) - E1(x, ex)) / (y - x);
}
```

(with analogous structure for $E_3$).

Either approach replaces the `min_P`-keyed threshold and ports cleanly to
the derivative-aware overloads (the formulas are smooth, so
`Derivative<double>` propagation through them is automatic).

## Migration

1. Replace `E1`, `Ee` with `expm1`-based forms; drop the threshold
   parameter from their signatures.
2. Replace `E2`, `E3` either with the Opitz form (cleanest, modest extra
   cost) or with a Stewart-recursion + relative-threshold Taylor branch
   (cheaper, more lines).
3. Update the call sites in `calc_Qdt*_eig` to remove the `t_min_P`
   argument passed into these helpers — the helpers no longer need a
   precision-tuning constant from the model.
4. The `min_P` model parameter is then used only by the Bayesian
   shrinkage in `bayesian_prior_regularization_of_Qdt.md` and is renamed
   `min_P_prior` accordingly.

## References

- **[VanLoan1978]** C. F. Van Loan (1978). "Computing integrals involving
  the matrix exponential." *IEEE Trans. Automatic Control* 23(3): 395–404.
- **[MolerVanLoan2003]** C. Moler & C. Van Loan (2003). "Nineteen dubious
  ways to compute the exponential of a matrix, twenty-five years later."
  *SIAM Review* 45(1): 3–49.
- **[Higham2008]** N. J. Higham (2008). *Functions of Matrices: Theory
  and Computation.* SIAM. (§4.5, §10.5–10.6, §11)
- **[Higham2002]** N. J. Higham (2002). *Accuracy and Stability of
  Numerical Algorithms*, 2nd ed. SIAM. (Ch. 1: catastrophic cancellation)
- **[Najfeld-Havel1995]** I. Najfeld & T. F. Havel (1995). "Derivatives
  of the matrix exponential and their computation." *Advances in Applied
  Mathematics* 16(3): 321–375. (Opitz form; sensitivity analysis)
- **[Opitz1964]** G. Opitz (1964). "Steigungsmatrizen." *ZAMM* 44: T52–T54.
  (Original Opitz representation)
- **[deBoor2005]** C. de Boor (2005). "Divided differences." *Surveys in
  Approximation Theory* 1: 46–69.
- **[McCurdy-Ng-Parlett1984]** A. C. McCurdy, K. C. Ng & B. N. Parlett
  (1984). "Accurate computation of divided differences of the exponential
  function." *Mathematics of Computation* 43(168): 501–528.
- **[Stewart1996]** G. W. Stewart (1996). "Afternotes on Numerical
  Analysis." SIAM. (Recursion with Taylor cutoff for clustered nodes)
- **[Ng1992]** K. C. Ng (1992). "Argument reduction for huge arguments:
  good to the last bit." Sun Microsystems technical report. (`expm1`
  background)
