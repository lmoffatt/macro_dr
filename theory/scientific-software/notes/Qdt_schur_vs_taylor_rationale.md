# Qdt via Schur vs via Taylor — rationale

## The question

The codebase has three `calc_Qdt*` paths: `_eig` (eigendecomposition),
`_taylor` (uniformization + composition-rule doubling on N×N moment
matrices), and `_schur` (Padé scaling-and-squaring on a 3N×3N Van Loan
augmented matrix). They co-exist; cost analyses sometimes claim large
performance gaps between them. This note records that those gaps are an
artifact of dense vs sparse implementation, not algorithmic differences.

## The algebraic identity

Both `_schur` and `_taylor` are evaluating the **Van Loan integral identity**.
For a generator $Q$ and conductance $G = \mathrm{diag}(\boldsymbol\gamma)$,
build the block-upper-triangular augmented matrix

$$
M \;=\;
\begin{pmatrix}
Q & G & 0 \\
0 & Q & G \\
0 & 0 & Q
\end{pmatrix}\cdot dt.
$$

Then [Van Loan 1978]:

$$
\exp(M) \;=\;
\begin{pmatrix}
P & F & B \\
0 & P & F \\
0 & 0 & P
\end{pmatrix},
$$

with

$$
P \;=\; e^{Q\,dt},
\quad
F \;=\; \int_0^{dt} e^{Q s}\,G\,e^{Q (dt-s)}\,ds,
\quad
B \;=\; \int_0^{dt}\!\!\int_0^s e^{Q r}\,G\,e^{Q (s-r)}\,G\,e^{Q(dt-s)}\,dr\,ds.
$$

`(P, F, B)` is the triple every Qdt-construction path must produce.
Algorithmically, both `_schur` and `_taylor` are computing this same
triple. Different storage conventions and different small-interval
approximations are the only substantive differences.

## Storage equivalence

| Implementation | What's stored | Total entries |
|---|---|---|
| Dense Schur (current) | The full 3N×3N matrix $M$ and its exponential | **9 N²** (with diagonal blocks $P$ stored 3× and $F$ blocks stored 2×) |
| Sparse-aware Schur | The 3 unique blocks $P$, $F$, $B$ | **3 N²** |
| Taylor (`Qn`) | `P`, `PG_n` (=$F$), `PGG_n` (=$B$) | **3 N²** |

Sparse Schur and Taylor have identical storage. The dense form just
duplicates blocks for the sake of working with one matrix object.

## Compute equivalence

The composition rule when concatenating two intervals of length $\tau$
into one interval of length $2\tau$:

$$
P(2\tau) = P(\tau)^2,
\qquad
F(2\tau) = P(\tau)\,F(\tau) + F(\tau)\,P(\tau),
$$

$$
B(2\tau) = P(\tau)\,B(\tau) + B(\tau)\,P(\tau) + F(\tau)\,F(\tau).
$$

This is exactly $\exp(M\,\tau)^2 = \exp(M\,2\tau)$ written out block-by-block.

| Implementation | Per-doubling cost |
|---|---|
| Dense Schur squaring | 27 N×N matmuls (much redundancy: 9 output blocks, most duplicates) |
| Sparse-aware Schur squaring | **6 N×N matmuls** ($P^2$, $P F$+$F P$, $P B$+$B P$+$F F$) |
| Taylor doubling (`sum_Qn`) | **6 N×N matmuls** (the same three formulas) |

Sparse Schur and Taylor have identical compute cost. The dense form
wastes ~4.5× on redundant block multiplications.

## The only substantive difference

The two paths differ in how they compute the **small-interval
expansion** (after scaling so $\|Q\cdot dt/2^n\|$ is small):

- **Schur**: Padé(13,13) rational approximant of $\exp$.
- **Taylor**: polynomial truncation of $\exp$ at order $k \approx 5\text{–}6$.

Both achieve double-precision accuracy in the small-norm regime. Padé
covers a wider $\|Q\cdot dt\|$ band per call (so fewer doublings
needed); Taylor is cheaper per call (so more doublings if `dt` is large).
For macro_dr's typical regime ($\|Q\cdot dt\| \ll 1$, with `0`–`2`
doublings), the two roughly tie.

## Derivative payload — where dense Schur loses

For derivative-aware computation, each $N\times N$ matmul costs
$(1 + 2p)\,N^3$ ($p$ = number of model parameters). The dense 3N×3N
Schur matrix has every block carry its own derivative payload, so
dense-Schur derivative cost is $27\,(1+2p)\,N^3$ per matmul: 27× over
the equivalent N×N op.

Sparse-aware Schur and Taylor only ever instantiate the three unique
N×N blocks, so the derivative payload is the same as a single
$(1+2p)\,N^3$ matmul. Identical to each other; 27× cheaper than dense
Schur.

This is why the current `calc_Qdt_schur` delegates to `calc_Qdt_taylor`
for derivative inputs ([qmodel.h:2479-2488](../../../legacy/qmodel.h)).
The 16×-ish speed-up that observation reflects isn't algorithmic — it's
dense-vs-sparse, exposed because the current Schur path is the dense
implementation.

## Conclusion

The three Qdt paths sit on a spectrum of *implementation choices for
the same algebraic identity*:

- `_eig`: diagonalize $Q$, evaluate `(P, F, B)` via divided differences of
  $\exp$ on the spectrum. Numerically distinct from `_schur` and
  `_taylor` because it factors through eigenvectors.
- `_schur` (dense, current): construct $M$ as one 3N×3N object, apply
  Padé scaling-and-squaring, harvest blocks. Algebraically equivalent
  to `_taylor`; computationally wasteful by ~4.5–27×.
- `_taylor`: build `(P, F, B)` block-wise from the outset, double via the
  composition rule. The natural sparse form of the Van Loan algorithm.

There is no architectural redesign to do here. `_taylor` is the right
workhorse for macro_dr's regime, `_eig` is a numerically distinct
cross-check, `_schur` provides a Padé/Van Loan reference for
plain-double cross-validation. The asymmetric delegation in
`calc_Qdt_schur` (Schur for plain double, Taylor for derivatives) is
the correct expression of this: dense Schur's only deficit is its
performance with derivatives, and the delegation routes around it.

## If you ever rewrite Schur

A sparse-aware, derivative-aware Schur implementation would:

- Store only the three unique blocks `(P, F, B)`.
- Express scaling-and-squaring as the doubling formulas above
  (`P^2`, `P F + F P`, `P B + B P + F F`) — identical operations to
  `sum_Qn`.
- Replace Taylor truncation of $\exp(Q\,dt/2^n)$ with a Padé(13,13)
  rational approximant of the small-interval $P$, then build the
  small-interval $F$, $B$ via Van Loan + Padé applied block-wise.

The result is operationally **the same algorithm as `_taylor`**, with
Padé instead of Taylor for the small-interval expansion. It would
remove the asymmetric delegation in `calc_Qdt_schur` and unify the
codebase around a single computational pattern. Whether the gain
justifies the work depends on whether anyone wants the Padé small-step
expansion (broader convergence band, classical algorithm with textbook
error analysis) as a default rather than a fallback.

## Reference

- C. F. Van Loan (1978). "Computing integrals involving the matrix
  exponential." *IEEE Trans. Automatic Control* 23(3): 395–404.
- N. J. Higham (2008). *Functions of Matrices: Theory and Computation.*
  SIAM. §10.4 (scaling-and-squaring).
- C. Moler & C. Van Loan (2003). "Nineteen dubious ways to compute the
  exponential of a matrix, twenty-five years later." *SIAM Review*
  45(1): 3–49.
