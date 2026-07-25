# macroir v1 — implementation plan for parallel agents

Status: draft. Audience: implementing agents and the maintainer.
Companion document: `flat_model_language.md` (the level-0 model declaration).

---

## 0. What this is

`macroir` is a small, self-contained C++ library that computes, for a kinetic
scheme and a current recording:

- the log-likelihood `logL`,
- its gradient `score`,
- the Gaussian Fisher information `FIM`,

using the MacroIR algorithm, plus forward simulation of recordings. It is
packaged for R and for Python.

It is **not** a reimplementation of `macro_dr` and it is **not** intended to
reproduce the eLife paper. It reimplements one algorithm, validated numerically
against `macro_dr`.

### In scope for v1

- Flat (level-0) runtime model declaration, single agonist.
- Algorithm family IR. Then R and NR as branches of the same recursion.
- LSE (classical non-linear least squares) as a separate, optional module.
- Simulation by uniformization.
- `logL`, `score`, `FIM`.
- R and Python bindings.

### Explicitly out of scope for v1

Do not build these. If a package seems to need one, stop and ask.

- Allosteric / conformational model generation (level 1). It gets its own paper.
- MR, VR, IRT, micro, and every Taylor variance-correction variant.
- The diagnostics battery: `Probit_statistics`, distortion measures, affine
  invariant distance, correlation distortion, autocorrelations, bootstrap.
- MCMC, parallel tempering, thermodynamic evidence, Bayes factors, priors.
- The `.macroir` DSL, a command manager, a script interpreter.
- Memoization / function tables. Measured on the reference cell: memoization is
  worth zero (0.980 s with, 0.951 s without, same RSS). Do not port it.
- An optimizer in C++. See P10/P11: Levenberg-Marquardt lives in R/Python,
  which is trivial once `logL`, `score` and `FIM` are available.
- OpenMP anywhere in the core. Parallelism belongs to the host language.

---

## 1. Ground rules for every agent

1. **Never modify `macro_dr`.** Read it freely. It is the reference. The only
   write into it is fixture generation, which is P0.2 and is done once.
2. **Never rebuild `macro_dr`.** A full build is ~8000 s of CPU and destroys
   the frozen binary's provenance. Read source, run the frozen binary.
3. **Own your files.** Each package below lists the files it owns. Do not write
   a file owned by another package. If you need something from another package,
   it is in a header that P0.3 already pinned; code against the header.
4. **Do not widen scope.** If your acceptance test passes and you are tempted to
   add a feature, stop. Write it in `NOTES.md` under your package instead.
5. **No dependencies** beyond Eigen and the standard library in `core/`.
   No Boost, no LAPACK link, no OpenMP, no fmt, no spdlog.
6. **C++17.** Not 20. Portability to older toolchains and to CRAN matters more
   than concepts do.
7. Every public function that can fail returns the error type pinned in P0.3.
   No exceptions across the library boundary, no `assert` as error handling.
8. When you port a formula from `macro_dr`, cite the source location in a
   comment: `// ported from legacy/qmodel.h:5629-5921`.

---

## 2. Repository skeleton

The maintainer creates the repo and this skeleton. Agents fill it in.

```
macroir/
  CMakeLists.txt
  core/
    include/macroir/          # header-only library
      error.hpp               # P0.3
      der.hpp                 # P0.3  (derivative bundle)
      model.hpp               # P0.3 / P1
      linalg.hpp              # P0.3 / P2
      matrixfun.hpp           # P0.3 / P2
      qdt.hpp                 # P0.3 / P6
      filter.hpp              # P0.3 / P7
      simulate.hpp            # P0.3 / P4
      lse.hpp                 # P0.3 / P8
      experiment.hpp          # P0.3 / P5
    tests/
      harness/                # P3
      model/                  # P1
      linalg/                 # P2
      qdt/                    # P6
      filter/                 # P7
      simulate/               # P4
  fixtures/                   # P0.2, generated, committed
    README.md                 # how each fixture was produced
  tools/
    gen_fixtures.sh           # P0.2
  bindings/
    r/                        # P10
    python/                   # P11
  docs/                       # P12
```

---

## 3. Sequencing

```
P0.1 ─ P0.2 ─┐
P0.3 ────────┼─→ wave 1:  P1  P2  P3  P4  P5   (parallel, no coupling)
             │
             └─→ wave 2:  P6 (needs P2)
                          P7 (needs P6, P1, P5)
                          P8 (needs P1, P5)
                          P9 (needs P7)
             └─→ wave 3:  P10  P11  P12
```

P0 is sequential and blocking. Everything in a wave runs in parallel.

---

## 4. Phase 0 — blocking, do first, in order

### P0.1 Preserve the reference binary

`build/gcc-release/macrodr_cli` is built from commit `0ffbda7`, which is the
provenance hash written into row 1 of the reference CSVs. `HEAD` is already
past it. Any rebuild destroys the only bit-exact reference and the only
generator of new reference cases.

- Copy it (and record `sha256`) to a location outside `build/`.
- Record: commit hash, compiler, flags, date.
- Write `fixtures/README.md` documenting how to invoke it.

Acceptance: the binary runs from its new location and reproduces one stored CSV.

### P0.2 Generate fixtures

Three tiers. All go in `fixtures/`, all committed.

**Tier A — interior fixture.** From
`projects/eLife_2025/figures/data/figure_1_likelihood_diagnostic_IR.csv`
(105 KB, 6 intervals, 2 states, `N_ch=20`). This is the only artefact that
exposes intermediate quantities: `P`, `gmean_i`, `gvar_i`, `gmean_ij`,
`gtotal_ij`, `d_GS`, and three sequential checkpoints of the filter. It is the
fixture the numeric packages develop against.

**Tier B — per-interval fixture.** From a short seeded run of
`figure_3_time` (`seed = 20260722`). Per-simulation seeds are drawn from a
master stream before the OpenMP loop, so a prefix of 4 recordings reproduces
the first 4 of the 1000 exactly (~0.78 s). Gives per-interval `score` and
`FIM`. **Filter `segment_index == 0`**: every `(sim, sample)` row is written
twice.

**Tier C — coverage fixtures. This is the one that does not exist yet and
matters most.** Every existing fixture is 2-state, so nothing exercises the
loss of precision in dividing by `(lambda_k - lambda_j)` as two eigenvalues
approach each other, nor the guard that handles it.

The frozen binary can only run schemes compiled into it, so tier C uses
**adversarial parameter values on the already-registered schemes**, not new
schemes. No recompilation; the 0ffbda7 binary stays the generator.

Two facts shape this:

- the eigendecomposition is of `Q`, not of `Q*dt`
  (`calc_eigen(calc_Qx(m, x))`, `legacy/qmodel.h:1099`), so eigenvalues carry
  units of inverse time and are of order 1e2 to 1e3 here;
- the guard at `legacy/parameters_derivative.h:1474` and `:1488` is
  **absolute**, `100*sqrt(eps) = 1.49e-6`. Against eigenvalues of order 1e3
  that is roughly 1e-9 in relative terms, so it fires only on numerical
  degeneracy and does not protect against the precision loss that begins
  several orders of magnitude earlier. The ladder therefore sweeps the
  **relative** gap and crosses the guard only at the bottom.

Generate, with the frozen binary:

- `scheme_CO`, `scheme_CCO`, `scheme_COC` and `scheme_1` at their reference
  parameter values (the baseline cases);
- the degeneracy ladder on `scheme_1`, five rungs. Family:
  `kon = 560.763`, `gating_off = 3*kon = 1682.289`, `[agonist] = 1`,
  `koff = gating_on = eps`. Two eigenvalues collide when the fastest binding
  rate meets the gating-off rate while the reverse rates are slow, and `eps`
  is the knob (the relative gap goes as `sqrt(eps)`):

  | rung | `eps` | absolute gap | relative gap | guard |
  |------|-------|--------------|--------------|-------|
  | C1 | 8.640953e+00 | 1.708e+01 | 1e-2 | – |
  | C2 | 8.413654e-02 | 1.682e-01 | 1e-4 | – |
  | C3 | 8.411467e-04 | 1.682e-03 | 1e-6 | – |
  | C4 | 8.411445e-06 | 1.682e-05 | 1e-8 | – |
  | C5 | 8.411442e-08 | 1.682e-07 | 1e-10 | **fires** |

  `tools/find_degenerate_params.py` in the `macroir` repo regenerates this
  table and prints the `log10` values the scripts need.

- one 3-state case from `scheme_CCO` as a sanity fixture. A 3-state chain has
  far less room: with rates in `[1, 1000]` the smallest reachable relative gap
  is about 2e-3, and even across `[1e-3, 1e6]` only about 9e-6. The deep rungs
  need `scheme_1`.

**Label rungs C4 and C5 in `fixtures/README.md` as numerical stress tests, not
physiology.** A rate of `8.4e-8` per second is not a channel anybody has
measured. They exist to probe the arithmetic.

For every fixture record: the model in level-0 form, the parameter vector, the
protocol, the recording, and every output quantity available.

Acceptance: `tools/gen_fixtures.sh` regenerates every fixture byte-identically
from the frozen binary.

### P0.3 Pin the interfaces

No implementations. Headers with signatures, types and doc comments only,
compiling against a stub. This is what makes wave 1 parallel, so it must be
complete before wave 1 starts.

Two decisions must be made here and cannot be deferred:

**D1 — derivative representation. DECIDED: mirror `macro_dr`.**

Carry the value plus one derivative object of the same shape per parameter,
which is what `macro_dr` already does with `Derivative<X, Parameters>`:

```cpp
template <class X>
struct Der {
    X value;
    std::vector<X> d;      // d.size() == n_params, one partial per parameter
};
```

Not a scalar dual number: that keeps matrix operations as matrix operations
(BLAS-friendly, no per-element allocation), and the existing codebase
deliberately templates each input separately because templating on a single
scalar `T` breaks derivative propagation.

### The mirroring rule

The derivation strategy in `macro_dr` works and is the author's. **Mirror the
mathematics and the numerical strategy; do not mirror the genericity
machinery.**

Mirror, structure for structure, name for name where it helps:

- the forward propagation architecture (`Derivative<X, Parameters>` and its
  operations),
- the eigendecomposition derivative, including the near-degeneracy guard,
- the divided differences `E1 / Ee / E2 / E3`,
- the Qdt assembly,
- the IR filter step and the accumulation of `logL`, `score`, `FIM`.

Do **not** mirror:

- `constexpr_Var_domain`, the 39-alternative `std::variant`,
  `merge_Maybe_variant`, and the `if constexpr` dispatch on algorithm flags.
  That machinery exists because `macro_dr` carries seven algorithm families;
  `macroir` carries one, plus `R` and `NR` as a runtime enum. Importing it
  would import exactly the complexity the rewrite exists to shed.
- memoization and function tables (measured: worth zero on the reference cell),
- the `Vector_Space` tag machinery for the diagnostic slots.

**Why mirroring is worth more than it looks.** If the propagation is
structurally the same, a disagreement can be localised by comparing intermediate
derivative objects one to one, instead of staring at a difference in the final
score. It converts validation into directed debugging.

**And what it costs.** A mirrored implementation inherits any bug in the
original, and the fixtures cannot see it, because the fixtures come from the
original. The finite-difference cross-check is therefore not optional and not a
nicety: it is the only test in the whole plan capable of catching a defect that
was copied faithfully. Every package that produces a derivative carries one.

**D2 — error type.** Recommendation: a `Maybe<T>` carrying either a value or a
message, mirroring `macro_dr`'s `Maybe_error`, with an explicit `and_then`.
Pin the exact signature; every package returns it.

Acceptance: all headers compile standalone; a stub program that calls every
public entry point links.

---

## 5. Wave 1 — parallel

### P1 — Model layer

**Owns:** `core/include/macroir/model.hpp`, `core/tests/model/`.

Implement `flat_model_language.md`: parse and validate a level-0 declaration,
and assemble from it, at runtime:

- `Q(a)` for a given agonist concentration `a`,
- `dQ/dtheta_p` for every parameter — a **constant sparse matrix per parameter**,
  computed once at load, with entries `ln(10) * exponent_p * q_ij` for
  `log10` parameters. There is no autodiff and no expression AST in this
  package. If you find yourself writing an expression evaluator, you have
  misread the spec.
- the conductance vector `g`, sign included,
- the initial distribution,
- parameter names, in order.

Also implement the reverse map (`Q, g -> parameter vector`), which in the
current code is a hand-written third copy and here must be **derived** from the
same declaration.

Implement the validation rules listed in `flat_model_language.md`.

**Reference:** `legacy/models_simple.h` and `legacy/models_MoffattHume_linear.h`
show the three current representations. Note the two known drifts: `scheme_CO`
declares `on`/`off` but its formula says `kon`; `scheme_1` puts the current's
sign only in the lambda (`p[4] * -1.0`), not in `g_formula`.

**Acceptance (exact, no tolerance):** for `scheme_CO`, `scheme_CCO`,
`scheme_COC` and `scheme_1`, transcribed to level 0:
same state count; same sparsity pattern of `Q0` and `Qa`; same parameter names
in the same order; same coefficients; same conductance vector including sign;
and `Q` assembled at the reference parameter values equal element-by-element to
the fixture. Plus: every validation rule has a test that trips it.

### P2 — Linear algebra and matrix functions

**Owns:** `core/include/macroir/linalg.hpp`, `matrixfun.hpp`,
`core/tests/linalg/`.

- Eigendecomposition of a real non-symmetric `Q` (Eigen `EigenSolver`; do not
  link LAPACK).
- `exp(Q * dt)` and the divided-difference quantities the Qdt assembly needs.
- The derivative of the above with respect to `Q`.

**Port, do not reinvent.** The derivation of the eigendecomposition derivative
is already solved in `macro_dr`: read
`legacy/parameters_derivative.h:1465-1520` and the divided differences at
`legacy/qmodel.h:763-870`, and port them. In particular understand the role of
the `100 * sqrt(eps)` guard before changing anything about it, and understand
whether `eig_enforce_q_mode` (defined at `parameters_derivative.h:1516` and
`matrix.h:2100`, currently never called) is needed here.

Also note `legacy/schur_parlett.h` (622 lines) exists as an alternative path.
v1 uses the eigen path only; record in `NOTES.md` any case where it is not
adequate.

**Acceptance:** against Tier A and Tier C fixtures, `exp(Q*dt)` and its
derivative to 1e-13 relative, floor 1e-15. Independently, the derivative agrees
with central finite differences on `Q` to 1e-6 relative. The Tier C
near-degenerate cases must pass, and a test must exist that fails if the
degenerate handling is removed.

### P3 — Fixture harness and comparison infrastructure

**Owns:** `core/tests/harness/`.

This package delivers before any numerics exist, and every other package
depends on its output format.

- Reader for each fixture tier.
- Comparator with a per-quantity tolerance policy. **Absolute floor 1e-15**:
  `macro_dr` writes with `setprecision(digits10+1)` = 16 digits, not 17, so
  doubles do not round-trip. Asking for more is chasing a ghost.
- Default tolerances: primitives 1e-13, derivatives 1e-11, aggregated `logL`
  and `score` 1e-12, `FIM` 1e-10, all relative.
- A report that localises a mismatch to (quantity, interval, parameter index)
  rather than printing a scalar difference.
- **A mutation test of the harness itself.** Inject a known error (e.g. drop
  the factor `2.0` at `legacy/qmodel.h:1740`, or perturb one fixture value by
  1e-9) and assert the harness reports it. Without this, "the battery passes"
  means nothing.

**Do not** use "sum of per-step values equals the total" as a correctness
check. It is tautological: `legacy/qmodel.h:6092` is literally
`logL = logL + t_logL` and the per-step values dumped **are** the summands. It
returns 1e-16 even when the formula for `y_var` is wrong. Keep it as a smoke
test, label it as such.

**Acceptance:** the harness reports green on the fixtures against the frozen
binary's own output, and red on every injected mutation.

### P4 — Simulation

**Owns:** `core/include/macroir/simulate.hpp`, `core/tests/simulate/`.

Exact CTMC realisation by uniformization: given a model, parameters and a
protocol, produce a current recording. Independent of the derivative machinery;
nothing here carries a `Der`.

**Decide the RNG contract now, not after v1.** Users in R expect `set.seed()`
to control it; users in Python expect a `Generator` or an explicit seed
argument. Changing this after release invalidates third parties' results. The
core takes an explicit seed or engine; the bindings adapt.

**Do not port the `seed = 0` sentinel.** In `macro_dr`, `seed = 0` means
`std::random_device` and the resolved seed is never logged, which is why
figures 4 and 5 of the paper cannot be reproduced. In `macroir` a seed is
always explicit and is always recorded in the output.

**Reference:** `src/core/simulate.cpp`, in particular the master-stream seed
derivation around lines 292-302.

**Acceptance:** distributional tests (mean and variance of the current against
the analytic `y_var = e + N*gSg + N*ms`); and, on a fixed seed, a recording
that is stable across runs, thread counts and platforms.

### P5 — Experiment, protocol and recording types

**Owns:** `core/include/macroir/experiment.hpp`.

Small but everything depends on it: agonist concentration as a function of
time (step protocol), sampling frequency, interval boundaries, sub-intervals,
and the recording (a vector of currents aligned to intervals). Handle missing
samples (the "predict, do not update" path) from the start; retrofitting it is
painful.

**Acceptance:** round-trips the protocol of every fixture; the interval
boundaries it computes match the fixture's `step_start` / `step_end` /
`n_step` columns exactly.

---

## 6. Wave 2

### P6 — Qdt assembly

**Owns:** `core/include/macroir/qdt.hpp`, `core/tests/qdt/`. **Needs P2.**

From `exp(Q*dt)` and the divided differences, assemble the per-interval
quantities the filter consumes: `P`, `gmean_i`, `gvar_i`, `gmean_ij`,
`gtotal_ij`, and their derivatives.

**Reference:** `legacy/qmodel.h:1581-1677`. Note that the current assembly
includes a Bayesian shrinkage step (commit `a3241c0`); understand what it is
for before deciding whether v1 needs it, and record the decision.

**Acceptance:** every quantity against Tier A column-by-column at 1e-13, and
against Tier C at the same tolerance. Derivatives cross-checked against central
differences at 1e-6.

### P7 — IR filter, logL, score, FIM

**Owns:** `core/include/macroir/filter.hpp`, `core/tests/filter/`.
**Needs P6, P1, P5.**

The predict/update recursion with the integrated measurement, the boundary-state
conditioning, and the accumulation of `logL`, `score` and `FIM`.

The FIM to reproduce is the analytic Gaussian one accumulated per step:

```
t_GFI = XXT(d_y_mean)/r_y_var + XXT(d_y_var)/(2 * r_y_var^2)
```

which is PSD by construction and is what an optimizer wants as curvature.
(`legacy/qmodel.h:6109-6119`.) The numerical FIM is not part of the core.

**Reference:** `legacy/qmodel.h:5629-5921` for the step; `:4124-4318` for the
trust coefficient and softmin. **Warning:** `legacy/qmodel.h:4160` contains a
second LogSumExp formulation of the trust coefficient marked
`NOT CURRENTLY WIRED IN`. Porting that branch instead of the active one yields
a different score, a different Fisher and different confidence intervals. Port
the active path (`calculate_trust_coefficient`).

**Acceptance:** Tier A checkpoints at 1e-13; Tier B per-interval `score` and
`FIM` at 1e-11 and 1e-10; totals at 1e-12; Tier C at the same tolerances. Score
cross-checked against central finite differences in `theta` at 1e-6 on at least
five schemes, which is the independent check that catches sign errors.

### P8 — LSE

**Owns:** `core/include/macroir/lse.hpp`. **Needs P1, P5.** Independent of
P6/P7 and can run alongside them.

Classical non-linear least squares with marginalised noise (the Moffatt & Hume
2007 JGP method). No filter. Same `logL`/`score`/`FIM` interface.

**Acceptance:** against fixtures generated with `family = nonlinearsqr` at the
same tolerances. Note the known tautology in the reference: `r̄²_std ≡ 1` for
LSE by construction; do not treat it as a passing check.

### P9 — R and NR branches

**Owns:** additions to `filter.hpp` behind an enum. **Needs P7.**

`R` and `NR` are variations of the same recursion (recursive versus
non-recursive), not separate algorithms. Add them as a runtime enum on the
filter, not as template parameters, and not as a variant type.

**Acceptance:** fixtures generated with `macro_R` and `macro_NR` at the same
tolerances as P7.

---

## 7. Wave 3

### P10 — R binding
### P11 — Python binding

**Own:** `bindings/r/`, `bindings/python/`. Neither depends on the other.

Surface, minimal and identical in both:

```
scheme  <- read/construct a level-0 model
sim     <- simulate(scheme, theta, protocol, seed)
ll      <- loglik(scheme, theta, protocol, recording)   -> logL, score, FIM
fit     <- mle(...)                                     -> in the host language
```

`mle` is Levenberg-Marquardt / Fisher scoring written **in R and in Python**,
roughly thirty lines each given `logL`, `score` and `FIM`. Do not write an
optimizer in C++.

Both packages vendor `core/include` at release time so they build standalone.
The core is header-only precisely so this is a file copy.

R specifics: `Rcpp` + `RcppEigen`. No writes to the working directory, ever
(the current code writes `scheme_N_model_description.txt` during static
initialisation; that pattern is forbidden here). Distribution via r-universe or
GitHub first; CRAN is a later decision.

Python specifics: `pybind11` or `nanobind`, Eigen vendored or as a submodule.

**Acceptance:** in both languages, the three fixtures of Tier C give `logL`,
`score` and `FIM` matching the C++ core bit-for-bit; and a worked example
defines a scheme from scratch, simulates, and fits.

### P12 — Documentation and examples

**Owns:** `docs/`.

The one thing `macro_dr` never had: a document that tells a user how to define
their own kinetic scheme and fit their own recording. Three worked examples:
2-state, the 5-state `scheme_1`, and one with missing samples.

---

## 8. Decisions still open for the maintainer

1. **D1 and D2** in P0.3 (derivative representation, error type). Blocking.
2. Repo name. `macroir` collides with the `.macroir` DSL extension and with the
   algorithm name.
3. Whether the Bayesian shrinkage in the Qdt assembly (P6) is in v1.
4. Whether `coef` is stored as an integer (recommended, see
   `flat_model_language.md`) or as a double.
5. Licence, and whether `macroir` is public from the first commit.

---

## 9. Notes on using an agent for each package

- Give the agent this document plus `flat_model_language.md`, plus its package
  section, plus read access to `macro_dr`.
- Tell it its acceptance test is the definition of done, and that passing it
  while having written extra features is a failure, not a bonus.
- The packages that most reward a careful agent are P2 (matrix function
  derivatives) and P7 (the filter). P1 and P3 are mostly mechanical but are
  blocking, so they should go first or in parallel with everything.
- P3 should be finished and red before P2, P6 and P7 start producing numbers,
  so that the first number any of them produces is already compared.
