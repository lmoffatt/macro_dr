# macroir v1 — implementation plan

Rewritten 2026-07-25, replacing a version that was wrong in structure. What it
got wrong is recorded in "How the first attempt failed" at the end, because the
failure is instructive and cheap to repeat.

Companion: `flat_model_language.md` (the level-0 model declaration).
Repo: `/home/lmoffatt/Code/macroir`, MIT. Reference: `macro_dr`, GPL-3, same
author.

---

## What macroir is

A small C++ library that computes, for a kinetic scheme and a current recording,
the log-likelihood, its gradient, and the Gaussian Fisher information, by the
MacroIR algorithm, exposed to R and Python.

Its first job is to **reproduce macro_dr's numbers**. Not to improve on them.

---

## The three rules

### R1 — Faithful. The default is always a transcription.

Transcribe macro_dr line for line, thresholds and branches included, with a
comment naming the file and lines. Where its formula has an odd constant, the odd
constant goes in, with a note and nothing more.

Anything that differs is a defect, whatever its merits, because a single
deliberate divergence anywhere in the chain makes every downstream comparison
uninterpretable: when logL disagrees you cannot tell your bug from your
improvement.

A variant that looks better is **not a code change, it is a research claim**.
Demonstrating that one likelihood algorithm beats another is precisely the
apparatus the eLife paper builds, and it is harder than writing the variant.
Do not write alternatives. Not behind a flag, not "for comparison". If you find
one, write a line in `NOTES.md` and keep going.

### R2 — A discrepancy is not a finding until it moves logL or the score.

Intermediate quantities are diagnostics for locating a disagreement once logL
already disagrees. They are not evidence of one.

Worked case: an intermediate conditional variance was found to lose all its
precision at small dt and even go negative, with relative errors reaching 1e9.
Irrelevant. `y_var = e + N gSg + N ms` with `e = Current_Noise/dt`, so the
residual term's share of `y_var` goes as dt². It matters at large dt, where the
computation is exact, and degrades at small dt, where it has vanished. Worst
measured impact on `y_var`: **2.4e-7**. Hours were spent on it and it was
reported twice as a finding.

Before writing the word "finding", compute how much it moves logL.

### R3 — Read the design before the formulas.

The single most costly mistake of the first attempt. Read `macro_dr` to learn
**how it is shaped**, then go back for the formulas. Specifically, before writing
anything, read:

- `legacy/derivative_operator.h` — the operators on `Derivative`
- `legacy/parameters_derivative.h:88-140` — the `d_d_` storage
- one full algorithm function in `legacy/qmodel.h`, start to finish

and see the property that organises the whole codebase, stated in R4.

### R4 — One body of code, two modes. This dictates the order of work.

macro_dr writes each algorithm **once**. Run it with `double` and you get the
value. Run the same function with `Derivative<double, Parameters_transformed>`
and you get the value and the gradient, because the operators propagate both:

```cpp
// legacy/derivative_operator.h:556
auto operator*(const T& x, const S& y) {
    return Derivative<F, X>(primitive(x) * primitive(y),
                            derivative(x)() * primitive(y) + primitive(x) * derivative(y)(),
                            get_dx_of_dfdx(x, y));
}
```

Three consequences, and they are why this is not a stylistic preference:

1. **You cannot have a sign error in the derivative that is not also in the
   value**, because there is no derivative code. Hand-writing the chain rule
   creates a second implementation that drifts from the first silently, with the
   value still correct.
2. **Refactoring the algorithm cannot break the derivative.** It follows.
3. **Branches are selected on `primitive(x)`** (`qmodel.h:816`, `:847-849`), so
   the derivative run takes the same branch as the plain run. Only possible
   because it is one body of code.

Each function templates **each argument separately** (`operator*(const T& x,
const S& y)`), not on a single `T`, so a derivative type can meet a plain one.

**Therefore the route is:**

> **First get logL working with plain doubles, end to end. Then template the same
> code and instantiate it with the derivative type. The score and the FIM come
> out of that, not out of new code.**

There is no work package called "derivatives". If one appears in your plan, you
have the architecture wrong.

---

## Route

Three stages. Each ends at a number that can be compared with macro_dr. Do not
start a stage before the previous one's checkpoint passes.

### Stage 1 — logL, plain doubles, end to end

Everything with `double`. No templates, no derivative type, no alternatives.

**Already in the repo and usable as is** (value-only, tested, 4 suites green):

| file | what | ported from |
|---|---|---|
| `matrix.hpp` | small dense row-major matrix | — |
| `error.hpp` | `Maybe<T>`, `Status` | `Maybe_error` |
| `eigen.hpp` | `dgeevx_` with macro_dr's gauge and sort | `lapack_headers.h:1278-1669` |
| `model.hpp` | level-0 scheme, `assemble_Q`, `assemble_g`, `validate` | `flat_model_language.md` |
| `schemes.hpp` | the four registered schemes, transcribed | `models_simple.h`, `models_MoffattHume_linear.h:12` |
| `faithful.hpp` | `Ee`, `E3` and their branches | `qmodel.h:814`, `:860` |
| `matrixfun.hpp` | `exp(Q dt)` from the decomposition | — |
| `qdt.hpp` | the whole interval assembly | `qmodel.h:1681-1780` |

`qdt.hpp` is complete: `to_transition_probability`, `kappa_F(V)`, the conjugate
shrinkage with its two priors, `gvar_ij`, the back-conversion, both range
canaries, the marginals.

**To write:**

1. **Experiment, protocol, recording.** Agonist concentration as a step function
   of time, sampling frequency, interval boundaries, the recorded current per
   interval. Handle missing samples (the predict-without-updating path) from the
   start; retrofitting it is painful. Reference: `qmodel.h:4483-4487` for the
   no-data branch.

2. **The IR filter step.** Read the whole function before writing any of it.
   - prediction, `y_mean`, `y_var`, `gS`, `gSg`, `ms`, `sigma_pre`:
     `qmodel.h:4488-4612`
   - the rank-1 update: `qmodel.h:5630-5702`
   - `qmodel.h:5703-5909` is diagnostic output, not part of the algorithm
   - the trust coefficient and its softmin: `qmodel.h:4124-4318`.
     **Only `calculate_trust_coefficient` is live.** `qmodel.h:4160` holds a
     second LogSumExp formulation marked `NOT CURRENTLY WIRED IN` and
     `calculate_psd_trust_coefficient` is computed and stored but never applied
     (`:5648`, `:5667`, `:5697`). Porting either gives a different score and
     different confidence intervals.
   - verified decomposition, useful as a check:
     `y_var = e + N gSg + N ms`, with `e = Current_Noise * fs / n_samples`
     (`qmodel.h:4539`), i.e. `Current_Noise / dt`.

3. **logL accumulation.** `l_t = -0.5 log(2 pi v_t) - 0.5 chi2_t` with
   `chi2_t = (y_t - mu_t)^2 / v_t` (`qmodel.h:4092`, assembled at `:4416-4417`).
   The Poisson branch at `:4095` only runs when `Proportional_Noise != 0`, which
   is 0 in the paper's models and drags in GSL. Skip it.

**Checkpoint 1.** Compute logL for one stored recording and compare against
macro_dr's. Nothing downstream starts until this matches.

### Stage 2 — the same code, with derivatives

1. **`Der<X>`**, mirroring `Derivative<X, Parameters_transformed>`: the value plus
   one object of the same shape per parameter (`parameters_derivative.h:94-99`
   for the scalar, `:131-136` for the matrix). Not a scalar dual number.

   macro_dr also carries a raw `Parameters_transformed const*` to check that two
   derivatives are with respect to the same theta, guarded by
   `MACRODR_DX_ASSERT`, which is compiled to a no-op in every build
   (`MACRODR_STRICT_DX_ASSERT` is defined nowhere). Leave the pointer out; it
   carries no behaviour.

2. **The operators**, from `derivative_operator.h:556-610`. Template each argument
   separately. Seeding is the identity in transformed space
   (`parameters_derivative.h:1068-1079`), and the log10 Jacobian is applied once
   explicitly (`:1028`, `parameters.h:164-166`).

3. **Template the stage 1 chain** on its input types and instantiate with `Der`.
   Do not write a second version of anything. Where stage 1 has `assemble_Q` and
   would want an `assemble_dQ`, there is one templated function.

   One thing must NOT carry a derivative: `kappa_F(V)`, the shrinkage
   pseudo-count, is computed on the primitive only (`qmodel.h:1076-1092`).
   Giving it one puts `d eps/d theta` into the shrinkage ratio and changes the
   score.

4. **The Fisher**, accumulated per step:
   `t_GFI = XXT(d_y_mean)/r_y_var + XXT(d_y_var)/(2 r_y_var^2)`
   (`qmodel.h:6109-6119`). PSD by construction. The score is not accumulated
   separately; it is the derivative of the accumulated logL.

**Checkpoint 2.** Score against central finite differences of the stage 1 logL,
and against macro_dr's score. FIM against macro_dr's.

### Stage 3 — the rest, in this order

Simulation by uniformization (`src/core/simulate.cpp`; decide the RNG contract
before v1, and do not port the `seed = 0` sentinel, which means `random_device`
with the resolved seed never logged). Then R and NR as a runtime enum on the
filter. Then LSE. Then the R and Python bindings, which are a few days each once
the core is header-only and the API is four functions. The optimizer is thirty
lines in R and in Python given logL, score and FIM; do not write one in C++.

---

## Out of scope for v1. Do not build these.

Allosteric model generation (it gets its own paper), MR, VR, IRT, micro, the
diagnostics battery (`Probit_statistics`, distortion, bootstrap), MCMC and
thermodynamic evidence, the `.macroir` DSL, memoization (measured on one cell:
0.980 s with, 0.951 s without), and OpenMP anywhere in the core.

---

## Testing

Only what a checkpoint needs. The first attempt built an elaborate degeneracy
ladder, quadrature references and branchless comparisons before there was a
likelihood, and none of it decided anything.

- Exact identities cost nothing and catch real errors: rows of `P` sum to one,
  `P(0) = I`, the semigroup property, and best of all a **constant conductance**,
  where `Abar = g` with zero variance whatever the trajectory, which exercises
  the whole divided-difference machinery against a closed-form answer.
- Finite differences are the only independent check on a derivative. Compare
  with a tolerance derived from the reference's own error budget
  (`eps ||f|| / h` for round-off plus `h^2` for truncation), not tuned until
  green.
- Compare against macro_dr at the checkpoints, not continuously.

Reference material: the frozen binary is at
`/home/lmoffatt/Code/macrodr-reference-0ffbda7/` with its provenance, including
the BLAS/LAPACK versions, which are part of it because the numbers depend on
which implementation is installed. It is the only generator of new reference
cases and cannot be rebuilt. `macrodr_cli --commit` self-identifies.

---

## Traps, all verified

- **Every `(sim, sample)` row of the dlik dumps is written twice**
  (`likelihood.cpp:2421`). Filter `segment_index == 0`.
- **"The sum of the per-step values equals the total" is tautological.**
  `qmodel.h:6092` is literally `logL = logL + t_logL` and the dumped per-step
  values are the summands. It returns 1e-16 with a wrong `y_var`.
- **`seed = 0` means `random_device`** and the resolved seed is never logged.
- **Four drifts between the symbolic model formulas and the running lambdas**,
  found by transcribing four schemes. Always follow the lambda:
  `scheme_CO` names `on`/`off` while its formula says `kon`; `scheme_1` carries
  the current's sign only in the lambda (`p[4] * -1.0`); `scheme_COC` conducted
  in state 1 in the lambda and state 2 in the formula, and declared `scheme_CCO`'s
  parameter names. The last two were fixed in `models_simple.h` on 2026-07-25.
- **`scheme_CO` gives `Current_Baseline` the value 0 with a log10 transform**,
  which has no logarithm. `scheme_1` overrides that parameter to linear.
- **The eigendecomposition is of `Q`, not of `Q dt`**
  (`calc_eigen(calc_Qx(m, x))`, `qmodel.h:1099`).
- **What IR consumes is `sum_j P_ij gvar_ij`**, not `gvar_i`, which is MR's.

---

## Open, and not the code's to decide

**What happens when a safeguard fires.** In macro_dr almost none is fatal:
`calc_Qdt_agonist_step` (`qmodel.h:3147-3166`) discards the error unread and
falls back to `calc_Qdt_taylor`, which carries a different pseudo-count
(`10 sqrt(N) eps` at `:3017` against `eps kappa_F(V)` at `:1694`). The same
pattern appears in `calc_Qdtg_agonist_step` and `calc_Qdtm_agonist_step`. So
production silently changes algorithm where the eigen path fails. macroir is
eigen-only and would stop. This matters at checkpoint 1: a macro_dr run that fell
back to Taylor is not comparable against a macroir run that refused.

**The production dt.** Whether `Simulation_n_sub_dt(100)` subdivides the interval
that reaches `calc_Qdt` or is only used by the simulation. Not traced. It decides
which regime the comparison exercises.

**The name.** `macroir` collides with the `.macroir` DSL extension and with the
algorithm's own name.

---

## How the first attempt failed

Recorded because it is cheap to repeat and expensive to discover.

**It built bottom-up from the pieces it understood.** Matrix, eigendecomposition,
divided differences, interval moments: 1.807 lines and 86 tests, all passing, and
no likelihood. The deliverable never got closer.

**It wrote value code and derivative code separately**, which made "derivatives
of the interval moments" look like a work package. Under R4 that package does not
exist. Discovering this meant undoing the architecture.

**It substituted its own judgement four times** — a different route for
`dP/dtheta`, a different third divided difference, and two safeguards macro_dr
does not have — each defensible in isolation and each a divergence that would
have made checkpoint 1 uninterpretable. All four were reverted.

**It reported two numerical differences as findings** that moved logL by 2.4e-7,
after hours of characterising them.

**It read macro_dr for formulas, not for shape.** Mining `qmodel.h` for `E3` and
`parameters_derivative.h` for `Omega` while never seeing that the codebase is
organised around one body of code running in two modes. That is the root of the
other four.
