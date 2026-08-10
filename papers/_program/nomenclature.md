# Nomenclature: naming the methods

> Updated: 2026-07-29. Shared across the two papers; cited, never restated. Settled items graduate
> to `decisions.md`. Which methods sit in the body: `program.md` §1.
>
> Scope: what we call every method in the program, what the letters mean, how to describe them so each
> name picks out exactly one method, and what the regions of the usage map are called.

## The constraint

`IR` / `MacroIR` is in print (Communications Biology 2025, P2X2). The acronym is fixed; nothing here
renames it. What was open was the *gloss*, and the name of the one new member, `VR`.

**`VR` keeps its letter** (Luciano, 2026-07-28; `decisions.md` §2). It earned the rung by behaving
distinctly: it was predicted to come out over-confident and more so than MR, and it did. The hazard
below still applies.

## The structure has two levels, not one lattice

The earlier "one object with two knobs" description was true of the five likelihoods and became false
once least squares joined the comparison. The honest structure is a **root question with a ladder
hanging from it** (`program.md` §1):

- **Root: do you model the gating fluctuations at all?** `LSE` (classical nonlinear least squares on
  the mean current) answers no. In the engine it is `family_approximation = 2`; it carries the same two
  knob settings as `INR` and is distinguished only by that third flag. **It has no
  compositional name**, because it is a different kind of object, named for what it is.
  **Corrected 2026-07-28:** this used to read "no rung and no gloss", which is now misleading in the
  one place it matters. `LSE` is the bottom rung of the cost ladder the paper walks and it is the
  comparison anchor (`decisions.md` §1). What it lacks is a *compositional* name, not a place.
- **Given yes, what does the Gaussian condition on, and how is its variance accounted?** That is the
  lattice below.

## The lattice, as implemented

Verified against `projects/eLife_2025/runs/run-20260418-182133/script.macroir`. `recursive_approximation`
says whether the occupancy covariance is propagated between intervals; `averaging_approximation` (`av`)
says how the single-channel conductance is treated within an interval.

| Label | recursive | av | interval variance |
|-------|-----------|-----|---|
| `NR`  | false     | 0   | — |
| `R`   | true      | 0   | — |
| `INR` | false     | 1   | total (**= the published `MacroINR`**; supplement member) |
| `NMR` | false     | 1   | **none — the defective implementation, see below** |
| `MR`  | true      | 1   | total |
| `VR`  | true      | 1   | residual |
| `IR`  | true      | 2   | residual (+ boundary gain) |

**The last column is the internal conductance-variance form, NOT the predicted observable variance.**
Reading it as the latter produces the wrong conclusion that `IR` predicts less than `MR`; it predicts
the same, at the same prior, because the boundary term `IR` restores in `gSg` is exactly the one the
residual form drops from `ms`. See the qualifier paragraph below.

**`INR` (false, 1, total) is the published `MacroINR`, the Comm Biol 2025 control.** Two readings of
the `I` are available and **they name the same algorithm**, which is why the bridge needs no
lawyering. Read compositionally it is **I**nterval (the averaged conductance, `av = 1`) +
**N**on-**R**ecursive, the conductance prefix on the non-recursive base. Read as `IR`'s
boundary-conditioned `I` it is the non-recursive member of the `IR` family. **Those coincide**: the
end state is unobservable without an update, so `gmean_ij` can only enter through its row marginal
`gmean_i` (`legacy/qmodel.h:1666`, `gmean_i = gtotal_ij · u`), and `(recursive=false, av=1)` and
`(recursive=false, av=2)` are the same code path, branching only on `averaging > 0`. **This closes
the "open item: the MacroINR bridge" at the foot of this file** on much stronger ground than the
2026-07-29 parse argument it replaces: it does not matter which `I` was meant. The Comm Biol
Introduction's "ignores time averaging" is loose prose for "does not do `IR`'s boundary-conditioned
interval treatment", not a claim about the `av` flag.

**`NMR` is a second-class citizen and it is kept on purpose: a badly implemented algorithm that
exists** (Luciano, 2026-07-31). It is `INR` minus the `N·ms` term, the variance of the interval-mean
conductance given the start state. The regression entered at `a3e0a89` (2025-12-02), which split the
algo-state computation into recursive and non-recursive copies; the copy re-implemented the variance
without `ms` and accepted-and-ignored the `variance` flag, although every dispatch script passes
`variance_approximation = 1`. Fixed at `1f7138b` (2026-07-31). So `NMR` names exactly what ran in that
window, and it is **not** `MacroINR`: the submitted Comm Biol source
(`macro_dr_submission` `b4a0e28`) adds the term in the same place, gated on `variance` alone.

**The data key `macro_NMR` is therefore not a legacy spelling, it is a different algorithm.** Every
file in `figures/data/` carrying `macro_NMR` was produced by the defective engine; `macro_INR` marks
the corrected one, and the `.Rmd` readers tell them apart by that label. Do not sweep `macro_NMR` out
of the data or the readers. (Line 1 of every produced CSV is the git hash, so the two are also
separable by provenance.)

**What `NMR` measures, and what it does not.** `NMR` is **numerically indistinguishable from `NR`** on
the freeze, and that now has a mechanism rather than being a coincidence: with `ms` gone, the only
things left separating them are `g → gmean_i` and `P_half → P`. **Do not carry that forward as
"the interval-mean conductance buys nothing without recursion"** — that sentence was measuring the
code as written.

**Answered 2026-08-01: `INR` is NOT indistinguishable from `NR`, in either moment.** The re-run is
`figures/data/1f7138b/` (17 cells at `nsim` 10⁴). At noise 0.1, N_ch 10⁴, on `k_off`, information
distortion total / sample / correlation: `NR` 78.9 / 47.2 / 21.2 against `INR` 22.2 / 1.00 / 22.1. And
median |bias| over every cell and interval, in `N_ch`: `NR` 0.099 against `INR` 0.003 log10. So the
`N·ms` term is what carries the first moment and the per-sample fidelity, and what it leaves behind is
pure correlation distortion, which only recursion removes. Recompute:
`papers/1_method/decisions/recompute/d3_interval_vs_recursion_2x2.py`; roster consequence: Q-5 in
`papers/1_method/decisions.md`.

The suffix (`NR` / `R`) is the occupancy axis: non-recursive or recursive. The prefix is the
conductance axis: none for the instantaneous conductance, `M` for the mean conductance, `I` for the
interval-conditioned (boundary) mean conductance, which above the non-recursive base collapses onto
the mean conductance.

> **CORRECTED 2026-08-06 (Luciano).** This line read "The suffix (`N` / `R`)", which does not parse:
> `NR` ends in `R`, so the one member the letter is there to mark as non-recursive comes out
> recursive. The split that works on all six names is `NR` against `R` (`NR`, `INR` end in the first;
> `R`, `MR`, `VR`, `IR` in the second). Caught while writing the naming rule into the manuscript,
> where it had never appeared: `02_theory.tex`, the paragraph after the two axis paragraphs, now
> states this rule in the body, so the two must be kept in step.

**`VR` opens a third axis: the form of the interval variance.** `MR` and `VR` share `(recursive, av) =
(true, 1)` and differ only in whether the interval variance is the total per-start-state form or the
residual (boundary-conditioned) form. Before `VR`, no name needed to distinguish the two variance
forms because only one was ever run; now one does.

## What the letters mean

- `M` = **mean conductance**. Averaged over the interval, conditioned on the state at the interval's
  *start*. The K-vector (K = number of states)

  (γ̄₀)ᵢ₀ = Σᵢₜ Pᵢ₀→ᵢₜ(t) · Γ̄ᵢ₀→ᵢₜ

  with Γ̄ᵢ₀→ᵢₜ the mean interval-averaged conductance of a channel from state i₀ to state iₜ, and
  P(t) = exp(Qt).
- `I` = **interval-conditioned (boundary) mean conductance**. The same average conditioned on the
  states at *both* ends, so the object stays K×K (the pair (i₀, iₜ), the **boundary state**).
- `V` (in `VR`) = the **residual (boundary-conditioned) interval variance**, `Σⱼ P_ij·gvar_ij`, as
  opposed to the total per-start-state variance `gsqr_i − gmean_i²` that `MR` carries. The two differ
  by `Var_j[gmean_ij|i]`, the spread of the interval-mean across end states.

Two concrete differences drive the whole ladder:

- `MR → VR` replaces the total variance with the residual variance.
- `VR → IR` adds the boundary cross-covariance term N·γᵀΣγ, which enters **both** the predictive
  variance and the **gain**.

**Neither step is confined to what it flips, and "changes nothing else" must not be written of
either** (2026-08-06; retired in `program.md` §78, `../1_method/decisions.md:161`,
`../1_method/CONTINUE_HERE.md:88` and `figures_build_plan.md:212`, and now in the manuscript). The
predictive variance divides the gain, so flipping the variance form moves the update at once. The
two steps are separable in the algebra without being separable in what a recording does with them.

**And the identity needs its qualifier.** At the SAME PRIOR, `MR` and `IR` assign the interval the
same total predictive variance: the boundary term enters `gSg` with a plus and `ms` with a minus and
cancels, and `MR`'s contraction against the full `P_Cov` supplies exactly the `γ̄ᵀdiag(μ)γ̄` the total
form subtracts. So at equal state the whole `MR`→`IR` difference is the gain. Along a RECORDING they
do not report the same variance, because the gain makes the priors diverge: on the figure-1 dumps
`MR = IR = 1.048475625791748` at every interval where they still share a prior, then `MR` runs +69%
to +81% above `IR`. `VR` is below both from any state, by `μᵀVar_j[gmean_ij|i]` (0.378085 at that
same interval), so **`VR` is not "`MR`'s mean with `IR`'s variance"**: `IR`'s *total* predicted
variance is `MR`'s, and the residual form is a component `IR` uses inside a different decomposition.
Verified 2026-08-06 against `legacy/qmodel.h:4568-4619` term by term, numerically on random fields
(difference 0.0, structural), and on the dumps; numbers and the fuller statement in
`../1_method/figures_build_plan.md:195-245`.

`av` literally counts the conditioned endpoints (0 instantaneous, 1 start, 2 boundary); the flag is
self-documenting. The variance axis is not in `av`, which is why it needs its own letter.

## `VR`: the name is settled, the hazard is not

**`VR` is settled** (2026-07-28, closing `program.md` §9). The second caution below is discharged: VR
was measured and it does have a distinct behaviour. The first still binds.

- **The `V` collides with the cut Taylor variance-correction variants** `MRV`, `IRV`, and with the
  engine flag `taylor_variance_correction`. `VR` is *not* a Taylor variant. If the name survives,
  Methods must say so in one sentence, and the March Taylor data must be deleted so the two `V`s never
  appear in the same tree.
- The name earns its place only if `VR` turns out to have a distinct behaviour worth a rung. If it does
  not, describe it as "MR with the residual variance" and spend no letter. Precedent for waiting: `MR`
  was fixed as "strawman" before it was measured, and the label outlived the data that contradicted it
  in two cells.

## Describing `IR` uniquely: use "boundary-conditioned"

`MR` is also interval-averaged (av = 1) and recursive, so "the interval-averaged recursive likelihood"
does not identify `IR`. The unique identifier is that `IR` conditions on **both** interval endpoints.

Two ways a bare "interval-conditioned" misfires: it fails to discriminate `IR` from `MR`, and it can be
read as conditioning on *everything inside* the interval (the full trajectory), which would be the
exact likelihood `IR` does not compute. Both are avoided by one word.

**In the manuscript body, call `IR` boundary-conditioned**, not interval-conditioned. The word says
which part of the interval is conditioned on and makes explicit that the interior is marginalized. In
an abstract, either name MacroIR without a gloss, or describe the family by its axes so the unique
identifier falls out.

**Retired phrasing.** Do not write "the sole survivor" or "the only member calibrated across the
practical regime" (`decisions.md` §6). The ladder still makes `IR` the top rung, and that structural
point stands; what is retired is stating a one-band result as a global verdict.

## `ILSE`: the window axis crosses the root question (2026-08-06)

**The least-squares rung splits in two.** The window axis is not the family's property: the recording
is an interval average, so *any* method predicting it chooses between the mean at an instant and the
mean over the window, least squares included. `averaging_approximation` is live on the `family = 2`
branch and is the only flag that separates them.

| Name | av | data key | dump token | mean model |
|---|---|---|---|---|
| `LSE`  | 0 | `nonlinearsqr_g` | `..._LSE_av0` | deterministic mean at the sampling instant |
| `ILSE` | 1 | `nonlinearsqr`   | `..._LSE`     | deterministic mean averaged over the interval |

**THE TOKENS ARE CROSSED.** The file written `..._LSE` is av = 1 and is therefore `ILSE`. Kept that
way on purpose so nothing already reading it is silently relabelled; the crossing is resolved in
`figure_1.Rmd` and `figure_3.Rmd` and nowhere else.

**Each arm shares its mean model EXACTLY with the macroscopic member at the same `av`**, and differs
only in the predictive variance. Measured from the figure-1 dumps 2026-08-06, 12 intervals:
max|LSE − NR| = 8.882e-16, max|ILSE − INR| = 2.220e-16, max|LSE − ILSE| = 5.469e-02.

**Where each is measured.** `ILSE` is the ONLY least-squares arm on the design grid: zero
`nonlinearsqr_g` files under `figures/data`, so every plane/map/cloud number the paper quotes as
"LSE" is `ILSE`. `LSE` proper is measured only at the single cell of figure 3 (both arms produced by
one run of `figure_3_time.macroir`, git `ccd26f9-dirty`, 2026-08-05) and drawn in figure 1.

**The av = 0 branch was BROKEN before 2026-08-05**: it propagated the occupancy by `P_half` and
evaluated the mean before advancing, so the member fell half a step further behind every interval.
Fixed at `legacy/qmodel.h:7490`. Any `LSE_av0` dump predating that is wrong.

**FIGURE LABELS THAT STILL LAG.** These print `LSE` for the av = 1 arm and want `ILSE`:
`figure_2.Rmd:84`, `figure_4_common.R:35`, `figure_4.Rmd:258`, `figure_6.Rmd`,
`figure_3_supplement_{1,2,3}.Rmd`. **Hazard**: in the figure-4 family the label string is used as a
key (the N_ch exclusion in `figure_4_common.R`), so a blind rename can silently un-drop a column that
must not be drawn. `figure_4_rays.Rmd:83-84` has the crossing BACKWARDS in prose (calls the on-disk
arm interval-naive); the dispatcher case list is authoritative.

**Count.** The implemented and measured set is now EIGHT (`LSE`, `ILSE`, `NR`, `INR`, `R`, `MR`,
`VR`, `IR`); figure 3 draws all eight, the grid runs seven. The abstract still says "seven-rung",
scoped to the sweep — open, see `00_abstract.tex` COUNT OPEN note.

## The ladder, for the Methods presentation

| Conditioned on | Members | Conductance model |
|---|---|---|
| the gating fluctuations are not modelled | `LSE`, `ILSE` | the deterministic mean current only, at the sampling instant (`LSE`) or averaged over the interval (`ILSE`) |
| no endpoints | `NR`, `R` | instantaneous; the averaging is ignored |
| one endpoint (the start) | `MR`, `VR` | interval-mean given the initial state (`VR` uses the residual variance) |
| endpoints not distinguishable | `INR` | interval-mean; with no update the two endpoint columns coincide |
| two endpoints (the boundary) | `IR` | interval-mean given both boundary states; interior marginalized |
| the full trajectory | (exact) | intractable; the stochastic simulation supplies it as ground truth |

`IR` is the top rung below intractability. **Superseded 2026-08-06**: the body no longer walks a
single chain `LSE → NR → R → IR`. Figure 1 walks the 3×2 grid, three levels of the gating description
crossed with the window setting: `LSE`/`ILSE`, `NR`/`INR`, `R`/`IR`. `MR` and `VR` split the R → IR
step and live in a supplement (`program.md` §1).

`INR` gets its own row rather than sitting with `MR` and `VR`, because it does not condition on one
endpoint *by choice*: without a gain the question does not apply to it, which is the degeneracy in the
lattice table above. `NMR` is not on this ladder at all — it is `INR` with the interval variance
dropped, so it does not model the conductance the row claims.

## The regions of the usage map

Named here so the figure, its caption and the text do not drift into three vocabularies. **The
numbering runs bottom to top in instrumental noise at fixed channel count**, which is the direction
the built figure uses (Luciano, 2026-07-28: "la numeración es con 0 donde falla macroir y 4 donde
fallan todas"). Any prose that counts them the other way is wrong.

| Region | Name | What it means |
|---|---|---|
| **0** | closure fails | the Gaussian closure itself is misspecified, so even `IR` reports a wrong error bar. The one boundary with a *negative* slope, because instrumental noise makes the emission more Gaussian |
| **1** | `i` + `N_ch` + honest CI | the likelihood delivers the unitary conductance, the channel count and a calibrated interval; least squares reports one several times too narrow |
| **2** | rates + honest CI | the amplitude pair is no longer separable and the amplitude needs a prior; the likelihood still buys a calibrated interval |
| **3** | LSE is par | least squares is calibrated too and buys the same answer more cheaply |
| **4** | nothing estimable | no useful information left |

Call it a **usage map**, and say in the caption that it is a **concept map, not a phase diagram**: its
boundaries are level sets of continuous diagnostics, so a looser criterion moves each of them by up to
a decade in noise without changing the layout. Four documents currently forbid drawing a hard line on
the design plane; that sentence is what reconciles them with the figure, and it should be quoted
rather than re-derived.

The boundary criterion is a distortion of **1.15**, read as a variance ratio, which is a 7% error on
the reported standard deviation. Thresholds are owned by `machinery.md` §8, which must adopt this
value rather than continue to propose 1.1.

## Scoping the term "boundary state"

A **boundary state** is the pair (i₀, iₜ) of the channel's states at the two ends of an acquisition
interval. Chaining intervals, the end state of one is the start state of the next, so there is a single
state variable per junction: the filter conditions on the interval boundaries and leaves the trajectory
between them free. This is static condensation, or the spatial Markov property, on the time axis.

It is **not** a transition state in the mechanistic sense (a short-lived conformational intermediate).
Avoid the word *transition* anywhere near it (`project_boundary_state_naming`).

## ~~Open item: the `MacroINR` bridge~~ CLOSED 2026-07-29, re-grounded 2026-07-31

The bridge `INR = MacroINR` is **correct**, and the argument no longer rests on a parse. Both readings
of the `I` name the same algorithm, because without an update the boundary conditioning is
unobservable; see the lattice table above. The suspicion recorded here until 2026-07-29 (that a
published name carrying an `I` would misname a start-conditioned method) assumed the two readings were
in conflict. They are not.

What the 2026-07-29 closure got wrong is a different thing: it identified the *implementation in this
repo* with `MacroINR`. That implementation was missing the `N·ms` term and is now called `NMR`. The
bridge holds for `INR`, the fixed one.
