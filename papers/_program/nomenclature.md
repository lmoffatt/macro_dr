# Nomenclature: naming the methods

> Updated: 2026-07-28. Shared across the two papers; cited, never restated. Settled items graduate
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
  knob settings as the dropped `NMR` and is distinguished only by that third flag. **It has no
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
| `MR`  | true      | 1   | total |
| `VR`  | true      | 1   | residual |
| `IR`  | true      | 2   | residual (+ boundary gain) |

`NMR` (false, 1, total) was the sixth row until 2026-07-28 and is **dropped from the program**
(`decisions.md` §2). Recomputed on the freeze it is numerically indistinguishable from `NR`, so it
occupied a lattice cell without measuring anything `NR` does not. The row is kept in
`../1_method/decisions/D-4_ranking_verdict.md` §5 as the evidence for the drop.

The suffix (`N` / `R`) is the occupancy axis: non-recursive or recursive. The prefix is the
conductance axis: none for the instantaneous conductance, `M` for the mean conductance, `I` for the
interval-conditioned (boundary) mean conductance.

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

- `MR → VR` replaces the total variance with the residual variance and changes nothing else.
- `VR → IR` adds the boundary cross-covariance term N·γᵀΣγ in the **gain** and changes nothing else.

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

## The ladder, for the Methods presentation

| Conditioned on | Members | Conductance model |
|---|---|---|
| the gating fluctuations are not modelled | `LSE` | the deterministic mean current only |
| no endpoints | `NR`, `R` | instantaneous; the averaging is ignored |
| one endpoint (the start) | `MR`, `VR` | interval-mean given the initial state (`VR` uses the residual variance) |
| two endpoints (the boundary) | `IR` | interval-mean given both boundary states; interior marginalized |
| the full trajectory | (exact) | intractable; the stochastic simulation supplies it as ground truth |

`IR` is the top rung below intractability. **The body walks `LSE → NR → R → IR`**, a monotone ladder
of cost; `MR` and `VR` split the R → IR step and live in a supplement (`program.md` §1).

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

## Open item: the `MacroINR` bridge

`decisions.md` records: "Published-name bridge: IR = MacroIR, NMR = MacroINR." That looks wrong: `NMR`
(≡ `MNR`) runs at av = 1, which is `M` (start-conditioned), not `I`. A published name carrying an `I`
would misname it, unless `MacroINR` in Communications Biology denoted something else. Check against the
Comm Biol text before the family table is written; if the bridge is wrong, fix it in `decisions.md`
first.
