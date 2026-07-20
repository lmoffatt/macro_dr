# Nomenclature: naming the methods

> Updated: 2026-07-20. Shared across the three papers; cited, never restated. Settled items graduate
> to `decisions.md`. Which paper uses which method: `program.md` §1.
>
> Scope: what we call every method in the program, what the letters mean, and how to describe them so
> each name picks out exactly one method.

## The constraint

`IR` / `MacroIR` is in print (Communications Biology 2025, P2X2). The acronym is fixed; nothing here
renames it. What is open is the *gloss*, and the name of the one new member, `VR`.

## The structure has two levels, not one lattice

The earlier "one object with two knobs" description was true of the five likelihoods and became false
once least squares joined the comparison. The honest structure is a **root question with a ladder
hanging from it** (`program.md` §1):

- **Root: do you model the gating fluctuations at all?** `LSE` (classical nonlinear least squares on
  the mean current) answers no. In the engine it is `family_approximation = 2`; it carries the same two
  knob settings as `NMR` and is distinguished only by that third flag. **It has no rung and no
  compositional name** — it is a different kind of object, named for what it is.
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
| `NMR` | false     | 1   | total |
| `MR`  | true      | 1   | total |
| `VR`  | true      | 1   | residual |
| `IR`  | true      | 2   | residual (+ boundary gain) |

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

## `VR`: a provisional name, and a hazard

**`VR` is provisional** (`program.md` §9). Two cautions:

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
| no endpoints | `NR`, `R` | instantaneous; the averaging is ignored |
| one endpoint (the start) | `NMR`, `MR`, `VR` | interval-mean given the initial state (`VR` uses the residual variance) |
| two endpoints (the boundary) | `IR` | interval-mean given both boundary states; interior marginalized |
| the full trajectory | (exact) | intractable; the stochastic simulation supplies it as ground truth |

`IR` is the top rung below intractability. Paper 1 walks the recursive spine `R → MR → VR → IR`; the
non-recursive members `NR`, `NMR` and the off-lattice `LSE` are paper 2's.

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
