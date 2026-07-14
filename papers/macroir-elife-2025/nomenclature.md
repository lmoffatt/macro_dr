# Nomenclature: naming the likelihood family

Working discussion, opened 2026-07-14. Settled items graduate to `02_decision_log.md`.

Scope: what we call the five macroscopic likelihood algorithms, what the letters
in their names mean, and how to describe them in the manuscript so that each name
picks out exactly one member.

## The constraint

`IR` / `MacroIR` is already in print (Communications Biology 2025, P2X2). The
acronym is fixed. Nothing below proposes renaming it. What is still open is the
*gloss*: the words we use in the manuscript to say what the letters stand for.

## The family, as implemented

Verified against `projects/eLife_2025/runs/run-20260418-182133/script.macroir`,
which builds all members from two flags. `recursive_approximation` says whether
the occupancy covariance is propagated between intervals;
`averaging_approximation` (below, `av`) says how the single-channel conductance
is treated within an acquisition interval.

| Label | recursive | av |
|-------|-----------|-----|
| `NR`  | false     | 0   |
| `R`   | true      | 0   |
| `MNR` | false     | 1   |
| `MR`  | true      | 1   |
| `IR`  | true      | 2   |

So the names are compositional. The suffix (`N` / `R`) is the occupancy axis:
non-recursive or recursive. The prefix is the conductance axis: no prefix for the
instantaneous conductance, `M` for the mean conductance, `I` for the
interval-conditioned mean conductance.

## What the letters mean

- `M` = **mean conductance**. The conductance is averaged over the acquisition
  interval, conditioned on the channel's state at the interval's *start*. This is
  the K-vector (K = number of microscopic states)

  (γ̄₀)ᵢ₀ = Σᵢₜ Pᵢ₀→ᵢₜ(t) · Γ̄ᵢ₀→ᵢₜ

  where Γ̄ᵢ₀→ᵢₜ is the mean interval-averaged conductance of a channel that starts
  in state i₀ and ends in state iₜ, and P(t) = exp(Qt) is the transition matrix.

- `I` = **interval-conditioned mean conductance**. The same interval average, but
  conditioned on the states at *both* ends of the interval, so the object stays
  K×K (the pair (i₀, iₜ), which we call the **boundary state**) instead of
  collapsing to a K-vector.

The single concrete difference: `MR` drops the boundary cross-covariance term
N·γᵀΣγ that `IR` keeps. That term is why `MR` overestimates the predicted
variance and `IR` does not.

Useful sanity check: `av` literally counts the conditioned endpoints. Zero
endpoints (instantaneous conductance), one endpoint (the start), two endpoints
(the boundary). The flag is self-documenting.

## Is the `I` defensible?

Yes. An interval is determined by its two endpoints, so "conditioned on the
interval" can legitimately be read as "conditioned on both of its endpoints",
in contrast to `M`, which conditions on the single instant at which the interval
begins. The notation is coherent, and it is not a rationalization after the fact.

But it is a convention the reader has to be handed, and there are two ways it
misfires if we leave it bare.

## Problem 1: the abstract's descriptor is not unique

The current abstract says:

> ... a family of Gaussian likelihood approximations (non-recursive, recursive,
> interval-averaged, and their combinations) ... Only the interval-averaged
> recursive likelihood (MacroIR) stays calibrated across the practical regime ...

`MR` is interval-averaged (av = 1) and recursive. So "the interval-averaged
recursive likelihood" describes `MR` just as well as it describes `IR`. The phrase
does not identify MacroIR. This is the naming problem showing up in the one place
where it costs the most.

The unique identifier for `IR` is not that it averages over the interval (three of
the five members do that), it is that it conditions on **both** interval
endpoints.

## Problem 2: "interval-conditioned" can be read as an overclaim

"Conditioned on the interval" has a second reading: conditioned on *everything
inside* the interval, that is, on the full trajectory. That would be the exact
likelihood. `IR` does not do this. It conditions on the two endpoints and
marginalizes the interior analytically with closed-form moments.

A reader who takes the second reading will hear a claim of exactness and then find
it is not there. This is the more damaging misreading of the two, and it is easy to
avoid.

## Resolution

Keep the acronyms (`IR` is fixed anyway, and the compositional scheme is sound).
Fix the prose:

1. In the manuscript body, describe `IR` as **boundary-conditioned**, not as
   "interval-conditioned". The word *boundary* does two things that *interval*
   does not: it says which part of the interval is conditioned on, and it makes
   explicit that the interior is left free and marginalized.

2. In the abstract, drop the descriptor that does not discriminate. Either name
   MacroIR without a gloss ("Only MacroIR stays calibrated across the practical
   regime"), or describe the family by its two axes so that the unique identifier
   falls out: MacroIR is the only member that conditions on both ends of the
   acquisition interval.

## Recommended presentation: the endpoint ladder

The family is a ladder in how much of the interval each member conditions on. This
is the clearest way to introduce it in Methods, and it makes the verdict feel
inevitable rather than empirical.

| Conditioned on | Members | Conductance model |
|---|---|---|
| no endpoints | `NR`, `R` | instantaneous; the averaging done by the acquisition is ignored |
| one endpoint (the start) | `MNR`, `MR` | interval-mean conductance given the initial state |
| two endpoints (the boundary) | `IR` | interval-mean conductance given both boundary states; interior marginalized |
| the full trajectory | (exact) | intractable; this is what the stochastic simulation supplies as ground truth |

MacroIR is the top rung below intractability. Stated this way, the reason it is the
sole survivor is structural, and the abstract needs no extra words.

## Scoping the term "boundary state"

Introduce it once, with a bounding sentence, because it has a near neighbour that
means something else.

A **boundary state** is the pair (i₀, iₜ) of the channel's states at the two ends
of an acquisition interval. Chaining intervals, the end state of one interval is
the start state of the next, so there is a single state variable per junction:
the filter conditions on the states at the interval boundaries and leaves the
trajectory between them free. This is the same device as static condensation, or
the spatial Markov property, applied on the time axis.

It is **not** a transition state in the mechanistic sense (a short-lived
conformational intermediate). Avoid the word *transition* anywhere near it.

## Open item: the `MacroINR` bridge

`02_decision_log.md` records: "Published-name bridge: IR = MacroIR, NMR =
MacroINR."

That looks wrong. `NMR` (equivalently `MNR`) runs at av = 1, which is `M` (mean
conductance, start-conditioned), not `I`. A published name carrying an `I` would
be a misnomer for it, unless `MacroINR` in Communications Biology denoted
something else, for instance a non-recursive member at av = 2, which does not
exist in the present grid.

To check against the Communications Biology text before the family table is
written. If the bridge is wrong, fix it in the decision log first.
