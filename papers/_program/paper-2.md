# Paper 2 — the usage map (stub)

> Updated: 2026-07-20. A stub, not a folder. It becomes `2_map/` when it starts drafting.
> Its place in the program: `program.md` §1. Nothing here is settled unless `decisions.md` says so.

## The question

**Do you need a likelihood at all, and how much do you need to pay?**

Prior to both Gaussian closures. Not "which approximation is adequate" but "is the gating signal above
the instrumental noise in the first place". That is why it is a separate paper and not a section of
paper 1.

## Roster

`LSE`, `NR`, `R`, `IR`. The cheap end (least squares on the mean, and a cheap non-recursive
likelihood) against the expensive end, with `R` and `IR` shared with paper 1 as the hinge that lets a
reader connect the two.

`NR` is here, not in paper 1. It was moved rather than dropped: a cheap non-recursive likelihood is
exactly this paper's subject, and it is the roster's link to Milescu 2005 / QuB.

## Axis

Instrumental noise as a **fraction of the total** noise, i.e. referred to the population gating scale;
N_ch from 10² to 10⁶.

**The 10² floor is load-bearing, not a grid parameter.** At 10² channels and above the macro closure
is comfortable, so this paper never touches the multinomial boundary and therefore does not depend on
paper 3. That is what makes the publication order 1 → 2 → 3 valid structurally and not merely
convenient (`program.md` §3). If the floor drops, the order must be reconsidered.

## The claim it must make, and the one it must not

The obvious half: if instrumental noise swamps the gating noise, you cannot recover kinetic rates.
Luciano's own words on it, 2026-07-20: *"es una cuestión absolutamente lógica… ya sabemos que tiene
que ser así"*. **A paper that opens on that opens weak.** It is a statement about the point estimate.

The non-obvious half, and the one this paper is for: **the reported uncertainty**. LSE can return a
serviceable point estimate together with a badly wrong confidence interval over a broad region, and
where that region begins and ends has never been measured for macroscopic currents. That is what the
validation machinery measures and nothing else in the literature does.

So: the claim is about the **CI**, not the point estimate. Written that way the paper stops being
obvious and becomes the one that serves the person using least squares today.

## Status

- **LSE is implemented end to end in ops**: `dispatch_figure_3_LSE.sh` (`family=2`),
  `figure_3_mle_LSE.macroir`, `figure_3_time_LSE.macroir`, `figure_lse_smoke.macroir`, and a
  `figure_1_plus_lse.macroir` that builds the sixth column. Only fig 1's diagnostic panel is blocked,
  waiting on the guard being routed for `family==2`.
- Design, seams and config constraints: `theory/macroir/notes/nonlinearsqr_lse_plan.md`. It is
  detailed and current; do not re-derive it.
- **Runs dispatched 2026-07-20** to dirac. Data lands in `projects/eLife_2025/figures/data/82b956f`.

## Open, and blocking the map

- **[!] The dispatch does not match this roster.** It ran `nonlinearsqr` plus `macro_{IR,R,MR,NMR}`:
  MR and NMR are in but do not belong to this paper, and NR is out but does. It was launched before
  the three-paper split existed. Reconcile before drawing anything.
- **[!] Two identical `dispatch_figure_3_G.sh` submissions** went out minutes apart with different
  `DEPEND` (117193, 117258) and the same grid. The filepath carries no job id, so both write the same
  16 directories. Check `squeue` and cancel one.
- The dispatched fill is `n_sims = 1000`; paper 1's band-A cells are 10000. **Do not pool them in one
  panel** — the distortion scalars carry a Jensen bias in n_sims and the difference will read as a
  regime effect (`provenance.md`).
- Two diagnostics do not mean for LSE what they mean for a likelihood, and both must be annotated
  wherever an LSE cell sits beside a likelihood cell: `r̄²_std ≡ 1` is a tautology for LSE
  (Σ r_std² = n identically), and `F = Var(score)` requires homoscedasticity, which LSE assumes and
  the exact simulator violates. The second is a **result**, not an artifact: the size of that break
  is the classical method's overconfidence, measured.
- **[Q]** Venue.
