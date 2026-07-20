# Paper 3 — the multinomial boundary (stub)

> Updated: 2026-07-20. A stub, not a folder. It becomes `3_micro/` when it starts drafting.
> Its place in the program: `program.md` §1. Nothing here is settled unless `decisions.md` says so.

## The question

**Where does the Gaussian occupancy closure break, and where does the exact solver stop being worth
its cost?**

This is the macro closure: the multinomial distribution of channel occupancies replaced by a
multivariate Gaussian, valid for large N_ch by the central limit theorem, degrading at few channels.
Paper 1 owns the *other* closure (the interval likelihood). The two are siblings, one paper each.

## Roster

`micro_R` and `micro_IR`, against macro `R` and `IR`. The micro members keep the exact multinomial
occupancy; the comparison isolates the closure.

## The control variable is probably not N_ch

If low P_open behaves like reduced N_ch (recorded as a planned experiment in the design notes), then
what decides the multinomial regime is the number of channels that actually **fluctuate**, of order
**N_ch·p(1−p)**, not N_ch. At the P_open = 0.5 fixed everywhere today the two coincide up to a factor
and the difference is invisible. The moment this paper moves off 0.5 they separate, and the sentence
"the Gaussian breaks below N_ch = X" becomes false as written, because X depends on P_open.

Pick the effective count from the start and the result generalizes. Pick bare N_ch and it needs
correcting later. **[Q], and it should be settled before the analysis, not after.**

## It must have the same shape as papers 1 and 2

"The Gaussian breaks below N_ch = X" is a number, not a paper. The paper is: **micro is expensive,
here is where it stops being worth paying for** — a cost frontier, the same shape as the other two.
That sameness is what makes the three a program rather than three loose papers, and it is what lets
the validation machinery be written once in paper 1 and cited here.

## Data on disk

`projects/eLife_2025/figures/data/87889e6`, **all at noise 0.1**:

| Method | N_ch (n_sims) |
|---|---|
| `micro_IR` | 5 (100 / 1000 / 10⁴), 10 (100 / 1000 / 10⁴), 20 (100) |
| `micro_R` | 5 (100 / 10⁴), 10 (100 / 10⁴), 20 (100), 100 (100) |

Two cautions carried from `provenance.md`. The grid sits at **one noise level only**. And the n_sims
column is ragged: **never pair a 100-sim cell with a 10⁴-sim cell**, in either direction — the
distortion scalars carry a Jensen bias in n_sims and the comparison manufactures a difference.

**Bookkeeping correction owed:** the old decision log records `87889e6` as "micro, out of scope". It
now holds both these micro runs and the macro D-0 fill (NR, NMR, R, MR at noise 0.1/1/10), and paper 1
will cite one cell from it. Multi-commit provenance is already accepted (each CSV self-stamps its
engine hash), so this is a bookkeeping fix, not a policy change.

## What it owes paper 1

**One attribution anchor, already on disk:** `micro_IR` at N_ch = 10, n_sims = 10⁴, noise 0.1. Paper 1
runs down to 10 channels and reports IR's own degradation there, but cannot attribute it: all four of
its methods share the macro closure, so it cannot tell whether that degradation is its own subject
(the interval closure) or this paper's (the occupancy closure). micro_IR keeps the exact occupancy
*and* the interval treatment, so one cell decides it.

In paper 1 this is **one or two annotated cells, not a column** — a full micro column would re-open
the roster question the three-way split just closed. Paper 1 forward-cites this paper as a companion
for the full boundary, which is comfortable given that papers 2 and 3 are planned back to back.

## Open

- **[Q]** The control variable (above).
- **[Q]** The N_ch ceiling, decided jointly with the other two papers (`program.md` §2). Proposal: ≤ 20.
- **[Q]** Whether the noise axis needs extending here at all, or whether noise 0.1 is the right single
  slice for a paper about the channel-number boundary.
- **[Q]** Venue.
