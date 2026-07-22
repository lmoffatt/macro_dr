# The program: three papers on one axis

> Updated: 2026-07-20. **Decided this day** (working session + audios 10:54–10:59). Supersedes the
> single-paper plan, whose §0 goal sentence now describes paper 1 only.
>
> Owns: the three-paper map, what each paper owns and excludes, the publication order and why it is
> safe. Owns no science. Every fact it needs is owned in `_program/` or in a paper's own folder.

## 1. The axis

The method makes **two Gaussian closures**, and there is **one question prior to both**. Each paper
owns exactly one of the three. That is the whole partition, and it is why these are three papers and
not three slices of one.

| | The question | What is at stake | Control variable | Methods |
|---|---|---|---|---|
| **1 · method** | Given that you will compute a likelihood, what must it condition on? | the **interval-likelihood closure** (conductance over one interval → Gaussian); degrades at short Δ | noise on the **single-channel** scale (ν, `axes.md`); N_ch from 10; interval 1 → 0.01 τ | `R`, `MR`, `VR`, `IR` |
| **2 · map** | Do you need a likelihood at all, and how much must you pay? | nothing — this is **prior to any closure**: is the gating signal above the instrumental noise? | instrumental noise as a **fraction of total** noise; N_ch 10² → 10⁶ | `LSE`, `NR`, `R`, `IR` |
| **3 · micro** | Where does the Gaussian occupancy closure break, and where is the exact solver worth its cost? | the **macro closure** (multinomial → Gaussian); degrades at few channels | the **effective fluctuating count** N_ch·p(1−p), not N_ch alone (§5) | `micro_R`, `micro_IR` vs macro `R`, `IR` |

Papers 1 and 3 are siblings, one per closure. Paper 2 sits a level above both.

**All three have the same shape**: a cost frontier saying which is the cheapest method that still
reports its own uncertainty honestly, in that regime. That sameness is what makes this a program
instead of three loose papers, and it is what lets the validation machinery be written once (paper 1)
and cited twice.

## 2. The channel-number axis, partitioned

Low → **paper 3**. Middle → **paper 1**. High → **paper 2**. Each paper lives where its question is
live and the others have nothing to say.

**Decide the three ranges together, once, here** — not per paper. Today: paper 1 runs 10…10⁴,
paper 2 is specified 10²…10⁶, paper 3 has cells at 5, 10, 20, 100. Left to drift, this produces
either a gap or an undecided overlap, and in a three-paper program that shows.

**[Q] The three ranges are not yet fixed.** Proposal: paper 3 ≤ 20, paper 1 10…10⁴, paper 2 ≥ 10².

## 3. Publication order: 1, then 2 and 3 back to back

Schedule says 1 → 2 → 3: paper 1 is nearly drafted, paper 2 has runs in flight, paper 3 needs its
analysis. Structure would say 1 → 3 → 2, keeping the two closure-siblings together and letting the
applied paper close by citing both boundaries.

**The conflict dissolves on one decision: paper 2's floor of 10² channels.** At 10² and above the
macro closure is comfortable, so paper 2 never touches the multinomial boundary and therefore does
not depend on paper 3. It depends only on paper 1's interval result. **That floor is load-bearing for
the publication order, not a grid parameter.** If it drops, the order must be reconsidered.

## 4. What each paper owns, excludes, and still needs

### Paper 1 — method
- **Owns:** the validation machinery (written here, cited by 2 and 3); the endpoint ladder; the
  interval closure's failure and its mechanism.
- **The mechanism result, which is new:** `VR` inserts a third rung between `MR` and `IR` that keeps
  MR's boundary-free update and swaps the total per-start-state conductance variance for the residual
  one, which turns an algebraic claim into a measured one. This is what earns the paper its novelty
  independently of the map. **Do not state it as "MR → VR changes only the variance, VR → IR only the
  gain"** (the wording used here until 2026-07-21): the predictive variance divides the gain, so the
  variance step moves the update as well, and `MR` and `IR` predict the *same* observable variance
  from the same state, so the whole MR-to-IR difference is the gain. Measured in
  `../1_method/figures_build_plan.md` §F1-2.
- **Excludes:** `LSE`, `NR`, `NMR`. NR is not dropped, it **moves to paper 2**, where a cheap
  non-recursive likelihood is exactly the point. NMR is dropped outright: no literature attribution,
  no mechanistic role.
- **Literature anchor after the exclusions:** `R` carries it (Moffatt 2007; Münch 2022, a published
  Bayesian Kalman filter). The frame is "the published recursive filters treat each sample as
  instantaneous; here is what conditioning on the interval buys, and which half of the mechanism
  does it." It does not need NR.
- **Still needs:** `VR` implemented and run (small, if the residual form is computable from existing
  Qdtm fields as claimed — verify against current code, the note claiming it has a known-inverted
  verdict); the micro attribution anchor (§6); the scope declaration in band terms (§7).

### Paper 2 — usage map
- **Owns:** the usage map across the gating-noise crossover; the LSE arm.
- **Status:** LSE implemented end to end in ops; runs dispatched 2026-07-20.
- **[!] The dispatch does not match this roster.** It ran `nonlinearsqr` plus
  `macro_{IR,R,MR,NMR}`, i.e. it includes MR and NMR (not in paper 2) and omits NR (in paper 2). It
  was launched before this split existed. Reconcile before the map is drawn.

### Paper 3 — multinomial boundary
- **Owns:** where the macro closure breaks, and the micro cost frontier.
- **Status:** runs on disk (§6).
- **Must have the same shape as 1 and 2.** "The Gaussian breaks below N_ch = X" is a number, not a
  paper. The paper is: micro is expensive, here is where it stops being worth paying for.

## 5. Paper 3's control variable is not N_ch

If low P_open behaves like reduced N_ch (recorded as a planned experiment in the design notes), then
what decides the multinomial regime is the number of channels that actually fluctuate, of order
**N_ch·p(1−p)**, not N_ch. At the P_open = 0.5 fixed everywhere today the two coincide up to a
factor and the difference is invisible. The moment paper 3 moves off 0.5 they separate, and
"the Gaussian breaks below N_ch = X" becomes false as stated, because X depends on P_open.

Choose the effective count as the axis from the start and paper 3's result generalizes; choose bare
N_ch and it needs correcting later.

## 6. Data on disk, per paper

| Paper | Where | What |
|---|---|---|
| 1 | `projects/eLife_2025/figures/data/1c2ae6f`, `433ed13`, `87889e6` | the band-A grid; `87889e6` also holds the D-0 macro fill (NR, NMR, R, MR at noise 0.1/1/10) |
| 2 | `.../82b956f` | the LSE runs, starting 2026-07-20 |
| 3 | `.../87889e6` | `micro_IR` at N_ch 5 (nsim 100/1000/10⁴), 10 (100/1000/10⁴), 20 (100); `micro_R` at 5 (100/10⁴), 10 (100/10⁴), 20 (100), 100 (100). **All at noise 0.1 only.** |

**The micro attribution anchor for paper 1 already exists:** `micro_IR`, N_ch = 10, nsim = 10000,
noise 0.1 — paper 1's floor, its canonical n_sims, its canonical noise. Its job is to attribute IR's
own low-N_ch degradation: micro_IR keeps the exact multinomial occupancy *and* the interval
treatment, so if it is calibrated at 10 channels where macro IR is not, the degradation belongs to
the occupancy closure and therefore to paper 3. **One or two annotated cells, not a column** —
a full micro column re-opens the roster question the split just closed.

Two cautions. The anchor exists at **one noise level only**; if paper 1 makes its few-channel claim
across the noise fan, more cells are needed. And **never pair a 100-sim micro cell with a 10⁴-sim
macro cell**: the distortion scalars carry a Jensen bias in n_sims and the comparison will manufacture
a difference. Only `micro_IR` at N_ch 10 / nsim 10⁴ pairs cleanly.

**Correction owed:** the decision log records `87889e6` as "micro, out of scope". It now holds both
the micro runs and the macro D-0 fill, and paper 1 will cite a cell from it. Multi-commit provenance
is already accepted (each CSV self-stamps its engine hash), so this is a bookkeeping fix, not a
policy change.

## 7. The scope declaration paper 1 must carry

The reason LSE and the extended noise axis were added at all was that the paper lived only in the
regime that favours the gating-aware likelihoods. **The split does not fix that; a stated scope
does.** Paper 1 must say, in its own words, that it characterizes the gating-dominated regime, define
that regime by the crossovers in `axes.md`, and name the companion paper for the rest. A stated scope
is defensible and can be elegant. An unstated one is the flank that started this whole revision.

## 8. Citation directionality (this is what keeps four folders from becoming four copies)

The chronic failure of the previous single pack was duplication: the scope call written out in full in
four documents, the evidence formulas in five, the diagnostics list in three, the ranking table in
three. Four layers multiply the surfaces on which the same fact can be copied.

"One topic, one owner" no longer suffices alone. Add the direction:

- a paper **may cite** `_program/`;
- `_program/` **never cites** a paper;
- papers **do not cite each other** in the planning layer (only as companion papers in the finished
  manuscripts).

**If `_program/` finds itself citing a paper, the fact was filed in the wrong place.** That is the
whole test, and it is mechanical.

## 9. Open

- **[Q]** The three N_ch ranges (§2).
- **[Q]** Paper 3's control variable: effective fluctuating count or bare N_ch (§5).
- **[Q]** `VR`'s name. The letter reopens a third axis in a naming scheme that is currently
  compositional (prefix = conductance, suffix = occupancy). Cost is now small — with NR and NMR gone
  the third axis distinguishes only MR from VR — but `V` collides with the cut Taylor variants
  `MRV`/`IRV` and with the engine flag `taylor_variance_correction`. Viable if the March Taylor data
  is deleted and Methods states plainly that this V is not that V.
- **[Q]** Venue per paper. The old folder name assumed eLife for a single paper; with three, the
  question is per paper and is open for all three.
