# The program: two papers, split on the macro/micro boundary

> Updated: 2026-07-28. The three-paper split decided 2026-07-20 lasted three days; papers 1 and 2 were
> **merged on 2026-07-23** and the reasoning is in `decisions.md` §1. This file is rewritten around
> the surviving structure. The old three-row table is in git history.
>
> Owns: the two-paper map, what each owns and excludes, and the publication order. Owns no science.

## 1. The axis

The method makes **two Gaussian closures**. Each paper owns one, and the question prior to both
belongs to the macro paper because it is the question its anchor method answers.

| | The question | What is at stake | Control variables | Methods |
|---|---|---|---|---|
| **macro** | Do you need a likelihood at all, and if so what must it condition on? | the **interval-likelihood closure** (conductance over one interval → Gaussian), plus the prior question of whether the gating signal is above the instrumental noise | N_ch and instrumental noise, **faceted by N_ch**; interval 1 → 0.01 τ | body: `LSE`, `NR`, `R`, `IR`; supplement: `MR`, `VR`, `NMR` |
| **micro** | Where does the Gaussian occupancy closure break, and where is the exact solver worth its cost? | the **macro closure** (exact discrete occupancy distribution → Gaussian); degrades at few channels | the **effective fluctuating count** N_ch·p(1−p), not N_ch alone (§5) | `micro_R`, `micro_IR` against macro `R`, `IR` |

**The 1|2 cut was the artificial one.** It separated mid from high N_ch and the single-channel noise
scale from the fraction-of-total, and both distinctions dissolve once the figures facet by N_ch: at
fixed N_ch the gating variance is a fixed number, so the two noise conventions are the same variable
relabelled. The macro/micro cut is the real one, because the closure genuinely changes.

**Both papers have the same shape**: a frontier saying which is the cheapest method that still
reports its own uncertainty honestly, in that regime. That is what lets the validation machinery be
written once, in the macro paper, and cited by the micro one.

**The ladder, as the macro paper walks it**, is a monotone ladder of cost, and the reader enters at
whichever rung they are on today:

`LSE` (fit the mean, discard the fluctuations) → `NR` (use the gating variance, no filter) →
`R` (recursive, instantaneous sample) → `IR` (condition on the interval).

Reading LSE against IR answers "do you need a likelihood?". Reading R against IR answers "what does
conditioning on the interval buy?". The supplement carries `MR` and `VR`, which split the R → IR step,
and `INR` (the published `MacroINR`). **What `INR` shows is open** (2026-07-31): the near-identity with
`NR` that used to be quoted here as the finding was measured on `NMR`, a build missing the `N·ms`
interval-variance term, so it says nothing about the method. The re-run dispatched 2026-07-31 decides
it (`nomenclature.md`, the `NMR` entry).

## 2. The channel-number axis, partitioned

Low → **micro**. Everything above it → **macro**. The macro paper now spans the whole range where the
Gaussian closure holds, which is what the merge bought.

Today: the macro figures run N_ch 5 to 10⁴ for IR and 10 to 10⁴ for the rest; micro has cells at
5, 10, 20, 100.

**[Q] The macro floor is not fixed** (Luciano, 2026-07-28: to be revisited, "probablemente el piso sea
1 o 2 para ver la frontera"). This is not a grid parameter: the region map's lowest boundary, where
the Gaussian closure itself fails, is currently **extrapolated** below N_ch 10 and the figure's own
source note asks for N_ch 2 and 5 at noise 0.1 to 10, ten cells, to pin it. Whether that vertex is a
measured result or a dashed extrapolation is decided here.

## 3. Publication order: macro, then micro

The macro paper is the eLife shot and is nearly drafted; the micro paper follows and cites it for the
machinery. Nothing in the macro paper depends on the micro one except the attribution of IR's own
few-channel degradation, which one or two annotated cells already supply (§6).

## 4. What each paper owns, excludes, and still needs

### The macro paper
- **Owns:** the validation machinery (written here, cited by the micro paper); the endpoint ladder;
  the interval closure's failure and its mechanism; and **the usage map across the whole noise axis**,
  which is the region map (`../1_method/decisions.md`, "The figure set").
- **The anchor is least squares** (`decisions.md` §1). The frame is that the field fits the
  deterministic mean and discards the fluctuations, that the methods which use them have almost no
  uptake, and that the reason is a missing validity criterion. Not "here is my algorithm".
- **Literature positioning:** `R` carries the recursive lineage (Moffatt 2007; Münch 2022, a published
  Bayesian Kalman filter in the target journal, so the abstract must position against it).
  **The gap claim is about validity, never about absence:** the temporal correlation has carried
  kinetics since 1973 and a covariance likelihood predates MacroR by three years
  (`decisions.md` §6; `docs/bibliography/temporal_correlation_and_AR_errors_2026-07-28.md`).
- **The VR mechanism is a supplement, not the headline** (2026-07-23). It still earns its rung: `VR`
  keeps MR's boundary-free update and swaps the total per-start-state conductance variance for the
  residual one, turning an algebraic claim into a measured one, and it came out over-confident and
  more so than MR, exactly as predicted. **Do not state it as "MR → VR changes only the variance,
  VR → IR only the gain"**: the predictive variance divides the gain, so the variance step moves the
  update too, and `MR` and `IR` predict the *same* observable variance from the same state, so the
  whole MR-to-IR difference is the gain. Measured in `../1_method/figures_build_plan.md` §F1-2.
- **Excludes:** micro as a subject. Nothing else is excluded; `INR`, `MR` and `VR` are supplement members (`decisions.md` §2).
- **Still needs:** LSE and NR re-run at n_sims 10⁴ so they can share a panel with the band-A cells;
  the region map written into Results; the recording-condition overlay lifted out of the notebook
  into a citable Methods table; the sign convention verified against the producer.

### The micro paper
- **Owns:** where the macro closure breaks, and the micro cost frontier.
- **Status:** runs on disk (§6).
- **Must have the same shape as the macro paper.** "The Gaussian breaks below N_ch = X" is a number, not a
  paper. The paper is: micro is expensive, here is where it stops being worth paying for.

## 5. The micro paper's control variable is not N_ch

If low P_open behaves like reduced N_ch (recorded as a planned experiment in the design notes), then
what decides the microscopic regime is the number of channels that actually fluctuate, of order
**N_ch·p(1−p)**, not N_ch. At the P_open = 0.5 fixed everywhere today the two coincide up to a
factor and the difference is invisible. The moment the micro paper moves off 0.5 they separate, and
"the Gaussian breaks below N_ch = X" becomes false as stated, because X depends on P_open.

Choose the effective count as the axis from the start and the micro result generalizes; choose bare
N_ch and it needs correcting later.

## 6. Data on disk, per paper

| Paper | Where | What |
|---|---|---|
| macro, likelihood arm | `projects/eLife_2025/figures/data/1c2ae6f`, `87889e6` | the freeze carries IR everywhere and R/MR/NR at noise 0.1 only (plus R at 100); noise 1 and 10 for those three are on `87889e6`, with the now-dropped NMR. `433ed13` is the numerical-Fisher demo and **no paper number is quoted from it** |
| macro, LSE arm | `.../82b956f` | the LSE runs, dispatched 2026-07-20, **at n_sims 1000**, so they cannot share a panel with the 10⁴ cells until re-run |
| micro | `.../87889e6` | `micro_IR` at N_ch 5 (nsim 100/1000/10⁴), 10 (100/1000/10⁴), 20 (100); `micro_R` at 5 (100/10⁴), 10 (100/10⁴), 20 (100), 100 (100). **All at noise 0.1 only.** |

**The micro attribution anchor for the macro paper already exists:** `micro_IR`, N_ch = 10, nsim = 10000,
noise 0.1 — the macro paper's floor, its canonical n_sims, its canonical noise. Its job is to attribute IR's
own low-N_ch degradation: micro_IR keeps the exact occupancy distribution *and* the interval
treatment, so if it is calibrated at 10 channels where macro IR is not, the degradation belongs to
the occupancy closure and therefore to the micro paper. **One or two annotated cells, not a column** —
a full micro column re-opens the roster question the split just closed.

Two cautions. The anchor exists at **one noise level only**; if the macro paper makes its few-channel claim
across the noise fan, more cells are needed. And **never pair a 100-sim micro cell with a 10⁴-sim
macro cell**: the distortion scalars carry a Jensen bias in n_sims and the comparison will manufacture
a difference. Only `micro_IR` at N_ch 10 / nsim 10⁴ pairs cleanly.

**Correction owed:** the decision log records `87889e6` as "micro, out of scope". It now holds both
the micro runs and the macro D-0 fill, and the macro paper will cite cells from it. Multi-commit provenance
is already accepted (each CSV self-stamps its engine hash), so this is a bookkeeping fix, not a
policy change.

## 7. The flank, and what now closes it

The original problem: the paper lived only in the regime that favours the gating-aware likelihoods,
so it read as though IR won by construction. The 2026-07-20 answer was a **stated scope** in band
terms, on the ground that a least-squares arm was somebody else's paper.

**That answer is retired.** The merge puts the least-squares arm in this paper and extends the noise
axis across the crossover, so the flank is closed by **the comparison itself and by the region map**,
which shows the regions where LSE is as good as IR and the region where IR itself fails. A scope
declaration is still owed, but it is now about what the model is (two states, one open probability,
one jump protocol, simulation) and not about which band the paper dares to enter.

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

- **[Q]** The macro paper's N_ch floor (§2), which decides whether the region map's lowest vertex is
  measured or extrapolated. Leaning: 1 or 2, to see the frontier.
- **[Q]** The micro paper's control variable: effective fluctuating count or bare N_ch (§5).
- ~~**[Q]** `VR`'s name.~~ **CLOSED 2026-07-28: `VR` keeps the letter** (`decisions.md` §2). The
  collision caveat survives the closure: `V` still clashes with the cut Taylor variants `MRV`/`IRV`
  and with the engine flag `taylor_variance_correction`, so the March Taylor data must be deleted and
  Methods must say in one sentence that this `V` is not that `V`.
- ~~**[Q]** Venue per paper.~~ **CLOSED: eLife for the macro paper**, Biophysical Journal as the
  fallback (`decisions.md` §1). The micro paper's venue is open and is not on any critical path.
- ~~**[Q]** Narrow the macro paper to MacroIR and its validation, with the rest as an addendum.~~
  **NEVER OPEN. Withdrawn by the author on 2026-08-11, in those terms: a temporary weakness, and he
  went back to the merged version.** Raised in the audio of 2026-08-10 at 14.32.13 (*"creo que tengo
  que reorientar el paper totalmente ... hasta podría ser solo de macro IR ... y lo otro lo pongo
  como una adenda"*, with *"no sé"* on either side of it) and abandoned the same night. **The merge
  of `decisions.md` §1 stands, unamended.**
  This tombstone exists because the audio does not carry its own retraction and someone extracting
  that batch cold will read a reorientation that was never adopted. It also records what the audio
  was actually reaching for, which is sound and is not the narrowing: the 2007 flip state was
  readable from a delay that did not depend on fitting a model, and the mechanisms at stake now do
  depend on model validity, so there is no intrinsic signature in the current to fall back on. That
  point is written, in the Introduction's opening paragraph on the P2X$_2$ case
  (`../1_method/docs/manuscript-drafts/sections/01_introduction.tex:254-263`, committed at 23:25 the
  same night) and again at `05_discussion.tex:229`. It landed inside the broad paper, which is why
  the broad paper never had to give anything up for it.
