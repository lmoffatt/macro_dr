# Paper 1 figures: the build plan

> Opened 2026-07-21. **Owns: the order of runs and edits, and the command for each.** It does not own
> what a figure *shows* (`results.md`), the eLife/visual rules (`../_program/figures_system.md`), or
> which run made which data (`../_program/provenance.md`). Cite those; do not restate them.
>
> Telegraphic on purpose. Every state claim below is dated and carries its file:line, because the
> whole point of this document is to stop the figures being rebuilt against a stale belief.

## 0. State on 2026-07-21 (verified, not declared)

- **VR's C++ has landed on `main`.** Flag type `uses_variance_form_aproximation` with
  `variance_total=0` / `variance_residual=1` (`legacy/qmodel_types.h:100-106`); selector implemented
  as the spec's option (A) OR-guard (`legacy/qmodel.h:4557-4558`), so IR cannot regress whatever its
  variance_form. `Qdtm` already carried `gtotal_ij`, `gmean_ij`, `gsqr_i`, so no new Qdt field was
  needed (`legacy/qmodel_types.h:1234-1235`).
- **A binary with VR already exists**: `build/gcc-release/macrodr_cli`, mtime 2026-07-21 14:29, built
  from `0ffbda7` (the last VR commit). Nothing below needs a fresh compile *for Figure 1*.
- **VR data, as of 19:05 on 2026-07-21:** the Figure-1 diagnostic dump
  (`figure_1_likelihood_diagnostic_VR.csv`, at the data root) and **three cloud files** in
  `0ffbda7/` (N_ch 100, 1000, 10000 at nsim 10000, noise 0.1). No `_battery_pool_G` or
  `_battery_sim_G` for VR anywhere yet, and no N_ch = 10 cell. See §4.
- **VR is not a DSL algorithm name.** It is the flag triple `recursive=true`, `averaging=1`,
  `variance_form=1` on `build_likelihood_function_with_family`
  (`src/cli/command_manager.cpp:1024-1032`). The token `macro_VR` exists only as (a) a shell case
  label and (b) the DSL axis label that lands verbatim in the CSV `algorithm` column and the filename
  (`projects/eLife_2025/ops/slurm/dispatch_figure_3_G.sh:141,169,172`).
- **Only one dispatcher can emit it**, the Gaussian-Fisher one. VR therefore arrives on the G anchor
  (`figure_3_G_*` filenames), not on the `433ed13` numeric basis that every current cross-algorithm
  panel reads.
- **The R side has no shared algorithm vocabulary.** The roster is duplicated verbatim in 11
  notebooks under `projects/eLife_2025/figures/paper/`, in two incompatible key forms: bare (`"NR"`,
  the dlik-fed notebooks: `figure_3.Rmd:51`, `figure_4.Rmd:39`, `figure_S4_S5.Rmd`) and prefixed
  (`"macro_NR"`, the battery-fed ones: `figure_2.Rmd:42`, `figure_S3.Rmd:53`, `figure_5_master.Rmd`).
  This is the debt `../_program/00_index.md:68` already books against `figures_system.md`.

### The three ways a VR row disappears without an error

1. **Label miss.** `algo = factor(ALGO_LAB[algorithm], levels = unname(ALGO_LAB))`
   (`figure_2.Rmd:67`): an unlisted key gives `NA`, ggplot drops the rows with a warning. This is the
   MNR/NMR failure class exactly.
2. **Glob miss.** `figure_5_master.Rmd:84-85` and `figure_6_precision.Rmd` read via `Sys.glob`: a key
   in `ALGOS` with no file on disk renders an empty column, silently.
3. **Mis-grouping.** `figure_3.Rmd:160` splits the y-axis with `naive <- c("NR","NMR")` and a
   membership test; an unlisted VR falls into the *recursive* group by default. The companion vector
   `recg` is declared and never read, so it cannot catch the mistake.

Seven other notebooks build the path directly and `read.csv`, which errors loudly on a missing file.
Those are the safe ones. **Figure 1 is in the safe group** (`figure_1.Rmd:49`).

## 1. Figure 1 — what it is, and what blocks VR

**Figure 1 is the cheapest object in the paper.** One local, dispatcher-free invocation of
`macrodr_cli` on `projects/eLife_2025/ops/local/figure_1.macroir` (run 1 of the provenance table,
ledger `projects/eLife_2025/runs/run-20260711-174017/`, whose `script.macroir` is byte-identical to
the current script). It writes 8 loose CSVs straight into `figures/data/` with no hash subdirectory
and no provenance row; `figure_1.Rmd` reads 6 of them and emits `Figure_1.pdf` and `Figure_S1.pdf`.
Content is one simulated 20-channel, 6-interval, 12 ms recording plus the filter steps. Wall clock is
of order tens of milliseconds (inferred from file timestamps; no log was captured).

**The blocker is one block in the config, not C++.** Both figure-1 scripts build the five macro
algorithms with the *bool* builder `build_likelihood_function(...)` (`figure_1.macroir:39,45,51,57,64`),
which has no `variance_form` argument, so VR cannot be requested from them as written. Both builders
return the same `likelihood_algorithm_type` variant that `calc_likelihood_diagnostic` consumes
(`src/cli/command_manager.cpp:1014-1032,1077-1089`), so adding VR is config-only, and
`figure_1_plus_lse.macroir` already calls the family builder for LSE.

**The panel decision is not the data decision**, which is what let the data land before the panel
question was settled: the figure-1 data directory has long carried
`figure_1_likelihood_diagnostic_{IRV,MRV,MNRV}.csv`, three algorithms no panel ever loaded. The panel
question itself was then answered on 2026-07-21: **Fig 1 is the four-column ladder R, MR, VR, IR**,
superseding the 2026-07-20 two-column call (`decisions.md`, `results.md`).

> **Trap.** Those `IRV`/`MRV`/`MNRV` files are the *old* Taylor-variant V-suffix convention
> (2025-12/2026-06, slated for cleanup), unrelated to `macro_VR`. VR's own internal suffix is `_res`
> (`legacy/qmodel.h:644-645`). Do not fold the two cleanups together
> (`theory/macroir/notes/vr_variance_form_plan.md:79-81`).

**`00_plan.md` §5 is correct again**: it always specified the four-column Fig 1 including VR, which
is what was built. It was the two-column call that was the outlier.

### F1-0 · There are two producers at the same paths, and the live one is not `figure_1.macroir`

`projects/eLife_2025/ops/local/figure_1_plus_lse.macroir` (edited 2026-07-21 13:24) writes to **the
same nine fixed paths** as `figure_1.macroir` plus an LSE one, and its own header says it "supersedes
figure_1.macroir once the LSE diagnostic arm exists". So whichever of the two runs last owns the data,
silently. **Add VR to the live script, not the old one**, or the next LSE run clobbers it.

The live script has already settled the questions VR would otherwise raise, and settled them the same
way the code forces:

- `seed = 0` is the `std::random_device` sentinel, so the recording is a fresh draw per run and the
  resolved value is never logged (`figure_1_plus_lse.macroir:6`, and `../_program/decisions.md` §4).
- Therefore **all algorithms must be produced in one process**, against one `simulate`, or their
  panels do not line up (`:4-8`). There is no load-the-simulation path; a VR-only second invocation
  would describe a different recording and the comparison VR exists to make would be invalid.
- Therefore **running it replaces the whole figure-1 data set** and the macro panels change. The
  script's own verdict on that: "expected and harmless — figure 1 is illustrative, not quantitative"
  (`:10-14`). Take that as the standing decision; do not reopen it for VR.
- It is **already on the family builder** for LSE (`:113-117`, `family_approximation = 2`,
  `variance_form_approximation = 0`), so the VR block is a sibling of one that already works here.

Worth one decision, not more: since the data is being regenerated anyway, an explicit non-zero seed
costs nothing and makes this the one figure in the paper that is bit-reproducible. Optional, and
independent of VR.

**Routing flag, not a correction:** that script builds a six-algorithm Figure 1 including LSE, while
paper 1's roster excludes LSE (paper 2) and its Fig-1 decision is two columns. Which paper the LSE arm
of this script serves is unrouted (`../_program/00_index.md` rule 2).

### F1-1 · Add one VR block to each script — DONE 2026-07-21

Landed in **both** producers, deliberately: `figure_1.macroir` (which runs today) and
`figure_1_plus_lse.macroir` (which does not — its LSE diagnostic call is guarded for
`family == nonlinearsqr`, `src/core/likelihood.cpp`, so the whole script dies there until that C++
lands; its own STATUS block says so). **So the VR acceptance ladder runs from `figure_1.macroir`**,
and the LSE script carries the same block ready for when its arm unblocks.

The five macro algorithms stay on the bool builder untouched; only VR needs the family builder,
because that is the only builder exposing the variance form. Argument names verified against
`src/cli/command_manager.cpp:1024-1032` and against the LSE block already written in the second
script:

```
macro_VR = build_likelihood_function_with_family(
   model = model,
   adaptive_approximation = false,
   recursive_approximation = true,
   averaging_approximation = 1,
   variance_approximation = true,
   taylor_variance_correction = false,
   family_approximation = 0,
   variance_form_approximation = 1,
   taylor_qdt_approximation = false)

likelihood_diagnostic_VR = calc_likelihood_diagnostic(likelihood_algorithm = macro_VR,
   parameters = par_tr, experiment = experiment, simulation = simulation)

diagnostic_file_VR = write_csv(experiment = experiment, simulation = simulation,
   likelihood = likelihood_diagnostic_VR,
   path = "figures/data/figure_1_likelihood_diagnostic_VR")
```

Note the family builder takes **`family_approximation`, not `micro_approximation`** (family 0 =
macro), and it has no `micro_approximation` argument at all. Placement: **between `macro_MR` and
`macro_IR`**, build and diagnostic both, so the file reads in ladder order R → MR → VR → IR.

**The verbatim-`macro_VR` rule does not apply here.** Figure 1's CSVs have **no `algorithm` column at
all**: identity lives in the filename, and `figure_1.Rmd:24-39,47-51` attaches it **positionally**
from two parallel vectors (`filenames[]`, `algorithms[]`). A filename typo hard-errors at `read.csv`,
which is safe; a *misordering* of the two vectors silently relabels an entire column, which is not.
So the token to keep consistent is the file suffix, and the pair must be edited together.

Naming: the figure-1 files use the *display* token, and they already carry one divergence (the
mean-non-recursive file is `_MNR` while the notebook tag is `NMR`, `figure_1.Rmd:24-32`, inverted
relative to run 2). Use `VR` for both; do not copy the divergence.

Run: the existing local invocation against the already-built binary, from
`projects/eLife_2025/`. Tens of milliseconds of compute.

Two incidental defects found in the same read, neither blocking: `macro_NR` is built and never
referenced (the NR diagnostic goes through a different overload with its own inline flags,
`figure_1.macroir:99-102`), and NR/R are built with `variance_approximation = 0` while MNR/MR/IR use
`1`. Neither is VR's business; both are worth a look before the figure is called final.

### F1-2 · The acceptance ladder — RUN 2026-07-21, and it changes a sentence in the paper

Measured on the regenerated figure-1 recording (6 intervals, N_ch = 20, one realization, so
**illustrative and not evidence** — the parameter-space test still waits for the grid run):

| | sample 0 (pre-agonist) | sample 1 | later samples |
|---|---|---|---|
| `y_var` MR | 0.05 | **1.04848** | 2.92 → 3.24 |
| `y_var` VR | 0.05 | **0.378085** | 2.25 → 2.39 |
| `y_var` IR | 0.05 | **1.04848** | 1.72 → 1.79 |

1. **VR is not a no-op.** `variance_form` does reach the diagnostic path: VR's predicted variance is
   strictly below MR's at every gated interval, by 26% at the median.
2. **MR's and IR's predicted observable variance are the same object, for any prior.** Not a
   coincidence at one sample: the two expressions collapse to the same thing. With g = `gmean_i`,
   p = `P_mean`, Σ = `p_P_Cov`, and `gtotal_ij` = P_ij·gmean_ij, MR gives gᵀΣg + p·gsqr − p·g², and
   IR gives gᵀ(Σ − diag(p))g + p·(gtotal∘gmean)u + p·gsqr − p·(gtotal∘gmean)u, which is the same
   expression once gᵀdiag(p)g = p·g². The boundary term IR adds to `gSg` is exactly the one it
   removes from `ms`. Confirmed numerically at sample 1, the last step where the three still share a
   prior: MR = IR = 1.048475625791748, VR = 0.378085. **VR's is strictly smaller than both**, by
   E_i[Var_j(gmean_ij | i)]. From sample 2 on the three states have diverged, so those columns
   compare different priors and test nothing.
3. **So the variance ladder is not monotone: MR = IR, with VR below both.** The consequence for the
   paper is not a weakening — it is sharper than the current wording. If MR and IR agree on the
   predicted observable variance at equal state, then **the entire MR→IR difference is the gain**,
   which is what `decisions.md` claims and now has a direct measurement behind it.
4. **But "MR → VR changes only the variance, VR → IR only the gain" does not survive as written.**
   The predictive variance divides the gain, so VR's smaller variance changes the update immediately:
   the posterior departs from MR's at once (~0.01 in P_open here) and the propagated mean follows one
   step later. `VR → IR` therefore moves *both* the variance (back up to MR/IR's value) and the gain.
   Rewrite the one-liner in `00_plan.md:132` and `decisions.md:22` accordingly; the mechanism is
   unchanged, the slogan is wrong.
5. Cumulative logL over the six intervals: MR −6.896, VR −6.075, IR −6.065. **VR lands on IR, not
   between MR and IR.** One recording; do not put it in a caption.
6. **`VR` is therefore not "MR's mean with IR's variance"**, which is how
   `theory/macroir/notes/vr_variance_form_plan.md:9-12` and the paper summaries describe it. IR's
   *total* predicted variance is MR's; the residual form is a component IR uses inside a different
   decomposition, not IR's answer. Say instead: VR keeps MR's mean and boundary-free gain and
   replaces the total per-start-state conductance variance with the residual one.

**Trap found while checking this.** The `gvar_i` column in the diagnostic dump is copied straight
from the Qdt object (`legacy/qmodel.h:5748-5749`) and is **byte-identical between MR and VR**. It is
not the variance-form-selected quantity: the selection happens in the `ms` lambda
(`legacy/qmodel.h:4548-4570`), which never writes back to that slot. Do not use dumped `gvar_i` as
evidence of which variance form ran — derive the effective per-step gating variance as
(`y_var` − e)/N instead. No notebook plots `gvar_i` today, so this is latent.

The original statement of the ladder, kept because it is what was checked:

`vr_variance_form_plan.md:141-149` defines it; on the Figure 1 dump it is an `awk` over one file, and
it is the cheapest possible falsification of the implementation:

- **(0)** VR's predicted mean equals MR's to machine precision. If not, VR is not what it claims.
- **(1)** VR's predicted `y_var` is **strictly below** MR's, by `N·Σ_i p_i Var_j[gmean_ij | i]`.
- **(2)** VR's gain equals MR's (no boundary cross-covariance: both `gSg` and the gain key on
  `averaging::value == 2` alone, `legacy/qmodel.h:4526,4595`).

**Read the variance claim before writing a caption about it.** Reading the two branches side by side
(`legacy/qmodel.h:4562-4563` residual against `:4566-4567` total, with `gtotal_ij = P_ij·gmean_ij`
from `:1653`), the term MR keeps inside `ms` is the same term IR keeps inside `gSg`. It cancels
between them: **given the same prior, MR's and IR's predicted `y_var` are algebraically identical**,
and VR's is strictly smaller than both. If that survives measurement, the ladder is not monotone in
predicted observable variance, and the one-liner "MR → VR changes only the variance, VR → IR only the
gain" (`00_plan.md:132`, `decisions.md:22`) is loose in exactly the panel Fig 6 is built to show.
The spec states it correctly (`vr_variance_form_plan.md:143-145`); the paper-side summaries do not.
**This is a measurement, not an edit: F1-2 settles it on real numbers.**

### F1-3 · Decide the panel (yours)

Under the 2026-07-20 decision Fig 1 keeps two columns and the VR dump is a Methods/SI fact plus the
acceptance evidence for the engine. Reopening that decision is a content call, not a build step; if it
reopens, note the cost is not the run but `figure_1.Rmd`, which does **not** drive its panels off
`ALGOS`: five hand-written filter blocks (`:63-90`), 20 hand-written `cell()` calls, `ncol = 5`
hardcoded at `:592` against a fixed 7.0in page.

## 2. What each figure actually needs

| Fig | Needs VR? | Also blocked on |
|---|---|---|
| **1** | no (panel), yes as smoke test (data) | nothing |
| **2** | only for the A-strict roster re-render | it reads `433ed13` (numeric); VR exists only on the G anchor, so re-point or migrate the non-G dispatcher |
| **3, 4** | yes, if the roster is re-rendered | producer `figure_3_time.macroir` uses the bool builder → same migration as F1-1; dumps are ~1 GB per algorithm |
| **5** | yes | **and** R/MR exist at one noise level only on the G anchor, so the roster map cannot be built without filling them; plus anchor split, panel selection, and three conflicting validity thresholds |
| **6** | **yes, hard** (paper 1's headline) | job 1 needs per-interval `y_var` for VR → inherits the `figure_3_time.macroir` migration; job 2 is five-algorithm on the numeric basis |
| **7** | no (IR-only) | a scope decision, not data |

**"Three of paper 1's figures depend on VR"** (`decisions.md:66`) is never enumerated, and it conflicts
with `results.md`, which implies every five-algorithm figure. Enumerate it or drop the count.

## 3. Order

1. **F1-1 + F1-2.** Config edit, local re-run, acceptance ladder. Hours, no cluster, no compile.
2. **Name VR.** ~~The spelled-out name is written nowhere.~~ **Provisionally `"Variance Recursive"`**,
   set in `figure_1.Rmd`'s `FULL` map on 2026-07-21. The reasoning, so it can be overruled cheaply:
   the roster's prefix letter names *what the conditioning is about* — M for Mean, I for Interval, so
   V for Variance — and the four existing values are Title-Case `"<X> Recursive"`. The spec's only
   written expansion is the lowercase "residual-variance recursive", which is more accurate and too
   long for a column header. **One string, one place, for now.** It becomes eleven places the moment
   the other notebooks gain VR, so settle it before step 3.
3. **The roster edits**, once VR data exists and not before: 11 notebooks, two key forms, four
   `ncol = 5 → 6` bumps. `figure_1.Rmd` is **done** (2026-07-21): VR in both positional vectors, a
   `d_VR` block, `fB_VR`/`fC_VR`/`fD_VR` copied from MR (VR differs from MR only in `y_var`, so the
   panels are MR's on `d_VR`), VR in `.algos`, in `dW`, in both logL range lists, a sixth `colhead`,
   four more `cell()` calls, `ncol = 6`. The remaining ten should be done against the debt already
   booked in `figures_system.md` rather than as ten one-off patches, or the next algorithm repeats
   this.
4. **The VR grid run on dirac** (cluster, `dispatch_figure_3_G.sh`). This is now the blocking item for
   Figs 2, 5 and 6; see §5 for the exact cells and command.
5. **`figure_3_time.macroir` migration + rerun** for Figs 3/4 and Fig 6 job 1.

Steps 4 and 5 are independent of each other and both depend on step 1 passing.

## 3b. Figure 2 — DONE 2026-07-22, and the thesis held

Rendered with the full ladder R, MR, VR, IR on the Gaussian anchor. VR resolved from
`figures/data/0ffbda7` and the other three from `figures/data/1c2ae6f`, all at nsim 10000, so no
n_sim mismatch inside a panel. No edit was needed when VR's battery landed: `DATA_DIRS` is a search
path and the algorithm joined by itself, which is what that design was for.

Empirical-over-Fisher ellipse-area ratios, recomputed from the `battery_pool_G` files:

| | kinetic (k_on, k_off) | corrected | amplitude (N_ch, i) | corrected |
|---|---|---|---|---|
| R | 1.32 | 0.96 | 1.09 | 0.94 |
| MR | 1.97 | 0.93 | 1.53 | 0.88 |
| **VR** | **2.18** | 0.98 | **1.77** | 0.94 |
| IR | 1.02 | 1.01 | 1.00 | 1.07 |

**VR is over-confident and worse than MR in both pairs**, which is the falsifier's prediction, so the
mechanism claim stands and the Results arc collapses to its confirming branch (`decisions.md`).
Caveat to carry into every sentence written from this: **one cell**. The claim is scoped to the
gating-dominated regime and is currently supported at N_ch 100 / noise 0.1 / Delta = 0.1 tau alone.

Also settled by the same render: the two anchors agree. R 1.32 / MR 1.97 / IR 1.02 on the Gaussian
anchor sit where the numeric-Fisher run put them (`decisions/D-4_ranking_verdict.md`: R 1.18-1.42,
MR 1.52-2.07, IR ~1), so moving the main text onto one anchor costs no result.

## 3c. Figures 3 and 4 — producer and notebook prepared 2026-07-22, run pending

**Producer.** `ops/local/figure_3_time.macroir` gained a `macro_VR` block between MR and IR, using
`build_likelihood_function_with_family` (the bool builder used by the other five cannot express the
variance form). Verified against the DSL registration: `calc_dlikelihood_predictions` takes
`likelihood_algorithm_type`, which both builders return, so no C++ and no rebuild.

**The seed was fixed** at the same time, `seed = 0` → `seed = 20260722`. Zero was the
`std::random_device` sentinel, so the ensemble was redrawn every run and never logged: no number in
Figures 3, 4, S4 or S5 could be reproduced, and no algorithm could ever be added later without moving
every other column. From this run on, adding a seventh algorithm leaves the other six untouched.

**NR and NMR stay in the dump.** Only the panel roster drops to the recursive ladder. This costs
~2.3 GB and buys two things: `figure_4.Rmd` and `figure_S4_S5.Rmd` keep working with no edit, and
Figure 4's argument survives — its ratio rows exist to show the per-step identity holds for *everyone*,
so with four recursive columns all sitting on zero there would be no contrast left.

**Notebook.** `figure_3.Rmd` is now roster-driven end to end: `ALGOS <- c("R","MR","VR","IR")`,
a `LEAD` anchor replacing the eleven hardcoded `"NR"` tests that carried the y title, the six row
titles and the panel tags, `naive`/`recg` derived by intersection (with an all-recursive roster the
two-block y scale degenerates to one scale, which is correct), and `ncol = length(ALGOS)`.
Set `ALGOS` back to the five to recover the old figure.

**What the run costs and changes.** About four minutes locally; VR adds ~1 GB to the ~5.3 GB already
there, against 61 GB free. It **overwrites the five existing dumps with a different ensemble**, so
every number in Figures 3, 4, S4 and S5 moves by Monte-Carlo noise. The old numbers are not lost:
they are quoted in `Figure_3_caption.md` and `Figure_4_caption.md`, which are in git. Two of them are
known to change qualitatively and the captions must be rewritten, not patched: Figure 3 panel A's
87-nat spread becomes ~10 nats once NR and NMR leave the panel, and the "e^87 times more probable"
sentence goes with it.

```
cd projects/eLife_2025 && ../../build/gcc-release/macrodr_cli ops/local/figure_3_time.macroir
```

## 3d. Figure 5 — BUILT 2026-07-22. Subject: where IR stops being faithful, in both moments

**New subject.** Figure 5 is no longer the within-family validity map. `00_plan.md` §1 promised that
paper 1 would locate IR's own failure and no figure had ever been assigned to it; this is that figure.
`figure_5.Rmd` supersedes `figure_5_master{,_affine,_frobenius}.Rmd`, moved intact with `git mv` to
`projects/eLife_2025/figures/archive/paper_superseded_20260722/`.

**Design, after five rounds of author critique.** `macro_IR` alone. Columns N_ch {10, 20, 50, 100,
1000, 10000}, chosen to sample the transition rather than the decade grid. Inside a panel x = Δ·k_off
over the 7 free intervals, y = the 6 noise levels. Three blocks:

- **A, grouped MLE distributions** at three design points (**1 distorted**, **2 biased**,
  **3 faithful**), on the (k_off, N_ch) pair because those two carry opposite signs. Square panels,
  legend below, windows symmetric about the truth, three ellipses (empirical / reported / corrected),
  three markers (truth, empirical mean, truth + predicted bias), and **per-axis variance multipliers**
  coloured by the map's own ramp, printed grey when their interval covers one.
- **B, distortion-induced bias**, first moment, from `battery_sim_G` (θ_sim; at θ_pool the bias is
  zero by construction). Units are log10, **not "dex"** — that is astronomy jargon this readership
  does not use.
- **C, information distortion**, second moment, from `battery_pool_G` (θ_pool).

First moment before second, per the natural order. The design points keep their numbering from the x
axis so a reader finds them where the number says.

**The result the figure now carries.** The two moments fail at **opposite ends of the interval axis**:
bias worst at coarse Δ (0.064 in log10 for N_ch at Δ·k_off = 1, N_ch 10, falling to 0.002 at 0.01),
distortion worst at short Δ and shaped like a **U** (k_off 1.46, 1.21, 1.27, 1.40, 1.47, 1.56, 1.57
from Δ = 1 down to 0.01, minimum at 0.5). No single interval optimises both, which is a design
statement the paper can hand a reader. The two parameter rows carry opposite signs, so any scalar
averaged over parameters cancels them and reports a faithful algorithm.

**The panel-A limitation, measured and not hand-waved.** The ellipse compares the covariance of an
estimate against the sandwich MARGINALISED over the other parameters, while C compares score against
Fisher unmarginalised and fits nothing, so the ellipse is systematically milder. Measured at N_ch 10,
Δ 0.01, noise 0.05: ellipse 1.34 / 0.83 against map 1.62 / 0.65. It gets worse with N_ch, and by
N_ch 200 the ellipse reads *calibrated* (0.93 / 1.01) where the map says 1.25 / 0.90 — the panel would
contradict its own map. At the current point 1 (N_ch 20) the N_ch direction reads 1.13 against the
map's 0.851 and **is not resolved by 100 fits**; the caption says so.

**Two traps that cost real time to find, both now written into the notebook:**

1. **group_size 10 does not merely add noise at few channels, it reverses the sign** of the N_ch
   direction (1.58 measured against 0.65 true, N_ch 10). Its finite-sample inflation is comparable to
   the true deflation. Group 100 is the smallest defensible choice there, at the cost of 100 fits and
   14% relative uncertainty on a variance.
2. **The MLE distribution is not normal at few channels.** Skewness of the N_ch marginal is 4.13 at
   N_ch 10 group 10 (53 SE), 0.36 at N_ch 50, 0.14 at 200. A 5% trimmed variance flips every plain
   ratio below one and onto the map's value, so the plain sample variance is dominated by a right
   tail. **This is not a defect of the likelihood**: normality is an asymptotic property of the
   estimator, so a correct likelihood can produce a skewed MLE. It limits what an ellipse can mean,
   not what the likelihood is worth. Empirical coverage of the nominal 95% interval is 0.80–0.93 at
   N_ch 20 and 0.94–0.98 at 200.

## 4. The VR grid run (dirac)

**State on 2026-07-22 11:40.** The run is proceeding. `figures/data/0ffbda7/` now has the **complete
five-file set** (cloud, `battery_pool_G`, `battery_sim_G`, `empirical_G`, `pool_runs`) for N_ch 100,
1000 and 10000 at noise 0.1, plus N_ch 100 and 1000 at noise 1; N_ch 10 has clouds only at both
noises, and N_ch 10000 at noise 1 has its cloud. Figure 2's cell was complete as of 23:28 on 07-21 and
the figure is built (§3b).

**A cloud alone cannot draw Figure 2**, which is why the column waited: the panel carries three
ellipses and a bias marker, and only the empirical one comes from the cloud. Fisher and corrected come
from `_battery_pool_G`, the bias from `_battery_sim_G`.

**Nothing needed editing when they arrived, and nothing will for the rest.** `figure_2.Rmd`'s `DATA_DIRS` is a search path and each
algorithm resolves to the first directory holding a *complete* set for it, so VR joins by itself.
Mixing directories in one panel is sanctioned by D-0: every CSV stamps its own engine hash in row 1,
so the folder is a scratch grouping and the file is the provenance record.

### The cells VR still owes

To stand as a column beside R and MR, VR needs what they have: N_ch {10, 100, 1000, 10000} × noise
{0.1, 1, 10} at nsim 10000, twelve cells. At the measured 181 CPU-h per Gaussian cell
(`decisions/D-0_freeze_and_rerun_scope.md`) that is **~2,170 CPU-h**, against the ~42k left on the
account. Figure 2 alone needs exactly one of them, N_ch 100 at noise 0.1.

`dispatch_figure_3_G.sh` takes `NCHS`, `N_SIMS` and `N_NOISE` as **parallel arrays paired by index**
(one SLURM job per entry); the seven-point interval axis and the group-size axis are swept inside each
job, so they cost nothing extra. `RUN_DIR` overrides the per-commit output folder, which is how a
newer binary writes beside an older run.

Noise 0.1 belongs with the rest of that column in `1c2ae6f`:

```
RUN_DIR=1c2ae6f NCHS="10 100 1000 10000" N_SIMS="10000 10000 10000 10000" \
N_NOISE="0.1 0.1 0.1 0.1" N_ALGO="macro_VR" \
  projects/eLife_2025/ops/slurm/dispatch_figure_3_G.sh dirac
```

Noise 1 and 10 belong where R and MR keep theirs, in `87889e6`:

```
RUN_DIR=87889e6 NCHS="10 100 1000 10000 10 100 1000 10000" \
N_SIMS="10000 10000 10000 10000 10000 10000 10000 10000" \
N_NOISE="1 1 1 1 10 10 10 10" N_ALGO="macro_VR" \
  projects/eLife_2025/ops/slurm/dispatch_figure_3_G.sh dirac
```

**Precondition:** the dirac binary must bake a commit that has VR (`0ffbda7` or later);
`projects/eLife_2025/ops/build_cluster.sh dirac` builds it, and the dispatcher defaults `BIN` to
`build/macrodr_cli-dirac-current`. Cluster access is the author's (`01_writing_plan.md` §6).

## 5. Open, and blocking

- Does Fig 1 get a VR column (§1, F1-3)? Default: no, per 2026-07-20.
- VR's spelled-out display name (§3.2).
- Which figures the "three" are (§2).
- Fig 5's anchor and threshold conflicts, already open in `results.md`.
- The arc cannot be captioned before VR runs: the thesis inverts if VR comes out calibrated
  (`decisions.md:69-73`).
