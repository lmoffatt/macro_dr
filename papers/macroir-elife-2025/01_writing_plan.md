# Writing plan: the paper as a task graph agents execute

> Opened 2026-07-14. Rewritten 2026-07-14 to be **executable by agents under supervision**.
>
> Owns: **the objective, the task graph, the checks, and the human budget.** It does not own what any
> section *says* (the section plans do), nor the engine work (`09_carve_plan.md`, Freeze
> preconditions), nor which decisions were already settled (`02_decision_log.md`).
>
> Design rules for this document: every task names its inputs, its single output file, and a
> **command** that decides whether it is done. Every human touchpoint is a **yes/no or a correction
> of a prefilled artifact**, never a research assignment. No checkboxes: status is *derived* by
> running `check.sh`, never declared. Replaces `00_master_plan_v2.md` §9 and `01_workboard.md`.


## 0. Done means exactly this

`papers/macroir-elife-2025/check.sh` returns green when:

1. `elife_paper.tex` compiles with `latexmk -pdf`, zero LaTeX errors.
2. Zero `\todo{}` remain.
3. Every section file under `sections/` is non-empty and contains no `TBD`/`XXX`/`FIXME`.
4. **Every numeric literal in prose carries a source comment** `% src: <file>[#<locator>]` (LINT-SRC).
5. `\cite` count > 0, every `\cite` key resolves in `biblio.bib`, and `biblio.bib` contains no
   uncited entry (it is today the P2X2 paper's 170 references, and the `.tex` has **zero** `\cite`).
6. Six figures are `\includegraphics`'d, each with a caption, each caption's numbers passing LINT-SRC.
7. Data availability, Author contributions, Competing interests, Funding are non-empty.
8. Word count is under the limit for the article type chosen in **D-1**.

Green plus the seven section ACCEPTs of §3 = submittable. Nothing else counts as progress.
The previous version of this plan aimed at *a draft with two holes*, which is not the objective.


## 1. Three invariants (the first two are kept from the previous plan; they are now lints, not virtues)

- **LINT-SRC — no number without a source.** Every numeral an agent writes into the `.tex` is
  followed on the same line by `% src: <file>[#<locator>]` naming the file it was read from.
  This is the anti-hallucination gate, and it is also what makes the draft **survive a re-run**:
  if D-0 lands on "re-run", every number in the paper changes by a little, and refreshing them is
  then a mechanical grep-and-replace by an agent instead of a re-read of the whole manuscript.
- **LINT-TODO — a hole is written as a hole.** Where the evidence does not exist yet, write the
  paragraph that will surround it and mark the claim `\todo{}`. Never write a claim you cannot
  support and plan to check later.
- **LINT-SCOPE — an agent writes only its Output file.** Anything else it wants to change, it
  reports. This is what makes parallel drafting safe.


## 2. The shape of a task

Every task below is a brief an agent can execute cold, with no memory of this project:

| Field | Meaning |
|---|---|
| **Inputs** | The only files it may read. If it needs a fact that is in none of them, it **stops**. |
| **Output** | The single file it may write. |
| **Constraints** | Must-state and may-not-claim, verbatim. These are load-bearing; they encode errors already made. |
| **Check** | A command. Pass/fail, no judgment. |
| **Gate** | What Luciano sees. Default: nothing. |

Standing rules for any agent taking a task here:

1. **Never write a number you did not read from a file named in Inputs.** Attach `% src:`.
2. **If a needed fact is absent, do not guess and do not go looking outside Inputs.** Stop and emit a
   decision brief in the §3 format (evidence + recommended answer). A stopped task is cheap; a
   fabricated number costs the paper.
3. **Never restate a section plan.** Cite it. One topic, one owner (`00_master_list.md`).
4. Section plans are *inputs to prose*, not text to paste. The `.tex` is written for a referee.


## 3. The entire human budget (~3.5 hours)

The rule that makes this plan help rather than assign: **an agent never brings Luciano a question
without (a) the evidence already gathered and (b) a recommended answer.** His action is to approve,
pick, or correct a prefilled artifact.

| # | Decision | What he receives | His action | Min |
|---|---|---|---|---|
| **D-0** | **Freeze commit and re-run scope** (§4). The one that can add a month. | One page: three options, CPU-hours costed, what each one breaks. | Pick one | 15 |
| **D-1** | Research Article or Tools & Resources? T&R obliges benchmarking against existing methods (new science) and public deposition. | One page: what T&R obliges, what it buys, the word limits. | Pick one | 5 |
| **D-2** | The units of the six parameters. Written nowhere in the repo: not in the model header, not in a `.macroir`, not in a dispatcher. **Only Luciano knows.** | A six-row table **prefilled** with the agent's inference from the code and priors, one evidence line per row. | Correct the wrong cells | 15 |
| **D-3** | The novelty sentence. `introduction_plan.md` still claims "no published likelihood integrates the observable over the acquisition window"; the prior-art map retracts it (single-channel HMMs do exactly that). The novelty is the **scaling to channel populations**, the exact CTMC realization, and the validation. | A drafted replacement paragraph from `docs/bibliography/MacroIR_prior_art_map.md`. | Approve or edit | 10 |
| **D-4** | The ranking verdict. The table exists in four copies and two cells disagree with the data: MR's variance direction, and IR's corner bound (≤1.3 stated, ~0.5 observed) — probably a category error between predicted *observable* variance and reported *parameter* covariance. | The table recomputed from the CSVs, the two contested cells resolved, the diff against each of the four copies. | Yes/no | 10 |
| **D-5…D-11** | Seven section ACCEPTs (Theory, Diagnostics, Results, Methods, Discussion, Introduction, Abstract). | The one `.tex` diff, `check.sh` green, and the number→source table for that section. | Accept, or list changes | 20 ea |

**Not a decision, and not a gate:** the overconfidence factor (`abstract_draft.md` says 10-16 on line
158 and 14-21 on line 162). An agent resolves it from `Figure_S3_caption.md`, which is where it is
computed, and records it in `02_decision_log.md`. It was never Luciano's to settle.


## 4. D-0: the decision the previous plan missed

The previous version said the two missing figures need "exactly four cells (NMR × N_ch {10, 100,
1000, 10000} at noise 0.1), about 1,830 CPU-hours, one to two days of queue", and called it cheap.
`09_carve_plan.md` (Freeze preconditions) says something that does not fit inside that estimate:

> The binary stamps its own git hash into the output directory name and into row 1 of every CSV, and
> the dispatcher takes the data folder from `$BIN --commit`. So **any engine change alters the
> provenance key of everything produced after it.** Land these first, tag once, build once, run
> everything from that build. Otherwise the deposited code is not the code that made the deposited data.

E-1 (log the resolved seed) and E-3 (stop double-writing every evolution row, which **changes the CSV
schema**) both change the hash. Today's figures were built at `433ed13` and `1c2ae6f`. So the real
question is not "which four cells", it is:

| Option | What it costs | What it breaks |
|---|---|---|
| **(a) Freeze at the current hash.** Drop E-1/E-3. | 1,830 CPU-h (the four cells only). | Data Availability cannot claim reproducible ensembles (`seed = 0` means *random* and the resolved value was never logged), and the row-duplication bug ships in the deposited engine. |
| **(b) Land E-1…E-5, tag, build once, re-run everything from that build.** | The full grid. Dirac has no CPU-hour quota; the number needs costing. | Every figure is regenerated, so **every number in the manuscript changes slightly** (the ensembles are stochastic). LINT-SRC makes that a mechanical refresh pass. |
| **(c) Land E-1…E-5, re-run only the cells the paper's figures consume.** | Between (a) and (b). | The paper's data is one hash; older exploratory data is not, and must not be mixed into a panel. |

An agent costs (b) and (c) in CPU-hours from `05_experiment_grid.md` and the existing run logs, and
brings the one page. **Nothing that produces a number may start before D-0 is answered**, because
(b) and (c) invalidate the numbers that a drafted Results would cite. Prose that carries no numbers
is unaffected and starts immediately.

Second, smaller, and independent: **the schema port.** Every production notebook hardcodes the
numerical run and its column names (`Likelihood_Information_Distortion`, …). The Gaussian run emits
*different* names (`Gaussian_Distortion_Induced_Bias`, …) **and still carries the old columns filled
with NaN**, so a notebook repointed at the Gaussian directory does not fail: it renders grey. This is
real work between the run and Figure 5, nobody has costed it, and it is the thing most likely to be
mistaken for a re-knit.


## 5. Day zero: three agent tasks, no human input, start now

| ID | Task | Output | Check |
|---|---|---|---|
| **T-1** | Split `elife_paper.tex` (90 lines, skeleton, section titles + fill-hints in comments) into `sections/*.tex` pulled in with `\input`. One file per section = one agent owner = no merge conflicts. Carry the fill-hint comments across verbatim. | `docs/manuscript-drafts/sections/*.tex` + a 30-line `elife_paper.tex` | `latexmk -pdf` compiles, PDF unchanged |
| **T-2** | Write `check.sh` (the definition of done, §0). Everything else in this plan is measured by it. | `papers/macroir-elife-2025/check.sh` | It runs and reports the current state honestly (it will be red on almost everything: that is correct) |
| **T-3** | **Build `biblio.bib` from zero.** The longest pole in the project, and it has no blockers. Sources already on disk: `docs/bibliography/MacroIR_prior_art_map.md`, `docs/bibliography/identifiability/`, the Comm Biol and Moffatt 2007 lineage, Munch 2022, the Kalman prior art. | `docs/manuscript-drafts/biblio.bib` | Every entry has a DOI or a stable URL; `bibtex` runs clean |

T-3 is the item most likely to be underestimated. It is also the one an agent can do end to end.


## 6. The task graph

Writing tasks. `Needs` are hard preconditions. Everything with no `Needs` can run today, in parallel,
on separate files (that is what T-1 buys).

| ID | Section | Needs | Inputs (read only these) | Constraints (load-bearing) |
|---|---|---|---|---|
| **W-1** | Theory | T-1 | `theory_plan.md`, `nomenclature.md`, `theory/macroir/docs/Macro_IR/macroir_derivation.tex`, `macroir_macroir_paper_section.md` | Assume `nomenclature.md`, do not re-derive it (an MNR/NMR flip later is a `sed`, not a rewrite). **May not claim** the recursive members' Gaussian Fisher is the exact information matrix — it is not, the likelihood is misspecified by construction. (`docs/corrected covariance justification.md` says otherwise and is wrong on this point.) Every symbol used later in the paper is defined here. |
| **W-2** | Diagnostics | T-1 | `diagnostics_plan.md`, `correction_idm_reconstruction.md`, `supplement_information_distortion_main.tex`, `Gaussian_Fisher_Distortion_Family.md`, `03_metrics_diagnostics.md` | The identity is `IDM = K·CDM·Kᵀ`, **not** `SDM^½·CDM·SDM^½` — the supplement and `likelihood.cpp:3501` both print the false one (E-2). Be explicit, per identity, about which is **exact algebra**, which is **first-order**, and which is **convention**. |
| **W-3** | Results §1–§4 | T-1, D-0, D-4 | the four figure captions, `results_plan.md`, `docs/figure_provenance.md` (which run made which figure, and the three naming traps that have already misled readers) | Every claim traceable to a panel, every number to the file that computed it (LINT-SRC). Each figure answers the objection raised by the previous one. **Must state once:** which Fisher anchors which panel — the definitive figures were promised on the Gaussian anchor and Figures 2, S2, S3 currently sit on the numerical one (`433ed13`). |
| **W-4** | Methods | T-1, D-2 | `methods_plan.md` M1–M11, `docs/figure_provenance.md` | eLife excludes Methods from the word count, so be complete: the claims *are* about algorithms, so the Methods **are** the result. **Must state plainly:** `seed = 0` meant *random*, the resolved seed was never logged, so the ensembles are statistically equivalent but **not bit-reproducible** (`methods_plan.md` M5 currently claims the opposite and must be corrected — unless D-0 lands on (b)/(c), which fixes E-1 and makes them reproducible; the Methods text follows D-0); and the ground truth is exact CTMC uniformization, which is what licenses using it as the reference. |
| **W-5** | Results §5–§6 | C-lane (below) | Figures 5 and 6 once they exist | Until the figures exist: write the surrounding prose with the claim as `\todo{}`. Do not guess the numbers. Forty-four candidates sit in `figures/in_progress/`, sorted by `figure_5_PLAN.md`; four **unselected** PDFs already squat on the names `Figure_5_*` / `Figure_6_*` in `figures/paper/` and must be cleared before promotion or they will be cited by mistake. |
| **W-6** | Discussion | W-3, D-3, D-4 | `discussion_plan.md`, `docs/bibliography/identifiability/` | **Done when it delivers a decision rule** (how many channels, what interval relative to the relaxation, how much noise, what the cheap approximations cost you in that corner) — not "MacroIR is best". The decision rule is what gets cited. **Must contain:** the Kalman concession (MacroIR ≡ integrated-measurement augmented Kalman filter, verified to ~1e-8: concede it loudly, the novelty is elsewhere); the identifiability perimeter; and IR's own miscalibration at few channels (it over-predicts the observable variance by ~6% — the hero's own failure belongs in the paper that is about failure). |
| **W-7** | Introduction | W-3, D-3 | `introduction_plan.md`, the prior-art map | Written late on purpose: it sells what the paper delivers, which you only know once the Results exist. Done when the reader reaches "we provide that test" already convinced that no such test existed for this problem and that the field is choosing likelihoods by habit. |
| **W-8** | Abstract + Impact | all of the above | `abstract_draft.md`, the finished sections | The 209-word abstract in `elife_paper.tex` is **not a placeholder**: it is in the validity-map frame and it reads well. **Reconcile** it with the finished paper; do not rewrite from zero. Also produces the Impact Statement (eLife requires one, ≤40 words). |
| **W-9** | Front/back matter | D-1 | `09_carve_plan.md`, the venue instructions | Data availability, CRediT, competing interests, funding, acknowledgements, Key Resources Table if the article type requires one. Unowned until now; it is check-item 7. |

Compute lane (gated on **D-0**, then owned by `09_carve_plan.md` and `05_experiment_grid.md`):
**C-1** land E-1…E-5 → **C-2** tag + build once → **C-3** run the cells D-0 selected → **C-4** the
schema port (§4) → **C-5** promote Figures 5/6 from `in_progress/` and caption them → **C-6** the
number-refresh pass over every `% src:` in the manuscript.


## 7. The critical path, and what it means

```
D-0 ─► C-1 ─► C-2 ─► C-3 ─► C-4 ─► C-5 ─► W-5 ─┐
                                               ├─► W-6 ─► W-7 ─► W-8 ─► done
T-3 (bibliography) ────────────────────────────┤
W-1, W-2, W-4 (no numbers or few) ─────────────┘   [float: start today, finish whenever]
```

The finish date is set by **compute and the bibliography**, not by writing. Theory, Diagnostics and
Methods are pure float: large, mechanical, and blocked on nothing but D-2. That is why D-0 and T-3
start today and the drafting is allowed to lag them.

The previous plan sequenced the *sections* and left the two long poles (the compute freeze, and a
bibliography that does not exist) in a section called "what is deliberately not here". It optimised
the order of the work that was not the constraint.


## 8. Status is derived, never declared

`01_workboard.md` died with all 25 of its checkboxes unticked while the work they named was being
done. The fix is not to record less state; it is to **stop declaring state and start computing it**.

```
./check.sh          # the definition of done in §0, evaluated against the repo
```

There is no ledger in this file to keep in sync, and no checkbox to forget. If you want to know where
the paper is, run the command. Settled decisions go to `02_decision_log.md` and nowhere else; the
other copies become citations.
