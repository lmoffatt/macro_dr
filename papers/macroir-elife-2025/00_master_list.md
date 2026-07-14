# Master list: who owns what

> Opened 2026-07-14. This is the routing document for the MacroIR / eLife paper. It settles nothing about the science; it settles **where each thing is settled**.
>
> Built from a full read of the 26 documents in this pack plus the theory, figure, bibliography and source layers they depend on.

## The rule

**One topic, one owner. Everyone else cites.**

If you are about to write a fact, a number, a scope call or a verdict into a document, first find the topic in §3. If another document owns it, do not restate it: reference it. A restatement is not a summary, it is a second copy that will drift, and drift is what makes the same argument get had twice.

Three consequences worth naming, because the pack currently violates all three:

- A **section plan** (`results_plan.md`, `methods_plan.md`, and the rest) owns its section's argument and prose. It does **not** own cross-cutting decisions (scope, naming, the ranking, figure numbering). It consumes them.
- A **decision** is settled in exactly one place, `02_decision_log.md`. A document that argues toward a decision (`nomenclature.md`, `title_options.md`) owns the *argument*; the *verdict* graduates to the log.
- A document that is superseded says so **in its first three lines**, names its successor, and is never read as current again.

## 2. The registry

Status vocabulary: **LIVE** = the owner, read it as current. **REWRITE** = still the owner, but it currently says something wrong (see §4). **POINTER** = kept only to redirect. **RETIRED** = tombstone. **CREATE** = does not exist yet.

`git` column: `untracked` means the file is not in version control and one `git clean` deletes it. This is the most urgent problem in the pack and it applies to most of the documents everything else now depends on.

### A. Governance

| Document | Owns | Status | Updated | git |
|---|---|---|---|---|
| `00_master_list.md` (this) | The topic→owner routing table; the document lifecycle | LIVE | 2026-07-14 | new |
| `00_master_plan_v2.md` | The thesis, the two-approximation frame, the three regimes, the open-decision register (D-1…D-7), the task list | REWRITE | 2026-07-13 | untracked |
| `02_decision_log.md` | Every **settled** cross-cutting decision: scope, algorithm roster, the anchor, the ranking verdict, the repo boundary | REWRITE | 2026-07-13 | modified |
| `08_sources_audio_notes.md` | Which audio/chat is the source of record, and which wins when two disagree | LIVE | 2026-07-13 | modified |
| `09_carve_plan.md` | The repo boundary, the freeze trigger, what moves to `macroir-validity`, code availability | LIVE | 2026-07-13 | untracked |
| `07_code_tasks.md` | The engine backlog owed to this paper (regression test, the CSV row-duplication writer bug, the IDM call site) | REWRITE | 2026-07-13 | modified |
| `README.md` | Nothing. It routes, and it routes to three tombstones | POINTER | 2026-07-14 | modified |
| `AGENTS.md` (repo root) | Where an agent starts. Must register this file | REWRITE | — | modified |
| `00_master_plan.md` | — superseded by `00_master_plan_v2.md` | RETIRED | — | modified |
| `01_workboard.md` | — task tracking moved to master plan §9 | RETIRED | — | modified |
| `06_repro_pipeline.md` | — moved to `09_carve_plan.md` and `docs/figure_provenance.md` | RETIRED | — | modified |

### B. Argument and scope

| Document | Owns | Status | Updated | git |
|---|---|---|---|---|
| `From molecular … PROGRAM.md` | This paper's position in the program (bridge 2; Comm Biol is the demonstration) | LIVE | 2026-07-13 | untracked |
| `nomenclature.md` | The five names, the letter semantics, the MNR/NMR spelling, the data-key↔display-label map | REWRITE | 2026-07-14 | untracked |
| `title_options.md` | The title argument and the project's claim/accuracy rules (the rules do not belong here, see §5) | REWRITE | 2026-07-14 | untracked |
| `docs/bibliography/MacroIR_prior_art_map.md` | The novelty position, the Kalman concession, every citation attribution | LIVE | 2026-07-14 | modified |
| `docs/bibliography/identifiability/README.md` | The identifiability perimeter (what macroscopic data cannot recover) | LIVE | 2026-07-13 | untracked |

### C. Science content

| Document | Owns | Status | Updated | git |
|---|---|---|---|---|
| `03_metrics_diagnostics.md` | The metric **registry**: each diagnostic's definition, its sign convention, its threshold, and which CSV field carries it | REWRITE | 2026-07-13 | modified |
| `05_experiment_grid.md` | The sweep design: control variables, the definition of τ, the per-cell protocol | REWRITE | 2026-07-13 | modified |
| `correction_idm_reconstruction.md` | The correct distortion-decomposition identity (`IDM = K·CDM·Kᵀ`) and what composes exactly | LIVE | 2026-07-14 | untracked |
| `theory/…/Gaussian_Fisher_Distortion_Family.md` | What the anchor **H** is, and the numerical-vs-Gaussian Fisher bridge | LIVE | 2026-07-05 | committed |
| `theory/…/score_martingale_argument.md` | Why recursion fixes the reported uncertainty without fixing the bias | LIVE | — | committed |
| `theory/…/Macro_IR/` | The algorithm derivation feeding Theory and the supplement | LIVE | — | committed |
| `analysis_figure_S1_score_mean.md` | What the score-mean figure shows (and only that) | LIVE | 2026-06-30 | committed |
| `figures/paper/sample_correlation_distortion_analysis.md` | What the sample/correlation scaling laws do and do **not** support | LIVE | 2026-07-09 | committed |

### D. Figures and data

| Document | Owns | Status | Updated | git |
|---|---|---|---|---|
| `results_plan.md` | **The figure arc**: which figures exist, in what order, what claim each carries, the numbering register, the do-not-say list | LIVE | 2026-07-14 | untracked |
| `docs/figure_provenance.md` | **The run manifest**: which script produced which data, what is built, what is not, what is reproducible | LIVE | 2026-07-14 | untracked |
| `figures/in_progress/figure_5_PLAN.md` | Panel selection for Figure 5 (main / supplementary / drop), and nothing else | LIVE | 2026-07-08 | committed |
| `figures/paper/*_caption.md` | The caption text of the figure each is named for | REWRITE | various | committed |
| `figures/instructions.md` (this pack) | Figure assembly and the cross-figure visual system (colours, fonts, panel letters) | REWRITE | 2026-03-25 | committed |
| `04_figures_storyboard.md` | — the 4-figure arc is dead; the panel sketches survive only as reference | POINTER | 2026-07-13 | modified |

### E. Manuscript sections

Each owns its section's job, structure, argument and prose. None owns a cross-cutting decision.

| Document | Owns | Status | Updated | git |
|---|---|---|---|---|
| `abstract_draft.md` | The abstract and the Impact Statement, while drafting | REWRITE | 2026-07-14 | untracked |
| `introduction_plan.md` | The Introduction | REWRITE | 2026-07-14 | untracked |
| `theory_plan.md` | The Theory section (the family as one object with two knobs) | LIVE | 2026-07-14 | untracked |
| `diagnostics_plan.md` | The Diagnostics section (the validation machinery as a procedure a sceptic could run) | LIVE | 2026-07-14 | untracked |
| `results_plan.md` | The Results (see D: it also owns the arc) | LIVE | 2026-07-14 | untracked |
| `discussion_plan.md` | The Discussion, and the decision rule the paper is cited for | LIVE | 2026-07-14 | untracked |
| `methods_plan.md` | Materials and Methods, the emission model, and every run value | REWRITE | 2026-07-14 | untracked |
| `docs/manuscript-drafts/elife_paper.tex` | Nothing. It is the **vessel**, not a second author. Every claim in it is owned upstream | LIVE | 2026-07-13 | modified |

### F. Production (does not exist yet, see §5)

| Document | Owns | Status |
|---|---|---|
| `10_submission_pack.md` | Data Availability, code availability, CRediT, competing interests, funding, MDAR, Key Resources Table, preprint, article type | CREATE |
| `10_claims_and_style.md` | The accuracy rules (currently homeless inside a title brainstorm and already forked) | CREATE |
| `docs/manuscript-drafts/macroir.bib` | The citation database. The paper does not have one | CREATE |

## 3. Topic index

Before writing, look the topic up here.

| Topic | Owner |
|---|---|
| Thesis, the two Gaussian approximations, the three regimes | `00_master_plan_v2.md` |
| What is in and out of scope | `02_decision_log.md` |
| The algorithm roster (why exactly five) | `02_decision_log.md` |
| The names of the five, and the spelling | `nomenclature.md` → verdict to `02_decision_log.md` |
| The ranking verdict (IR sole survivor, MR strawman, …) | `02_decision_log.md` |
| The distortion anchor (Gaussian vs numerical Fisher) | `theory/…/Gaussian_Fisher_Distortion_Family.md`; the *choice* in `02_decision_log.md` |
| Diagnostic definitions, thresholds, CSV fields | `03_metrics_diagnostics.md` |
| The decomposition identity | `correction_idm_reconstruction.md` |
| The sweep axes and the definition of τ | `05_experiment_grid.md` |
| Which figures exist and in what order | `results_plan.md` |
| Which data commit each figure sits on; what is reproducible | `docs/figure_provenance.md` |
| Figure 5's panels | `figures/in_progress/figure_5_PLAN.md` |
| Novelty, prior art, the Kalman concession, citations | `docs/bibliography/MacroIR_prior_art_map.md` |
| The title | `title_options.md` → verdict to `02_decision_log.md` |
| The abstract | `abstract_draft.md` |
| Run values, the emission model, units | `methods_plan.md` |
| The repo boundary, the freeze, code availability | `09_carve_plan.md` |
| Which audio is the source of record | `08_sources_audio_notes.md` |
| Everything a journal asks for on submission day | `10_submission_pack.md` (to create) |

## 4. Live collisions

Places where two documents currently say different things about the same object. Each needs one decision, recorded once. Ordered by how much damage it does if it reaches print.

| # | The disagreement | Who says what | Owner who decides |
|---|---|---|---|
| C-1 | **The regime-map axes do not exist.** `02_decision_log.md`, `00_master_plan_v2.md` §6 and `05_experiment_grid.md` all commit to maps over **N_ch × K_off**. `figure_3_mle_G.macroir:46` hard-fixes `off = 100`, and every built and candidate panel is over **N_ch × interval** (and noise). The settled decision was never implemented. | decision log vs the pipeline | `05_experiment_grid.md` |
| C-2 | **The title.** `elife_paper.tex:13` = *"Validity limits of time-averaged likelihoods…"*. `title_options.md` CHOSEN = *"Information distortion in likelihood approximations…"*, and its own later revision prefers a third. The prior-art map warns that "information distortion" is taken in neural coding. | tex vs title_options | `title_options.md` |
| C-3 | **The anchor split blocks Results.** Master plan and decision log: *"definitive figures anchor on the Gaussian Fisher"*. `results_plan.md` replies that they currently cannot: the Gaussian run `1c2ae6f` contains **no NMR at all**, so every cross-algorithm panel must come from the numerical-Fisher run `433ed13`. | plan vs data on disk | `02_decision_log.md` (it is a compute decision) |
| C-4 | **MNR vs NMR.** Decision log mandates NMR. `nomenclature.md` derives MNR and rejects the "NMR = MacroINR" bridge as wrong. The `.tex` writes NMR; every built caption displays MNR; every data key is `macro_NMR`. | four-way | `nomenclature.md` |
| C-5 | **The seed.** `methods_plan.md` M5 records `seed = 0` as though it fixed the random stream. It is the sentinel for *random*: `calc_seed(0)` draws from `std::random_device`, and the resolved seed is logged nowhere. **The simulated ensembles behind every figure are not reproducible.** | methods_plan vs the source | `docs/figure_provenance.md` |
| C-6 | **The headline overconfidence factor**, 10–16 or 14–21. `abstract_draft.md` contradicts itself four lines apart. | abstract vs itself vs S3 caption | `Figure_S3_caption.md` (it is computed there) |
| C-7 | **The novelty sentence.** `introduction_plan.md` keeps *"no published likelihood integrates the observable over the acquisition window"*. The prior-art map retracts exactly that claim (single-channel HMMs do it; the novelty is the scaling to channel populations). | intro vs prior-art map | `docs/bibliography/MacroIR_prior_art_map.md` |
| C-8 | **The decomposition identity** printed in the supplement and computed in `likelihood.cpp:3501` is the one `correction_idm_reconstruction.md` proves false. | supplement + code vs the correction | `correction_idm_reconstruction.md` |
| C-9 | **The figure arc**, in six versions (storyboard, master plan §5, results_plan, figure_5_PLAN, figure_5_master_STATUS, and what is on disk). `figure_provenance.md` §8 already rules that `results_plan.md` wins. | six-way | `results_plan.md` |
| C-10 | **Figure numbering.** Two captions are numbered S4 while `Figure_S4.pdf` (bias) and `Figure_S5.pdf` (acf) both exist; `Figure_S1.pdf` has no caption and the one sentence describing it is wrong. | captions vs disk | `results_plan.md` |
| C-11 | **The validity threshold**: 1.1 (`03_metrics_diagnostics.md`), 1.15/1.5 (`figure_5_master_STATUS.md`), 1.15 (`figure_7_validity_map.Rmd`). | three-way | `03_metrics_diagnostics.md` (closes D-4) |
| C-12 | **The direction of MR's variance error.** The master plan says it overestimates; the S3 caption reports it as overconfident by 1.5–2×. Probably a category error (observable variance vs parameter covariance), which is exactly why it must be written down once. | plan vs caption | `results_plan.md` |

Duplication without contradiction, which is where the *next* collision comes from: the scope call is written out in full in four documents, the evidence-correction formulas (½ log det C, α⋆ = p / tr C) in five, the diagnostics list in three, the ranking table in three. Each copy is a place a future correction can fail to land.

## 5. Missing owners

Nothing in the pack owns these, and the paper cannot be submitted without them.

**Blocking**

1. **There is no bibliography.** `elife_paper.tex` has zero `\cite` commands and points at `biblio.bib`, which is the P2X2 paper's 170 references. The ~60 statistics and filtering works this paper must cite exist only as prose in the prior-art map. → create `docs/manuscript-drafts/macroir.bib`.
2. **Figures 5 and 6 do not exist.** No notebook selected, no PDF, no caption, and C-3 is unmade underneath them. `results_plan.md` calls Figure 5 the critical path.
3. **Data availability.** eLife makes it mandatory and rejects "available on request". The CSVs are gitignored, tens of GB, no repository and no DOI named anywhere, and per C-5 the ensembles cannot be regenerated bit for bit.
4. **The units of the six parameters.** `methods_plan.md` calls this the one gap only the author can close.
5. **The run manifest** for `433ed13` and `1c2ae6f`: both invocations are reconstructed from filenames, the dispatchers' committed defaults contradict the data, and which cluster ran what is written nowhere.
6. **The n_sims-uniformity audit.** `433ed13` mixes cells at n_sims = 10000, 1000 and 200, and every scalar summary of the distortion matrix is biased in n_sims (Jensen). `diagnostics_plan.md` calls this the likeliest way for a wrong result to reach print.

**Before submission**

7. Author list, CRediT contributions, competing interests, funding with grant IDs, MDAR checklist, Key Resources Table, running title, bioRxiv preprint (eLife only reviews papers already posted). → create `10_submission_pack.md`.
8. Article type: Research Article or Tools and Resources (D-6). Not cosmetic: T&R obliges benchmarking against existing methods and public deposition.
9. The cross-figure visual system (one algorithm→colour map, one font). Open in `figures/instructions.md` since March; every figure chose its own palette.
10. The generative-AI disclosure eLife requires in Methods.
11. Where the identifiability caveat lands (`discussion_plan.md`, a perimeter paragraph).

## 6. Keeping it in sync

The date column is not decoration. It is how you tell whether a document was written before or after the decision that invalidated it.

Three rules:

- **Stamp on edit.** Every document carries `Updated: YYYY-MM-DD` in its header. When you change it, change the stamp.
- **Cite, do not restate.** If a fact is owned elsewhere, write the pointer. A number copied into a second document is a number that will be wrong in one of them.
- **Propagate on change.** These four are load-bearing. When one changes, the documents listed after it must be re-read the same day:
  - `02_decision_log.md` → the master plan, every section plan, the `.tex`
  - `results_plan.md` → `figure_provenance.md`, every caption, the `.tex`
  - `docs/bibliography/MacroIR_prior_art_map.md` → `introduction_plan.md`, `abstract_draft.md`, `discussion_plan.md`
  - `nomenclature.md` → everything (a rename touches every figure, caption and script)

**Commit the pack.** Eleven of the documents this list depends on are untracked: `results_plan`, `methods_plan`, `figure_provenance`, `nomenclature`, `diagnostics_plan`, `theory_plan`, `discussion_plan`, `introduction_plan`, `abstract_draft`, `title_options`, `correction_idm_reconstruction`, `00_master_plan_v2`, `09_carve_plan`, `PROGRAM.md`, and the identifiability folder. They are one `git clean` from gone, and they are the newest and most load-bearing things in the project.

## 7. Needs the author

- **C-1**: build the K_off axis and re-run, or restate the maps on the axes the data actually have (N_ch × interval, N_ch × noise) and correct three documents. Recommended: restate. The interval axis is the one the paper's argument is about.
- **C-3**: spend compute to fill the missing Gaussian-Fisher cells (NMR at minimum), or ship the two-anchor split and disclose it. Recommended: fill them, if it costs less than a few thousand CPU-hours. Every headline number currently sits on the anchor the paper itself calls the less trustworthy one.
- **C-2** (title) and **C-4** (MNR vs NMR): both are one-line calls that unblock a lot of downstream text.
- **D-6**: Research Article or Tools and Resources.
- **Comm Biol / gvar_i**: does this paper say anything about the fact that the published P2X2 validation used the variance formula since corrected? A paper whose whole machinery exists to catch that class of error is the worst place for it to be unmentioned.
- **D-7**: this pack is in English and you work in Spanish. Keep it, or flip the planning layer.
