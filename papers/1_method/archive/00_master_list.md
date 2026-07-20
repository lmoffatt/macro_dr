> **RETIRED 2026-07-20. Superseded by `../../_program/00_index.md`.**
> It routed a single pack; the program now has four layers (_program + three papers).
> Do not read as current and do not edit. Kept for its content; the live owner is the successor.

---

# Master list: who owns what

> Opened 2026-07-14. Updated 2026-07-20. This is the routing document for the MacroIR / eLife paper. It settles nothing about the science; it settles **where each thing is settled**.
>
> Built from a full read of the 26 documents in this pack plus the theory, figure, bibliography and source layers they depend on.

> ## ⚠ Reframe in flight (2026-07-20): the usage map
>
> The paper was reframed on 2026-07-20 (audios 10:54–10:59). Two changes, both settled in `02_decision_log.md`:
>
> 1. **Least squares on the mean current enters as method zero** (`nonlinearsqr`, display `LSE`). The roster is six methods on two levels, not five on one lattice.
> 2. **The noise axis is swept through the gating-noise crossover.** Everything measured before that date sat in the band that favours the gating-aware likelihoods by construction.
>
> Consequently the deliverable is a **usage map** (cheapest adequate method per regime), not a ranking. `00_master_plan_v2.md` §0, §1, §1a and §4 carry the new frame.
>
> **Until the propagation in §4a is done, treat as stale on sight:** any sentence saying "five algorithms", "the family as one object with two knobs", "IR sole survivor", "only MacroIR stays calibrated", or "MR strawman". These are not errors of fact about band A; they are claims scoped to one band and stated globally.

## The rule

**One topic, one owner. Everyone else cites.**

If you are about to write a fact, a number, a scope call or a verdict into a document, first find the topic in §3. If another document owns it, do not restate it: reference it. A restatement is not a summary, it is a second copy that will drift, and drift is what makes the same argument get had twice.

Three consequences worth naming, because the pack currently violates all three:

- A **section plan** (`results_plan.md`, `methods_plan.md`, and the rest) owns its section's argument and prose. It does **not** own cross-cutting decisions (scope, naming, the ranking, figure numbering). It consumes them.
- A **decision** is settled in exactly one place, `02_decision_log.md`. A document that argues toward a decision (`nomenclature.md`, `title_options.md`) owns the *argument*; the *verdict* graduates to the log.
- A document that is superseded says so **in its first three lines**, names its successor, and is never read as current again.

**The completeness guarantee (this is what lets an agent trust the index).** §2 is **complete over the pack**: every file appears there with a status. If a file is not in §2, it is not part of the pack. §3 is **exhaustive over cross-cutting decisions**: if a topic is not in §3, it has no owner yet, which means **it is not decided**. So "I looked and it is not in the index" is a reliable negative, not a failed search. `check.sh` enforces §2 completeness mechanically (every `.md` in the pack must be registered), so the guarantee is a check, not a promise.

**On the `components/` experiment.** The atomic component-file scheme (`components/_SPEC.md`) is a **PARKED** proof-of-concept, not a live owner. It earns propagation only when `components/_LOG.md` records that this index failed (see the D-0 rerun refresh test). Until then, the owners are the documents below; `components/*` fact nodes are not canonical.

## 2. The registry

Status vocabulary: **LIVE** = the owner, read it as current. **REWRITE** = still the owner, but it currently says something wrong (see §4). **POINTER** = kept only to redirect. **RETIRED** = tombstone. **CREATE** = does not exist yet.

`git` column: everything listed here is in version control as of 2026-07-14. Keep it that way. Until that day, most of the documents everything else depends on were untracked, one `git clean` from gone.

### A. Governance

| Document | Owns | Status | Updated | git |
|---|---|---|---|---|
| `00_master_list.md` (this) | The topic→owner routing table; the document lifecycle | LIVE | 2026-07-14 | committed |
| `01_writing_plan.md` | **The order the manuscript gets drafted**: the gates upstream of prose, what each section needs before its first sentence, and when a section is done | LIVE | 2026-07-14 | committed |
| `00_master_plan_v2.md` | The thesis, the two-approximation frame, the three regimes, the open-decision register (D-1…D-7). **No longer the task list** (§9 is now a pointer) | REWRITE | 2026-07-14 | committed |
| `02_decision_log.md` | Every **settled** cross-cutting decision: scope, algorithm roster, the anchor, the ranking verdict, the repo boundary | REWRITE | 2026-07-13 | committed |
| `08_sources_audio_notes.md` | Which audio/chat is the source of record, and which wins when two disagree | LIVE | 2026-07-13 | committed |
| `09_carve_plan.md` | The repo boundary, the freeze trigger, what moves to `macroir-validity`, code availability, and **the engine work that must land in the frozen commit** | LIVE | 2026-07-14 | committed |
| `README.md` | Nothing. It routes | POINTER | 2026-07-14 | committed |
| `check.sh` | The mechanical done-oracle + index-completeness check (enforces §2) | LIVE | 2026-07-15 | uncommitted |
| `decisions/D-0…D-4*.md` | The four decision briefs (freeze/rerun, units, novelty, ranking), each researched with a recommendation, awaiting Luciano | LIVE | 2026-07-15 | uncommitted |
| `components/_SPEC.md`, `_LOG.md`, `_MAP.md` | The parked atomic-component scheme (methodology + does-it-work log). **PARKED POC, not a fact-owner** | PARKED | 2026-07-15 | uncommitted |
| `archive/components_poc/` | The three demo nodes (MR/sandwich/fact) that proved the scheme; archived, not canonical | RETIRED | 2026-07-15 | uncommitted |
| `AGENTS.md` (repo root) | Where an agent starts | LIVE | 2026-07-14 | committed |
| `00_master_plan.md` | — superseded by `00_master_plan_v2.md` | RETIRED | — | committed |
| `06_repro_pipeline.md` | — moved to `09_carve_plan.md` and `docs/figure_provenance.md` | RETIRED | — | committed |
| `archive/07_code_tasks.md` | — archived; its live items are the freeze preconditions in `09_carve_plan.md`. Its rule "we do not touch code" made it the owner of the code that had to be touched | RETIRED | — | committed |
| ~~`01_workboard.md`~~ | — deleted 2026-07-14. Born March with 25 checkboxes, retired July with every one still unticked, including the work that had actually been done. In git history | GONE | — | — |

### B. Argument and scope

| Document | Owns | Status | Updated | git |
|---|---|---|---|---|
| `From molecular … PROGRAM.md` | This paper's position in the program (bridge 2; Comm Biol is the demonstration) | LIVE | 2026-07-13 | committed |
| `nomenclature.md` | The five names, the letter semantics, the MNR/NMR spelling, the data-key↔display-label map | REWRITE | 2026-07-14 | committed |
| `title_options.md` | The title argument and the project's claim/accuracy rules (the rules do not belong here, see §5) | REWRITE | 2026-07-14 | committed |
| `docs/bibliography/MacroIR_prior_art_map.md` | The novelty position, the Kalman concession, every citation attribution | LIVE | 2026-07-14 | committed |
| `docs/bibliography/identifiability/README.md` | The identifiability perimeter (what macroscopic data cannot recover) | LIVE | 2026-07-13 | committed |

### C. Science content

| Document | Owns | Status | Updated | git |
|---|---|---|---|---|
| `03_metrics_diagnostics.md` | The metric **registry**: each diagnostic's definition, its sign convention, its threshold, and which CSV field carries it | REWRITE | 2026-07-13 | committed |
| `05_experiment_grid.md` | The sweep design: control variables, the definition of τ, the per-cell protocol | REWRITE | 2026-07-13 | committed |
| `correction_idm_reconstruction.md` | The correct distortion-decomposition identity (`IDM = K·CDM·Kᵀ`) and what composes exactly | LIVE | 2026-07-14 | committed |
| `theory/…/Gaussian_Fisher_Distortion_Family.md` | What the anchor **H** is, and the numerical-vs-Gaussian Fisher bridge | LIVE | 2026-07-05 | committed |
| `theory/…/score_martingale_argument.md` | Why recursion fixes the reported uncertainty without fixing the bias | LIVE | — | committed |
| `theory/…/Macro_IR/` | The algorithm derivation feeding Theory and the supplement | LIVE | — | committed |
| `analysis_figure_S1_score_mean.md` | What the score-mean figure shows (and only that) | LIVE | 2026-06-30 | committed |
| `figures/paper/sample_correlation_distortion_analysis.md` | What the sample/correlation scaling laws do and do **not** support | LIVE | 2026-07-09 | committed |

### D. Figures and data

| Document | Owns | Status | Updated | git |
|---|---|---|---|---|
| `results_plan.md` | **The figure arc**: which figures exist, in what order, what claim each carries, the numbering register, the do-not-say list | LIVE | 2026-07-14 | committed |
| `docs/figure_provenance.md` | **The run manifest**: which script produced which data, what is built, what is not, what is reproducible | LIVE | 2026-07-14 | committed |
| `figures/in_progress/figure_5_PLAN.md` | Panel selection for Figure 5 (main / supplementary / drop), and nothing else | LIVE | 2026-07-08 | committed |
| `figures/paper/*_caption.md` | The caption text of the figure each is named for | REWRITE | various | committed |
| `figures/instructions.md` (this pack) | Figure assembly and the cross-figure visual system (colours, fonts, panel letters) | REWRITE | 2026-03-25 | committed |
| `04_figures_storyboard.md` | — the 4-figure arc is dead; the panel sketches survive only as reference | POINTER | 2026-07-13 | committed |

### E. Manuscript sections

Each owns its section's job, structure, argument and prose. None owns a cross-cutting decision.

| Document | Owns | Status | Updated | git |
|---|---|---|---|---|
| `abstract_draft.md` | The abstract and the Impact Statement, while drafting | REWRITE | 2026-07-14 | committed |
| `introduction_plan.md` | The Introduction | REWRITE | 2026-07-14 | committed |
| `theory_plan.md` | The Theory section (the family as one object with two knobs) | LIVE | 2026-07-14 | committed |
| `diagnostics_plan.md` | The Diagnostics section (the validation machinery as a procedure a sceptic could run) | LIVE | 2026-07-14 | committed |
| `results_plan.md` | The Results (see D: it also owns the arc) | LIVE | 2026-07-14 | committed |
| `discussion_plan.md` | The Discussion, and the decision rule the paper is cited for | LIVE | 2026-07-14 | committed |
| `methods_plan.md` | Materials and Methods, the emission model, and every run value | REWRITE | 2026-07-14 | committed |
| `docs/manuscript-drafts/elife_paper.tex` | Nothing. It is the **vessel**, not a second author. Every claim in it is owned upstream | LIVE | 2026-07-13 | committed |

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
| **The usage-map frame: the root question + the ladder** | `00_master_plan_v2.md` §1 |
| **The noise bands A/B/C and where the two crossovers fall** | `00_master_plan_v2.md` §1a; the sweep in `05_experiment_grid.md` §2a |
| **LSE / `nonlinearsqr`: what it is, its seams, its config constraints** | `theory/macroir/notes/nonlinearsqr_lse_plan.md`; the *decision to include it* in `02_decision_log.md` |
| **Which diagnostics are meaningless or tautological for LSE** | `05_experiment_grid.md` §6 |
| What is in and out of scope | `02_decision_log.md` |
| The method roster (six, on two levels) | `02_decision_log.md` |
| The names of the methods, and the spelling | `nomenclature.md` → verdict to `02_decision_log.md` |
| The band-A column of the map (the former ranking) | `00_master_plan_v2.md` §4; the two contested cells in `decisions/D-4_ranking_verdict.md` |
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
| The engine work that must land before the frozen build | `09_carve_plan.md` (Freeze preconditions) |
| The order the manuscript gets written, and what gates each section | `01_writing_plan.md` |
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

## 4a. Propagation debt from the 2026-07-20 reframe

Done so far: `00_master_plan_v2.md`, `02_decision_log.md`, `05_experiment_grid.md`, this file. Everything below still says the old thing. Ordered by blast radius; each line is a site, not a summary.

**The structural rewrites (the "two knobs" framing is now literally false: LSE carries NMR's two knob settings and differs only by the third flag `family_approximation=2`).**

| Document | Sites |
|---|---|
| `nomenclature.md` | `:5` "the five macroscopic likelihood algorithms"; `:15` "## The family, as implemented" + the closed 5-row table `:23-29`; `:84` "the five members do that"; `:114-125` the endpoint ladder, which has **no rung for LSE** and needs the root question above it; `:127-128` "sole survivor" |
| `theory_plan.md` | `:9` "the five algorithms **one object with two knobs**"; `:15` "the members of the family differ in *what they condition that Gaussian on*" (LSE is a counterexample); `:44` "T3 — the family as a coordinate system"; `:55` "the five members are exactly the implemented grid"; `:59` "MacroIR is the top rung below intractability" (survives, but now as the top rung of the *lower* level); `:85` still calls the noise parameter "a misnomer", which `decisions/D-2_parameter_units.md:41-44` explicitly overturned |
| `methods_plan.md` | `:79` "## M6 — The five algorithms, as flags"; `:81` "one code with two flags"; `:91` "the family is a genuine grid, not a selection"; `:61-67` the noise-axis "naming trap" framing, also overturned by D-2; `:134`, `:136` "all five algorithms" |
| `sections/02_theory.tex` | `:2`, `:23`, `:45`, `:73`, `:75`, `:110`, `:115`, `:199` — the whole section is built on the two-knob coordinate system |

**The verdict rewrites.**

| Document | Sites |
|---|---|
| `02_decision_log.md` | done |
| `discussion_plan.md` | `:22-38` D1 "the verdict, as a decision rule" — closest to the new frame already, since it asks for a decision rule rather than a winner; `:74` "survives every approximation in the family" |
| `decisions/D-4_ranking_verdict.md` | `:28-51` the canonical recomputed table; rescope to band A rather than rewrite |
| `sections/05_discussion.tex` | `:2` the comment-form ranking summary |
| `sections/04_results.tex` | `:5`, `:11`, `:26`, `:41`, `:45`, `:49`, `:51` — "the five approximations", "all five members", the band-A-only verdict sentence |
| `components/_MAP.md` | `:12-16` the 5-bullet roster, `:16` "Sole calibrated survivor" |
| `From molecular … PROGRAM.md` | `:62` "this paper uses NR, NMR, R, MR, IR" |

**The front-matter rewrites (LSE moves from cited background to measured arm — it is currently the *premise* of the abstract's first sentence).**

| Document | Sites |
|---|---|
| `abstract_draft.md` | `:79` the 360-word draft, `:114` version A, `:118` version B — all three open on least squares as the thing the field does wrong, and must now open on it as the thing the paper measures; the Impact Statement drafts `:139-143`; the whole "what stays out" list needs the band-C result added |
| `introduction_plan.md` | `:56` "the algorithm family is the literature, not an invention" — **this is where LSE belongs**, and it has the strongest literature attribution of the lot (Clerx 2019, IonBench, Del Core & Mirams 2025); `:87` "for the five members of the family"; `:110` "five algorithms look like a coordinate system" |
| `title_options.md` | the CHOSEN title and its two rivals all describe a distortion study, not a usage map. Re-open (already collision C-2) |
| `elife_paper.tex:23` | the live abstract carries "Only the interval-averaged recursive likelihood (MacroIR) stays calibrated across the practical regime" |
| `results_plan.md` | `:10-17` the six-figure through-line; `:99` the headline numbers row, which is band A only; `:12-16` per-figure claims |

**The caption rewrites.** Every caption in `projects/eLife_2025/figures/paper/` enumerates the closed five: `Figure_1_caption.md:3,5,7`, `Figure_2_caption.md:3,5,7`, `Figure_3_caption.md:3,5,7`, `Figure_4_caption.md:3,5,11`, `Figure_4_alt_caption.md:5`, `Figure_S2_caption.md:5`, `Figure_S4_acf_caption.md:5`, `Figure_S4_bias_caption.md:5`. Two need arithmetic, not just wording: `Figure_S3_caption.md:3,5` says "the four approximations other than IR", and `Figure_3_caption.md:7` partitions the five exhaustively as "the recursive trio" and "the naive pair".

**Not a document, but it will bite:** `01_writing_plan.md:25` hard-gates the definition of done at exactly six figures, while `00_master_plan_v2.md` §8 D-2 leaves the count open. Resolve toward D-2.

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

- **Commit on write.** The whole pack went into git on 2026-07-14 (commits `0419c74`…`d2eb71e`). Before that, fourteen of the documents this list depends on were untracked, including every section plan, the run manifest and the master plan itself. A planning document that is not committed does not exist.

## 7. Needs the author

This section routes; it does not sequence. **The order in which any of it happens is `01_writing_plan.md`.**

Each live collision in §4 names the document that decides it. Beyond those, three calls are the author's alone and appear nowhere else:

- **The units of the six parameters.** They exist nowhere in the repo and there is no substitute path. It is an afternoon, and it blocks Methods. (`01_writing_plan.md` G-1.)
- **Research Article or Tools and Resources** (D-6). The only open decision that *adds* work rather than reordering it: T&R obliges benchmarking against existing methods. (`01_writing_plan.md` §5.)
- **Comm Biol / gvar_i.** Does this paper say anything about the fact that the published P2X2 validation used the variance formula since corrected? A paper whose whole machinery exists to catch that class of error is the worst place for it to go unmentioned.

Still open and unrouted: **D-7**, the language of the planning layer. The pack is in English and you work in Spanish.

Two of the recommendations in the first version of this section have since been settled by measurement, and are recorded here so they are not re-argued:

- **C-1 (map axes).** Recommendation stands: restate the maps on the axes the data actually have. No `K_off` axis exists in the scripts (`figure_3_mle_G.macroir:46` hard-fixes `off = 100`), the output path has no `K_off` component so cells would silently overwrite each other, and building it would cost roughly 35,000 CPU-hours for one supplementary figure.
- **C-3 (the anchor).** No longer an expensive decision. The gap is **four cells** (NMR × N_ch {10, 100, 1000, 10000} at noise 0.1), about 1,830 CPU-hours, one to two days of queue on a cluster with no CPU-hour quota. Fill them. But read `09_carve_plan.md` (Freeze preconditions) first: the engine fixes change the commit hash, and the hash is the data directory's name.
