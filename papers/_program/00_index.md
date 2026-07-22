# Index: who owns what, across four layers

> Updated: 2026-07-20. Successor to `macroir-elife-2025/00_master_list.md`, which routed a single
> pack. This routes the shared layer plus three papers. It settles nothing about the science; it
> settles **where each thing is settled**.
>
> **State (2026-07-20): the moves are done.** `papers/macroir-elife-2025/` is now `papers/1_method/`;
> the shared documents have been moved into `_program/` with `git mv`, so history follows each file.
> What remains is `check.sh` (step 4) and the five content rewrites (step 5), listed in §5.
> Rows marked **REWRITE** are at their final path but still say something scoped to the old
> single-paper frame.

## The rules

1. **One topic, one owner. Everyone else cites.** A restatement is not a summary, it is a second copy
   that will drift.
2. **Citation runs one way.** A paper may cite `_program/`; `_program/` never cites a paper; papers do
   not cite each other in the planning layer. If `_program/` needs to cite a paper, the fact is filed
   in the wrong place. (`program.md` §8 — this is the rule that keeps four folders from becoming four
   copies of the same drift.)
3. **A decision is settled in exactly one place.** Cross-paper decisions in `_program/decisions.md`;
   a paper's own decisions in its folder. A document that argues toward a decision owns the
   *argument*; the *verdict* graduates to the log.
4. **A superseded document says so in its first three lines**, names its successor, and is never read
   as current again.
5. **Stamp on edit.** Every document carries `Updated: YYYY-MM-DD`. The date is how you tell whether a
   document was written before or after the decision that invalidated it.
6. **Commit on write.** Before 2026-07-14 fourteen of these documents were untracked, one `git clean`
   from gone. A planning document that is not committed does not exist.

**The completeness guarantee.** §2 is complete over the four layers: every working `.md` appears with
a status. If a file is not here, it is not part of the pack. §3 is exhaustive over cross-cutting
topics: if a topic is not there, it has no owner, which means **it is not decided**. So "I looked and
it is not in the index" is a reliable negative, not a failed search. `check.sh` enforces §2
mechanically, which is what makes this a check rather than a promise.

## 1. The four layers

| Layer | Holds | Rule |
|---|---|---|
| `_program/` | everything the three papers share: the machinery, the model, the axes, provenance, nomenclature, bibliography, submission mechanics | never cites a paper |
| `1_method/` | paper 1: R, MR, VR, IR — what a likelihood must condition on | nearly drafted |
| paper 2 | the usage map: LSE, NR, R, IR — do you need a likelihood at all | **a stub in `_program/paper-2.md` until it starts drafting** |
| paper 3 | the multinomial boundary: micro_R, micro_IR | **a stub in `_program/paper-3.md` until it starts drafting** |

**Papers 2 and 3 get folders when they get manuscripts, not before.** Structure built ahead of content
rots: `01_workboard.md` was born with 25 checkboxes and retired with all 25 unticked, including the
ones whose work had actually been done.

## 2. The registry

Status: **LIVE** = the owner, read as current. **REWRITE** = still the owner, currently says something
wrong. **TO MOVE** = lives at the old path, belongs at the new one. **REWRITE+MOVE** = both.
**CREATE** = does not exist. **RETIRED** = tombstone.

### `_program/` — the shared layer

| Document | Owns | Status | Note |
|---|---|---|---|
| `00_index.md` (this) | the topic→owner routing across four layers | LIVE | — |
| `program.md` | the three-paper map, the N_ch partition, publication order, citation directionality | LIVE | — |
| `research_program.md` | the **wider** scientific program this sits in (theory → model → likelihood → evidence → ranking), and Comm Biol's place as the demonstration | LIVE | was `From molecular … PROGRAM.md`. `program.md` is the three-paper map; this is the context around it |
| `decisions.md` | every settled **cross-paper** decision | LIVE | — |
| `machinery.md` | the validation machinery: definitions, sign conventions, thresholds, the decomposition identity `IDM = K·CDM·Kᵀ`, the hazards | LIVE | merged from `03_metrics_diagnostics.md` + `correction_idm_reconstruction.md` + the definitional half of `diagnostics_plan.md` |
| `axes.md` | N_ch, interval, noise; the two crossovers and the three bands; the label→value maps and their three divergent copies in ops; the per-cell protocol | REWRITE | was `05_experiment_grid.md`. Rewritten 2026-07-20 for the bands; still needs its §3–7 split, with paper-1-specific cells going to `1_method/grid.md` |
| `nomenclature.md` | the letters, their semantics, the data-key↔display-label map | REWRITE | the endpoint ladder needs the root question above it, and `VR` needs deciding |
| `provenance.md` | which run made which data, what is reproducible, the seed sentinel, the n_sims/Jensen hazard | LIVE | was `docs/figure_provenance.md` |
| `figures_system.md` | the cross-figure visual system: one algorithm→colour map, fonts, panel letters | REWRITE | was `figures/instructions.md`; open since March, every figure chose its own palette |
| `elife-author-instructions.md` | the verified eLife author requirements | LIVE | — |
| `notation_map.md` | symbol↔code↔CSV notation | LIVE | — |
| `sources.md` | which audio/chat is the source of record, and which wins when two disagree | LIVE | was `08_sources_audio_notes.md` |
| `carve_plan.md` | the repo boundary, the freeze trigger, code availability, engine work owed | REWRITE | "one paper = one repo" is reopened by the split (`decisions.md` §5) |
| `paper-2.md`, `paper-3.md` | the two stubs | LIVE | promote to folders when they start drafting |
| `model_and_sim.md` | `scheme_CO`, the emission model, exact CTMC uniformization, the six parameters and their units | CREATE | merge the model half of `1_method/methods.md` + `1_method/decisions/D-2_parameter_units.md` |
| `submission.md` | front/back matter, CRediT, data availability, MDAR, article type | CREATE | from `elife-author-instructions.md` + the missing-owner list |
| `check.sh` | the done-oracle + index completeness, **taking the paper as an argument** | LIVE | done 2026-07-20: it sits at `_program/check.sh`, takes the paper folder as its argument, and item 9 checks this index. Do **not** triplicate it |

**Not in `_program/`, but shared and cited by it:** `docs/bibliography/` lives at the **repo root**, not
in the papers tree, because the engine and theory layers cite it too. It owns the novelty position,
the Kalman concession and every citation attribution. Do not move it here.

### `1_method/` — paper 1

All at their final path. **REWRITE** means the file is where it belongs and still says something
scoped to the old single-paper frame.

| Document | Owns | Status | Note |
|---|---|---|---|
| `00_plan.md` | paper 1's thesis, scope, roster, the band-A results, its open decisions | REWRITE | was `00_master_plan_v2.md`; reframed 2026-07-20, but its §0 and §1a still speak for the whole program — those parts now belong to `_program/` |
| `decisions.md` | paper 1's own settled decisions | REWRITE | was `02_decision_log.md`; the cross-paper half has been copied to `_program/decisions.md` and must now be **deleted from here**, not left as a second copy |
| `01_writing_plan.md` | the drafting task graph, the lints, the human budget | REWRITE | hard-gates the definition of done at six figures while the count is open |
| `title.md` | the title argument | REWRITE | three live conflicting titles, and all describe a distortion study rather than paper 1's scope |
| `abstract.md` | the abstract and Impact Statement | REWRITE | **opens on least squares**, which is now paper 2's subject |
| `introduction.md` | the Introduction | REWRITE | the roster paragraph and the literature attribution both change; `R` now carries the anchor |
| `theory.md` | the Theory section | REWRITE | "one object with two knobs" is now false; `:85` still calls the noise parameter a misnomer, which D-2 overturned |
| `diagnostics.md` | the Diagnostics **section**: how it is written and argued | REWRITE | the definitions themselves now live in `_program/machinery.md`; this keeps the prose strategy and the source lift-readiness table |
| `results.md` | the Results and **paper 1's figure arc** | REWRITE | its headline numbers row is band A only |
| `discussion.md` | the Discussion and the decision rule the paper is cited for | REWRITE | closest to the new frame already, since it asks for a decision rule rather than a winner |
| `methods.md` | Materials and Methods | REWRITE | the model and units halves move out to `_program/model_and_sim.md` |
| `README.md` | routes into paper 1; points at the program layer | LIVE | — |
| `00_master_plan.md`, `01_writing_plan.md`, `04_figures_storyboard.md`, `06_repro_pipeline.md` | retired/pointer stubs carried over from the pack | RETIRED/REWRITE | `04_*` and `06_*` are dead arcs; `01_writing_plan` is REWRITE (six-figure gate); `00_master_plan` is a tombstone |
| `analysis_figure_S1_score_mean.md` | what the score-mean figure shows | LIVE | — |
| `figures_build_plan.md` | the **order** of figure runs and edits, and the command for each | LIVE | opened 2026-07-21 for the VR re-runs; owns build order only — the arc is `results.md`, the visual system `_program/figures_system.md`, the run manifest `_program/provenance.md` |
| `decisions/D-0, D-3, D-4` | freeze scope, the novelty claim, the ranking verdict | LIVE | D-4 rescopes to band A rather than being rewritten |
| `docs/manuscript-drafts/` | the vessel. Owns nothing; every claim in it is owned upstream | LIVE | — |
| `grid.md` | paper 1's cells: which N_ch, which noise, which methods | CREATE | from `_program/axes.md` §3–7 once the N_ch ranges and `VR` are settled |
| `components/`, `archive/` | parked proof-of-concept and tombstones | RETIRED | never fact owners |

### Retired, do not read as current

`00_master_plan.md`; `06_repro_pipeline.md`; `archive/07_code_tasks.md`; `04_figures_storyboard.md`
(the 4-figure arc is dead); `components/` (a **parked** proof-of-concept, never a fact owner);
`archive/components_poc/`; `00_master_list.md` (superseded by this file once the moves land).

## 3. Topic index

Look the topic up here before writing.

| Topic | Owner |
|---|---|
| The three papers, what each owns, publication order | `_program/program.md` |
| The N_ch partition across papers | `_program/program.md` §2 |
| The noise bands and the two crossovers | `_program/axes.md` |
| Diagnostic definitions, thresholds, sign conventions, CSV fields | `_program/machinery.md` |
| The decomposition identity | `_program/machinery.md` |
| The model, the simulator, the parameter units | `_program/model_and_sim.md` |
| The method names and the letter semantics | `_program/nomenclature.md` |
| Which data commit holds what; what is reproducible | `_program/provenance.md` |
| Novelty, prior art, the Kalman concession, citations | `docs/bibliography/` (repo root, not under `_program/`) |
| The repo boundary, the freeze, code availability | `_program/carve_plan.md` |
| Everything a journal asks for on submission day | `_program/submission.md` |
| Cross-paper settled decisions | `_program/decisions.md` |
| Paper 1's thesis, scope, roster | `1_method/00_plan.md` |
| Paper 1's figure arc | `1_method/results.md` |
| Paper 1's cells | `1_method/grid.md` |
| What paper 2 is, until it has a folder | `_program/paper-2.md` |
| What paper 3 is, until it has a folder | `_program/paper-3.md` |

## 4. Carried over from the single pack

These were live collisions or missing owners before the split and none of them was closed by it. They
move with the topic.

**To `_program/`:** the seed sentinel (`seed = 0` means *random*, the resolved value was never logged,
so the ensembles are statistically equivalent but not bit-reproducible) → `provenance.md`; the
n_sims-uniformity hazard, now sharper with the 2026-07-20 fill at n_sims 1000 beside band-A cells at
10000 → `provenance.md`; the three divergent noise label→value maps in ops → `axes.md`; the units of
the six parameters, the one gap only the author can close → `model_and_sim.md`; the bibliography,
which does not exist (the `.tex` has zero `\cite` and `biblio.bib` is the P2X2 paper's 170 references)
→ `bibliography/`; data availability, which eLife makes mandatory and which has no repository or DOI
named anywhere → `submission.md`; the cross-figure visual system, open since March →
`figures_system.md`.

**To `1_method/`:** the direction convention of the distortion ratio, asserted in three places with
two different signs, on which every verdict depends; the ranking table's two contested cells (D-4);
figure numbering against what is on disk; the title, in three live conflicting versions.

**New, from the split:** the 2026-07-20 dispatch does not match paper 2's roster (`program.md` §4);
`87889e6` is recorded as "micro, out of scope" but holds both the micro runs and the macro fill
(`program.md` §6); `01_writing_plan.md` hard-gates the definition of done at six figures while the
figure count is formally open.

## 5. Order of operations

1. ~~`program.md` + this index.~~ **Done 2026-07-20.**
2. ~~The stubs and merges: `paper-2.md`, `paper-3.md`, `decisions.md`, `machinery.md`.~~ **Done.**
3. ~~The moves, with `git mv` so history follows the file.~~ **Done.** `macroir-elife-2025` →
   `1_method`; twelve shared documents into `_program/`; ten section plans renamed; three merged or
   superseded documents tombstoned into `1_method/archive/`.
4. ~~`check.sh`, retargeted to four layers and taking the paper as an argument.~~ **Done.** It sits at
   `_program/check.sh`, takes the paper folder as its argument, and item 9 enumerates the shared layer
   plus that paper.
5. The five real rewrites, one at a time, each touching a thesis: `00_plan.md` and `1_method/decisions.md`
   (both must shed the parts that moved to `_program/`, or the split has just manufactured two copies
   of everything it separated), `theory.md`, `nomenclature.md`, `abstract.md`.

**The immediate risk, named so it is not forgotten.** The cross-paper content was **copied** into
`_program/decisions.md`, not moved. Until step 5 deletes it from `1_method/decisions.md` and from
`00_plan.md`, the program has exactly the duplication this restructure exists to prevent. Do step 5
before anything else builds on either file.

**Not now:** `projects/eLife_2025/`. It is also named for a venue and a year, but every `.Rmd` and
every dispatcher hardcodes paths into it and there are jobs in flight. After the freeze.
