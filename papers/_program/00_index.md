# Index: who owns what, across four layers

> Updated: 2026-07-20. Successor to `macroir-elife-2025/00_master_list.md`, which routed a single
> pack. This routes the shared layer plus three papers. It settles nothing about the science; it
> settles **where each thing is settled**.
>
> **State: the layout below is decided, the moves are not done.** Only this file and `program.md`
> exist in `_program/` so far. Everything else still sits in `papers/macroir-elife-2025/`. Rows
> marked **TO MOVE** name the source. Until a row is moved, read it at its old path.

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

| Document | Owns | Status | Source if moving |
|---|---|---|---|
| `00_index.md` (this) | the topic→owner routing across four layers | LIVE | — |
| `program.md` | the three-paper map, the N_ch partition, publication order, citation directionality | LIVE | new; absorbs `From molecular … PROGRAM.md` |
| `decisions.md` | every settled **cross-paper** decision | CREATE | split from `02_decision_log.md` |
| `nomenclature.md` | the letters, their semantics, the data-key↔display-label map | REWRITE+MOVE | `macroir-elife-2025/nomenclature.md` — the endpoint ladder needs the root question above it, and `VR` needs deciding |
| `machinery.md` | the validation machinery: diagnostic definitions, sign conventions, thresholds, CSV fields, the decomposition identity `IDM = K·CDM·Kᵀ` | CREATE | merge `03_metrics_diagnostics.md` + `correction_idm_reconstruction.md` + the non-section half of `diagnostics_plan.md` |
| `model_and_sim.md` | `scheme_CO`, the emission model, exact CTMC uniformization, the six parameters and their units | CREATE | the model half of `methods_plan.md` + `decisions/D-2_parameter_units.md` |
| `axes.md` | N_ch, interval, noise; the two crossovers and the three bands; the label→value maps and their three divergent copies in ops | REWRITE+MOVE | `05_experiment_grid.md` §1–2b (already rewritten 2026-07-20) |
| `provenance.md` | which run made which data, what is reproducible, the seed sentinel, the n_sims/Jensen hazard | TO MOVE | `docs/figure_provenance.md` |
| `figures_system.md` | the cross-figure visual system: one algorithm→colour map, fonts, panel letters | REWRITE+MOVE | `figures/instructions.md` — open since March, every figure chose its own palette |
| `submission.md` | eLife mechanics, figure guidelines, front/back matter, CRediT, data availability | CREATE | `docs/elife-author-instructions.md` + the missing-owner list |
| `bibliography/` | the novelty position, the Kalman concession, every citation attribution | TO MOVE | `docs/bibliography/` |
| `sources.md` | which audio/chat is the source of record, and which wins when two disagree | TO MOVE | `08_sources_audio_notes.md` |
| `carve_plan.md` | the repo boundary, the freeze trigger, code availability, engine work owed | TO MOVE | `09_carve_plan.md` |
| `paper-2.md` | paper 2's stub: question, roster, axis, what is on disk, what it needs | CREATE | — |
| `paper-3.md` | paper 3's stub: same | CREATE | — |
| `check.sh` | the done-oracle + index completeness, **taking the paper as an argument** | REWRITE+MOVE | `macroir-elife-2025/check.sh`. Do **not** triplicate it |

### `1_method/` — paper 1

| Document | Owns | Status | Source |
|---|---|---|---|
| `00_plan.md` | paper 1's thesis, scope, roster, the band-A results, its open decisions | REWRITE+MOVE | `00_master_plan_v2.md` (rewritten 2026-07-20; the program parts lift out to `program.md`) |
| `01_writing_plan.md` | the drafting task graph, the lints, the human budget | TO MOVE | same name |
| `decisions.md` | paper 1's own settled decisions | CREATE | the paper-specific half of `02_decision_log.md` |
| `grid.md` | paper 1's cells: which N_ch, which noise, which methods | CREATE | `05_experiment_grid.md` §3–7 |
| `title.md`, `abstract.md` | its front matter | REWRITE+MOVE | `title_options.md`, `abstract_draft.md` — the abstract currently **opens on least squares**, which is now paper 2 |
| `introduction.md` | the Introduction | REWRITE+MOVE | `introduction_plan.md` — the roster paragraph and the literature attribution both change |
| `theory.md` | the Theory section | REWRITE+MOVE | `theory_plan.md` — "one object with two knobs" is now false, and `:85` still calls the noise parameter a misnomer, which D-2 overturned |
| `diagnostics.md` | the Diagnostics **section** (the machinery itself moves to `_program/machinery.md`) | REWRITE+MOVE | `diagnostics_plan.md` |
| `results.md` | the Results and **the figure arc** | REWRITE+MOVE | `results_plan.md` |
| `discussion.md` | the Discussion and the decision rule the paper is cited for | REWRITE+MOVE | `discussion_plan.md` |
| `methods.md` | Materials and Methods (minus the model, which moves to `_program/`) | REWRITE+MOVE | `methods_plan.md` |
| `docs/manuscript-drafts/` | the vessel. Owns nothing; every claim in it is owned upstream | TO MOVE | same |

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
| Novelty, prior art, the Kalman concession, citations | `_program/bibliography/` |
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

1. `program.md` + this index. **← you are here**
2. The four stubs and merges that have no source to move: `paper-2.md`, `paper-3.md`, `decisions.md`,
   `machinery.md`.
3. The moves, with `git mv` so history follows the file, **in batches by topic, not all at once**.
4. `check.sh`, retargeted to four layers and taking the paper as an argument.
5. The five real rewrites, one at a time, each touching a thesis: `00_plan.md`, `decisions.md`,
   `theory.md`, `nomenclature.md`, `abstract.md`.

**Not now:** `projects/eLife_2025/`. It is also named for a venue and a year, but every `.Rmd` and
every dispatcher hardcodes paths into it and there are jobs in flight. After the freeze.
