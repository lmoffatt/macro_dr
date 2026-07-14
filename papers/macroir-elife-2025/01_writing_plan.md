# Writing plan: the order the manuscript gets drafted

> Opened 2026-07-14. Owns **the order in which the manuscript is written, what each section needs before its first sentence, and when a section is done.** It does not own what any section *says*: that belongs to the section plans (`theory_plan.md`, `results_plan.md`, and the rest). It owns the sequence and the gates.
>
> Replaces `00_master_plan_v2.md` §9 and `01_workboard.md`. Engine work lives in `09_carve_plan.md` (Freeze preconditions); submission paperwork is out of scope (see §6).

## Why this is not a checklist

`01_workboard.md` was born in March with about twenty-five checkboxes and retired in July with **every one of them still unticked**, including the ones for work that had actually been done: the figures were built, the grids were run on Dirac, the score-mean analysis was written in June. Two commits in its whole life and no edit in between. It did not die because tracking work is a bad idea. It died because it recorded **state**, and the only thing that kept the state true was somebody remembering to open a file that nothing forced them to open.

So this document records what does not rot: **preconditions and definitions of done.** "Results §5 cannot be written until Figure 5 is selected" stays true whatever anybody does. "[ ] knit figure_2.Rmd" is false the moment somebody knits it and forgets to say so.

If you want to know what is done, look at the manuscript.

## 1. The strategy

Draft the whole paper against what exists now, leaving the two known holes with a **defined shape**, and polish afterwards. This is viable because the section plans are already written and good, and because four of the six figures are built and captioned. The draft is mostly transcription with judgment, not invention.

Two rules that make "draft now, polish later" safe rather than reckless:

- **No number enters the draft without a source.** Every figure, factor and threshold cites the file it was computed in. A number typed from memory into the Abstract propagates into the Introduction, the Results and the Discussion, and by then nobody knows which copy is the original. `results_plan.md` carries the do-not-say list; respect it.
- **A hole is written as a hole.** Where a figure does not exist, write the paragraph that will surround it and mark the claim `\todo{}`. Do not write a claim you cannot yet support and plan to check later. That is how the wrong sentence gets published.

## 2. The gates

Four things are upstream of prose. Each costs hours, and each poisons a section if you write before settling it.

| Gate | What it blocks | Who can settle it | Cost |
|---|---|---|---|
| **G-1. The units of the six parameters.** They are written nowhere in the repo: not in the model header, not in a `.macroir` script, not in a dispatcher (`methods_plan.md` M-units). | Methods, and the model table | **Only Luciano.** There is no substitute path. | an afternoon |
| **G-2. The novelty sentence.** `introduction_plan.md` still carries "no published likelihood for macroscopic currents integrates the observable over the acquisition window". The prior-art map retracts it: single-channel HMMs do exactly that. The novelty is the **scaling to channel populations**, plus the exact CTMC realization and the validation. | The Introduction, the Abstract, the Discussion's novelty paragraph | `docs/bibliography/MacroIR_prior_art_map.md` already has the replacement framing | an afternoon of reframing |
| **G-3. The overconfidence factor.** `abstract_draft.md` says 10-16 on line 158 and 14-21 on line 162. One of them is wrong. | The Abstract, the Results headline, the Discussion | Read `Figure_S3_caption.md`, which is where it is computed | an hour |
| **G-4. The ranking table.** It exists in four copies (master plan, decision log, `discussion_plan.md`, and the `.tex` fill-hint), and two cells disagree with the data: MR's variance direction, and IR's corner bound (≤1.3 stated, ~0.5 observed). Probably a category error between the predicted *observable* variance and the reported *parameter* covariance. | Results' verdict, the whole Discussion, Table 1 | `results_plan.md` decides; the verdict goes to `02_decision_log.md` and the other three copies become citations | a day, including the check against the data |

Settle these four before drafting the sections they feed. Nothing else is upstream of prose.

## 3. The order

Sections are drafted in dependency order, not manuscript order. Each entry gives what it needs, what it draws on, and what "done" means.

### 1st. Theory (`theory_plan.md`)

The only section with no empirical dependency, and every later section refers to its objects. Write it first so the vocabulary is fixed.

- **Needs:** nothing. `nomenclature.md` must be assumed, not re-derived (if the MNR/NMR spelling flips later it is a `sed`, not a rewrite).
- **Draws on:** `theory/macroir/docs/Macro_IR/macroir_derivation.tex` and `macroir_macroir_paper_section.md`.
- **Done when:** a reader can say "there are two Gaussian approximations, and the family differs only in how much of the acquisition interval the hidden state is conditioned on", and every symbol used later in the paper has been defined here.
- **May not claim:** that the recursive members' Gaussian Fisher is the exact information matrix. It is not; the likelihood is misspecified by construction. (`docs/corrected covariance justification.md` says otherwise and is wrong on this point.)

### 2nd. Diagnostics (`diagnostics_plan.md`)

The paper's genuinely new contribution, and the section most exposed to a methods referee.

- **Needs:** E-2 landed, or at least the correct identity known. The supplement and `likelihood.cpp:3501` both print `IDM = SDM^½ · CDM · SDM^½`, which is **false**. The exact form is `IDM = K·CDM·Kᵀ` (`correction_idm_reconstruction.md`). Verified 2026-07-14: fixing it **invalidates no figure already built**, so this is a correction to an equation, not a rebuild.
- **Draws on:** `theory/macroir/docs/Likelihood_Information_Distortion/supplement_information_distortion_main.tex`, `Gaussian_Fisher_Distortion_Family.md` (which owns what the anchor **H** is), `03_metrics_diagnostics.md`.
- **Done when:** a sceptic could run the three tests on their own likelihood tomorrow, and the text is explicit about which identities are exact algebra, which are first-order, and which are convention.

### 3rd. Results §1 to §4 (`results_plan.md`)

Figures 1 to 4 are built, captioned and on disk. This is the half of the Results that is fully writable today.

- **Needs:** G-3 and G-4.
- **Draws on:** the four captions, `docs/figure_provenance.md` (which run produced which figure, and the three naming traps that have already misled readers).
- **Done when:** every claim is traceable to a figure panel and every number to the file that computed it; each figure answers the objection raised by the previous one.
- **Must state, once:** which Fisher anchors which panel. The definitive figures were promised on the Gaussian anchor and Figures 2, S2 and S3 currently sit on the numerical one (`433ed13`).

### 4th. Methods (`methods_plan.md`)

eLife excludes Methods from the word count, so it can afford to be complete, and it should be: the paper's claims are about algorithms, so the Methods *are* the result. A vague Methods section here does not inconvenience the reader, it makes the paper unfalsifiable.

- **Needs:** G-1 (the units).
- **Draws on:** `methods_plan.md` M1-M11 and `docs/figure_provenance.md` (the computational path).
- **Done when:** someone in 2029 could re-run this. Two things must be said plainly rather than hidden: that `seed = 0` meant *random* and the resolved seed was never logged, so the ensembles are statistically equivalent but **not bit-reproducible** (`methods_plan.md` M5 currently claims the opposite and must be corrected); and that the ground-truth simulator is exact CTMC uniformization, which is what licenses using it as the reference.

### 5th. Discussion (`discussion_plan.md`)

- **Needs:** G-2 and G-4, and the Results drafted.
- **Done when:** it delivers a **decision rule** (how many channels, what interval relative to the relaxation, how much noise, and what the cheap approximations cost you in that corner). Not "MacroIR is best". The decision rule is what gets cited.
- **Must contain:** the Kalman concession (MacroIR is equivalent to an integrated-measurement augmented Kalman filter, verified to ~1e-8: concede it loudly, the novelty is elsewhere), the identifiability perimeter (`docs/bibliography/identifiability/`), and IR's own miscalibration at few channels (it over-predicts the observable variance by about 6%: the hero's own failure belongs in the paper that is about failure).

### 6th. Introduction (`introduction_plan.md`)

Written late on purpose. The Introduction sells what the paper actually delivers, and you only know that once the Results exist.

- **Needs:** G-2, and the Results drafted.
- **Done when:** the reader arrives at "we provide that test" already convinced that no such test existed for this problem and that the field is choosing likelihoods by habit.

### 7th. Abstract (`abstract_draft.md`)

A 209-word abstract already exists in `elife_paper.tex`. It is not a placeholder; it is in the validity-map frame and it reads well. Treat it as a draft to be **reconciled** with the finished paper, not rewritten from zero.

- **Needs:** G-3, and everything else drafted.
- **Also produces:** the Impact Statement (eLife requires one, 40 words or fewer).

## 4. The two holes

**Results §5 (the regime maps) and §6 (the distortion decomposition).** No figure is selected. Forty-four candidates sit in `projects/eLife_2025/figures/in_progress/`, already sorted into main / supplementary / drop by `figure_5_PLAN.md`, but none is promoted, none has a caption, and four unselected PDFs already occupy the numbers `Figure_5_*` and `Figure_6_*` in `figures/paper/`, which will cause confusion if not cleaned.

Two things gate them, and the good news is that neither is expensive:

- **The anchor.** The paper promises Gaussian-Fisher-anchored figures, and today no cross-algorithm panel can be: `1c2ae6f` is missing **NMR entirely**. The fill is exactly **four cells** (NMR × N_ch {10, 100, 1000, 10000} at noise 0.1), about 1,830 CPU-hours, one to two days of queue on a cluster with no CPU-hour quota. This decision is cheap and should be taken, not deferred. See `09_carve_plan.md` (Freeze preconditions) first: landing the seed fix changes the commit hash, so the new cells write to a new data directory.
- **The schema port, which nobody has costed.** Every production notebook hardcodes the numerical run and its column names (`Likelihood_Information_Distortion`, `Likelihood_Sample_Distortion`, ...). The Gaussian run emits *different* names (`Gaussian_Distortion_Induced_Bias`, `Likelihood_Gaussian_Information_Distortion`, ...) **and still carries the old columns filled with NaN**. So a notebook repointed at the Gaussian directory does not fail. It renders grey. This is real work between the fill and Figure 5, and it is the thing most likely to be mistaken for a re-knit.

Write §5 and §6 as prose around a `\todo{}` claim. Do not guess the numbers.

## 5. The one decision that adds work

Every other open decision **reorders** work. This one **adds** it:

**Research Article or Tools and Resources?** *Tools and Resources* obliges benchmarking against existing methods, which is new science, and public deposition of code and major datasets. It is the only branch in the plan that could add a month. `diagnostics_plan.md` argues the T&R case; `abstract_draft.md` defers it "to submission". Deferring it is the one deferral that costs something.

Everything else (the title, the MNR/NMR spelling, the validity threshold, the figure count, the Comm Biol disclosure) is consumed by prose that is not written yet, and can be settled at any point before the sections that use it. Sequence them freely.

## 6. What is deliberately not here

- **Engine work** (the seed, the IDM call site, the row duplication, the regression test, the valgrind read): `09_carve_plan.md`, Freeze preconditions. It gates the *build*, not the prose.
- **Compute** (which cells to run, what they cost): `09_carve_plan.md` and `05_experiment_grid.md`.
- **Submission paperwork** (bibliography, Data Availability, CRediT, competing interests, funding, MDAR, Key Resources Table, the preprint, the publication fee): unowned today. It needs a `10_submission_pack.md`, and the bibliography is the item to start now, because the paper has **no citation database at all** (`biblio.bib` is the P2X2 paper's 170 references and the `.tex` contains zero `\cite`), it has no blockers, and it takes longer than anyone expects.
- **What each section says:** the section plans. This document must never restate them.
