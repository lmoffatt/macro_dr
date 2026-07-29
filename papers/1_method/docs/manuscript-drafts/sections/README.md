# Rules for every manuscript section

> **Updated 2026-07-29.** One brief per section lives beside its `.tex` in this directory
> (`00_abstract.md` … `07_backmatter.md`). This file holds only what governs **all** of them, so a
> section brief links here rather than repeating it. Together they replace `SPINE.md`, which
> replaced the eight standalone section plans; the originals are in
> `../../../archive/section-plans-20260729/`.
>
> **The rule that keeps these files alive: a section brief holds no prose and no numbers.** The
> plans that rotted did so because they carried drafts the manuscript then overtook. Policy does not
> rot the same way.
>
> | Section | prose | brief |
> |---|---|---|
> | Abstract | `00_abstract.tex` | `00_abstract.md` |
> | Introduction | `01_introduction.tex` | `01_introduction.md` |
> | Theory | `02_theory.tex` | `02_theory.md` |
> | Diagnostics | `03_diagnostics.tex` | `03_diagnostics.md` |
> | Results | `04_results.tex` | `04_results.md` |
> | Discussion | `05_discussion.tex` | `05_discussion.md` |
> | Methods | `06_methods.tex` | `06_methods.md` |
> | Back matter | `07_backmatter.tex` | `07_backmatter.md` |

Accuracy constraints, not taste. They came out of a four-review plus nine-judge convergence on the
title and they bind the whole manuscript.

1. **Distortion, not loss.** The effect is bidirectional: some approximations over-state the
   information the data carry, others under-state it, and `IR` itself runs to both sides in the
   few-channel corner. Any phrasing built on "information loss" or "preserving information" is
   factually wrong.
2. **The approximation distorts, not the time averaging.** Time averaging is the physical reality that
   the boundary-conditioned likelihood handles correctly. "Time averaging degrades the information"
   inverts the causation and contradicts the thesis.
3. **Continuous, not a phase transition.** Avoid "breaks down at", "fails beyond". **[SETTLED SINCE]**
   This once forbade the region map outright; it does not. The reconciliation is the map's own
   sentence, and it should be quoted rather than re-derived: *the map is a concept map, not a phase
   diagram; its boundaries are level sets of continuous diagnostics, so a looser criterion moves each
   of them by up to a decade in noise without changing the layout.*
4. **No "only", no universal claims.** The study is two-state. "Only X stays calibrated" invites a
   scope objection that the data cannot answer.
5. **No method promotion.** The paper characterises; it does not re-present MacroIR. That bridge is
   published (Comm Biol 2025). Every sentence beginning "MacroIR is…" is a candidate for deletion.
6. **Do not claim the test.** The diagnostic is classical: Huber 1967, White 1982, and the generalized
   form of Golden, Henley, White & Kashner. "There was no way to test whether a likelihood is faithful"
   is false and a statistically literate referee will say so. What is new is the *measurement* of it
   for this class of likelihood, enabled by a process that can be simulated exactly.
7. **The gap is validity, never absence.** The temporal correlation of macroscopic currents has carried
   kinetics since 1973. Do not write that it was unused. See `../../../decisions.md` and
   `docs/bibliography/temporal_correlation_and_AR_errors_2026-07-28.md`.
8. **Register.** Plain scientific prose, the register of the 2007 and 2025 papers. No aphorisms, no
   rule-of-three triads, no em-dashes.

**Two readers, and every section has to hold both.** The *electrophysiologist* who fits schemes to
macroscopic currents, has never computed a score, and must not need to know what a sandwich estimator
is; and the *inference-methods reader* (the Mirams / Münch / Del Core axis), for whom the words that
buy attention are misspecification, score, and information matrix equality. One sentence each, no more,
and gloss each in physical terms in the same sentence.

---

## Title

**Chosen, for now:** *Information distortion in likelihood approximations for macroscopic ion-channel
currents.* "for", not "of": what is approximated is the likelihood, not the currents.

**Currently in the manuscript** (`elife_paper.tex:32`): *A validity criterion and a usage map for
likelihoods of macroscopic ion-channel currents.* **[OPEN]** The two disagree and the `.tex` one is
post-merge, so it is probably the survivor; decide and record which, because rule 5 above and the
abstract both key off the title.

eLife rule that binds: if the biological system is not in the title, it must be in the abstract. It is
currently in both.

---

## Standing blockers, which belong to no single section

1. **The direction convention.** `../../../../_program/machinery.md` §4 asserts it three times with two different
   signs and it has never been checked against the producer. Every verdict and every boundary of the
   usage map inherits it, and the manuscript carries eight `% TODO-SIGN` markers waiting on it. It is
   about an hour of reading code and it should be done before any result sentence is finalised.
2. **`Figure_2.pdf` is not rendered in `paper_both`.** The notebook and the caption are there and the
   caption was updated 2026-07-23, but the PDF exists only in `paper_1` from before the merge. It is the
   only one of the manuscript's 23 referenced graphics that is missing, and the `.tex` already tracks it
   as `PENDING-FIG2`. One R render.
3. **LSE at n_sims 1000** against the likelihood arm at 10⁴. Mixing them in a panel manufactures a
   regime effect through the Jensen bias. Re-run, or restrict every LSE panel to a uniform n_sims and
   say so in the caption.
4. **Telegraph noise.** The Discussion asserts MacroIR fails there. Theory names the regime, but no
   simulator capability and no run exists on any freeze commit. Either measure it or demote it to an
   argument with a citation, and say which.
5. **`check.sh` tests `N_CAP >= 5`.** The figure set is six. Until that is raised the done-oracle stays
   green on a manuscript short one body figure.
