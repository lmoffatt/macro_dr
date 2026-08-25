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

**SETTLED 2026-08-04**, live in `elife_paper.tex` and nowhere else:

> *Likelihood approximations distort the ion channel kinetic information in macroscopic currents.*

The declarative finding, which is eLife's dominant shape and the grammar of Luciano's most-cited
paper (evidence in `../../../archive/section-plans-20260729/title.md`), reordered to the grammar of
`moffatt2007estimation` so that "ion channel" is open and attached to the kinetic information while
"macroscopic currents" stays unmodified. **No hyphen in "ion channel" anywhere in the manuscript.**

Retired in the same pass, with reasons, so neither comes back. *A validity criterion and a usage map
for likelihoods of macroscopic ion-channel currents*, which sat in the `.tex` and which nobody had
chosen: "usage map" promises the recommendations the Discussion's demarcation explicitly declines.
"Map" is not the problem and stays available in the body; "usage" is. And *Information distortion in
likelihood approximations for macroscopic ion-channel currents*, a noun phrase, which is eLife's
secondary shape. Gerund forms ("Measuring how...", "Quantifying how...") were considered and dropped:
tertiary shape, and they read as a procedure rather than a finding. The novelty scoping a gerund
would buy is already done, and done better, by the narrow-novelty paragraph in the Introduction.

eLife rule that binds: if the biological system is not in the title, it must be in the abstract. It is
currently in both.

---

## Standing blockers, which belong to no single section

1. ~~**The direction convention.**~~ **CLOSED 2026-08-04.** Read against the producer and confirmed as
   already written: `Lapack_PSD_Normalized_Congruence_Matrix` (`legacy/lapack_headers.h:2519`) computes
   H^(−1/2) J H^(−1/2) with H first, and `legacy/distributions.h:470-478` names the reported Fisher
   information as that H. So C > 1 is over-confidence, which is what `machinery.md` §4,
   `03_diagnostics.tex` and the theory supplement all already said. The warning that it was asserted
   with two different signs was stale; no dissenting copy exists. The six `% TODO-SIGN` markers are
   cleared and the verification is recorded at the top of `03_diagnostics.tex`. The seventh, in
   `00_abstract.tex`, is prose recording that the marker was retired from that section and stays.
   **Canonical phrasing, from rule 1 above:** an approximation *over-states* or *under-states the
   information the data carry*. Interval-width wording is a permitted gloss in the same sentence, never
   a substitute, so that the two directions read the same way everywhere.
2. ~~**`Figure_2.pdf` is not rendered in `paper_both`.**~~ **CLOSED, stale since 2026-08-03.** It is
   there: `figures/paper_both/Figure_2.pdf`, 1.08 MB, rendered 2026-08-03 13:21, which is newer than
   the `paper_1` copy of 2026-07-22. Drop any `PENDING-FIG2` marker that survives in the `.tex`.

   **The "real graphics blocker" that stood here is CLOSED, stale, verified 2026-08-25.** The
   supplement set was renumbered per parent figure (`Figure_2_S1`, `Figure_3_S1..S3`,
   `Figure_4_S1..S4`, `Figure_5_S1..S3`) and every one of the eighteen graphics the `.tex`
   references exists in `paper_both`, rendered between 2026-08-05 and 2026-08-12 (checked file by
   file against the `\includegraphics` set). The pre-merge `Figure_4_supplement_N` naming this
   entry described no longer appears in any section file. If a `% src:` trail to an old
   `figure_4_supplement_N.html` resurfaces, that single trail is what needs repair, not this set.
3. ~~**LSE at n_sims 1000** against the likelihood arm at 10⁴.~~ **CLOSED 2026-08-04, the blocker was
   stale.** Checked on disk: `figures/data/0ffbda7/` holds 99 least-squares and `nonlinearsqr` files
   and every one of them is `nsim_10000`. The only `nsim_1000` least-squares cells are eight files in
   `82b956f/`, which is not one of the three freeze directories the manuscript quotes from
   (`1c2ae6f + 87889e6 + 0ffbda7`), so no quoted number comes from them. The arm is uniform at 10⁴ and
   Luciano is extending it at the same n_sims. The magnitude argument closes it independently: the
   measured anisotropy of `ILSE` is 3.17 to 3.22 across four channel counts, while the low-n_sims
   inflation of that statistic is about 1.15 against an estimation floor of 1.043, so the artifact is
   negligible next to the signal even where it applies.
4. ~~**Telegraph noise.**~~ **CLOSED as stale, 2026-08-25:** the word and the assertion no longer
   appear anywhere in `05_discussion.tex` (grep over the non-comment text). The claim was removed
   rather than measured; if it is ever wanted back, it re-enters as an argument with a citation,
   never as an assertion, and this entry does not reopen.
5. ~~**`check.sh` tests `N_CAP >= 5`.**~~ **CLOSED 2026-08-25.** The script already tested `>= 6`
   when this entry still said 5; the body set is seven (`Figure_1`..`Figure_7`) and the threshold
   was raised to `>= 7` in the same pass, so the done-oracle no longer stays green one figure short.
6. **The 2026-08-25 coherence pass, so its load-bearing points are findable from here.** The
   abstract was regenerated (note 1n in `00_abstract.tex` lists every change and its reason); the
   Discussion carries the P2X2 verdict sentence and the regime clause; the Introduction's Comm Biol
   paragraph was rewritten against a full read of the 2025 paper; Results' 93/70/31 sentence moved
   to the shared sub-grid (97/78/37) and its cell-set defect closed; item (f) closed with the
   fisher_only run (`figures/data/a202e03`), its table in `supplementary_file_1.tex` S1.5.
   **Still open after the pass:** the h-step robustness cell for NR's differenced-Fisher number
   (`jobs_to_be_run.md`, Job 2, aviso 5 — the number is now printed, so the cell is owed), and the
   abstract's two-state objection (note 1n's closing line).
