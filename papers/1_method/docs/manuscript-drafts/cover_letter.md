# Cover letter and submission-form fields

Canonical text of the eLife cover letter. v6, 2026-08-30, rewritten after the outside review of
2026-08-30 and then revised on the findings of four independent readers (fact-check against the
.tex sources, referee-style refutation, style lint, editor triage). Finding first, then the stake,
the method, the fit to eLife, and availability. v5 (official 2026-08-28) is archived at
`archives/cover_letter_v5_SUPERSEDED_20260830.md`; v1 and the drafts v2 to v4 are beside it. The
PDF is `cover_letter.pdf`, compiled from `cover_letter.tex`, which mirrors this file and must be
edited with it; the PDF was NOT recompiled on 2026-08-30.

Addressed to "Dear Editors": eLife assigns the Senior Editor after submission, and any suggestion
of editors goes in the submission form's own field, not in the letter.

Voice: "I" for the submitting author and for the measurement; "we" and "our" only for the 2025
paper, which has a coauthor.

One thing to settle before submission: the bracketed related-manuscripts sentence, which eLife
asks for and which only Luciano can confirm (the Introduction says the evidence correction "is derived elsewhere" and Results name "a
companion paper"; if either is under consideration anywhere, say so instead).

Dropped or changed on the readers' findings (two passes), each with the reason: "and was later
confirmed" (the manuscript no longer says it; the Jiang 2012 sentence left the Discussion on
2026-08-27, so restoring it there is the only way to put it back here); "rested on the evidence
over nine schemes" (the 2025 abstract adds atomistic simulations); "had never been measured"
without "how much and in which directions" and "for macroscopic currents" (Münch 2022 did count
coverage; the single-channel literature exists); "are now chosen" (the Introduction says current
practice is deterministic fits; now "beginning to be chosen ... as we did"); "eight likelihoods in
use" (MR and VR are our own partial corrections, tried and rejected); "the filter is not the
novelty ... known since 1988" (contradicted Appendix 5; now: form is the 1988 filter, what is new
is the score and Fisher information and the measurement); "the most recent calibration" for Münch
(the letter's own first paragraph says it had never been measured; now "the closest prior check");
"stated as a conjecture" for the nine-scheme transfer (the Discussion's hierarchy files it under
suggested, not shown); "the same test runs unchanged on any scheme from the deposited library" (no
scheme library is deposited, and generality is suggested, not shown); "the exact simulator" moved
from the library to the engine deposit (the backmatter lists it for the engine only; if the
bindings do expose simulate(), add it to the backmatter first); "under either criterion" attached
to the factor of a hundred (false: the loose gap is a factor of about three); "arm", "recursive
likelihood", "decades", "unitary current", "one pass", "evidences", "the reference condition",
"accumulates" (glossed or replaced for an editor outside the field); "I know the stake from the
inside" (carried no fact).

Numbers, each with its source: "about twenty times ... at one representative condition" = INR's
accumulated ratio at the Figure 3 cell (Results, 20.8; N_ch = 100, S~ = 0.01, Delta~ = 0.1, at the
truth, 1,000 recordings); "560 ... 210" = counted 2026-08-30 over the battery_pool_G files of all
data directories at n_sim = 10^4: IR 80 (N_ch, noise) combinations x 7 intervals = 560, and
exactly 30 combinations shared by all eight members = 210 (S1.1 was stale on this and was
corrected the same day); "an order of magnitude, two- to threefold on the error bar" = the abstract's (re-synced 2026-08-30 with the approved 188-word abstract; "up to" dropped, it turned a typical value into a ceiling)
wording (maxima are 57, 33 and 74 in Results; if the abstract changes, change this too); "a factor
of a hundred ... a factor of three" = 1.7 to 2.0 decades strict, 0.4 to 0.6 loose (Discussion,
Figure 7 caption, measured 2026-08-30); "1988" = Zadrozny 1988, form only, per Appendix 5; the
three DOIs = biblio.bib (moffatt2026code, moffatt2026data, moffatt2026macroirlib), to be checked
as public before sending. "Whether its subunits couple symmetrically or in sequence" is not in the
letter any more (the 2025 abstract says "sequential, asymmetric coupling mechanism" against
"symmetric, concerted gating", verified in the PDF, if the phrase is wanted back).


---

Dear Editors,

I am submitting "Likelihood approximations distort the ion channel kinetic information in
macroscopic currents" for consideration as a Research Article.

Kinetic mechanisms of ion channels are beginning to be chosen by comparing the Bayesian evidence
for rival schemes fitted to macroscopic currents, as we did for nine P2X2 schemes in 2025. How
much the likelihoods behind such a comparison misreport the information a recording holds, and in
which directions, has never been measured for macroscopic currents. This manuscript measures both
for a ladder of eight likelihoods, from least squares to a filter that uses the channel state at
both ends of each sampling interval, against ten thousand exact simulations at every combination
of channel number, noise and sampling interval tested.
Every likelihood that ignores the correlation between successive samples misreports its
information by an order of magnitude, two- to threefold on the error bar, wherever that correlation
is present in the record, and no single rescaling repairs it, because the distortion differs
between parameter directions. The filter conditioning on both ends reports the uncertainty it
delivers, which I call calibrated, except where single openings become resolvable in the record,
and there its departure follows one measured power law. No measured condition gives both a
calibrated least-squares error bar and information about the single-channel current; the two lie
about a factor of a hundred apart in noise, and still a factor of three apart under the paper's
looser criterion.

In our 2025 ranking of nine P2X2 schemes the winning mechanism changed when the likelihood's
recursion was removed, overturning a factor of six in evidence, and we published that
sensitivity. In the two-state scheme this manuscript shows that the likelihood without the
recursion, our control, reports at one representative condition about twenty times the
information the record holds, while the likelihood we argued for is calibrated there; whether
that carries to nine schemes is suggested and not shown, and the evidence is not recomputed here.

The measurement needs no fit to a real recording. At known parameters a correct likelihood
satisfies two identities: its score, the gradient of the log-likelihood, averages to zero, and the
covariance of the score equals the Fisher information the likelihood itself reports. The
calibrated filter was tested over 560 such combinations, all eight likelihoods over 210 of them.
The scheme has two states at an open probability of one half, the smallest that isolates the
approximation's own error; the failure comes from a correlation the likelihood does not carry,
which richer schemes share but are not tested here.

The closest prior check on a filter of this kind appeared in eLife, by Münch and colleagues in
2022. They counted how often the true value fell inside the reported interval over repeated fits,
a check that costs a fit per simulation and so covers a line of conditions. The two identities
cost one likelihood evaluation per simulation and return a matrix over a plane: which parameter
is misreported, and whether at each sample or across the record. The filter itself was published
with the P2X2 analysis and has the form of a 1988 integrated-measurement Kalman filter; what is
new is its score and Fisher information, and the measurement built on them.

The engine with its exact simulator, the configurations and
the notebooks are at 10.5281/zenodo.22168409, the simulation output at 10.5281/zenodo.22167744,
and the likelihood with its score, Fisher information and diagnostic as a small library with R
and Python bindings at 10.5281/zenodo.22168263, so a reader can apply the same two identities to
their own scheme. [No related manuscript is under
consideration elsewhere.]

Yours sincerely,

Luciano Moffatt

---

## Impact statement (submission form, 15 to 30 words, third person)

A likelihood for macroscopic ion-channel currents whose reported information is measured against
the exact process, with the apparatus to repeat the measurement on any scheme, available from R
and Python.

(30 words, Luciano's decision of 2026-08-25.)

Finding-forward alternative, trimmed to 30 words on 2026-08-28 (the v1 note called its longer form
29 words; it was 36):

In a published comparison of nine gating schemes the ranking moved with the likelihood; measuring
score and information against exact simulation shows which reported uncertainties can be believed,
and where.

---

## Other form fields (2026-08-30)

**Article type.** Research Article, decided 2026-08-30 against the review's Tools and Resources
suggestion: same editors and the same assessment vocabulary for every type, and the seven eLife
papers in the bibliography, Münch 2022 included, are all Research Articles; T&R would move the
evaluation onto the artefact, the least mature part (K = 2, no experimental data, bindings checked
on one cell).

**Subject areas.** Structural Biology and Molecular Biophysics; second area Physics of Living
Systems or Computational and Systems Biology (Luciano to choose).

**Reviewing Editor suggestion.** Marcel P. Goldschen-Ohm (University of Texas at Austin) was the
Reviewing Editor of Münch et al. 2022 (first page of the PDF). Experience with the competitor's
paper, not a conflict; Luciano's decision.

**Funding.** No grant funded this work. Funder: CONICET (the author is a member of its Carrera del
Investigador Científico); grant number: none. Matches the manuscript's Funding section.

**Preprint.** Not in the letter: eLife posts the preprint itself (Luciano, 2026-08-30), so no DOI
is quoted and nothing is posted beforehand.

**Related manuscripts.** eLife asks. The bracketed sentence in the letter is a placeholder for
Luciano's answer.

**Suggested and excluded reviewers.** Form fields, not the letter. To decide.

---

## Notes that still apply

**Why the atomistic simulations are not in the letter.** v1 named and bounded them; v5 and v6 do not
mention them, so the 2025 abstract's "supported by atomistic simulations" is not contradicted, only
left out. The reasoning that governed v1 is kept for the referee stage: The first sentence of the 2025 abstract
says the kinetics are "supported by atomistic simulations", so a letter describing that claim as
resting on the recordings alone is contradicted by a document the editor can open in one click, and
it drops the coauthor's contribution on the way. They are named for that reason and scoped for
another: they are on a zebrafish P2X4 homologue, they corroborate the asymmetry and fix its
direction, which the 2025 introduction states electrophysiology alone cannot resolve, and no
atomistic result orders kinetic schemes. What the evidence ratios produced on their own is the
choice among the nine, and that is the quantity this manuscript's diagnostic bears on. The sentence
about the record is therefore scoped to the record and says nothing about the conclusion.


**The answer to have written before submitting, and probably never to use.** A referee who opens the
cited 2025 paper reaches its supplement and finds two things: Table S4, whose validation against
simulated data leaves the unitary current and the channel number outside their 90% intervals with
the product preserved, and Table S1, which shows two of the non-recursive runs at effective sample
sizes of 3 and 4 with R-hat 1.21 and 1.16 against a declared standard of 1800 and 1.01. The question
is why either case should be believed. The answer, in order: the manuscript's claim is that the
ordering moved between two likelihood configurations differing in the recursion, and poor
convergence in two of the non-recursive runs weakens that ordering further without rescuing the
other, which is exactly why this paper declines to say which to believe; the displacement of the two
amplitude parameters in that validation table is a property of that run, its cause is identified,
and it is being handled at the journal that published it; the member measured here is the one
specified in Theory and its calibration is what this paper reports, point by point, against exact
simulation. On the engine itself, `sections/05_discussion.tex` already carries the bound: the
diagnostic compares two independent implementations, the likelihood and the exact simulator,
everything downstream of where they separate is tested by their disagreement, and the only thing
they share upstream is the two-state specification, which is written out in Methods to be read by
eye.
