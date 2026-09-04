# Cover letter and submission-form fields

Canonical text of the eLife cover letter. Made official 2026-08-28 (v5); the PDF is
`cover_letter.pdf`, compiled from `cover_letter.tex`, which mirrors this file and must be edited
with it. Earlier versions and the notes on the restrictions withdrawn on 2026-08-10 are in
`archives/cover_letter_v1_SUPERSEDED_20260828.md`; drafts v2 to v4 of 2026-08-28 are beside it.

Addressed to "Dear Editors": eLife assigns the Senior Editor after submission, and any suggestion
of editors goes in the submission form's own field, not in the letter.

"How subunits couple and whether they do so symmetrically" is a gloss of the 2025 question and is
the one phrase that has not been checked against that paper's own wording.

---

Dear Professor [Senior Editor],

I am submitting "Likelihood approximations distort the ion channel kinetic information in
macroscopic currents" for consideration as a Research Article.

We have reached the limit of what a macroscopic current shows by eye. The flip state we proposed
for P2X2 in 2007 was still on the visible side, a delay before the current rises, later confirmed.
The questions now asked of these receptors, how subunits couple and whether they do so
symmetrically, leave no such mark. They are settled by ratios of evidence between kinetic schemes,
and whether the likelihoods behind those ratios report the information a recording actually holds
had never been measured. This manuscript measures it, and finds that most of them do not.

We know the stake from the inside. Our 2025 report that P2X2 activates asymmetrically rested on
the evidence over nine schemes, and under a control likelihood differing only in its recursion a
different mechanism came first and a factor of six in evidence disappeared. We published the
sensitivity and argued for one arm. On a recording there is no truth to measure a likelihood
against, and fitting well does not supply one: a recursive likelihood follows the record whether
its assumptions hold or not.

So we went outside the fit. Channel gating simulates exactly where its likelihood cannot be
evaluated exactly, and at known parameters a correct likelihood must satisfy two identities: its
score averages to zero, and the covariance of that score equals the Fisher information it reports.
We measured both for eight likelihoods, from least squares to a filter conditioning each interval
average on both endpoints, over 210 combinations of channel number, noise and sampling interval,
ten thousand simulated recordings at each, in a two-state scheme so that the error measured is the
approximation's own.

Every likelihood that leaves the gating correlation unmodelled misreports its information, by up
to an order of magnitude, and no rescaling repairs it because the distortion has directions. Our
published control is one of them: across a record it accumulates about twenty times the
information it holds, and the recursion is what brings that factor back to one. We do not
recompute the P2X2 evidences. We show that the arm favouring the alternative mechanism is
distorted in the simplest model there is and the calibrated arm is not, and we audit that arm on the same grid. It departs where the record begins to resolve single-channel openings, the frontier with the microscopic regime, whichever of channel number, noise or sampling interval carries it there, and we say by how much.

Münch and colleagues asked this of their own filter in eLife in 2022 by counting coverage, which
tests a posterior and costs a fit per replicate. The identities cost one pass and return a matrix:
which parameter is misreported, and whether the error is committed at each sample or accumulated
across the record. We are not offering another filter; ours agrees to one part in 10^8 with an
integrated-measurement Kalman filter known since 1988. The measurement is what is new.

For an experimentalist the map has one line that holds however the thresholds are drawn: no
region gives both a trustworthy least-squares error bar and information about the unitary current,
the two boundaries one to three decades of noise apart. The recordings rich enough to determine
the amplitudes are the ones whose classical error bars are the wrong width. The calibrated
likelihood, its score and Fisher information, the exact simulator and the diagnostic are a small
library with R and Python bindings, so a reader can put the same two identities on their own
scheme.

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

## Notes that still apply

**Why the atomistic simulations are not in the letter.** v1 named and bounded them; v5 does not
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
