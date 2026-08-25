# Cover letter and submission-form fields

Drafted 2026-08-07 for the eLife submission. Not part of the compiled PDF.

The axis wording here matches the manuscript and is accurate: the two runs behind Figures 1d and 1e
of Moffatt & Pierdominici-Sottile (2025) shared the interval averaging and the gating variance and
differed in the recursion alone. That stays.

The naming restriction that stood here is **withdrawn, 2026-08-10**. It read: "Never name them as
MacroIR and MacroINR, and never call either the boundary-conditioned member: the object that ran
there is not the member this manuscript specifies in Theory." The published names may be used; see
`sections/00_abstract.tex` note 1f(i) ("NAMING, settled 2026-08-10"), and the bridges
`IR = MacroIR` and `INR = MacroINR` are closed at `papers/_program/nomenclature.md:286-296`.

Why it was wrong, so it is not rebuilt: two objects carried one name. The `macro_NMR` in this
repository's freeze omits the N·ms interval-variance term (refactor regression at a3e0a89) and is
genuinely not the member Theory specifies; the **published** MacroINR is a different object and does
carry the term, confirmed 2026-07-31 in that same nomenclature section and again 2026-08-10 against
the run's own `_likelihood_model.csv`. The restriction generalised a fact about a local build to a
peer-reviewed algorithm, on a premise that had already been refuted in writing eight days earlier.
Whether this letter should use the names is now a free choice; the prohibition is what is dead.

One claim below is deliberately weaker than it could be, for the same reason it is weaker in the
manuscript: nothing here says which of the two orderings is right, which this paper does not
measure.

The ranking restriction that stood here is **withdrawn, 2026-08-10**. It required "the two rankings
were not the same" rather than "put different mechanisms first", on the ground that the paper
numbers nine schemes I–IX in its main text and eleven I–XI in Table S1 with no published mapping.
That premise is true and irrelevant: Figure 1 is in the main text, uses main-text numbering, and
draws the schemes in panels a–c, so d and e give both rankings without touching the SI. IX first
under the recursive likelihood, VI first under the non-recursive one, and the main text supplies the
direction ("the control method systematically underestimated evidence for schemes with
conformational intermediates") and the margin it erases (a Bayes factor of 6.4 for IX over VI). See
`sections/00_abstract.tex` note 1f(i).

---

## Cover letter

Dear Professor [Senior Editor],

I am submitting "Likelihood approximations distort the ion channel kinetic information in
macroscopic currents" for consideration as a Research Article.

A method that cannot resolve a mechanism cannot mislead you about one. The likelihoods now applied
to macroscopic ion-channel currents can resolve mechanisms that leave no visible trace in the
record, and they return posteriors tens of times narrower than the priors that went in. That
resolution is what makes their calibration a substantive question. Parameter intervals somewhat too
narrow are a nuisance. A mechanism ranked wrongly is a conclusion.

We have been on both sides of that line. In 2007 we proposed that P2X2 passes through a flip state
before opening. That claim had a signature you could point to in the record, a delay before current
appears, and other groups later confirmed it. Last year we reported that the same receptor
activates asymmetrically (Communications Biology, 2025). Nothing in those recordings shows asymmetry
the way the delay showed flip. Atomistic simulations of a zebrafish P2X4 homologue corroborate it
and supply a direction the electrophysiology cannot resolve, and that support comes from outside the
record. What it does not settle is which of nine kinetic schemes the recordings prefer, and that
came from ratios of evidence alone. When we ranked the nine twice, under two likelihoods that shared
the interval averaging and the gating variance and differed in whether each interval is conditioned
on the data already observed, the two orderings put a different mechanism first: the conformational
asymmetric scheme under one, the subunit-specific allosteric alternative under the other. We
reported the sensitivity, and the direction we gave for it rested on argument: on recordings there
is no ground truth to measure a likelihood against.

Nothing in ordinary practice adjudicates it. A recursive likelihood follows the record whether its
assumptions hold or not, so its residual cannot be read the way a least-squares residual can, and
fitting well says nothing about whether the information a likelihood reports is the information it
has. The error we were afraid of is the kind that passes every check we know how to run, which is
why we went outside the fit for one.

The manuscript measures both classical identities, the score mean and the information equality,
along a ladder of macroscopic likelihoods from least squares to a filter conditioned on both ends of
each acquisition interval, against exact simulations of the process they approximate. The measured
object is a matrix, and it is the same matrix any correction to an evidence ratio has to start from,
through a volume term and an effective-sample rescaling that the Discussion sets out. That is why it
had to be measured first. We do not recompute the P2X2 evidences here, and the Discussion says so.

eLife published Münch and colleagues on this problem class in 2022, and they asked our question of
their own filter, counting how often the true rate matrix falls inside a credibility volume of given
mass. What their count tests is a posterior, so the prior and the sampler are in the verdict with the
likelihood, and it costs a fit per replicate at every design point. The identities we use are
evaluated on the likelihood alone at parameters that are known, they cost one pass, and what comes
back is a matrix that says by how much and in which parameter the reported uncertainty departs, and
whether the departure is committed at each sample or accumulated across them. That is what puts a
plane of 560 design cells within reach. We do not offer another filter. Our best member coincides
with an integrated-measurement Kalman filter known since 1988 to about one part in 10^8, and we say
so in the Discussion. What is new is the measurement.

The measurement does adjudicate the two instruments of the P2X2 comparison, and that is the one
statement about our own record the manuscript now makes: the control that reordered the schemes is
the non-recursive member of the ladder, measured to misreport its information about twentyfold at
macroscopic channel counts, which is exactly the factor the recursion restores. For the job of
reporting uncertainty, a likelihood whose identities have never been measured supplies estimates
whose stated precision cannot be acted on; in that operational sense the field did not yet have a
calibrated instrument for macroscopic currents. It has one now, with its own measured limits at few
channels and low noise, and it has the apparatus that certified it, which runs on any scheme an
exact simulator can reach.

The study is simulation on a two-state scheme sampled through a uniform acquisition window, and no
experimental recordings are analysed. The scheme is deliberately small so that what the diagnostic
reports is the approximation's own error, with no misspecified mechanism mixed into it, and the
manuscript states those limits beside the results. The software carries the calibrated member with
its score and its Fisher information together with the exact simulator, and the measuring apparatus
itself: the distortion matrix, its eigenvalues and the first-order bias come back from one call in R
or Python, so a reader can run the same two identities on their own scheme rather than extrapolate
from ours.

One reading of the map does not depend on where a level set is drawn. No region of the measured
plane gives both information about the unitary current and a trustworthy least-squares error bar,
and the two boundaries sit one to three decades of instrumental noise apart everywhere they were
measured. The recordings rich enough to determine the amplitudes are the ones whose classical
intervals are the wrong width.

Yours sincerely,

Luciano Moffatt

---

## Impact statement (submission form, 15–30 words, third person)

A likelihood for macroscopic ion-channel currents whose reported information is measured against
the exact process, with the apparatus to repeat the measurement on any scheme, available from R
and Python.

(30 words. Replaced 2026-08-25, Luciano's decision: the tool-forward form, the home found for the
deliverable trio the abstract panel evicted from the abstract's closer. The previous
finding-forward version, kept in case the venue prefers a finding: "In a published comparison of
nine gating schemes for a ligand-gated receptor the ranking moved with the likelihood; measuring
each likelihood's score and information against exact simulation shows which reported
uncertainties can be believed, and where." 29 words.)

---

## Notes on choices made here

**Why the scope paragraph sits before the result and not at the end.** The last thing the editor
reads should be the finding. Putting "two-state scheme, no experimental data" last leaves the
submission on its own weakest sentence.

**Why the personal register appears twice and no more.** "We have been on both sides of that line"
and "the error we were afraid of" are the only two, and they carry the motive that makes
self-validation something other than circular. A third would turn the letter into a narrative.

**What is deliberately absent.** No mention of Supplementary Table S4 of the 2025 paper, and no
claim about which of the two orderings is correct. The first belongs to the record where the record
lives; raising it here asks the editor to decide whether this submission is an article or a
correction. The second is not measured anywhere in this manuscript.

**Why the atomistic simulations are named and then bounded.** The first sentence of the 2025 abstract
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
