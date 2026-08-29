# Cover letter, draft v2

Drafted 2026-08-28 from the v1 letter (cover_letter.md, untouched), the dictated version of
2026-08-28 15:46, and the current sections. Not part of the compiled PDF.

What changed against v1, so the corrections can be aimed:
1. 560 design cells is gone. The live plane is 210 cells per member (04_results.tex:127,
   12_appendix_priorart.tex:114); 560 was retired 2026-08-27 and the results note forbids quoting it.
2. The twentyfold is anchored. It is the information accumulated across a record, which is the
   supplement's reading, and the letter now says so, since the abstract's "order of magnitude" is
   the body figure's and Results forbids matching the two without the anchor.
3. Third beat added. v1 said "the measurement does adjudicate" and stopped. It now delivers the
   verdict at the strength the Discussion already carries: indirect support for the published
   reading, no evidences recomputed.
4. Grafted from the audio: the lens, the code verification, the closing on sharing.
5. Two states now carries both halves, the constraint and the design reason.
6. "The field did not yet have a calibrated instrument" is out. The non-stationary fluctuation
   analysis concession the Discussion makes is in.
7. The margin of six is stated.
8. The scope paragraph moved up so the letter has one ending.
9. DROPPED, needs checking before it goes back: "posteriors tens of times narrower than the priors
   that went in". It appears nowhere in the manuscript and nowhere else in the repo. If the 2025
   paper supports it, it is a good sentence and belongs in the first paragraph.

---

## Cover letter

Dear Professor [Senior Editor],

I am submitting "Likelihood approximations distort the ion channel kinetic information in
macroscopic currents" for consideration as a Research Article.

A likelihood that cannot resolve a mechanism cannot mislead you about one. The likelihoods now
applied to macroscopic ion-channel currents do resolve mechanisms that leave no visible trace in the
record, and they do it by amplifying small differences into ratios of evidence. That makes them a
lens, and a lens has to be characterised before its images are read. Parameter intervals somewhat
too narrow are a nuisance. A mechanism ranked wrongly is a conclusion.

We have been on both sides of that line. In 2007 we proposed that P2X2 passes through a flip state
before opening. The claim had a signature you could point to in the record, a delay before the
current appears, and other groups later confirmed it. Last year we reported that the same receptor
activates asymmetrically (Communications Biology, 2025). Nothing in those recordings shows asymmetry
the way the delay showed flip. Atomistic simulations of a zebrafish P2X4 homologue corroborate it
and supply a direction the electrophysiology cannot resolve, and that support comes from outside the
record. Which of nine kinetic schemes the recordings prefer came from ratios of evidence alone. We
ranked the nine twice, under two likelihoods that shared the interval averaging and the gating
variance and differed in whether each interval is conditioned on the data already observed, and the
two orderings put a different mechanism first, erasing a margin of six. We reported the sensitivity.
The direction we gave for it rested on argument, because on recordings there is no ground truth to
measure a likelihood against.

Nothing in ordinary practice adjudicates it either. A recursive likelihood follows the record
whether its assumptions hold or not, so its residual cannot be read the way a least-squares residual
can, and fitting well says nothing about whether the information a likelihood reports is the
information it has. The error we were afraid of is the kind that passes every check we know how to
run, which is why we went outside the fit for one.

This class of model allows it. A Markov ensemble simulates exactly where its likelihood cannot be
evaluated exactly, so the question can be put from outside, over recordings drawn from parameters
that are known. The manuscript measures both classical identities against those simulations, the
score mean and the information equality, along a ladder of eight macroscopic likelihoods running
from least squares on the mean current to a filter that conditions each interval average on both of
its endpoints. Ten thousand recordings at each of 210 design cells per member, over channel number,
instrumental noise and acquisition interval. Two states is a limit and also a choice: three would
have been more than we could carry while still understanding what we were seeing, and two isolates
the error the approximation itself commits, with no misspecified mechanism mixed into it. An
implementation error shows up as a departure like any other, so the same pass tests the code, and
these are not small codes.

The measured object is a matrix, and it is the same matrix any correction to an evidence ratio has
to start from, through a volume term and an effective-sample rescaling that the Discussion sets out.
That is why it had to be measured first. Members that leave the correlated fluctuation unmodelled
misreport their information by up to an order of magnitude, threefold on the error bar, wherever the
residual retains memory, and no single rescaling repairs it, because the distortion has directions
in parameter space. Conditioning on one endpoint of each interval is worse than ignoring the
averaging altogether, which we had not expected and can now explain. Our own calibrated member is
audited on the same grid rather than exempted from it: in the few-channel corner its distortion runs
from 0.645 to 1.701, it errs in both directions, and the manuscript says where.

That measurement reaches the two instruments of our own comparison. The likelihood used as the
published control is one of the members that fail, accumulating about twenty times the information
it holds across a record at macroscopic channel counts, and the recursion is exactly what returns
that factor to one. We do not recompute the P2X2 evidences here, and the Discussion says so. What
the measurement does say about them is that the arm which favoured the alternative mechanism is
distorted in the simplest model there is, so the calibrated instrument supports the reading we
published rather than the control that displaced it.

eLife published Münch and colleagues on this problem class in 2022, and they asked our question of
their own filter, counting how often the true rate matrix falls inside a credibility volume of given
mass. What their count tests is a posterior, so the prior and the sampler are in the verdict with the
likelihood, and it costs a fit per replicate at every design point. The identities we use are
evaluated on the likelihood alone at parameters that are known, they cost one pass, and what comes
back is a matrix that says by how much and in which parameter the reported uncertainty departs, and
whether the departure is committed at each sample or accumulated across them. That is what puts 210
cells per member within reach. We do not offer another filter. Our best member coincides with an
integrated-measurement Kalman filter known since 1988 to about one part in 10^8, and we say so in
the Discussion. What is new is the measurement.

The study is simulation on a two-state scheme at an open probability of one half, driven by a single
concentration jump, and no experimental recordings are analysed. Where a preparation yields many
exchangeable sweeps the classical route remains strong and the manuscript says so: non-stationary
fluctuation analysis returns the channel number and the unitary current with no likelihood at all,
and a sweep-level bootstrap returns honest intervals with no model. The case for a filter is
strongest where that route is weakest, at few sweeps, under rundown, on a single precious recording,
and whenever kinetics and amplitudes are wanted from the same record.

One reading of the map does not depend on where a level set is drawn. No region of the measured
plane gives both information about the unitary current and a trustworthy least-squares error bar,
and the two boundaries sit one to three decades of instrumental noise apart everywhere they were
measured. The recordings rich enough to determine the amplitudes are the ones whose classical
intervals are the wrong width.

Praising a method and withholding it comes to very little, so the calibrated member, its score, its
Fisher information and the exact simulator are in a small library with R and Python bindings, where
the distortion matrix, its eigenvalues and the first-order bias come back from one call. A reader
can put the same two identities on their own scheme instead of extrapolating from ours.

Yours sincerely,

Luciano Moffatt

---

## Impact statement (submission form, 15 to 30 words, third person)

A likelihood for macroscopic ion-channel currents whose reported information is measured against
the exact process, with the apparatus to repeat the measurement on any scheme, available from R
and Python.

(30 words, unchanged from v1.)

Finding-forward alternative, if the venue prefers a finding. The v1 note called its longer form 29
words; it was 36 and would have been rejected by the form. Trimmed to 30:

In a published comparison of nine gating schemes the ranking moved with the likelihood; measuring
score and information against exact simulation shows which reported uncertainties can be believed,
and where.

---

## Notes carried over from v1 that still apply

The naming and ranking restrictions withdrawn on 2026-08-10 stay withdrawn. The published names may
be used; the bridges IR = MacroIR and INR = MacroINR are closed at
papers/_program/nomenclature.md:286-296.

No mention of Supplementary Table S4 of the 2025 paper. The prepared answer for a referee who opens
that supplement and finds Table S4 and the two non-recursive runs at effective sample sizes of 3 and
4 is at the foot of cover_letter.md and is unchanged. Adding the third beat above raises the chance
that question is asked, which is the price of having an ending.

The personal register appears twice, "we have been on both sides of that line" and "the error we
were afraid of", as in v1. The closing sentence about withholding a method is a third voice of a
different kind, borrowed from the dictated version, and it is the one place to check whether it
reads as warmth or as a flourish.
