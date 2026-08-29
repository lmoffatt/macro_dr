# Cover letter, draft v3

2026-08-28. Same content decisions as v2, written in the voice of the Introduction: short sentences,
open on what can be seen, no defensive clauses. v2 can be deleted. Numbers checked against the
current sections (210 cells per member, the twenty as accumulated information, the margin of six).
Still dropped and worth restoring if the 2025 paper supports it: "posteriors tens of times narrower
than the priors that went in".

---

Dear Professor [Senior Editor],

I am submitting "Likelihood approximations distort the ion channel kinetic information in
macroscopic currents" for consideration as a Research Article.

Some mechanisms leave a mark you can see without fitting anything. The flip state we proposed for
P2X2 in 2007 was one, a delay before the current rises, and other groups confirmed it later. Last
year we reported that the same receptor activates asymmetrically. That claim leaves no such mark. It
rests on ratios of evidence between nine kinetic schemes, and when we ranked them a second time
under a control likelihood, one that shares the interval averaging and the gating variance and
differs only in whether each interval is conditioned on the data already observed, a different
mechanism came first and a margin of six disappeared. We published the sensitivity and argued for
one of the two arms. Argument is what it was. On a recording there is no truth to measure a
likelihood against.

Nothing inside a fit settles this. A recursive likelihood follows the record whether its assumptions
hold or not, so its residual cannot be read the way a least-squares residual can, and fitting well
says nothing about whether the information a likelihood reports is the information it has. The error
we were afraid of is the kind that passes every check we know how to run.

So we went outside the fit. Channel gating is a Markov ensemble, and it simulates exactly even where
its likelihood cannot be evaluated exactly. That gives ten thousand recordings from parameters we
chose, and two classical identities to hold each likelihood to: at the true parameters the score
averages to zero, and its covariance equals the Fisher information the likelihood reports. We
measured both, for eight likelihoods, from least squares on the mean current to a filter that
conditions each interval average on both of its endpoints, across 210 design cells of channel
number, noise and sampling interval. It is a two-state scheme and simulated data throughout. Three
states would have been more than we could hold and still understand what we were seeing, and two
isolates the error the approximation itself commits. An implementation error shows up as a departure
like any other, so the same pass tests the code.

The results are blunt. Every member that leaves the correlated fluctuation unmodelled misreports its
information, by up to an order of magnitude, threefold on the error bar, and no single rescaling
repairs it, because the distortion points in particular directions in parameter space. Conditioning
on one endpoint of the interval turns out worse than ignoring the averaging altogether, which we did
not expect. Our own member is on the same grid as the rest, and in the few-channel corner it departs
in both directions, from 0.645 to 1.701. We say where.

One of the members that fail is the control from our own comparison. Across a record it accumulates
about twenty times the information it holds, and the recursion is what brings that factor back to
one. We do not recompute the P2X2 evidences here. What we can say is that the arm which favoured the
alternative mechanism is distorted in the simplest model there is, and the calibrated one is not.

eLife published Münch and colleagues on this problem in 2022, and they asked our question of their
own filter by counting coverage. A coverage count tests a posterior, so the prior and the sampler
are in the verdict with the likelihood, and it costs a fit per replicate at every design point. The
identities cost one pass at parameters that are known, which is what makes 210 cells affordable, and
what comes back is a matrix that says in which parameter the reported uncertainty is wrong and
whether the error is committed at each sample or accumulated across the record. We are not offering
another filter. Our best member agrees to one part in 10^8 with an integrated-measurement Kalman
filter known since 1988, and we say so. The measurement is the new thing.

One line of the map holds however the level sets are drawn. No region gives both a trustworthy
least-squares error bar and any information about the unitary current, and the two boundaries sit
one to three decades of noise apart. The recordings rich enough to determine the amplitudes are the
ones whose classical error bars are the wrong width.

Where a preparation gives many exchangeable sweeps, fluctuation analysis still answers the amplitude
question with no likelihood at all, and the manuscript says so. The case for a filter is strongest
on the single precious recording. Praising a method and not handing it over is worth little, so the
calibrated member, its score, its Fisher information and the exact simulator are in a small library
with R and Python bindings. The distortion matrix, its eigenvalues and the first-order bias come
back from one call, and a reader can put these two identities on their own scheme instead of taking
ours on trust.

Yours sincerely,

Luciano Moffatt

---

## Impact statement (submission form, 15 to 30 words, third person)

A likelihood for macroscopic ion-channel currents whose reported information is measured against
the exact process, with the apparatus to repeat the measurement on any scheme, available from R
and Python.

(30 words.)

Finding-forward alternative, trimmed to 30 (the version noted in v1 as 29 words was 36):

In a published comparison of nine gating schemes the ranking moved with the likelihood; measuring
score and information against exact simulation shows which reported uncertainties can be believed,
and where.
