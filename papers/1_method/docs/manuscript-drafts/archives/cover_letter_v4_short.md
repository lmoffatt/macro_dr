# Cover letter, draft v4, one page

2026-08-28. Half the length of v3. Everything the abstract already says is out; what stays is what
only a letter can say: the motive, the verdict on our own control, and the map line.

---

Dear Professor [Senior Editor],

I am submitting "Likelihood approximations distort the ion channel kinetic information in
macroscopic currents" for consideration as a Research Article.

We have reached the limit of what can be read off a macroscopic current by eye. The flip state we
proposed for P2X2 in 2007 sat within it: a delay before the current rises, later confirmed. The
questions now asked of these receptors, allosteric coupling between subunits and the asymmetries in
it, leave no such mark, and answering them takes an inference that can tell schemes apart when the
record alone cannot. We reported last year that P2X2 activates asymmetrically. That conclusion rests
on ratios of evidence between nine schemes, and under a control likelihood that differs only in the
recursion, a different mechanism came first and a margin of six disappeared. We published the
sensitivity and argued for one arm. Argument is what it was: on a recording there is no truth to
measure a likelihood against.

Nothing inside a fit settles this. A recursive likelihood follows the record whether its
assumptions hold or not, and fitting well says nothing about whether the information it reports is
the information it has. The error we feared is the kind that passes every check we know how to run.

So we went outside the fit. Channel gating simulates exactly even where its likelihood cannot be
evaluated exactly, which gives ten thousand recordings at known parameters and two classical
identities to hold each likelihood to: the score averages to zero, and its covariance equals the
Fisher information reported. We measured both for eight likelihoods, from least squares to a filter
conditioning each interval average on both endpoints, over 210 cells of channel number, noise and
sampling interval, in a two-state scheme so that what we measure is the approximation's own error.
The same pass tests the code.

The result reaches our own comparison. The control that reordered the schemes accumulates about
twenty times the information it holds across a record, and the recursion is what brings that factor
back to one. We do not recompute the P2X2 evidences; the arm that favoured the alternative is
distorted in the simplest model there is, and the calibrated arm is not. That member is on the same
grid as the rest and departs in the few-channel corner, in both directions, and we say where.

Münch and colleagues asked our question of their own filter in eLife in 2022 by counting coverage,
which tests a posterior and costs a fit per replicate. The identities cost one pass and return a
matrix that says which parameter is misreported and whether the error is committed at each sample
or accumulated across the record. We are not offering another filter; ours agrees to one part in
10^8 with an integrated-measurement Kalman filter known since 1988. The measurement is the new thing.

One line of the map holds however the level sets are drawn: no region gives both a trustworthy
least-squares error bar and information about the unitary current, and the two boundaries sit one
to three decades of noise apart. The recordings rich enough to determine the amplitudes are the
ones whose classical error bars are the wrong width.

The calibrated likelihood, its score and Fisher information, the exact simulator and the
distortion diagnostic are in a small library with R and Python bindings, so a reader can put the
same two identities on their own scheme rather than take ours on trust.

Yours sincerely,

Luciano Moffatt
