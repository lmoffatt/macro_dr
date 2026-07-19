# Context, the gap we fill, and the semiseparability = state-space = Kalman equivalence

> Discussion note, 2026-07-17. Synthesizes the narrative reframing (context → gap → contribution →
> payoff), the confidence-interval-validity gap, the bootstrap history, and the cross-field
> equivalence (ion-channel recursive filters = astronomy semiseparable Gaussian processes =
> Kalman). Feeds the Introduction (P1-P4), the Discussion (feasibility + Kalman concession), and D-3.
> Most papers cited here have a PDF in `docs/bibliography/` and an entry in `manuscript-drafts/biblio.bib`
> (exceptions: the classic statistics references Godambe 1960 and White 1982 have bib entries but no PDF).

---

## 1. The narrative spine (context → gap → fill → payoff)

**Context.** A macroscopic ion-channel current carries the kinetic information twice: in the
*deterministic* relaxation of the mean, and in the *stochastic* fluctuations around it (the number of
channels, the unitary current, and the temporal correlation that encodes the rates). Using both at
once requires a likelihood that models the mean and the fluctuation structure jointly. Such
likelihoods exist, a small lineage of Gaussian approximations (Milescu 2005; Moffatt 2007 BJ;
Celentano & Hawkes 2004; Stepanyuk 2011/2014; Münch 2022; MacroIR / Moffatt & Pierdominici-Sottile
2025). They differ along two axes: whether they filter (recursion) and how they treat the acquisition
window. Milescu 2005 sits at one end of the first axis: it keeps only the diagonal of the fluctuation
covariance and discards the temporal correlation, trading that information for speed.

**The gap.** These are approximations of an intractable exact likelihood, and no one has established
**when they are faithful**: when the reported uncertainty is the real one and when it distorts the
inference. Goodness of fit does not answer it (an approximation can fit the mean beautifully and report
an information matrix off by more than an order of magnitude). So the field either reverts to
deterministic least-squares on the mean, discarding the fluctuation information, or uses the
approximations on faith.

**What we fill.** Because the forward process is exactly simulable, we measure, against that ground
truth, how faithfully each approximation reports what the data know. We map where each holds and breaks
across channel number, acquisition interval, and instrumental noise, and identify the one that stays
calibrated.

**Payoff.** The kinetic information in the fluctuations stops being something the field discards on
distrust and becomes something it can use with a validated tool and a known regime of trust.

---

## 2. Two problems with least-squares on the mean

Least squares on the deterministic mean current (the dominant practice: Clerx et al. 2019; IonBench,
Owen & Mirams 2025; Wang et al. 2012 opts out of the likelihood explicitly on cost grounds) has two
distinct problems:

1. **It throws away the fluctuation information** — the unitary current, the channel number, and the
   temporal correlation that identifies the rates. Del Core & Mirams 2025 states the consequence
   directly: "standard methods for deterministic models do not distinguish between stochastic channel
   gating and measurement error noise, resulting in biased estimates."
2. **It has no validated estimator of the uncertainty on the recovered parameters.** This is the sharper
   half, and it is the paper's real target. See §3.

Non-stationary fluctuation analysis (NSFA) is the standard objection to problem 1, and it must be
scoped precisely: NSFA uses the second moment, but as a two-moment regression that returns the unitary
current, the channel number, and peak open probability, **not** the temporal correlation structure in a
likelihood that identifies kinetics. Stepanyuk et al. (2014) say it themselves: "the unitary current is
virtually the only parameter that can be reliably obtained from this type of analysis," and "to the best
of our knowledge, kinetic rates have never been estimated for any synaptic receptors in their intrinsic
environment."

---

## 3. How the field quantifies parameter uncertainty, and why it is not resolved

The field DOES report confidence intervals on recovered kinetic parameters. What is **not** established
is whether those intervals are valid. The landscape:

| Method | How it reports parameter uncertainty |
|---|---|
| Moffatt 2007 (BJ, the methods paper) | inverse **Fisher** information (variance of the score), used because the Hessian was not always positive-definite; inverse Hessian as fallback; FIM ellipses cross-checked against Monte-Carlo scatter |
| Celentano & Hawkes 2004 | empirical standard errors from the scatter of the MLE over simulated datasets (mean ± s/√N); estimates shown to be approximately normal, as ML theory predicts |
| Milescu et al. 2005 | empirical scatter of the MLE over many simulated datasets (Monte Carlo), noting the estimate distribution is non-Gaussian |
| **Moffatt & Hume 2007 (JGP, the P2X2 application)** | **bootstrap** (balanced resampling, 100 replications, percentile intervals) |
| Stepanyuk 2011 / 2014 | **bootstrap** (resample traces with replacement, refit each replicate), cross-checked vs the inverse Hessian (observed information) |
| Münch et al. 2022 | full Bayesian posterior sampled by Hamiltonian Monte Carlo / NUTS |
| Moffatt & Pierdominici-Sottile 2025 (Comm Biol) | full Bayesian posterior via parallel-tempered MCMC + thermodynamic integration |

**The gap, precisely.** The reported model-based CI is the model's own Fisher (H). Under the Gaussian
approximation (misspecified by construction) the correct frequentist covariance is the sandwich
`H⁻¹ J H⁻¹` (Godambe 1960; White 1982), and the two differ (the information distortion). **No one in the
ion-channel field has formed the sandwich `H⁻¹ J H⁻¹`, nor used the discrepancy between the model
information `H` and the score covariance `J` as a misspecification / CI-validity diagnostic against exact
simulation** (Moffatt 2007 already computes the score, the Fisher, and the model-based score covariance;
what is new is estimating `J` against ground-truth simulation and forming the sandwich). The closest,
Münch 2022, gives an `N_ch` rule of thumb (be careful for both algorithms when `N_ch` is roughly in
[10, 100]) but does not compute any of these, and explicitly notes the Fisher information
matrix is singular (breaking Cramér-Rao) and discusses bootstrap as the route an ML method would need
while arguing it fails under bias. Stepanyuk compared bootstrap vs the inverse Hessian at one operating
point and found agreement, which is exactly what happens in a benign regime; the paper shows that agreement is the
exception, and maps where the model-based CI becomes overconfident (up to ~16x in the failure corners).

So the contribution is not "we invented CI validation" (the sandwich/White/Godambe machinery is old
statistics; bootstrap and Fisher both existed in-domain). It is the first characterization, against exact
simulation, of **where the reported confidence intervals of these likelihoods are trustworthy and where
they are overconfident** — plainly: the field reports error bars on kinetic rates and never validated
them.

---

## 4. Bootstrap history (corrected)

An earlier automated search concluded "Moffatt 2007 does not use bootstrap." That was reading the wrong
2007 paper. There are **two** Moffatt 2007 papers:
- **Moffatt 2007, Biophys. J. 93(1):74-91** — the methods paper. Reports the **inverse Fisher
  information** (the Hessian was not always positive-definite; inverse Hessian as fallback).
- **Moffatt & Hume 2007, J. Gen. Physiol. 130(2):183-201** — the P2X2 application (same receptor as the
  Comm Biol 2025 demonstration). Uses **bootstrap**: "To get estimates of the standard errors in the
  parameter determinations we used a bootstrap approach [...]. A balanced resampling was obtained by
  concatenating 100 copies of the original set of traces ... 100 bootstrap replications ... Each bootstrap
  replica was refitted ... Errors were estimated from percentile intervals."

So bootstrap for kinetic-parameter CIs is not unheard of in-domain (Moffatt & Hume 2007; Stepanyuk
2011/2014), but it is rare, it is expensive (refit per replicate), and it does not tell you *why* or
*where* the cheaper model-based CI fails. That localization is what the information-distortion measurement
provides.

---

## 5. Computational scaling, stated honestly

Celentano & Hawkes 2004 is the full-covariance likelihood: it forms and solves the `n x n`
between-time-point covariance matrix. The authors report the binding cost as scaling with the **square**
of the number of **time points** `n` (the `O(n^2)` covariance-matrix generation; a dense solve is
`O(n^3)`). They test their method (CVF) extensively, on simulated models and on real GABA_A traces, but
must subsample each recording to 100-200 time points rather than use the full trace, which is why the
tractable methods exist.

Two routes make it tractable, and they compute the same likelihood:
- the **recursive filter** (Moffatt 2007; a Kalman filter), and
- the **quasiseparable covariance** (Stepanyuk 2011/2014; called semiseparable in celerite/astronomy).

Both remove the **time-point** cubic (near-linear in `n`). Both also pay a one-time `O(k^3)` for the
**propagator** `exp(Q dt)` (the eigendecomposition `Q = U Λ U⁻¹`) per distinct rate matrix `Q`, since `Q`
depends on the parameters. Stepanyuk's "linear in the number of states" is the leading per-sample term,
and their `O(k^3)` eigendecomposition is disclosed in their own operation count (amortized over the few
perturbations), not hidden. Their linear-vs-cubic contrast with Kalman is real at the per-sample level:
Stepanyuk's own count gives Moffatt's Kalman recursion `(2k^3 + 5k^2 + 3k·n_O)·N·N_T` operations (cubic
in states **per sample**) against their near-linear-per-sample recursion, a measured ~58x fewer operations
on a 7-state, 200-trace benchmark.

**Honest statement for the paper:** the tractability comes from removing the time-point cubic
(recursion/filtering); the residual state cost is the `O(k^3)` eigendecomposition of the propagator,
incurred **once per distinct rate matrix** (per perturbation), with `O(k^2)` matrix-vector propagation
per interval thereafter. Report **"O(k^3) per distinct rate matrix, O(k^2) per interval, linear in the
number of intervals, independent of channel number"** and do not claim linear-in-states. (Whether the
current MacroIR code re-forms `exp(Q dt)` at full `O(k^3)` each interval or reuses the eigenbasis across
intervals that share `Q` is an implementation detail to read from the code, not assume.) Handling
arbitrary interval spacings (the MacroIR advantage) is what lets a
recording be summarized by intervals spaced logarithmically after each perturbation; this is a structural
property, not a benchmarked speed result (cost benchmarking is deferred to the program's efficiency
component).

---

## 6. semiseparability = state-space = Kalman: one class, three representations

The recursive filter, the semiseparable-covariance Gaussian likelihood, and the state-space Gaussian
process are three equivalent representations of a single class of Gaussian process: a generative model
(the linear stochastic differential equation / state-space form), an algorithm that evaluates its
likelihood (the Kalman filter), and a structural property of its covariance (semiseparability). A
Gaussian process whose kernel has a rational power spectral density is exactly the stationary solution of
a linear time-invariant stochastic differential equation, whose likelihood is evaluated by a Kalman
filter in `O(n)`; the covariance of such a process is semiseparable, which is what makes the direct
covariance algebra `O(n)` too (Hartikainen & Särkkä 2010; Särkkä & Solin 2019; Foreman-Mackey et al.
2017). The three communities reached it by three routes:
- ion channels via the recursive filter (Moffatt 2007);
- ion channels via quasiseparable covariance (Stepanyuk 2011/2014);
- astronomy via fast Gaussian-process regression (celerite).

External corroboration in our own domain: **Stepanyuk et al. 2011 explicitly identify Moffatt 2007 as a
Kalman filter** — "the problems mentioned above are overcome in a recursive algorithm, which utilizes
Kalman filter for the maximum likelihood estimation of kinetic parameters [ref: Moffatt 2007]." So
conceding that MacroIR is an (integrated-measurement) Kalman filter is not self-criticism; it agrees with
what an independent group already stated in print.

---

## 7. The astronomy connection (why the technique is proven, and where the barrier really is)

The semiseparable / state-space Gaussian-process machinery **took off**, spectacularly, in astronomy
(stellar and exoplanet time-series analysis), not in ion channels:

- **celerite** (Foreman-Mackey, Agol, Ambikasaran & Angus 2017) factorizes a semiseparable covariance in
  `O(N J^2)` rather than `O(N^3)` and is a widely adopted tool for stellar variability, exoplanet transit
  and radial-velocity time series, and Gaussian-process regression on large astronomical datasets.
- The integrated / augmented-measurement construction that is the core of MacroIR (augment the state with
  the running integral, propagate by a van Loan block matrix exponential, condition the endpoint on the
  observed integral) is an active area in astronomy: Rubenzahl et al. 2026 (*Astron. J.* 171:259,
  arXiv:2601.02527) uses exactly this van-Loan-on-augmented-integral pattern for exposure-time-averaged
  and overlapping radial-velocity measurements.

**The point for the paper.** The same technique that never became standard in ion channels is a proven,
widely used tool one field over. Its failure to become the field's standard tool is plausibly a matter of
habit (least-squares on the mean) and the unvalidated status of these likelihoods, rather than any failure
of the mathematics (proven in astronomy). Moffatt, Stepanyuk, and Münch did adopt the machinery in-domain;
it simply never displaced the deterministic-mean default. The paper removes the "unvalidated" barrier,
which is the one thing standing between the field and a device the rest of quantitative time-series science
already trusts.

---

## References (PDFs in `docs/bibliography/`, entries in `manuscript-drafts/biblio.bib`)

**Ion-channel kinetic likelihoods and uncertainty**
- Moffatt, L. (2007). Estimation of ion channel kinetics from fluctuations of macroscopic currents.
  *Biophys. J.* 93(1):74-91. doi:10.1529/biophysj.106.101212. [inverse Fisher information; Hessian fallback]
- Moffatt, L. & Hume, R. I. (2007). Responses of rat P2X2 receptors to ultrashort pulses of ATP...
  *J. Gen. Physiol.* 130(2):183-201. doi:10.1085/jgp.200709779. [**bootstrap**; PDF added 2026-07-17]
- Celentano, J. J. & Hawkes, A. G. (2004). *Biophys. J.* 87(1):276-294. doi:10.1529/biophysj.103.036632.
- Milescu, L. S., Akk, G. & Sachs, F. (2005). *Biophys. J.* 88(4):2494-2515. doi:10.1529/biophysj.104.053256.
- Stepanyuk, A. R., Borisyuk, A. L. & Belan, P. V. (2011). Efficient maximum likelihood estimation of
  kinetic rate constants from macroscopic currents. *PLoS ONE* 6(12):e29731. doi:10.1371/journal.pone.0029731.
  [**calls Moffatt 2007 a Kalman filter**; bootstrap; PDF added 2026-07-17]
- Stepanyuk, A., et al. (2014). *Front. Cell. Neurosci.* 8:303. doi:10.3389/fncel.2014.00303.
  [source of the two NSFA quotes in §2; PDF on disk]
- Münch, J. L., Paul, F., Schmauder, R. & Benndorf, K. (2022). *eLife* 11:e62714. doi:10.7554/eLife.62714.
- Moffatt, L. & Pierdominici-Sottile, G. (2025). *Communications Biology* (P2X2 demonstration).

**Misspecification / sandwich machinery**
- Godambe, V. P. (1960). *Ann. Math. Stat.* 31:1208. White, H. (1982). *Econometrica* 50(1):1-25.

**Astronomy: semiseparable = state-space GP = Kalman**
- Foreman-Mackey, D., Agol, E., Ambikasaran, S. & Angus, R. (2017). Fast and scalable Gaussian process
  modeling with applications to astronomical time series (**celerite**). *Astron. J.* 154:220.
  doi:10.3847/1538-3881/aa9332. arXiv:1703.09710. [PDF added 2026-07-17]
- Hartikainen, J. & Särkkä, S. (2010). Kalman filtering and smoothing solutions to temporal Gaussian
  process regression models. *2010 IEEE Int. Workshop on Machine Learning for Signal Processing (MLSP)*,
  pp. 379-384. doi:10.1109/MLSP.2010.5589113. [GP regression ≡ Kalman/state-space; PDF added 2026-07-17]
- Särkkä, S. & Solin, A. (2019). *Applied Stochastic Differential Equations*. IMS Textbooks 10,
  Cambridge UP. doi:10.1017/9781108186735. [PDF added 2026-07-17]
- Rubenzahl, R. A., Hattori, S., Särkkä, S., Farr, W. M., Luhn, J. K., Bedell, M. & Foreman-Mackey, D.
  (2026). Scalable Gaussian Processes for Integrated and Overlapping Measurements via Augmented
  State-space Models. *Astron. J.* 171:259. doi:10.3847/1538-3881/ae4d0b. arXiv:2601.02527.
  [integrated-measurement augmented state-space, = MacroIR's device; PDF on disk. §7 prose de-duplicated
  2026-07-17 to cite this once (it had appeared as both "arXiv:2601.02527" and "Rubenzahl et al. 2026").]

*[Verified 2026-07-17: Hartikainen & Särkkä 2010 = IEEE MLSP pp. 379-384, doi:10.1109/MLSP.2010.5589113
(from the paper + the Särkkä 2013 reference list). arXiv:2601.02527 = Rubenzahl et al. 2026, published
*Astron. J.* 171:259, doi:10.3847/1538-3881/ae4d0b (one paper, not two). PDFs for celerite, Hartikainen &
Särkkä 2010, Särkkä & Solin 2019 now in docs/bibliography/ with biblio.bib entries. No §6/§7 citation
remains unverified.]*
