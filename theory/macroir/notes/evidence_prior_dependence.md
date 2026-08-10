---
date: 2026-08-08
status: working note, not a supplement yet
scope: marginal likelihood (evidence), its exact decompositions, its dependence
       on the prior, and where a correction under approximate likelihood would
       be a contribution
feeds: the planned evidence-verification paper; §6 of this note is the honest
       novelty ledger and should be read before claiming anything
---

# Evidence, prior dependence, and the information-gain matrix

This note consolidates a discussion held on 2026-08-08. It has two purposes. The
first is to write down, once, the exact identities relating the marginal
likelihood to the Kullback–Leibler (KL) divergence between prior and posterior,
so they stop being re-derived. The second, and the more important one, is to
record **which of these are classical and which are not**, because the honest
answer is that almost all of the mathematics here is between 1950 and 1998, and
the contribution lies elsewhere. Section 6 is that ledger.

Notation throughout: $\theta$ the parameter vector ($p$ components), $y$ the
data, $L(\theta) = p(y \mid \theta)$ the likelihood, $p(\theta)$ a **proper,
normalized** prior, and

$$Z = p(y) = \int L(\theta)\,p(\theta)\,d\theta$$

the marginal likelihood, also called the evidence. $\mathrm{KL}(q \| r) = \int q
\log(q/r)$. We write $\mathrm{KL}_{+} = \mathrm{KL}(\text{post} \| \text{prior})$
and $\mathrm{KL}_{-} = \mathrm{KL}(\text{prior} \| \text{post})$.

---

## 1. Exact identities

Nothing in this section assumes normality, large samples, or any approximation.
The only requirement is a proper prior, so that $Z_{\beta=0} = 1$.

### 1.1 The two KL identities

Starting from $p(\theta \mid y) = L(\theta) p(\theta) / Z$, so that
$p(\theta \mid y)/p(\theta) = L(\theta)/Z$:

$$\mathrm{KL}_{+} = \int p(\theta\mid y) \log \frac{L(\theta)}{Z} \, d\theta
= \mathbb{E}_{\text{post}}[\log L] - \log Z$$

$$\mathrm{KL}_{-} = \int p(\theta) \log \frac{Z}{L(\theta)} \, d\theta
= \log Z - \mathbb{E}_{\text{prior}}[\log L]$$

Rearranged, these are the two decompositions of the evidence:

$$\boxed{\ \log Z = \mathbb{E}_{\text{post}}[\log L] - \mathrm{KL}_{+}
\qquad\text{and}\qquad
\log Z = \mathbb{E}_{\text{prior}}[\log L] + \mathrm{KL}_{-}\ }$$

The first is the "fit minus complexity" or Occam decomposition; it is the exact
statement of which the Bayesian Information Criterion is an asymptotic limit.
The second measures the same quantity from the other end.

**A trap worth naming.** The subtrahend in the first identity is $\log Z =
\log \mathbb{E}_{\text{prior}}[L]$, the *log of the expectation*. It is not
$\mathbb{E}_{\text{prior}}[\log L]$, the expectation of the log. By Jensen the
two differ by a non-negative gap, so

$$\mathbb{E}_{\text{post}}[\log L] - \mathbb{E}_{\text{prior}}[\log L]
= \mathrm{KL}_{+} + \big(\log Z - \mathbb{E}_{\text{prior}}[\log L]\big)
= \mathrm{KL}_{+} + \mathrm{KL}_{-}$$

Using the difference of expected log-likelihoods as if it were $\mathrm{KL}_{+}$
overestimates the information gain by exactly $\mathrm{KL}_{-}$. The error is
silent: it is positive, of the right order, and grows with sample size just as
the KL does. See §7.

### 1.2 The variational identity, which contains both

For **any** distribution $q$ over $\theta$:

$$\log Z = \mathbb{E}_{q}[\log L] - \mathrm{KL}(q \| \text{prior})
+ \mathrm{KL}(q \| \text{post})$$

The first two terms are the evidence lower bound (ELBO). Setting $q =$ posterior
kills the third term and gives the first identity of §1.1; setting $q =$ prior
kills the second and gives the other. So the two identities are the two
endpoints of one family, and every other $q$ gives a strict lower bound. This is
the standard variational / free-energy decomposition.

### 1.3 The tempering path

Define the geometric path $p_\beta(\theta) \propto p(\theta)\,L(\theta)^\beta$,
which interpolates from the prior at $\beta = 0$ to the posterior at $\beta = 1$,
with normalizer $Z_\beta$. Let

$$U(\beta) = \mathbb{E}_{\beta}[\log L]$$

Two exact identities follow by differentiating $\log Z_\beta$:

$$\frac{d}{d\beta}\log Z_\beta = U(\beta)
\qquad\Longrightarrow\qquad
\log Z = \int_0^1 U(\beta)\, d\beta$$

$$\frac{dU}{d\beta} = \mathrm{Var}_{\beta}[\log L] \ \ge 0$$

The first is thermodynamic integration. The second is the fluctuation–dissipation
relation of the path: in the physical analogy $U$ is the internal energy and
$\mathrm{Var}_\beta[\log L]$ the heat capacity. It also says $U$ is
non-decreasing, so $U(0) \le \log Z \le U(1)$.

This is Fisher information in the literal sense, not by analogy: the tempered
family is exponential in $\beta$, and its Fisher information with respect to
$\beta$ is exactly $\mathrm{Var}_\beta[\log L]$.

### 1.4 The Jeffreys divergence is the whole span

Adding the two identities of §1.1:

$$\boxed{\ \mathbb{E}_{\beta=1}[\log L] - \mathbb{E}_{\beta=0}[\log L]
= \mathrm{KL}_{+} + \mathrm{KL}_{-} \ }$$

that is, $U(1) - U(0)$ is the **symmetrized KL (Jeffreys divergence)** between
prior and posterior. Combined with $U(0) \le \log Z \le U(1)$:

$$\log Z - U(0) = \mathrm{KL}_{-}, \qquad U(1) - \log Z = \mathrm{KL}_{+}$$

so the evidence is a point that cuts the segment $[U(0), U(1)]$ into the two KL
divergences. **The symmetrized KL is not a component of the evidence; it is the
length of the segment.** Any attempt to write $\log Z = A - \mathrm{Jeffreys}$
forces the antisymmetric combination $\mathrm{KL}_{+} - \mathrm{KL}_{-}$ into
the expression, which is not sign-definite and has no useful reading. The
components are the two endpoint fits, $\mathbb{E}_{\text{prior}}[\log L]$ and
$\mathbb{E}_{\text{post}}[\log L]$, depending on which end you measure from.

### 1.5 Everything as a weighted integral of one variance

Integrating $dU/d\beta = \mathrm{Var}_\beta$ and exchanging the order of
integration gives all three quantities as the same variance weighted differently:

$$\mathrm{KL}_{+} = \int_0^1 \beta\,\mathrm{Var}_{\beta}[\log L]\, d\beta$$
$$\mathrm{KL}_{-} = \int_0^1 (1-\beta)\,\mathrm{Var}_{\beta}[\log L]\, d\beta$$
$$\mathrm{KL}_{+} + \mathrm{KL}_{-} = \int_0^1 \mathrm{Var}_{\beta}[\log L]\, d\beta$$

with $\log Z = \int_0^1 U(\beta) d\beta$ as the fourth member. All exact, all
free of any normality assumption.

### 1.6 A diagnostic that comes for free

Where $\log Z$ falls inside $[U(0), U(1)]$ is decided by the shape of
$\mathrm{Var}_\beta$ along the path. If the variance concentrates near
$\beta = 0$, information arrives early, the prior is overwhelmed quickly, and
$\mathrm{KL}_{-}$ dominates. If it concentrates near $\beta = 1$, the prior
holds until the end and $\mathrm{KL}_{+}$ dominates. **The ratio
$\mathrm{KL}_{+}/\mathrm{KL}_{-}$ therefore reports where along the path the
learning happens**, and the curve $\mathrm{Var}_\beta$ against $\beta$ is a
diagnostic of the temperature ladder itself.

This costs nothing new in a program that already computes the evidence by
thermodynamic integration with parallel tempering: the samples at each $\beta$
already exist and $U(\beta)$ is already being averaged. The variance at each
$\beta$ comes from the same samples.

### 1.7 Other exact decompositions

**Chib's pointwise identity.** For any $\theta^{*}$, with no integral at all:

$$\log Z = \log L(\theta^{*}) + \log p(\theta^{*}) - \log p(\theta^{*} \mid y)$$

Exact; the whole difficulty moves into evaluating the posterior density at a
point. This is the basis of Chib's method.

**Path sampling.** Thermodynamic integration is one case of Gelman & Meng's path
sampling: any continuous path between prior and posterior works. Changing the
path changes the estimator variance, not the value.

**The chain rule (prequential decomposition).**

$$\log Z = \sum_t \log p(y_t \mid y_{1:t-1})$$

The evidence is the sum of one-step-ahead predictive log-densities. Exact, no
assumptions.

**Averaged over data.** $\mathbb{E}_y[\mathrm{KL}_{+}] = I(\theta ; Y)$, the
mutual information between parameters and data. This is Lindley's expected
information gain and is why Bayesian optimal experimental design maximizes mutual
information. And $\mathbb{E}_y[\log Z] = -H(Y)$, minus the entropy of the prior
predictive.

### 1.8 The chain rule and this algorithm

MacroIR computes, at fixed $\theta$,

$$\log L(\theta) = \sum_t \log p(y_t \mid y_{1:t-1}, \theta)$$

marginalizing the hidden channel state but not $\theta$. The evidence is the same
telescoping structure one level up, marginalizing $\theta$ as well. The two are
not the same object and the note should not conflate them.

What follows from putting them side by side is the useful part: **each term of
the sum is a predictive density, and the per-interval distortion diagnostic of
the main line of work is a statement about whether those predictive densities are
calibrated.** If they are not, the error propagates term by term into
$\log L(\theta)$ and from there into the evidence. The evidence correction
derived in the Posterior Information Distortion supplement is exactly that
propagation accounted for.

---

## 2. What the evidence is not

### 2.1 $Z$ carries units

$Z$ is a density in data space. Rescaling $y$ rescales $Z$. Its numerical value
is therefore not interpretable, and in particular **$Z = 1$ is not a meaningful
condition**: it can be arranged by a change of units. Only ratios between models
on the same data (Bayes factors) are meaningful.

### 2.2 When is the posterior equal to the prior

$p(\theta \mid y) = p(\theta)$ if and only if $L(\theta) = Z$ for all $\theta$ in
the support, that is, the likelihood is **constant** in $\theta$. In that case
$Z$ equals that constant, which is not in general 1. In the picture of §1.4 this
is the degenerate case: the segment collapses, $U(0) = U(1) = \log Z$, both KLs
vanish, and $\mathrm{Var}_\beta[\log L] \equiv 0$.

The scale-free criterion for "nothing was learned" is $\mathrm{KL}_{+} = 0$,
which is invariant under reparametrization of both $y$ and $\theta$. The evidence
has neither invariance.

### 2.3 Improper priors

With an improper prior, $Z$ is defined only up to an arbitrary multiplicative
constant, and that constant does not cancel in a Bayes factor unless the two
models share exactly the same improper prior on exactly the same parameters. This
is why the Bayes factor is undefined under improper priors while posterior
estimation is usually unaffected: in estimation the constant cancels between
numerator and denominator of Bayes' rule.

Jaynes' position on this is often misreported. He insisted that an improper
density has meaning **only** as the limit of a well-defined sequence of proper
ones, and argued at length that the marginalization paradoxes come from skipping
that limit. That position is not the blind spot. The asymmetry that matters here
is older and is Jeffreys': improper priors are acceptable for estimation and not
for testing, which is why Jeffreys insisted on a proper (Cauchy) prior on the
parameter under test.

---

## 3. The Gaussian layer, in terms of $\mathbf{G}$

Everything above is exact. This section is the Laplace approximation, and it is
where the connection to the existing diagnostic machinery is made.

Let the prior be $\mathcal{N}(\mu_0, \mathbf{H}_0^{-1})$, let $\hat\theta$ be the
maximum a posteriori estimate, $\mathbf{H}_L$ the likelihood precision (observed
or Gaussian Fisher) and $\mathbf{H}_{\text{post}} = \mathbf{H}_L + \mathbf{H}_0$.
Define the **information-gain matrix** in prior-whitened coordinates,

$$\mathbf{G} = \mathbf{H}_0^{-1/2}\,\mathbf{H}_{\text{post}}\,\mathbf{H}_0^{-1/2}$$

which is the object already defined in
`docs/Posterior_Information_Distortion/supplement_information_gain.tex`.

### 3.1 Three-term decomposition of the evidence

$$\log Z \approx \log L(\hat\theta)
\;-\; \tfrac{1}{2}\,\|\hat\theta - \mu_0\|^2_{\mathbf{H}_0}
\;-\; \tfrac{1}{2}\log\det \mathbf{G}$$

with $\|x\|^2_{\mathbf{H}_0} = x^{\mathsf T}\mathbf{H}_0 x$. The three terms read
as: best attainable fit; prior misfit, i.e. how far the data dragged the estimate
from the prior mean measured in prior units; and volume contraction from prior to
posterior, which is the information gain in nats.

### 3.2 The KL splits along the same seam

$$\mathrm{KL}_{+} \approx \tfrac{1}{2}\Big[\mathrm{tr}(\mathbf{G}^{-1}) - p
+ \log\det\mathbf{G}\Big]
+ \tfrac{1}{2}\|\hat\theta - \mu_0\|^2_{\mathbf{H}_0}$$

The bracket depends only on the spectrum of $\mathbf{G}$; in eigenvalues it is
$\tfrac{1}{2}\sum_i (1/\lambda_i - 1 + \log\lambda_i)$, which vanishes at
$\mathbf{G} = \mathbf{I}$ and is non-negative. The second term is the same prior
misfit as above.

**This separates the two axes that a prior-sensitivity study needs.** The
*width* of the prior enters through the spectrum of $\mathbf{G}$; the *distance*
between prior and truth enters through $\|\hat\theta - \mu_0\|^2_{\mathbf{H}_0}$.
They do not mix.

Consistency check: $\mathbb{E}_{\text{post}}[\log L] \approx \log L(\hat\theta) -
\tfrac{1}{2}\mathrm{tr}(\mathbf{H}_L\mathbf{H}_{\text{post}}^{-1})$, and since
$\mathbf{H}_L = \mathbf{H}_{\text{post}} - \mathbf{H}_0$ this equals
$\log L(\hat\theta) - \tfrac{1}{2}(p - \mathrm{tr}\,\mathbf{G}^{-1})$.
Subtracting the $\mathrm{KL}_{+}$ above returns §3.1 exactly.

### 3.3 Prior-induced bias and variance, in closed form

$$\hat\theta_{\text{MAP}} - \hat\theta_{\text{ML}}
= -\,\mathbf{H}_{\text{post}}^{-1}\mathbf{H}_0\,(\hat\theta_{\text{ML}} - \mu_0)$$

In prior-whitened coordinates $u = \mathbf{H}_0^{1/2}(\theta - \mu_0)$ this is
simply

$$\text{bias} = -\,\mathbf{G}^{-1} u_0$$

with $u_0$ the prior offset in prior units. The prior-induced bias is linear in
how off-centre the prior is, and the coefficient is the inverse information-gain
matrix. Likewise the whitened MAP covariance is $\mathbf{G}^{-1}$ and the
whitened maximum-likelihood covariance is $(\mathbf{G} - \mathbf{I})^{-1}$.

### 3.4 Four readings of one eigenvalue

Per eigendirection $i$ with eigenvalue $\lambda_i$ of $\mathbf{G}$:

| quantity | value |
|---|---|
| bias shrinkage factor | $1/\lambda_i$ |
| variance reduction vs ML | $(\lambda_i - 1)/\lambda_i$ |
| information gain (nats) | $\tfrac{1}{2}\log\lambda_i$ |
| effective parameters | $1 - 1/\lambda_i$ |

$\lambda_i \gg 1$: the direction is data-informed, the prior barely moves it and
barely shrinks it. $\lambda_i \to 1$: the direction is uninformed, the prior owns
it entirely. The effective-parameter column is the local-Gaussian $p_D$ of the
Deviance Information Criterion; the supplement's $p_{\text{eff}}$ is that with
$\mathbf{H}_{\text{post}}$ replaced by the sandwich-corrected precision.

### 3.5 Why this makes a prior sweep cheap

In this project the parameters are fitted in base-10 logarithmic coordinates, so
a log-normal prior on the natural parameters is a normal prior here and
$\mathbf{H}_0$ is diagonal with entries $1/\sigma_{0i}^2$. Sweeping the prior
width is therefore sweeping $\lambda_i$ directly, and the full dependence of
bias, variance, information gain and effective parameter count on the prior's
width and centre is available in closed form from the $\mathbf{H}_L$ that is
already computed. No new runs. What does need checking at a few points is whether
the Gaussian approximation holds, which is precisely what the main line of work
characterizes.

---

## 4. Prior dependence and robustness

### 4.1 The asymmetry between estimation and comparison

In estimation, prior details wash out as data accumulate. In model comparison
they never do: the prior width enters the Bayes factor as an offset of order the
log prior volume, which does not decay with sample size. This is a fact about the
arithmetic, not a philosophical position, and it is the reason the same prior can
be innocuous in one use and decisive in the other.

### 4.2 Bartlett's paradox

As the prior variance on the parameter under test grows, the Bayes factor tends
to favour the null regardless of the data. The conclusion is then governed by the
tail of the prior, a region where there is typically neither knowledge nor data.
Jeffreys–Lindley is the companion phenomenon driven by sample size instead of
prior diffuseness.

### 4.3 The Jaynesian reply, and what survives it

The reply is that the prior encodes a state of knowledge, so dependence on it is
correct behaviour rather than a defect. Against a genuinely elicited prior this
is right, and the "prior sensitivity" objection is weak. Three things survive it:

1. The sensitivity is concentrated in regions where no knowledge is claimed, and
   the maximum-entropy / transformation-group programme does not deliver a unique
   prior for most real problems, so convention re-enters.
2. The estimation/comparison asymmetry of §4.1 is arithmetic, not doctrine.
3. The strongest objection is not about the prior at all. If no model in the list
   is true, the posterior model probability converges on whichever is closest in
   KL, which is a statement about the list rather than about the world. A perfect
   prior does not repair that.

Point 3 loses much of its force if the list is treated as explicitly provisional
and the Bayes factor is read as a relative standing within the current list, to
be revised by expanding the list. That is a coherent position and is close to
what Gelman & Shalizi themselves advocate.

### 4.4 to 4.6 The robustness programme

The productive answer to prior dependence is to report a range or a threshold
instead of a number. Three named lines:

**Robust Bayesian analysis** (Berger). Specify a class $\Gamma$ of priors and
report $[\inf_\Gamma \mathrm{BF}, \sup_\Gamma \mathrm{BF}]$. If the ordering of
two models survives the whole class, the conclusion is robust to everything the
class represents.

**Bounds on Bayes factors** (Edwards, Lindman & Savage 1963; Berger & Sellke
1987; Sellke, Bayarri & Berger 2001; Held & Ott 2018). Take the prior in the
class most favourable to the alternative and report that extreme.

**Reverse-Bayes** (Good 1950). Invert the problem: fix the posterior conclusion
that would convince you and derive the prior that would be needed to reach it,
then ask whether that prior is plausible. Matthews' analysis of credibility and
the sceptical-prior literature are modern instances. A 2025 paper works out the
critical prior-variance threshold at which two analysts reach opposite
conclusions on the same data ("Bayes factor reversal").

**This is exactly the "map which region of prior space flips the ordering"
idea**, and it is a 75-year-old programme rather than a gap.

### 4.7 Priors as hyperparameters

Treating the prior's width as a hyperparameter with its own prior and integrating
is hierarchical modelling, and is the continuous model expansion that Gelman
prescribes instead of discrete comparison. The robustness map and the
hierarchical model are two readings of the same surface: the map displays the
sensitivity, the hierarchical model integrates over it with a weight. The
disagreement is doctrinal, not computational. The robust camp declines to
integrate on the grounds that there is no credible prior over priors.

---

## 5. Why Stan does not compute the evidence

Two independent reasons, worth keeping separate.

**Technical.** Hamiltonian Monte Carlo samples one target distribution and yields
no normalizing constant; obtaining $Z$ requires a temperature ladder, bridge
sampling or nested sampling, which are different algorithms. More specifically,
Stan's `~` statement drops normalization constants for speed, so in an ordinarily
written Stan program the log density is defined only up to an unknown constant
and $Z$ is not even well defined. Bridge sampling against a Stan fit requires
rewriting the model with `target +=` and full `lpdf` calls. This is the same
units-and-normalization point as §2.1, embedded in a tool.

**Doctrinal.** Gelman et al. (BDA3 §7.4) recommend against Bayes factors because
the marginal likelihood is highly sensitive to aspects of the model that are
typically assigned arbitrarily and are untestable from data. Above that sits the
Gelman & Shalizi rejection of the discrete-model-space framing. The offered
replacement is predictive: leave-one-out cross-validation and the expected log
predictive density.

**LOO is not immune, and this matters here.** The `loo` machinery is built from
the matrix of **pointwise** log-likelihood terms. Those are the same objects the
distortion diagnostic measures. If the per-interval predictive densities are
miscalibrated, the expected log predictive density inherits the error exactly as
the evidence does, with the difference that nobody is currently checking.

The consequence for framing: **both sides of the Bayes-factor-versus-LOO argument
assume the likelihood is correctly computed.** The validity of the per-interval
likelihood is a precondition for both, and it is not part of either literature.

### 5.1 Simulation-based calibration is a real competitor and must be cited

Simulation-based calibration (Cook, Gelman & Rubin 2006; Talts et al. 2018) draws
$\theta$ from the prior, simulates $y$ from the generative model, fits, and checks
that the rank of the true $\theta$ among posterior draws is uniform. If the
implemented likelihood is an approximation of the generative model, **SBC will
detect it**. Any claim that the Bayesian workflow has no tool for this is wrong
and a referee will say so.

Four differences survive, and they are the argument:

1. **SBC returns a verdict, not a correction.** A rank histogram says "something
   is off". It does not say that $k_{\text{off}}$ under-reports by 1.62 while
   $N_{\text{ch}}$ over-reports by 0.65. Here the magnitude *is* the correction.
2. **Cost.** SBC needs a prior and a full posterior fit per replicate, i.e.
   thousands of posterior samplings. The Bartlett check needs neither prior nor
   sampler (§7).
3. **SBC averages over the prior.** It reports calibration on prior average. The
   quantity a practitioner needs is local at $\hat\theta$, where the standard
   error will be reported.
4. **Marginal calibration can hold while joint calibration fails.** LOO-PIT and
   per-observation calibration checks see the per-sample component of the
   distortion. The correlation component, which is the dominant one for the
   recursive members, passes a marginal check while the joint distribution is
   wrong. This is precisely the split the distortion matrix decomposes into.

The defensible sentence for the discussion is therefore not "the predictive
workflow has no tool", but: *predictive checks detect the marginal component and
average over the prior; the correlation component and the local correction
require the geometry of the likelihood.*

---

## 6. Where the problem lives, and the vocabulary that already exists

### 6.1 The criterion is not "time"

The failure condition is **not** "the model evolves in time". It is that the
likelihood being evaluated is not the exact likelihood of the generative process.
Temporal autocorrelation is the most common *mechanism* by which that gap turns
into variance inflation, because the per-observation errors correlate and the
naive information adds them as if independent, but the criterion is the gap.

### 6.2 This has a name: information-biased composite likelihood

In the composite-likelihood literature the second Bartlett identity fails by
construction, and the machinery is standard. Following Lindsay (1982), a
composite likelihood is **information-unbiased** when $H(\theta) = J(\theta)$ and
**information-biased** otherwise, and the Fisher information is replaced by the
**Godambe information** $G = H J^{-1} H$, which is the sandwich. See Varin, Reid
& Firth (2011) for the review.

**The distortion matrix of this project is the information bias of that
literature, measured and resolved by direction.** Using their vocabulary plugs
the work into a community that will recognize the object immediately, and it is
also the honest thing to do.

### 6.3 The map of model classes

| class | status |
|---|---|
| exact likelihood, computed exactly (linear-Gaussian state space via Kalman, iid, fully observed Markov chains) | Bartlett holds; cheap and worth showing as a control |
| composite / pseudo-likelihood, GEE, block likelihood | failure known since Lindsay 1982; the sandwich is routine. Showing it fails is not news |
| **approximate filters and moment closures**: EKF, UKF, linear noise approximation, moment closure, Euler–Maruyama discretization of SDEs | **the unattended class.** The likelihood is presented as the model's own, and nobody checks. This is the bibliographic neighbourhood of this project (LNA filtering, aggregated-trajectory inference, stochastic compartmental models) |

### 6.4 What a systematic study would have to claim

Not "does the Bartlett identity fail" — that is settled case by case, and in
composite likelihood it is the founding fact. The defensible framing is
**magnitude and anisotropy**: how much it fails, in which regime, along which
directions, and the demonstration that the existing scalar repairs (effective
sample size, learning rate) are insufficient when the distortion is anisotropic.
The composite-likelihood community computes the Godambe sandwich as a correction
but rarely maps when and how much. The map is the contribution.

---

## 7. Cost, and why the check belongs in the library

### 7.1 A pointwise check is cheap

At a single $\theta_0$ the check needs $N$ simulated recordings and, for each, one
likelihood pass carrying score and Fisher. With $N = 10^4$ the relative error on
each variance component is about $\sqrt{2/N} \approx 1.4\%$, enough to resolve
distortions from 1.1 upwards. The number is not arbitrary; it is what the target
resolution requires.

**The cost is dominated by the simulator, not by the diagnostic.** Exact CTMC
simulation with a thousand sub-steps per measurement interval is roughly three
orders of magnitude more expensive per replicate than the likelihood pass. That
expense is buying the independence between simulator and likelihood, which is the
condition that makes the test valid (§7.3).

### 7.2 The expensive part of the campaign was neither

What made the production campaign expensive was the maximum-likelihood stage and
the sweep, not the diagnostic. The repository keeps them separate: there is an
older diagnostic battery that never calls the MLE stage. The optimization was
forced later and for a specific reason: for an algorithm **with bias**, the
distortion measured at the simulation point is contaminated by that bias, so one
must go to the maximum first.

So the precise statement is: **the pointwise check is cheap for an unbiased
algorithm and stops being cheap for a biased one.** For MacroIR, which is the use
case, it is cheap.

### 7.3 The load-bearing condition is independence, not availability

The check does not require simulating the truth; it requires simulating **the
model**. It is a self-consistency test, which is much weaker than knowing the
truth, and it is why the test is feasible at all: simulation is usually easy where
the likelihood is hard, which is the premise of the whole simulation-based
inference field.

What it does require is that **the simulator be independent of the approximation
under test**. Here it is: the ground truth comes from exact uniformization of the
continuous-time chain and shares nothing with the Gaussian closure. If one
simulated with the same closure, the test would pass identically and mean nothing.
This must be stated in Methods; it matters more than the existence of a simulator.

Where it cannot be done: models defined only through an unnormalized density
(Gibbs random fields, doubly-intractable problems), where sampling is as hard as
the likelihood; black-box models with no generative story; and the dangerous case,
where the only available simulator reuses the approximation being tested.

### 7.4 The honest limitation: it cannot be run on the experiment you care about

On real data there is no known truth and no ensemble. What can be done is to
validate a region of design space by simulation and then locate the real
experiment inside it. **That is what the region map is for**, and stated this way
it stops being a closing figure and becomes the mechanism by which the result
transfers to a real experiment.

There is a stronger and cheaper variant worth implementing: fit real data, obtain
$\hat\theta$, then simulate $N$ replicates **at $\hat\theta$** and check the
identity there. That is a local validation at the point that matters, for a
fraction of the fitting cost, and it does not rely on extrapolating from a generic
map.

### 7.5 Therefore it belongs in the library, as a first-class function

The stated reason nobody uses these estimators is that their validity is
uncharacterized and there is no usable implementation. Shipping the R/Python
`macroir` without the self-check reproduces exactly that situation. The check
should be a first-class entry point returning the distortion matrix and a verdict,
not a vignette. This is also, verbatim, what the project asked for on 2025-08-27:
the test declared next to the function, callable from the command line.

One deliberate decision: a shipped check will tell users that MacroIR itself has
distortion above 1 in some regimes. That is the selling point, not the problem,
but it should be a decision rather than an accident.

---

## 8. Automatic differentiation: why reverse mode does not serve

The Gaussian Fisher is assembled from the per-interval Jacobians of the predictive
moments,

$$I(\theta) = \sum_t \left[ \frac{\partial_\theta \mu_t\, \partial_\theta \mu_t^{\mathsf T}}{\sigma_t^2}
+ \frac{\partial_\theta \sigma_t^2\, \partial_\theta \sigma_t^{2\,\mathsf T}}{2\sigma_t^4} \right]$$

so what is needed is $\partial_\theta \mu_t$ and $\partial_\theta \sigma_t^2$ for
**every** $t$, not the gradient of the total scalar.

Reverse-mode automatic differentiation (backpropagation) returns the gradient of
one scalar output in a single backward pass. Applied to $\log L = \sum_t \ell_t$
it gives $\sum_t \nabla \ell_t$, the total score, and nothing per interval.
Recovering the $T$ individual Jacobians would need $T$ backward passes. Forward
mode costs $p$ forward passes and delivers all of them at once. With $p = 6$ and
$T$ large, forward-mode propagation is the right choice for the **analytic**
Gaussian Fisher, and it is what the engine implements.

**But the availability argument is weak, and it should not be made.** The Fisher
can always be obtained numerically, by finite-differencing the gradient, at a cost
of order $p$ gradient evaluations — the same order as forward mode. That is
exactly the route this project ended up taking as well: the numerical Fisher was
introduced in April 2026 and is still required in some cases. So the choice of AD
mode does not block anything, and claiming it does would be wrong.

What actually differs is not availability but **properties**, and the project paid
to learn each one:

- The analytic Gaussian Fisher is positive semidefinite **by construction**. The
  numerical Hessian evaluated away from the maximum can be indefinite, which is
  what forced the maximum-likelihood stage into a design that had deliberately
  avoided optimization (May–June 2026), and is why the distortion anchor was moved
  to the Gaussian form in July 2026.
- The numerical route inherits step-size sensitivity and is vulnerable to
  non-smoothness in the objective. The most expensive bug of the project was
  exactly a finite-difference discontinuity, produced by a branch on the sign of a
  displacement inside a safety factor.
- And the analytic Gaussian form is not universally valid: for the non-linear
  least-squares member it disagrees with the numerical one, so the numerical route
  is required there. Both are needed; neither dominates.

So the correct statement is: **forward propagation is the efficient route to the
analytic form, the numerical form is always available at comparable cost, and the
two are not interchangeable because they fail differently.** That last clause is
itself one of the results.

**On the 2025-10-17 decision.** Backpropagation was considered and discarded then,
with the stated reason that it prevented the individual gradients needed for a
Levenberg–Marquardt Hessian — an argument that became moot when that optimizer was
dropped. In hindsight the decision was still right, but for a smaller reason than
one might claim: forward propagation is the natural fit for the per-interval
Jacobians that the Gaussian Fisher would later need. It was not a necessary
condition, only a convenient one.

**And Stan — but only half of it.** Stan is reverse-mode over a single scalar
target, so per-observation scores are not exposed. Yet Stan users routinely
compute the pointwise `log_lik` vector already, because `loo` requires it, and
finite-differencing those terms would deliver the per-observation scores and hence
the variability matrix $J = \sum_t s_t s_t^{\mathsf T}$. **That half is easy, and
it is easy for everyone**: $J$ is a sum of outer products, so it is positive
semidefinite by construction and cannot go indefinite. The ingredients are already
being computed for another purpose and nobody assembles them.

The other half does not come for free, and this is the point that matters. The
sandwich is $H^{-1} J H^{-1}$ and also needs $H$, the curvature. **$H$ is not
positive definite by design.** A numerical Hessian evaluated away from the maximum
can and does go indefinite, which is exactly the failure that forced this project
into an optimization stage it had deliberately planned to avoid. So the cheap
route gives the well-behaved factor and leaves the ill-behaved one untouched.

What the analytic Gaussian Fisher supplies is precisely a **positive-semidefinite
by construction stand-in for $H$**: being itself a sum of outer products of
per-interval moment Jacobians weighted by positive scalars, it cannot go
indefinite, at any parameter value, without needing to sit at a maximum. That is
model-specific — it has to be derived for the observation model at hand — and it
is not something a general-purpose probabilistic programming language can hand
you.

So the correct division is: the variability matrix is universally available and
well behaved; a usable curvature matrix is the hard part, is model-specific, and
is what this work provides for this class of models. Stating it that way is both
more accurate and a stronger claim than an architectural impossibility argument.

Note also that this is a *different* mechanism from the evidence problem of §5,
which is about normalizing constants and the `~` statement. Both follow from the
same design goal, efficient posterior sampling, but they are separate obstructions
and should not be conflated in writing.

---

## 8bis. A hard-failure diagnostic: how often is the numerical Fisher not positive definite

### 8bis.1 Why it is worth reporting

Every other diagnostic in this programme reports a **graded** failure: your error
bar is off by a factor. This one reports a **hard** failure: you cannot form an
error bar at all, because the curvature is not usable. That is a qualitatively
different statement and deserves to be visible.

It is also cheap. At the simulation point it needs no optimization and no
bootstrap for the statistic itself, only a minimum-eigenvalue test per replicate.

And it was already observed. On 2026-06-10: individual numerical Hessians carry
negative eigenvalues "for all the algorithms, even for macro IR and even at
10,000 channels", and **on averaging they almost all disappear**, "except in the
algorithms that do not work" — macro R at $\Delta k_{\text{off}} = 1$, macro NR,
and macro MR at 100 channels. The observation discriminated between members and
was never turned into a reported quantity.

**The statistic is the mean, not the per-replicate fraction. This corrects an
earlier version of this note, which had it backwards.**

Per-replicate non-PD at $\theta_{\text{sim}}$ happens everywhere, including for
the good algorithm at 10,000 channels. It is finite-sample fluctuation of the
curvature and it measures how poorly determined the curvature is with that much
data — close to a restatement of "small effective sample" or "design near
non-identifiability". It is not a property of the algorithm and it does not
discriminate.

What discriminates is the definiteness of the **mean**. Under a correct
likelihood at the true parameter, $\mathbb{E}[H_{\text{obs}}] = I(\theta_0)$,
which is positive semidefinite always. So a negative eigenvalue in the averaged
numerical Fisher is not a graded measure, it is a **certificate of failure**: no
small-sample regime explains it, because the average converges to the Fisher
information. The hard-failure framing survives; only the estimator changes.

**Two consequences.** First, the statistic is available from the bootstrap
summaries already on disk (mean and percentiles of the minimum eigenvalue), so no
new emission and no re-run are needed. Second, the earlier claim in this note —
that the non-PD fraction is "the probability that a user cannot form an error
bar" — is **wrong and must not be used**. A user fits their data, converges to an
interior maximum, and there the Hessian is positive semidefinite by construction;
a non-PD result there is a convergence or numerical artifact, not a statistical
event.

**Which reframes the justification of the optimization stage, with less ambition
and better support.** It is not that the user needs to optimize in order to have
an error bar; the user already optimizes. It is that *this measurement* needs to
evaluate curvature away from the maximum, at $\theta_{\text{sim}}$, and that is
where indefiniteness lives. The optimization is a necessity of the measurement
design, not a recommendation of practice.

### 8bis.6 Three places to evaluate definiteness, three meanings

| evaluated at | meaning |
|---|---|
| $\theta_{\text{sim}}$, per replicate | finite-sample noise in the curvature; not diagnostic |
| $\theta_{\text{sim}}$, averaged | certificate of failure if not PSD; the Bartlett identity is violated |
| $\hat\theta$ (converged interior maximum) | PSD by construction; says nothing |
| over posterior samples | **validity of the quadratic approximation over the region that matters** |

The last row is the one with a theorem behind it. Bernstein–von Mises says that
under regularity the posterior converges to a Gaussian centred at the maximum with
the inverse Fisher as covariance. So non-PD over posterior draws is the failure of
those regularity conditions, or of being far from the regime where the theorem has
taken effect. It is not a numerical nuisance: it measures whether the theorem
underwriting every reported error bar is in force yet.

Caveats for that row: the causes are several (few data, non-linear
parametrization, multimodality), and it is **not invariant under
reparametrization**, unlike the KL. It answers "does the quadratic approximation
hold *in the coordinates I chose*", which is operationally what one wants but is
not a property of the model alone. Fitting in $\log_{10}$ coordinates is already
a choice of this kind.

**Its strongest use is narrative.** At present the maximum-likelihood stage enters
the paper as an unexplained methodological complication; the project deliberately
planned not to optimize and was forced to. If the non-PD fraction is reported,
the optimization stops being an implementation choice and becomes the consequence
of a measured result: in region X the curvature at the simulation point is
unusable in Y% of replicates, therefore the analysis must be anchored at the
maximum.

### 8bis.2 The absolute level is not interpretable

The fraction has three contributions and only one is the object of interest:

1. **Genuine non-identifiability of the design.** A flat direction gives a
   singular or indefinite curvature regardless of the algorithm. Identified in
   this project on 2026-02-13: measurements before the channels open carry no
   kinetic information, so the matrix is singular there.
2. **Finite-sample noise.** At the true parameter the *expectation* of the
   observed Hessian is the Fisher information, which is positive definite under
   identifiability; but the observed Hessian of a single replicate is a random
   matrix fluctuating around it, and in directions where the information is small
   relative to its own dispersion the eigenvalue changes sign with appreciable
   probability. **With few channels or short recordings this happens with nothing
   wrong at all.**
3. **Algorithmic distortion**, the quantity of interest.

Therefore: report members **side by side within a cell**, or relative to the exact
member, never as an absolute level. The cell fixes design, channel count, record
length and parameter count, so contributions 1 and 2 are shared and the difference
between algorithms isolates 3.

This also improves the narrative rather than weakening it: if part of the failure
is a property of the regime rather than of the algorithm, the optimization stage
becomes a consequence of working at few channels instead of an admission about
MacroIR.

### 8bis.3 Optimizing avoids the problem, it does not solve it

At an interior maximum the Hessian is positive semidefinite by definition, so
evaluating there guarantees the property by construction. What is gained in
definiteness is paid in anchoring: the curvature is now measured at $\hat\theta$,
which is random and correlated with the data, instead of at the fixed
$\theta_{\text{sim}}$. This is exactly the sim-versus-pool anchor distinction
already present in the code and in the `_sim` variants of the Figure 4
supplements, and it should be stated when the optimization stage is justified.

### 8bis.4 What the data on disk can and cannot support

The June numerical battery (`433ed13`) does emit
`Likelihood_Numerical_Fisher_Information`, together with eigenvalue spectra and
minimum eigenvalues for several matrices. But **every component in that file is a
`Probit_statistics` bootstrap summary** — mean and empirical percentiles — not
per-replicate matrices.

Consequence: the *fraction of replicates* with a non-PD numerical Fisher is **not
recoverable from what is on disk**. What is recoverable today is a diluted
version: whether a low percentile of the minimum eigenvalue crosses zero. That
gives a continuous map and a boundary, but not the sentence that speaks to a user
("in what fraction of your recordings will you be unable to form an error bar").

Getting the real statistic needs one extra scalar emitted per replicate — the
minimum eigenvalue of the numerical Fisher, or a boolean — which is a one-line
addition to the emitter since the matrix is already computed, plus a re-run.

**And the member where it matters most has no data at all.** The numerical Fisher
was never computed in production for the non-linear least-squares member: the LSE
scripts omit that stage deliberately and the dispatcher declares it. So if a
re-run is done, the one that pays best is not the five macro members of June but
LSE with the numerical Fisher — the variant script has existed since 2026-08-01
and has never been executed.

### 8bis.5 Note on the Figure 6 LSE boundary

This is the one place where a body figure rests on something the text contradicts,
and it should be resolved or disclosed.

Figure 6 draws a **dotted** boundary whose own legend in `figure_6.Rmd` reads
"LSE's reported error bar becomes honest (distortion < 1.15)", fitted as a power
law with slopes 1.01 in $k_{\text{off}}$ and 1.03 in $N_{\text{ch}}$. That
distortion is computed with the **Gaussian** Fisher, because the numerical column
in the LSE files is empty. Yet the project's own conclusion of 2026-08-02 is that
for least squares the Gaussian and numerical Fisher disagree and the numerical is
the appropriate one — a finding important enough that it is listed as one of the
two incidental results worth reporting, tied to the "deceptive simplicity" of
least squares.

The boundary is not necessarily wrong, and it may move very little. The point is
that there is currently no way to know.

There is a second asymmetry to disclose regardless. The neighbouring **dashed**
boundary is IR's, and states the converse about the same threshold: IR's error bar
*stops* being honest when the distortion exceeds 1.15. Both use the same numeric
criterion, but the Gaussian Fisher is the appropriate estimator for IR and,
by the project's own finding, is not for LSE. If the two boundaries end up drawn
with different estimators, the caption must say so; a reasonable reader will
assume they are comparable.

---

## 8ter. Multimodality, bifurcation, and cheap detection — parked for the evidence paper

This section is for the evidence-verification paper, not the current one. It is
recorded so it is not lost.

### 8ter.1 Bifurcation is the right word, phase transition is the metaphor

Phase transition is vocabulary for observable behaviour: first order is a jump in
the quantity, second order is continuity with a break or divergence in its
derivative (the susceptibility, which is a variance). Bifurcation is vocabulary
for the landscape: critical points appearing, disappearing or exchanging
stability. In mean field they are the same object seen from two sides.

The technical link that matters here: **a bifurcation requires a zero eigenvalue
of the Hessian**. A degenerate critical point is exactly $\det H = 0$. So the
non-PD diagnostic is not merely reporting that the Gaussian approximation is poor;
it is detecting that the posterior support includes the neighbourhood of a
degenerate critical point, that is, proximity to a bifurcation of the likelihood
surface. In catastrophe-theory vocabulary the codimension-1 degeneracy is the
fold and the codimension-2 one is the cusp, and the cusp is the classic picture of
a continuous parameter change producing a jump.

This separates the two mechanisms discussed earlier:

- **Multimodality is a genuine bifurcation.** Two local maxima; as data
  accumulate one disappears or loses the competition.
- **The identifiability threshold is not.** Below it a direction is flat: there is
  no second mode, only absent curvature, which then grows continuously. It may
  look abrupt because information appears suddenly with the design — in this
  system the mechanism is the number of transitions observed within an interval —
  but structurally it is a **crossover**, and calling it a phase transition in
  writing would be wrong.

### 8ter.2 A small theorem

If $f = \log L$ (one dataset) has two local maxima $\theta_1, \theta_2$, then along
any continuous path joining them the observed information $J = -\nabla^2 f$ fails
to be positive semidefinite somewhere. Proof in one line: along the path $f$
decreases then increases, so at some point the second derivative in the path
direction is positive, i.e. $v^{\mathsf T} J v < 0$ for $v$ the path direction.
(The mountain-pass theorem gives more — an actual saddle — but is not needed.)

This is about the **observed** information for a single dataset, never the
expectation, which is the Fisher information and is PSD under correct
specification. That is why the average sees nothing.

Three caveats on using it:

1. **The converse fails.** Non-PD in posterior draws does not imply two maxima; a
   curved ridge with a single mode produces it too. Necessary, not sufficient.
2. **The bad point may not be in the typical set.** If the typical set is a level
   set above the saddle, it has two disconnected components, one per mode, and the
   valley is never sampled: both modes are visited and everything seen is PD.
3. **Tempering hides it further.** With replica exchange, mode switches happen by
   temperature swap rather than by crossing the barrier at $\beta = 1$. So read
   the non-PD fraction **as a function of $\beta$**: zero at $\beta = 1$ and
   non-zero at intermediate $\beta$ is the signature of modes connected only
   through high temperature — which is itself a measure of how necessary the
   ladder was. It reads alongside $\mathrm{Var}_\beta[\log L]$ on the same axis:
   one says where the information arrives, the other where the landscape stops
   having a single mode.

### 8ter.3 Cheap detection, and why the two-point gradient test fails

The natural cheap test — is the gradient monotone along a segment — is a valid
certificate of non-concavity:
$(\nabla f(\theta_1) - \nabla f(\theta_2))\cdot(\theta_1-\theta_2) > 0$ implies
$f$ is not concave somewhere on the segment, and it costs no extra gradient
evaluations if consecutive sampler states are reused.

**But it is systematically blind to the case of interest.** The quantity
telescopes:
$d\cdot g = \int_0^1 d^{\mathsf T}\nabla^2 f\, dt = d\cdot\nabla f(\theta_1) -
d\cdot\nabla f(\theta_2)$, so it depends only on the endpoint gradients. At two
exact maxima both gradients vanish and the test returns exactly zero. In 1-D with
$f = -(x^2-1)^2$, the pair $x = \pm 1$ gives $d\cdot g = 0$ while the curvature at
$x = 0$ has the wrong sign.

**The right cheap test needs three collinear points and no gradients at all.**
With positions $t_1 < t_2 < t_3$ and values $f_1, f_2, f_3$, the sign of
$\frac{f_3-f_2}{t_3-t_2} - \frac{f_2-f_1}{t_2-t_1}$ is the sign of the average
directional curvature. It looks at the middle, which is where the information is.

Better still for basin membership: **evaluate the log density along the chord** and
look for a dip below the lower endpoint. That detects a barrier and, more useful,
**measures its height** — which is what decides whether tempering was needed and at
what temperature.

### 8ter.4 The ensemble sampler already computes this and throws it away

The Goodman–Weare stretch move proposes $X_k' = X_j + Z(X_k - X_j)$, so **every
proposal is a point on the chord between two walkers, and its log density is
already evaluated** to accept or reject. Two diagnostics are therefore free:

- **Acceptance resolved by $Z$ and by pair.** Small $Z$ returns a near-identity
  move and should almost always be accepted; mid-range $Z$ lands between the
  walkers. Acceptance high at small $Z$ and low at mid $Z$, for a given pair, is a
  barrier.
- **The chord profile.** Recording $\log p$ of the proposal together with $Z$ gives
  $\log p$ as a function of position along the chord, accumulated over thousands of
  proposals, and the depth of the dip is the barrier height.

Status of prior art, checked: the **coarse** version is established folklore in
the emcee community — a low acceptance fraction is read as evidence of
multimodality with wide low-probability valleys, and the standard advice is to
switch to tempering. The **resolved** version (by $Z$, by pair, the chord profile,
the barrier height) was not found. And the general statistics diagnostic for
multimodality is R-hat with overdispersed starts, which detects disagreement
between chains but gives neither the location nor the height of the barrier.

Why it is not done, in order of weight: rejected proposals are discarded by design,
so a large amount of already-evaluated density is thrown away; the diagnostic
culture is convergence-oriented rather than geometry-oriented; the coarse version
already exists so the perceived marginal value is low; and it is not obviously
actionable, since the answer to multimodality is tempering, which one would adopt
anyway.

**It is actionable here specifically**, because the ladder already exists and the
barrier height tells you whether it was well placed. Note also the cross-field
observation, of the same shape as the White/Gelman split of §5: measuring a density
profile along a coordinate to obtain barrier heights is routine in computational
chemistry (umbrella sampling, metadynamics, the string method), where it is the
objective rather than a diagnostic, and it has not crossed into statistics.

**Measured in this project**: at $\beta = 1$ the stretch acceptance was very low
and mixing happened only through temperature exchange.

### 8ter.5 Sources of multimodality, and how to tell them apart

Three distinct sources, requiring opposite responses:

1. **Label symmetry.** The number of equivalent maxima is the order of the
   automorphism group of the kinetic scheme, restricted to permutations preserving
   the conductance classes. This is computable exactly and is typically far smaller
   than the state count: a cyclic $n$-mer gives $n$ or $2n$, not $2^n$. **Do not
   confuse the size of the state space with the mode multiplicity.** The fix is a
   canonical ordering, which removes the modes by construction; it is not a
   sampling problem.
2. **Aggregation non-identifiability.** Genuinely distinct $Q$ matrices producing
   the same observable aggregated process, not related by any relabelling
   (Fredkin–Rice, Kienker; the bound is in `docs/bibliography/identifiability/`).
   **Not fixable by ordering.** What breaks it is the experimental design;
   non-stationary protocols do, which is what this project's protocol already
   exploits.
3. **Sloppy ridges.** Long curved near-flat valleys, not discrete modes. These
   defeat the stretch move for a different reason: the proposal is along the
   straight line between two walkers, and a *curved* ridge is left by that line.
   Affine invariance protects against linear correlation, not against a banana.

**They are distinguished by peak heights**: symmetry copies are exactly equal;
aggregation-degenerate solutions are also exactly equal but unrelated by
permutation; a ridge gives similar but unequal values.

**The cheap mapping step**: multi-start Gauss–Newton, which already works in this
codebase, launched from a few hundred dispersed points. It returns how many
distinct optima there are, their heights, and whether the count matches the
theoretical $|\mathrm{Aut}|$. Any excess over the symmetry prediction is
aggregation degeneracy or ridge. Choosing a sampler before this is done is blind.

### 8ter.6 The item that bears on the published results

**If two compared schemes have different mode multiplicities, the Bayes factor is
biased by the ratio of those multiplicities.** A scheme with six equivalent
solutions carries six times the posterior volume of one with a single solution,
all else equal, and the evidence rewards it for that. This is the same correction
that requires dividing by $k!$ in mixture models, where it is standard practice.

Concrete question for the Communications Biology results: **did the compared
schemes have the same symmetry order?** If yes, the factors cancel and there is
nothing to do. If not, the Bayes factor is off by a known integer ratio and is
**correctable after the fact without re-running anything**, since it only requires
counting automorphisms of each scheme. This belongs on the same list as the
`gvar_i` question of §10, and it is independent of it.

### 8ter.7 Method options, if a different sampler becomes necessary

Only relevant once §8ter.5 has established which of the three problems is present.

- **HMC alone does not serve**: it gives no $Z$ and does not tunnel between modes.
  If the diagnosis is ridge rather than modes, HMC with a Fisher-based mass matrix
  is probably the best move, keeping thermodynamic integration for $Z$.
- **Nested sampling** was designed for this request: the evidence is the primary
  output and multimodality is handled natively by clustering live points. It
  scales poorly with dimension; PolyChord tolerates more than MultiNest. Verify the
  dimension limit against the parameter count before investing.
- **Sequential Monte Carlo with tempering** is the least disruptive: it reuses the
  existing temperature structure, returns $Z$ as the product of normalizing-constant
  increments with a variance estimate, admits any kernel within each temperature
  including HMC, and parallelizes almost perfectly.
- **Levenberg–Marquardt / Gauss–Newton is an optimizer, not a sampler.** It was
  slower per step because it needs the Jacobian; it does not compete with MCMC. Its
  use here is the multi-start mapping of §8ter.5, where per-step cost is irrelevant
  because the starts are few and independent.

---

## 9. Novelty ledger

---

## 6. Novelty ledger

Read this before claiming anything.

### Classical, and not ours

The KL identities (§1.1), the variational decomposition (§1.2), thermodynamic
integration and path sampling (§1.3, Gelman & Meng 1998), Chib's identity,
the chain rule, Lindley information and its identification with mutual
information (Lindley 1956), the Laplace and BIC expansions, the shrinkage
formulas of §3.3 (textbook), $p_D$ and DIC (Spiegelhalter et al.), Bartlett's
paradox, Jeffreys–Lindley, robust Bayes (Berger), Bayes factor bounds, and
reverse-Bayes (Good 1950). None of this is a contribution. The derivations in
§§1 and 3 are recorded here for use, not for credit.

### The nearest prior art for a correction

Generalized Bayesian inference with a **learning rate**: SafeBayes
(Grünwald & van Ommen 2017 and following). Under misspecification the posterior
concentrates too fast; the repair is to temper the likelihood by an exponent
$\eta$ and to learn $\eta$ from the data. This is the closest existing tool to
"correct the evidence when the likelihood is not trustworthy", and it must be
cited.

### Where the difference is, and it is measurable

1. **The learning rate is a scalar; the distortion here is a matrix.** SafeBayes
   tempers the whole likelihood uniformly. The measured distortion in this
   project is anisotropic: in the worst cell of the design plane, $k_{\text{off}}$
   under-reports its standard error by a factor 1.62 while $N_{\text{ch}}$
   over-reports by 0.65 **in the same cell**. A scalar $\eta$ repairs one
   direction and worsens the other. The anisotropy is not a technicality, it is
   the reason the existing repair does not apply.
2. **The origin of the error is different.** SafeBayes addresses a model that is
   wrong relative to nature. Here the model may be right and the *computation* of
   the likelihood is approximate, through a Gaussian closure over the state space.
   That is approximate-likelihood territory (synthetic likelihood, pseudo-marginal
   methods) rather than classical misspecification, and the intersection of the
   two literatures is thin.
3. **The distortion is measurable per interval, and has been measured.** For this
   class of recursive algorithms nobody had done it. That is empirical, and
   empirical results are results.
4. **The propagation to the evidence is computed rather than assumed**, via the
   correction in the Posterior Information Distortion supplement.

### What is a unification and should be called one

That Lindley information, the $p_D$ of DIC, the sandwich correction and
robustness to the prior are four readings of the spectrum of a single matrix
$\mathbf{G}$ (§3.4) is a useful reorganization of known material. It should be
presented as "here it is all together and computable", not as discovery. A method
paper is allowed to do that, provided it says so.

### And the library

The stated reason nobody uses the sophisticated estimators is that their validity
is uncharacterized and there is no usable implementation. If that argument is
accepted for MacroR, it applies here too: a validated implementation is a
contribution in this field, not a consolation prize.

---

## 7. What to check in the code before any of this is used

1. **Which quantity is actually implemented** where the evidence appears: is the
   subtrahend $\log Z$ or $\mathbb{E}_{\text{prior}}[\log L]$? Per §1.1 these
   differ by $\mathrm{KL}_{-}$, and the error would be silent: positive, of the
   right magnitude, and growing with sample size like the correct quantity. This
   is the same failure profile as the bugs that have cost the most in this
   project.
2. **Is the prior normalized** in the code path that computes the evidence. §2.3.
3. **$\mathrm{Var}_\beta[\log L]$ from the existing tempering runs.** The samples
   already exist. This yields $\mathrm{KL}_{+}$, $\mathrm{KL}_{-}$, their ratio
   as a "where does learning happen" diagnostic, and a check on the temperature
   ladder, at no additional compute.
4. **The published Bayes factors.** It remains open whether correcting the
   `gvar_i` defect changes the Bayes factors reported in the Communications
   Biology paper. The right way to answer is a robustness map in the sense of
   §4.4: not a single corrected number, but the region of prior space over which
   the ordering of the models is preserved.

---

## 8. Open

- An anisotropic analogue of the learning rate. If a scalar $\eta$ is
  insufficient, what is the matrix-valued object that plays its role, and does it
  admit a variational justification comparable to the one SafeBayes has?
- Whether the evidence correction can be stated without the Gaussian layer of §3.
  Sections 1 and 2 are exact; §3 is not, and the correction currently lives in
  §3.
- A validity criterion computable from a **single** sample rather than from
  replicates, which is the standing gap in the diagnostic programme.

---

## References to obtain

Not yet in `docs/bibliography/`. Fetch before citing.

- Lindley (1956), On a measure of the information provided by an experiment,
  *Ann. Math. Statist.* 27(4):986–1005.
- Good (1950), *Probability and the Weighing of Evidence* (reverse-Bayes).
- Bartlett (1957), A comment on D. V. Lindley's statistical paradox.
- Berger (1994), An overview of robust Bayesian analysis, *TEST*.
- Berger & Sellke (1987), *JASA*.
- Gelman & Meng (1998), Simulating normalizing constants: from importance
  sampling to bridge sampling to path sampling, *Statist. Sci.*
- Spiegelhalter et al. (2002), DIC, *JRSS-B*.
- Grünwald & van Ommen (2017), Inconsistency of Bayesian inference for
  misspecified linear models, and a proposal for repairing it, *Bayesian Anal.*
  12(4).
- Vehtari, Gelman & Gabry (2017), Practical Bayesian model evaluation using
  leave-one-out cross-validation and WAIC, *Stat. Comput.*
- "The Bayes factor reversal paradox" (2025), arXiv:2511.22152.
