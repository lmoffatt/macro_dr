# CCO vs COC: what can be told apart, and with what

Status 2026-09-01. Everything marked VERIFIED was checked numerically; the
scripts live in `tmp/` (`check_cco_coc_equivalence.py`,
`twin_coc_from_cco.py`, `cco_coc_design.py`, `cco_coc_freefit.py`,
`cco_coc_bounded.py`). Trigger: the MacroIRNext audios of 2026-09-01.

Topologies, as they stand in `legacy/models_simple.h`:

- `scheme_CCO`: C1 --kon·x--> C2, C2 --koff--> C1, C2 --gating_on--> O,
  O --gating_off--> C2. Open = state 2.
- `scheme_COC`: C1 --on·x--> O, O --off--> C1, O --inactivating_on--> C2,
  C2 --inactivating_off--> O. Open = state 1.

## 1. The twin map (closed-form, exact), VERIFIED

Given a CCO with k1 = kon·x, km1 = koff, k2 = gating_on, km2 = gating_off,
the equivalent COC **at that concentration x** comes out with no solver:

    θ₁, θ₂ = roots of  s² − (k1 + km1 + k2)·s + k1·k2 = 0
    a₁ = (θ₂ − k2)/(θ₂ − θ₁)          (mixture weight on θ₁)

    on·x = θ₁          inactivating_off = θ₂
    off  = a₁·km2      inactivating_on  = km2 − off

Both root assignments (θ₁ or θ₂ to the agonist branch) give valid twins;
they only swap the branch labels.

For `scheme_CCO_par.csv` (kon 6.73, koff 166, gating_on 743, gating_off
45.3) at 10 µM, the twin is in `tmp/scheme_COC_twin_par.csv`:
on 5.423005, off 9.347195, inactivating_on 35.952805,
inactivating_off 922.069949.

The inversion from observables (λ₁, λ₂, P_o, w) also works, but the
correct anchor is

    K_out = (1 − P_o)·(w·λ₁ + (1 − w)·λ₂)

The (1 − P_o) factor is mandatory: it follows from C(0) = i²P_o(1−P_o)
and C′(0) = −i²P_o·K_out. Without it the inversion falls off the feasible
region (negative rates and negative discriminant with tame observables).
Minimal control case: C⇌O has λ = α+β but K_out = α.

## 2. At equilibrium at one concentration: indistinguishable, always, VERIFIED

With a single open state the aggregated process is alternating renewal:
every closure leaves from O and every opening enters O, so the process
regenerates at each opening. What remains is iid open dwells exp(K_out)
and iid closed dwells with a two-exponential mixture density. Both
topologies generate the **same family** of those laws:

- The CCO mixture weights are proper iff θ₁ < k2 < θ₂, which is exactly
  km1 > 0: always true.
- 20000 random CCOs → all have a COC twin with positive rates, and vice
  versa. Closed dwell densities agree to 1e-15.
- P_O(t) trajectories of the twin pair started at equilibrium: agree to
  8.9e-16 over 2000 samples.

Hard consequence: **the joint law of the macroscopic current is
identical** for any N, any noise, any record length. No equilibrium
statistic separates them. It is not a limitation of MacroIR nor of noise
analysis; it is a theorem. The 4 degrees of freedom of the law
(λ₁, λ₂, P_o, w) saturate the 2·n_O·n_C = 2·1·2 = 4 bound from the
identifiability folder: the two topologies are minimal realizations of
the same thing.

This retro-explains 2007: it was not lack of power, it was impossible.

## 3. With an agonist jump: they separate, but through the MEAN, VERIFIED

The repo protocol (pre 0 µM / 10 µM / post 0 µM, 50 kHz) breaks the
degeneracy. The signature is structural and classical:

    CCO:  dP_O/dt(0⁺) = 0        (two steps from C1 to open)
    COC:  dP_O/dt(0⁺) = on·x > 0

That is, CCO has a sigmoidal delay on the rising front and COC does not.
Sweep: 24000/24000 random twin pairs differ in the transient amplitude,
always with the same sign. No twin reproduces the jump.

And against a **free** COC (not the twin), fitting its four rates to the
CCO mean trajectory over the whole protocol, in the N→∞ limit which is
the most favorable case for the mimic:

- The optimum runs to the degenerate limit (inactivating_on/off →
  0.1 / 1e4 in the bounded box, or 1e11 / 1e17 unbounded, both with the
  same residual): the best COC mimic **is not a three-state COC**, it
  collapses to an effective C1⇌O.
- Residual: RMS 7.46e-3, peak 3.32e-2 in P_O units (4.4% of the plateau),
  **concentrated in the first 2 ms** after the jump. Bounded and free
  fits agree, so the optimum is genuine.

## 3b. The EXACT recursive filter on data, VERIFIED, decisive

Luciano's objection (audio and chat 2026-09-01): being stochastically out
of equilibrium the system relaxes back, the kinetics show up in that
relaxation, and the recursive algorithm should detect it.

Tested directly, with no approximation at all: EXACT hidden chain over
occupancy compositions (n_C1, n_C2, n_O), 66 states for N=10 channels,
Gaussian emission, exact forward filter. That is, the likelihood MacroIR
approximates, computed without approximating. Script
`tmp/cco_coc_exact_filter.py`.

Equilibrium data at 10 µM simulated from CCO, scored with CCO and with
the twin. Large genuine spontaneous excursions (n_O sweeping 5→10 out of
10 channels, mean 8.25):

    record   excursion   logL_CCO     logL_COC     difference
       0      5 → 10     -1209.3206   -1209.3206    0.0
       1      6 → 10     -1190.4913   -1190.4913    4.5e-13
       ...
    prefixes:  n=500 5.7e-14 | n=2000 0.0 | n=8000 1.8e-12

Δlog-likelihood = 0 at machine precision, and **it does not grow with
record length**. The fluctuations are there, they are large, the exact
recursive filter follows them, and it distinguishes nothing.

Why the intuition fails here and where it holds. The mean relaxation from
a spontaneous fluctuation IS the autocorrelation function (Onsager
regression), which the twin matches by construction. In general that
would not be enough: matching the autocorrelation does not imply matching
the full law, and there the recursive filter would win, since it sees the
higher-order statistics of the relaxation trajectories, not just their
average. What kills the CCO/COC case is that with ONE open state the full
law is pinned by the dwell distributions, which were matched too. Nothing
of higher order is left to see.

**The intuition is correct and holds as soon as there are two open
states**, where the equivalence class does not saturate and the
autocorrelation no longer determines the law. That is the model pair to
choose if the goal is to show that the recursion adds something over
noise analysis. (See the correction in 3quater for the precise
mechanism.)

## 3bis. Following the POSTERIOR, not isolated points, VERIFIED

Refined objection: one must follow the posterior sequentially as MacroIR
does, not evaluate correlations at isolated times. That is exactly what
was done (`tmp/cco_coc_posterior_track.py`): predict (a ← a·T) / update
(a ← a·emission) recursion accumulating the normalizer, with the EXACT
posterior over the 66 composition states instead of the moment
approximation.

On the same equilibrium record, step by step:

    INTERNAL posterior, E[n_C1] (channels in C1, out of 10)
      t=0    CCO 0.598372   COC 0.685820   diff −0.087
      t=3    CCO 1.252771   COC 1.493904   diff −0.241
      t=203  CCO 1.665315   COC 1.838723   diff −0.173
      max |diff| over the record: 4.39e-01     ← NOT zero

    one-step prediction of the OBSERVABLE
      max |diff E[n_O]|                 7.11e-15
      max |diff Var[n_O]|               6.40e-14
      max |predictive density ratio − 1|    4.55e-15
      accumulated logL: −232.9895356481 in both, difference 0.000e+00

That is: the two models **genuinely disagree** about where the closed
channels sit (up to 0.44 channels out of 10), the filter faithfully
follows that different posterior, and yet **every prediction it makes is
identical**. The discrepancy lives entirely in the internal coordinate
the observable does not see, and cancels exactly in every prediction.

## 3ter. Why: the factorization at every order, VERIFIED

With ONE open state, conditioning on "open" at an intermediate time pins
THE state, so the Markov property factorizes the whole hierarchy:

    E[1_O(t₁)·1_O(t₂)·…·1_O(t_n)] = π_O · P_OO(t₂−t₁) · … · P_OO(t_n−t_{n−1})

Verified to 1e-16 on CCO, and CCO against the twin agree to **1.4e-15 at
orders 2, 3 and 4** over all lag combinations tried.

All multi-time correlations, at every order, are generated by the single
function P_OO(t), which is exactly what the autocorrelation measures.
Matching it matches every fluctuation statistic there is, including
anything MacroIR integrates over the sampling interval. There is no
higher-order residue. That is the deep reason for the zero, and it does
not depend on the point sampling of test 3b.

And it breaks where it should: with TWO open states the same
factorization fails (relative difference up to 6.2e-4 at third order,
verified on a C⇌O1⇌O2 model). That gap between third order and the
autocorrelation is precisely what a full-likelihood recursive filter can
exploit and a spectrum cannot.

## 3c. Information per record, and the protocol that maximizes it

Same exact filter, jump protocol, truth = CCO, competitor = twin. KL per
record (2200 samples in every case), `tmp/cco_coc_protocol_scan.py`:

    protocol                            N=10              N=20
    one long application (40 ms)     0.179 ± 0.026     0.370 ± 0.058
    train of 11 pulses 2 ms on/off   8.392 ± 0.321    13.127 ± 0.605

**The pulse train yields 47× more with exactly the same data.** Direct
consequence of the discrimination living in the first ~2 ms of the rising
front: one long application buys a single front and then pays 38 ms of
plateau that contributes zero. Pulse.

CORRECTION to section 4: the independent-sample Gaussian bound
**overestimates by ~100×** (it predicted ~15 nats at N=10 against the
best COC; the real value against the twin, which is more favorable, is
0.18). Use the numbers of this section, not section 4's. The scaling with
N also turns out sublinear (2× from N=10 to N=20), not ∝N.

With 8.4 nats per record at N=10, a decisive Bayes factor (log BF ≈ 5)
comes from **a single** 44 ms pulsed record. The experiment is cheap.

## 3quater. THE PAIR THAT WORKS: two open states, VERIFIED, this is the figure

Built and measured (`tmp/two_open_demo.py`). This is the experiment that
demonstrates what CCO/COC cannot.

    Model A (TWO open states, same conductance):  C ⇌ O1 ⇌ O2
      a = 385.530 (C→O1)   b = 2234.368 (O1→C)
      c =  49.786 (O1→O2)  d =  103.159 (O2→O1)

    Model B (ONE open state):                     C1 ⇌ C2 ⇌ O
      k1 =   52.985  km1 =   99.959
      k2 = 1112.851  km2 = 1507.048

B is built by inverting A's observables with the closed-form map of
section 1. All at equilibrium at fixed concentration: **no jumps, no
protocol**, spontaneous fluctuations only.

What noise analysis sees, i.e. everything it measures:

              λ1            λ2          P_o          w          variance
      A  108.64825523  2664.19500182  0.20370646  0.30193812  0.16221014
      B  108.64825523  2664.19500182  0.20370646  0.30193812  0.16221014
      max relative difference: 2.0e-15

And the full autocorrelation, lag by lag (0.05 to 10 ms): identical to
3e-15. Mean, variance and spectrum are the same object. No noise
analysis, however good and with however much data, can separate them.

Third order, where they part:

      t1[ms] t2[ms]        A                B          rel. diff.
       0.460  0.460   0.087037990096  0.072244973272    17.0%
       1.841  1.841   0.055864147997  0.033359629127    40.3%
       4.602  4.602   0.037465552910  0.024888218808    33.6%

Exact recursive filter, N=10 channels, equilibrium data:

      samples   Δlog-likelihood       per 1000 samples
          250      0.598 ± 0.256            2.39
         1000      3.348 ± 0.684            3.35
         2000      7.776 ± 0.854            3.89
         4000     18.189 ± 1.447            4.55
         6000     27.849 ± 2.227            4.64
      reverse direction (truth = B): 26.912 ± 1.977 at 6000, 4.49 per 1000

**Δlog-likelihood grows linearly with the record**, ~4.6 nats per 1000
samples (20 ms at 50 kHz), symmetric in both directions. With N=10
channels a decisive Bayes factor (log BF ≈ 5) comes from ~22 ms of
recording.

### CORRECTION to the mechanism (Luciano's objection, 2026-09-01)

Objection: open and closed are interchangeable (subtract from the maximum
current), so "2 open + 1 closed" should not differ from "2 closed +
1 open". It is correct, and my original framing was WRONG on two counts.

ERROR 1. I said A is not a renewal process. **It is renewal.** A enters
the open class always through O1 (C→O1 is the only way in) and the closed
class always at C (single state): deterministic entries, iid dwells. B
too. BOTH are alternating renewal, and the law of each is exactly its
pair of dwell densities (f_conducting, f_non_conducting).

ERROR 2. The contrast is not "two open states break the factorization" as
a state count. Verified (`tmp/two_open_whats_really_going_on.py`):

    mean dwell times (ms)      conducting   non-conducting     P_o
      A                          0.663550      2.593832     0.203706
      B                          0.663549      2.593836     0.203706

    are they exponential? (ratio against the equal-mean exponential)
       t[ms]   A conduct.   B conduct.   A non-cond.   B non-cond.
        0.00     1.48146      1.00000      1.00000      2.88417
        3.00     0.25134      1.00000      1.00000      0.26783
        6.00     7.38548      1.00000      1.00000      0.10445

**A has NON-exponential conducting dwells and exponential non-conducting
ones; B exactly the other way around. With the same means on both
sides.** The multi-exponential structure simply MOVED from one side of
the cycle to the other.

The true mechanism, then: in a 3-state model the autocovariance has 2
exponentials, i.e. exactly 4 numbers (λ1, λ2, P_o, w), and the model has
4 rate constants. With ONE conducting state those 4 numbers determine the
whole law (the factorization of 3ter, which forces exponential conducting
dwells). With TWO they do not. So B is the unique one-conductor model
with those 4 numbers, and A is a different model with the same 4. They
differ in the full law and agree on everything second-order.

**Why Luciano's symmetry does not collapse the pair.** Inverting the
current sends A to a model A' with one conductor and two non-conductors,
structurally like B. But A' is NOT B: A' has P_o = 0.796294 and B has
0.203706. The inversion maps the discrimination problem to an equivalent
problem; it does not collapse one model onto the other. The question that
survives the relabeling is **whether the extra kinetic component sits in
the open times or in the closed times**, and noise analysis cannot answer
it.

That is the paper's question, and it is more classical and sharper than
"one or two open conformations". The numbers above (identity to 2e-15 in
everything second-order, Δlog-likelihood growing linearly at 4.6
nats/1000 samples) do not change: they are measured and did not depend on
the wrong framing.

OPEN: the scaling with N_ch is unmeasured; sweep it. And it must be
confirmed that MacroIR (the moment approximation) recovers this
separation, which is distinct from the separation existing: the exact
filter sees it, how much the approximation keeps is unknown. That is the
paper's second result.

## 3quinquies. The COMPLETE classification of 3-state chains

The 3quater pair is **CCO against COO**: the same linear chain, with the
extra state on the closed side or on the open side. Named that way,
current inversion (y → N·i − y, a bijection on the data, so it maps
equivalences to equivalences) sorts the four topologies into two classes:

    CCO   C1 − C2 − O     one conductor, at the end
    COC   C1 − O  − C2    one conductor, in the middle
    COO   C  − O1 − O2    two conductors   = mirror of CCO
    OCO   O1 − C  − O2    two conductors   = mirror of COC

From CCO ≡ COC (verified, section 2) the mirror PREDICTS COO ≡ OCO with
no new algebra. Verified (`tmp/three_state_classification.py`), building
the OCO twin of 3quater's COO by mirroring the CCO→COC map:

    COO  a=385.530000 b=2234.368000 c=49.786000 d=103.159000
    OCO  O1→C=100.806706  C→O1=9.196541  C→O2=376.333459  O2→C=2286.506294

    second order: λ1, λ2, P_o, w agree to 1.4e-14
    THIRD order:  agree to 1.9e-14  (where CCO vs COO parted by 17-40%)
    exact filter, N=10: Δlog-lik. −5.3e-13 ± 2.4e-13 at 6000 samples

**Result: two equivalence classes, and nothing else.**

    class I  (one conductor):  {CCO, COC}   exponential conducting dwells
    class II (two conductors): {COO, OCO}   non-exponential conducting dwells

Within each class, exactly indistinguishable, forever, with any amount of
equilibrium data. Between classes, 4.6 nats per 1000 samples at N=10
channels.

So the only thing macroscopic equilibrium data can determine about a
3-state chain is **on which side of the cycle the extra kinetic component
sits**, and that they determine well. The position of the extra state
within its side (end or middle) is inaccessible in principle, not for
lack of data.

That is the paper: a classification theorem with the two classes
exhibited explicitly by closed-form maps, the interior of each class
proven inaccessible, and the boundary between classes proven accessible
and measured with the full likelihood exactly where noise analysis is
blind.

## 4. Where the fluctuation information lives, upper bound, see 3c

Decomposing the expected separation into mean part and variance part
against the best COC:

    N_ch      mean       variance    variance share
       10     14.7        29.2         66.6%
       30     44.0        29.2         39.9%
      100    146.6        29.2         16.6%
     1000   1466.3        29.2          1.9%

The structure is clean and worth more than the numbers: the mean
information scales **∝ N** (the mean difference is N·i·ΔP_O against
binomial variance N·i²·P(1−P)), while the fluctuation information is
**independent of N** (both variances scale with N, their ratio does not).
Hence the crossover near N ≈ 20.

Translation: at large N the model is decided by the mean and no MacroIR
or anything like it is needed, a deterministic fit suffices. The window
where the autocorrelation decides is the **few-channel** one, which is
exactly the N·S̃ = N/v² axis of the macro/micro boundary in the eLife
paper.

CAVEAT: these numbers are an upper bound. They come from a Gaussian
independent-sample approximation, and the correlation time here is
~16 ms ≈ 800 samples, so both columns are inflated. The mean also has its
residual concentrated in 2 ms, shorter than the correlation time, so it
inflates differently from the variance. **The real ratio must be measured
by running the likelihood, not estimated.** Second caveat: with realistic
current noise (here 1e-3, negligible) the variance term degrades, since
v_B → noise² and the variance ratio → 1.

## 5. The two experiments

**A. Null test / evidence calibration.** Equilibrium data at 10 µM, CCO
against its exact twin. The log-likelihood ratio is identically 0 for
every record, so:

- Whatever MacroIR reports there is **approximation error**, measurable
  against known truth.
- Δlog-evidence must stay O(1) and **not grow with record length**. An
  estimator that "detects" the true model is miscalibrated.
- It serves to compare adaptive-beta MCI, Levenberg-Marquardt and
  whatever comes next, which is the evidence validation needed before
  P2X2.

Design TRAP: the current protocol's agonist segment (2000 samples =
40 ms) **does not equilibrate**: τ_slow ≈ 16 ms, and it reaches 91.5% of
the equilibrium P_o (0.755 against 0.8255). For the null test start at
equilibrium and run ≥10τ, or the entry transient contaminates the zero.

**B. Discrimination, and whether the autocorrelation adds anything.**
**Pulsed** protocol (see 3c: 11 pulses of 2 ms yield 47× one long
application), both models **free** (the twin is not the competitor here,
it is the null-test construct). The quantitative question is not whether
they separate, they do; it is how much the fluctuations add over the
mean. Design: same data, two likelihoods, (a) full MacroIR and (b) the
mean-only limit, sweeping N_ch. The gap is the fluctuation contribution.

## 6. What the win is NOT

The sigmoidal delay of the rising front is a classical, well-known
discriminator. A paper saying "with jumps I separate CCO from COC" adds
nothing. What is defensible, and follows from the above:

1. An impossibility theorem at equilibrium, with the explicit twin map as
   a constructive tool (sections 1 and 2).
2. Bayesian evidence correctly reporting that impossibility: a
   calibration test with known truth (experiment A).
3. The discrimination extending to the few-channel regime, where the
   rising front is buried in noise and classical analysis would need
   ensemble averaging, while the full likelihood extracts it from a
   single record (experiment B).

The paper's axis is (3), and its charm is that the regime is quantified
by the same N·S̃ that the eLife paper already characterized.

## 7. Open

- Actually run A and B; the section 4 numbers are an upper bound.
- Matching CCO/COC at **several** concentrations at once is impossible
  under mass action (the twin's rates depend on x nonlinearly). So even
  classical noise analysis at several [A] breaks the degeneracy. What
  remains is to quantify how much data that route needs, to avoid
  overselling the Bayesian advantage.
- All of this is for ONE open state and two closed ones. In richer models
  the equivalence classes get richer. Caution from the 3quater
  correction: two open states alone do not create adjacent-dwell
  correlations (COO is still renewal, since both class entries are
  deterministic); correlations between adjacent dwells require
  multi-state classes with multiple entry states on both sides (e.g.
  2 open + 2 closed with multiple connections). That is where the
  "algebra of distinguishable models" from the audio has genuinely
  nontrivial content.
