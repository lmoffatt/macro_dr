# MacroIR-Bessel: implementation plan

2026-08-31. Written after a design pass in which three candidate constructions were developed independently and each was then attacked by an adversarial reviewer. This note explains the winning construction and the plan around it, and is meant to be readable on its own. Nothing is implemented yet. Companion voice notes: program/source-notes/audios/audios/MacroIRNext/ (2026-08-31).

> ADDENDUM 2026-08-31 (later, same day): phase 0 implemented and gated.
> `legacy/acquisition_filter.h` (pole tables by recurrence + root finder, no
> transcribed constants; complex dd1/dd2; real-pair dd1; noise closed forms) +
> `tests/math/test_acquisition_filter.cpp` (Catch2; includes the appendix
> bound-table row 61.5%/0.641 as a gate) + `theory/macroir/notes/
> bessel_reference/` (oracle: exact filtered-truth simulator + the augmented
> recursion for native/grouped/box reads from one generic code path; block
> equations in `augmented_recursion_blocks.md`). Run `python3 bessel_oracle.py`
> for the current verdict. IRT is NOT inherited by this member (decided in the
> same session: only the Taylor branch consumes pair-resolved seconds). The
> HEAD command surfaces and the two Qdt flavors are mapped in the session of
> this date. LATER SAME DAY (restriction lifted by Luciano): the mathematical
> body of calc_Qdtf_eig and of the recursion sibling is implemented in
> `legacy/qdtf_engine.h` (self-contained, gated by
> `tests/math/test_qdtf_engine.cpp`, cross-checked digit-exact against the
> oracle). THIRD PASS (same day, "dale adelante"): the member is wired into
> macro_dr — `legacy/qdtf_member.h` (the likelihood loop over a recording:
> filter warm-up, per-sub-step composition through the carried state, mixed
> native/grouped reads, box config as anchor), `include/macrodr/cmd/
> qdtf_likelihood.h`, the DSL command `calc_qdtf_likelihood(model,
> parameters, experiment, data, n_poles, cutoff)`, and the parity test
> `tests/macroir/test_qdtf_member.cpp` (box vs the av=2 member). All TUs pass
> g++ -fsyntax-only with forced instantiation; the parity RUN needs the full
> build. Not yet wired: derivative mode (GN fitting), the simulator twin
> (M1), and the MacroR2-variant integration for the generic MLE/evidence
> commands (the direct entry point covers evaluation and samplers that need
> only plain logL).

## 1. The problem

The patch-clamp amplifier low-pass filters the current before the digitizer samples it. In the Moffatt & Hume 2007 recordings (Methods of the JGP paper, verified in the PDF) the filter is a 4-pole Bessel with cutoff f_c = 10 kHz, and the sampling rate is fs = 50 kHz.

A low-pass filter has memory. What the digitizer records at time t is a weighted average of the recent past of the current i(t):

    y(t) = ∫ h(s) · i(t − s) ds  +  filtered instrument noise,

where h(s) is the filter's impulse response. For a 4-pole Bessel at 10 kHz, h is appreciable for about 100 µs, which spans five sampling intervals at 50 kHz.

MacroIR instead models each recorded value as the plain average of the current over its own interval: a "box" weighting, confined to that one interval. That is a good model when the interval is long compared to the filter's memory and a poor one when it is short. The dimensionless group that controls the error is f_c · Δ, where Δ is the interval length. The manuscript's own bound table (08_appendix_derivation.tex:1010-1020) gives, at f_c·Δ = 0.2, a variance deficit of 61.5% and a correlation of 64.1% between neighbouring recorded values, both of which the box model ignores.

On the real 2007 record this is the dominant regime exactly where the information is. The record has 1078 intervals with log-spaced binning: the single-sample intervals around the ATP pulses are 20 µs long, so they sit at f_c·Δ = 0.2, the worst column of the table. Roughly a quarter to a third of the intervals sit below f_c·Δ ≈ 0.5 and they are the ones that carry the kinetics; 742 of 1078 sit above f_c·Δ = 10, where the filter has fully settled within the interval and the box model is fine. The 187 µs test pulses are only about 2.7 times the filter's fixed smear of 0.69/f_c = 69 µs. The remedy the manuscript itself offers (block-average the data until f_c·Δ = 5) would decimate the record 25-fold and discard the interval containing the concentration jump, which is precisely the data the ultrashort-pulse experiment exists to produce. A Bessel-aware likelihood buys that interval back.

Two separate things go wrong when the weighting spills across interval boundaries:

1. A recorded value depends on what the channels did during earlier intervals, and the current recursion never conditions on that.
2. Consecutive recorded values share signal history and share filtered noise, so they are correlated even once the occupancies at the interval boundaries are known. The recursion treats them as conditionally independent.

The theory is already half prepared for the fix. The derivation note (theory/macroir/docs/Macro_IR/macroir_derivation.tex:148) declares the per-channel boundary-conditioned mean and variance of the recorded quantity to be pluggable inputs ("e.g. via a Bessel filter or similar"), and every recursion equation downstream of them is agnostic to the weighting. What no plug can fix on its own is the two breakages above, and the manuscript's appendix names the obstruction explicitly: in the augmented formulation, the equivalence holds because an accumulator is reset to zero at each window (08_appendix_derivation.tex:954-956), and a filter that remembers across windows forbids that reset.

## 2. The construction: carry the amplifier's state inside the recursion

A k-pole filter is itself a small physical system. It has k internal variables (think of the voltages on the capacitors inside the filter), collected in a vector φ, which evolve linearly, driven by whatever current enters:

    dφ/dt = A_f φ + b_f · i(t),        recorded output = c_fᵀ φ,

with A_f, b_f, c_f fixed known matrices determined by the pole positions and the cutoff. For the 4-pole Bessel, k = 4.

The channel occupancies n(t) (the vector counting how many of the N channels sit in each of the K states) form a Markov process with generator Q and conductance vector γ, so the current is i(t) = nᵀγ + baseline + ξ(t), with ξ white instrument noise of power spectral density S₀. Since φ is driven by i, the joint pair (n, φ) is Markov: its present value determines the statistics of its whole future. Everything the filter remembers about the past current is stored in φ.

Now recall what MacroIR's central trick was: condition each channel on its occupancy at both ends of the interval; then the interval average becomes a function of states, and the chain of intervals stays Markov. The Bessel kernel broke this because the recorded value remembers more about the past than the endpoint occupancies do. The repair is to enlarge what the recursion tracks so that the extra memory is part of the tracked state. MacroR currently carries a mean vector μ and covariance matrix Σ_nn for the occupancies. The Bessel member carries, in addition, the filter mean m = E[φ], the cross-covariance Σ_nφ (how uncertainty about the occupancies is tied to uncertainty about the filter's internal state), and the filter covariance Σ_φφ. Between intervals nothing is reset. The obstruction of section 1 disappears by construction: the quantity that could not be reset is now simply carried.

The observation step then becomes almost trivial, with one case split:

- Native-rate interval (one raw sample, the 20 µs intervals): the recorded value is the filter output at the sampling instant, a linear readout c_fᵀφ(Δ). The Bayes update is ordinary Gaussian conditioning of the joint (n, φ) moments on that scalar.
- Binned interval (the log-spaced averages of the 2007 record): the recorded value is the average of the filter output over the interval. Add one more variable, a cumulator u with du/dt = c_fᵀφ, reset u to zero at the start of each recorded interval (this reset is legitimate: the recorder's averaging genuinely restarts), and read z = u(Δ)/Δ.

The member switches between the two reads per interval based on n_samples. Applying the cumulator read to a native-rate sample would be wrong; it would impose an extra averaging the hardware never performed.

This construction is one Luciano had already derived, twice, before this plan (found 2026-08-31):

- In writing, in the first draft of the eLife manuscript (macro_dr-46a60dd/docs/eLife 2025/version inicial/elife-macroir.tex:71, Scope): "Incorporating explicit analog low-pass filtering (e.g. Bessel) is possible by augmenting the system with the filter state-space (Kalman/Bucy style) or by convolving the spectral kernels with the filter poles; both routes are substantially more expensive computationally." Those two routes are precisely Designs A and B below; the review process of this plan picked A and killed B.
- In full, in the Gemini conversation of 2025-11-11, archived at ~/Projects/gemini/archive/threads/2025/2025-11-11_1218_what-is-the-final-property-of-the-basement-filte.md ("basement filter" is the speech-to-text mangling of "Bessel filter"; follow-up on colored noise in the same folder, ..._1549_...). The thread contains the complete augmentation skeleton, in Kalman language: the 4-pole Bessel as a linear system dx_f/dt = A_f·x_f + B_f·u built from the reverse Bessel polynomial θ₄(s) = s⁴ + 10s³ + 45s² + 105s + 105 (companion A_f, B_f = [0,0,0,105]ᵀ, readout C = [1,0,0,0]); input u = gᵀp (the channel current); the augmented generator Q′ = [[Q, 0], [B_f·gᵀ, A_f]] with observation row g′ᵀ = [0ᵀ, C]; discrete propagation e^(Q′Δt); and the resolution that killed the old impossibility belief: a serial LTI filter ADDS its states (K+4), while colored noise modeled as Markov states MULTIPLIES the state space (that had been the confusion with Kim's construction, a parallel process). The 2025-11-12 voice note (program/source-notes/audios/audios/Chat de WhatsApp con MacroIR 8/Gemmini transcripto 8.md:49) retells it. What the thread leaves open, and this plan supplies: the ensemble second-moment structure (multinomial boundary counts, Γ̄/V̄ interval moments), the binned-interval read (the cumulator), the colouring of the instrument noise by the same filter, and the endpoint-conditioned pole-shifted integrals. One algebraic bonus already visible in the thread's Q′: it is block-triangular, so spec(Q′) = spec(Q) ∪ {filter poles}, which is exactly why the new integrals are the existing divided differences evaluated at pole-shifted arguments.

Two useful things fall out of the same construction:

- The instrument-noise colouring is handled exactly. The filter that colours the signal also colours the white noise ξ, and its contribution to the covariance of φ over an interval has a closed form, Σ_ξ(Δ) = S₀ ∫₀^Δ e^{A s} b bᵀ e^{Aᵀ s} ds. This replaces the current white-noise term e = Current_Noise·fs/n_samples (legacy/qmodel.h:3740 and siblings), of which it is the generalization.
- The box member is recovered as an exact special case. Delete the φ rows and keep only the cumulator u (so A = 0, b = 1): the equations collapse, term by term, to today's MacroIR. This is the main regression anchor: the new code, run in that configuration, must reproduce the existing member to machine precision, and the existing member is validated against the frozen reference binary. (The limit f_c → ∞ also recovers MacroIR, but only approximately at finite f_c, so it serves as a secondary check.)

## 3. Why the new integrals need no new machinery

To advance one interval, the recursion needs expectations of the form: given that a channel started the interval in state a and ended it in state b, what are the mean and variance of its weighted conductance integral? For the box weighting, macro_dr computes exactly these in calc_Qdt via the eigendecomposition Q = VΛW: the answers come out as "divided difference" functions of the eigenvalues, the kernels named E2 (for first moments) and E3 (for second moments) in legacy/qmodel.h:770-866, evaluated at arguments λ_a·Δ, λ_b·Δ.

The Bessel impulse response is a sum of decaying exponentials, one term e^{−ν_p s} per pole ν_p (the poles come in complex-conjugate pairs, and everything recombines to real quantities at the end). The new integrals weight the conductance by these exponentials:

    W_p = ∫₀^Δ e^{−ν_p (Δ−s)} γ(X_s) ds,

and an exponential weight does nothing except shift arguments inside the very same divided differences: where the box case evaluates E2(λ_a·Δ, λ_b·Δ), the Bessel case evaluates E2(λ_a·Δ, (λ_b − ν_p)·Δ), and similarly for E3 with pairs of poles. No new class of integral appears; the same code shapes run with shifted inputs. Setting ν_p = 0 gives back today's formulas identically, which is the algebraic face of the "box as special case" statement above.

The count of tables grows, but less than a naive count suggests (refinement 2026-08-31, from Luciano's established box-member result that gvar_ij is never consumed: gvar_i plus gmean_ij suffice). The same structure holds per pole, with an identity that makes it cheap: end-marginalizing a pole-discounted leg collapses it to a scalar, e^((Q−νI)t)·1 = e^(−νt)·1. Consequences: (a) FIRST moments are needed pair-resolved, one E2-class table per mode (~5 vs 1 today) — they feed Σ_nφ, the macro variance part and the update vector; (b) SECOND moments are consumed only start-marginalized, and after the collapse each pole pair costs a matrix-vector chain, O(K²), not a full E3 table — ~15 vector-level integrals, not 15 tables; (c) the sub-interval composition monoid inherits the rule (the discount of already-accumulated reward through a later sub-interval is the state-independent scalar e^(−νΔ_B)), so it too needs pair-resolved firsts and marginal seconds only. The 10-40× cost estimate in section 8 is therefore an upper bookkeeping; ~5-15× is the likely range on short intervals. TO PROVE as an M0 lemma (chat derivation, not yet verified): the no-pair-resolved-seconds claim through every consumer, including whether the variance-correction branches (MacroIRT, qmodel.h:4680+) consume anything pair-resolved that the base member does not, if the Bessel member inherits them. An independent cross-check evaluator exists for all integrals either way: Carbonell, Jiménez & Pedroso 2008 generalize Van Loan's block-matrix-exponential method to nested integrals of exactly this shape, and the macroir rewrite already records the Van Loan identity for the box case as its decision D6.

## 4. What stays approximate

The adversarial review verified by independent re-derivation that the moment propagation itself is exact: every expectation the recursion needs is an affine or bilinear function of the moments already carried, so no information leaks that the construction silently drops. What remains approximate is exactly what is approximate in MacroIR today: (a) the Bayes update treats the joint distribution as Gaussian, and (b) between intervals only two moments are carried. Both closures now act on a larger object (the joint of occupancies and filter state), so their error is of the same kind but need not be of the same size; the plan measures it (milestone M2) rather than assuming it. The paper must say "extended closure, measured", never "exact likelihood".

## 5. The two alternatives, and why they lost

Design B ("generalized kernel"): keep the recorded value as a weighted path integral and handle the multi-interval memory by conditioning on the occupancies at the last m interval boundaries instead of 2. Killed by the review: two overlapping filter windows share the actual within-interval fluctuations of the path, and no stack of boundary occupancies records those, so the correlation between consecutive recorded values is systematically missed. The filter state φ is precisely the variable that does record them. B was also the most expensive of the three in schedule and had the shakiest cost accounting.

Design C ("FIR ladder"): chop time to the raw 20 µs grid, treat each raw sample as a small box-average observation, and write each recorded value as a discrete weighted sum (a finite impulse response, FIR) of the last L of them, updating only when a value was actually recorded. The review found it sound. It is kept in two roles: its first rungs need no C++ at all and become Phase 0 of this plan (a truth generator plus a measurement of how much damage the box member suffers on filtered data), and the member itself is the fallback if Design A's cost in derivative mode turns out unacceptable. It is deliberately never built as a second member alongside A, because it is exact only with respect to a discretized stand-in for the filter, and building both would duplicate the same physics.

## 6. Where to build it

Recommendation: develop in the rewrite repo (Code/macroir), and port to macro_dr HEAD only at the end. Reasons: the rewrite builds in seconds against macro_dr's ~25 minutes per translation unit; its "one body of code, two modes" automatic differentiation (der.hpp) yields the score and the Fisher information from the same source with no hand-written derivative code; and its validation apparatus (per-header test suites checked against independent numerical oracles, decision D9, plus the calibrate/distortion machinery) is exactly what a new member needs.

There is a governance catch that must be resolved first. The rewrite's charter rule R1 says every line must transcribe macro_dr, because a deliberate divergence anywhere makes parity comparisons uninterpretable. A Bessel member has no reference implementation anywhere to transcribe, so under R1 it is a research claim, and admitting it requires an explicit owner decision recorded in docs/decisions.md, in the mold of the existing precedent D12 (the one prior deliberate divergence). Validation then follows the D9 pattern (quadrature of the defining integrals, exact limits, finite differences) plus the exact box subcase of section 2, which chains the new member back to the frozen reference binary.

macro_dr HEAD enters at the last milestone, because two deliverables live only there: thermodynamic evidence (model comparison) and the DSL/cluster pipeline for the real P2X2 runs.

## 7. Plan of work

Each milestone has a gate: a check that must pass before moving on.

- M0, foundations (1 week). Write the derivation note for the new endpoint-conditioned integrals; record the charter decision; implement the pole tables (4- and 8-pole Bessel, normalized to the −3 dB cutoff) and the pole-shifted E2/E3 evaluations. Derive here, deliberately early, the composition rule for intervals containing an agonist concentration step (the analogue of the existing sub-interval composition in calc_Qdt): an error in it would corrupt exactly the pulse intervals that carry the kinetics. Gates: every new integral agrees with adaptive numerical quadrature of its definition; the Carbonell block-matrix route agrees with the spectral route; setting the poles to zero reproduces gtotal_ij/gtotal_sqr_ij to machine precision; near-coincident shifted arguments (the numerically dangerous case of divided differences) are stress-tested.
- M1, truth generator (1 week). Extend the simulator so it produces Bessel-filtered records exactly: over each constant-conductance piece of a channel's simulated path the filter equation has a closed-form solution, so the convolution is computed piecewise exactly, with the noise filtered by the same system, and the poles, cutoff and seed logged in the output. Gates: agreement with an independent high-rate digital filter built by matrix-exponential discretization (a generic off-the-shelf digital filter is a biased oracle and is not accepted); with the filter disabled, bit-compatibility with the existing validated simulator. In the same week, before any likelihood work: measure how badly the current box member does on filtered truth (interval coverage, lag-1 residual autocorrelation): these are the paper's motivation numbers. Also run one cheap positioning check: whether a filter-augmented linear noise approximation coincides numerically with this construction, as the box versions did at the 1e-8 level; if it does, the paper is positioned on cost and on independence from channel number, and no accuracy claim is drafted.
- M2, the member on the two-state C-O model (2-3 weeks). Implement the augmented recursion (the three MacroR blocks extended with m, Σ_nφ, Σ_φφ). Inside this milestone, on the critical path: a proof that the assembled joint covariance stays positive semidefinite; a new canary type for the filter blocks (following the existing canary-family convention); and an asymptotic branch for intervals with f_c·Δ > 10, where the filter has settled and the computation must gracefully revert to box cost (without this branch the 60 s intervals of the record produce numerical garbage). Gates: the deleted-rows subcase equals MacroIR to machine precision; predicted per-interval mean and variance match ensemble moments over 10⁴ simulated truths, including the interval containing the concentration jump; the lag-1 autocorrelation of standardized residuals is near zero where the box member on the same data shows the predicted 0.075/(f_c·Δ); the score matches finite differences.
- M3, one more state (1 week). Three-state C-C-O with a fast rate at or above f_c; run the calibrate() coverage battery at f_c·Δ in {0.2, 1, 5}. Gate: the box member fails the factor-of-two coverage criterion where the Bessel member passes. (The 61.5% figure from the appendix is an analytic bound on the white-noise component alone, so it is not used as a pass number.)
- M4, evidence on small models (1 week). Compute evidence for Bessel member vs box member on Bessel-filtered synthetic data, by direct quadrature or importance sampling through the R bindings (2-4 free parameters make this cheap). Gate: the evidence selects the Bessel member on Bessel truth, and the size of the box member's evidence bias is measured.
- M5, port and P2X2 (2-3 weeks plus cluster wall-clock). Port to macro_dr HEAD at the Qdt production seam (calc_Qdt_agonist_step and the Calc_Qdt_step memoizer, now keyed also on the pole set), one new recursion branch, and one kernel flag on set_Likelihood_algorithm plus its simulate twin. Gate: cross-repo parity on shared cells. Then P2X2 in three steps: regenerate the scheme_10 Bessel synthetic with the M1 generator and refit it; fit the real 2007 record at native rate, with the filter fixed at 10 kHz and 4 poles (never fitted); then the unliganded-gating models.

Total: about 8 to 10 focused weeks, with M2 the schedule risk. Division of labor as established: Luciano writes the hot production paths and compiles; Claude derives, transcribes and audits.

## 8. Cost and risks

- Runtime. Plain (no-derivative) mode: roughly 10-40 times MacroIR per distinct interval, from the table count of section 3 plus complex arithmetic, halved by conjugate symmetry. Two mitigations are verified against the actual record: only 103 distinct interval lengths exist, so memoizing on (Δ, agonist, pole set) removes most of the work; and 742 of 1078 intervals have the filter fully settled, reverting to box cost. The binding constraint is thermodynamic-evidence wall-clock on the cluster.
- Derivative mode is the open cost question. The pole-shifted spectral route needs scalar arithmetic on complex-conjugate pole pairs composed with the derivative type Der<double>, which does not exist yet; the plan commits to implementing it as real 2×2 blocks in der.hpp (decision 2 below). The fallback (block matrix exponentials) is correct but costs 100 to 1000 times baseline, which would sink fitting and evidence.
- Numerical stiffness. ν_p·Δ ranges from about 1 (20 µs intervals) to about 4·10⁶ (60 s intervals); the settled-filter branch of M2 is load-bearing, and without it the long intervals produce NaNs.
- Identifiability. The filter parameters (cutoff, pole count) are fixed from the rig and never fitted: they trade almost degenerately against fast rates and against Current_Noise. Rates faster than f_c remain weakly identified because the filter genuinely destroyed that information; wider posteriors there are the honest outcome, and a defect would be narrow ones. Pink_Noise (1/f) has no finite-dimensional linear realization and stays outside the augmentation as an additive per-observation floor. Proportional_Noise must either be declared zero for this member or given an explicit place in the update; inheriting it silently would mis-scale it.

## 9. Decisions that are Luciano's

1. Blocks all code: the repo and charter decision of section 6 (develop in the rewrite under an explicit R1 extension, D12-style).
2. Blocks M0: commit to the real 2×2 pole-pair arithmetic over Der<double>, so that derivative mode runs at the 10-40× cost rather than the 100-1000× fallback.
3. RESOLVED 2026-08-31: the generator of Sim_scheme_10_bessel.txt is ~/Code/python/pythonProject/bessel.ipynb (outside the repo, which is why no in-repo search found it). Recipe: read a full-rate 50 kHz scheme_10_inact simulation from the build dir (seed 836713249886511131, binary 0f15336), fill the gap rows with 0, design scipy.signal.bessel(order 4, cutoff 10 kHz, digital), filter with filtfilt (ZERO-PHASE, i.e. non-causal), block-average to the 1078-interval idealized grid, restore the real record's NaN gaps, write the file. Consequences: the old file was filtered acausally (no amplifier does that), its gaps were zeroed before filtering, and the noise was filtered together with the signal; it can serve as a rough sanity target, never as truth for a causal-filter likelihood. M1's regeneration stands, with this recipe documented as the contrast. Remaining action: copy bessel.ipynb into the repo for provenance.
4. The Current_Noise convention under a filter. σ²/fs and σ²/B_filter differ by a factor of 5 on this rig (docs/bibliography/recording_configurations/SOURCES.md:138-145), and the two scheme_10 posterior files already disagree by 4.24× in exactly this parameter. The Bessel member reads Current_Noise as the pre-filter power spectral density S₀; priors must move accordingly, and any evidence comparison across members must normalize the convention first.
5. Whether M4 runs in the rewrite (recommended: cheap, catches design errors before the expensive port) or waits for HEAD.

## 10. How the eventual paper is positioned

Concede early, as the MacroIR paper conceded the Kalman device: filters inside single-channel likelihoods exist (Michalek 1999; Venkataramanan, Kuc & Sigworth 2000; Qin, Auerbach & Sachs 2000; all discrete-time, FIR-truncated, with cost exponential in the filter memory per channel); absorbing a linear filter into an augmented state is textbook (Jazwinski 1970); the integral primitive is Carbonell 2008; overlapping-window augmented state-space models exist (Rubenzahl 2026). The defensible claim: the endpoint-conditioned mean, variance and cross-sample covariance of a Bessel-filtered current for N independent channels, computed from the generator's spectral decomposition at a cost polynomial in state count and pole count and independent of N, exact for the analog kernel where the published treatments are exact only for a sampled truncation; and the payoff of fitting the real P2X2 record at its native rate. Never write "nobody modeled the filter", and write no accuracy-over-LNA claim unless the M1 check leaves room for one.

## 11. Loose ends inherited from the survey

- Qin, Auerbach & Sachs 2000 is still unread (paywalled); the honest cost comparison against the FIR route (prior-art map item E11) is unmeasured.
- Jalali & Hawkes, "more realistic filters" (cited as unpublished results in Colquhoun & Sigworth 1995, p. 530): search for whether it ever appeared, before any first-ness sentence.
- Sigworth 1980 and Heinemann & Conti 1992 PDFs are still missing from docs/bibliography (standing rule: every paper looked at gets its PDF and a biblio entry).
