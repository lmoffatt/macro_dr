# Figure 5: is there a predicted breakdown boundary, and what may the figure claim?

> Updated: 2026-08-03. Owns **what Figure 5 is allowed to say**: whether a theoretical breakdown
> boundary can be drawn on it, which claims survive a circularity audit, and what the figure should
> carry instead. It does not own the figure's build (`../../projects/eLife_2025/figures/paper_both/`),
> the diagnostic definitions (`../_program/machinery.md`), or the axes and bands
> (`../_program/axes.md`), though §9 lists corrections owed to that file.
>
> **The verdict, up front: no predicted line goes on Figure 5.** Five candidate criteria were derived
> and every one was refuted, for one shared structural reason rather than five independent ones. The
> figure carries a decomposition instead, and the negative result is the stronger claim. §7 and §8 are
> what may be written; §3 to §5 are why nothing else may.
>
> PROVENANCE. A twelve-agent workflow of 2026-08-03 (`wf_cd2f528f-2c3`): five independent derivations,
> one measurement of the observed boundary from `data/1c2ae6f` over a dense ten-point N_ch axis (10,
> 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000 x 6 noise labels x 7 intervals x 2 anchors, 840
> cells), an adversarial verification per criterion, and a synthesis. Zero criteria survived. Scratch
> in `tmp/fig5_synth/`. Numbers below are from that run unless a file:line says otherwise; anything
> load-bearing should be re-measured before it goes into the manuscript with a `% src:`.

## 1. The question

Figure 5 shows where MacroIR stops being calibrated. The intent was to draw a breakdown boundary
computed WITHOUT reference to the code under test, so that "it breaks where theory says it must"
becomes an external validation rather than a self-report. The candidates were: occupancy
non-normality, the emission mixture and the gating-to-instrumental variance ratio, the interval
closure, the implementation's own guards, and the established LNA / moment-closure criteria.

## 2. Where the two Gaussians actually enter

Established by reading the code, and it corrects a loose statement that was circulating in the
planning documents.

**The state Gaussian is NOT in the propagation.** The covariance recursion
`P_Cov <- P^T (P_Cov - diag(P_mean)) P + diag(P_mean * P)` (qmodel.h:4389) is EXACT for the
multinomial; two-state occupancy stays exactly binomial under it. Propagating moments approximates
nothing.

It enters at the **measurement update**, where the gain `gS = Cov(X_end, y)` with `chi = dy/y_var`
(qmodel.h:~4617) is the linear projection of the occupancy onto the observed current. The linear
projection equals E[n | y] only under joint normality. So the state Gaussian is precisely the
assumption that **the occupancy conditioned on the measured current is Gaussian**, i.e. that the
conditional mean is linear in y. At ten channels and two states the occupancy takes eleven values and
its conditional mean is a saturating nonlinear function; the linear gain over- or undershoots.

**The measurement Gaussian** is a marginal statement about one sample. Conditioned on the endpoint
states, the interval-averaged current is a Markov-bridge functional: bounded, carrying atoms (the
path may have no transition), generally skewed. MacroIR matches its first two moments and adds the
exactly-Gaussian instrumental noise.

**The mapping that organises everything.** For a correctly specified sequential likelihood the score
increments are a martingale difference sequence, so Var(sum s_t) = sum Var(s_t). If the conditional
mean is wrong the innovations carry residual dependence, the increments correlate, and that equality
breaks. Therefore:

  state Gaussian   -> fails as the CROSS-TIME term (CDM)
  measurement Gaussian -> fails as the PER-SAMPLE term (GSDM)

All five candidate criteria are criteria on per-sample cumulants. Every one of them predicts GSDM and
none of them predicts CDM. That is the shared reason they fail.

## 3. The three refutations

### 3.1 Same gating fraction, boundary 30x apart

The decisive one, and it fits nothing. Two points on the measured anisotropy = 1.15 boundary, both at
noise label 0.05:

  interval 0.01  ->  N_ch = 370
  interval 0.5   ->  N_ch = 12.4

At those two points the gating fraction R/(1+R) is 0.9946 and 0.9957, identical to 0.1%. P_open is
pinned at 0.5 in all 294 cells. The interval-averaging factor differs by at most 1.35. Every criterion
in the set is a function of exactly those quantities, so none can move its boundary by more than
~1.35x between two points where the data moves it 30x. On the old four-point N_ch axis the same pair
reads 455 and 12, i.e. 38x.

Same inputs, thirty-fold different output: the mechanism is not in those inputs. This kills the family
at once rather than each member separately.

### 3.2 The N_ch exponent matches no cumulant mechanism

A cumulant-driven excess scales as a fixed power of N_ch: skewness-driven -0.5, kurtosis-driven -1.0.
Measured, floor-corrected, over the ten-point axis: **-0.237 (sim), -0.246 (pool)**. Neither, and no
standard cumulant mechanism produces a quarter power.

It is not even a power law. The exponent drifts from -0.20 at interval 0.01 to -0.42 at interval 1,
and from -0.32 at label 0.05 to -0.11 at label 10. Floor-free version: any exact cumulant argument
gives a^2 = a0^2 + A^2/N_ch, so A^2 from consecutive N_ch pairs must be constant. It runs 0.35, 0.51,
2.04, 2.47, 4.00, 4.30, 1.63, 14.2, 7.88.

### 3.3 The per-sample factor is not the factor that breaks

k_off diagonal at N_ch = 10, noise label 0.05, across the interval axis:

```
                 0.01   0.02   0.05   0.1    0.2    0.5    1
per-sample      1.058  1.105  1.202  1.240  1.174  1.014  1.013
cross-time      1.496  1.426  1.232  1.154  1.095  1.114  1.193
total           1.597  1.593  1.490  1.423  1.265  1.112  1.208
```

The per-sample part is a ridge peaking mid-axis at 1.24. The total peaks at the shortest interval at
1.60, where the per-sample part is 1.06. Over all 1470 macro_IR sim rows the cross-time factor exceeds
the per-sample factor in 66%, median |log m| 0.0116 against 0.0053. The criteria predict the middle
row; the boundary is set by the bottom one.

## 4. Two reasons a line is the wrong deliverable regardless

**The boundary is not one curve.** At the same cell the strict boundary is N_ch ~ 370 read by
anisotropy and N_ch ~ 4800 read by lambda_max of the same matrix: a 13x spread from the choice of
summary scalar. And lambda_max is still 1.165 at N_ch = 10000 in that corner against an estimation
floor of 1.063, so by that reading the algorithm never clears the strict criterion there.

**The inversion is ill-conditioned.** Distortion scales as roughly N_ch^-0.25, so solving for N*
amplifies error by four in the exponent: a 10% error in predicted distortion moves N* by 1.5x, a
factor-2 error by 16x. The MEASURED boundary already moved 20% (455 -> 370) just by densifying the
N_ch axis from four points to ten. Any threshold obtained by inverting a weak power law is a badly
determined target, whatever theory produces it.

## 5. What actually breaks

The state Gaussian at the update, with two measurable signatures.

The cross-time term, by the martingale argument of §2.

And r̄²_std, which is not contaminated by the score algebra. At N_ch = 10, label 0.05 it falls
monotonically across the interval axis: 0.9988, 0.9970, 0.9932, 0.9838, 0.9721, 0.9423, 0.9364. At
Δt = τ with ten channels MacroIR **over-predicts its own residual variance by 6.4%**, i.e. the linear
gain is leaving information on the table. Flat within 0.7% at every larger channel count.

## 6. The circularity ledger

**N_ch.** An absolute integer, not rescalable. But an absolute axis does not make a PREDICTED position
on it absolute: every candidate's N* formula contains the noise label and/or the interval, both of
which carry k_off, so the predicted crossing moves 1:1 with a mis-specified rate. **Claimable:**
direction, the measured exponent with the explicit caveat that it drifts and is therefore not a clean
power law, and the one-sided exclusion of §7. **Not claimable:** a predicted crossing point.

**Noise.** The label is 10·ν with ν = Current_Noise·k_off/g², normalised by k_off and the conductance.
Circular. One escape: the ratio R = gating variance / instrumental variance is a ratio of two variances
of the same record, both in current², whose only k_off dependence is the within-interval attenuation,
spanning 0.74 to 0.99 over the whole grid. R is therefore near-absolute and measurable by
non-stationary fluctuation analysis with no fit. It is the only non-circular combined quantity here,
and §3.1 is precisely the data rejecting it as the controlling variable (R² = 0.17 in the boundary
fits). **Claimable:** monotone improvement with instrumental noise, as a measurement plus the
cumulant-dilution derivation. **Not claimable:** a predicted noise threshold.

**Sampling interval.** interval_in_tau = Δt·k_off, fully normalised, and additionally confounded with
observation count: the recording is 10 τ long at every interval, so it holds 1000 observations at
d = 0.01 and 10 at d = 1, for six parameters. Nothing in the run separates coarser resolution of the
transient from fewer observations. **Claimable:** only relative statements in which k_off cancels
between the two axes. Everything else on this axis is descriptive and must be labelled descriptive.

## 7. What Figure 5 may claim: the five that survive

Each is derived first and confirmed after, which is what makes it a prediction rather than a
description.

**(1) Instrumental noise repairs, with the mechanism.** Gaussian noise is independent and contributes
nothing above the second cumulant, so with f the gating share of the predictive variance the
standardized cumulant of order k falls as f^(k/2): skewness as f^(3/2), excess kurtosis as f². A
second route runs the same way: higher noise weakens the gain, the filter tracks less, the predictive
gating variance rises toward its marginal value, and non-Gaussianity falls as its inverse square root.
There is no regime in which added instrumental noise hurts. Confirmed with no exception: at N_ch = 10,
interval 0.01 the anisotropy runs 1.358, 1.280, 1.225, 1.132, 1.121, 1.041 across labels 0.05 to 10.

**(2) The one-sided exclusion.** Occupancy non-normality alone cannot produce a factor-2 information
distortion anywhere on this grid. Derived from the binomial cumulants; confirmed by the grid maxima
(anisotropy 1.42, magnitude 1.077, lambda_max 1.78). Reaching factor 2 by that route needs fewer than
about one resolved channel event per informative sample. One-sided bounds are the most defensible
class here because they do not depend on the exponent being right.

**(3) The square-root-of-noise relation, the only non-circular quantitative prediction in the set.**
The interval at which per-sample non-Gaussianity is worst moves as (noise label)^(1/2). It follows on
paper: the instrumental variance per averaged sample goes as 1/Δt while the gating variance of the
average falls with Δt through the averaging factor, so holding the ratio fixed gives Δt* ∝ √noise.
Measured exponent **+0.44 ± 0.04**, consistent with +1/2 at 1.6σ. k_off cancels between the two axes,
so it is genuinely external. It belongs on the per-sample component, never on the total.

**(4) The two Gaussians have DIFFERENT interval signatures, and the data separates them.** Predicted
before measuring: the measurement Gaussian is marginal and must be a ridge with an interior maximum,
because short intervals are masked by an instrumental variance growing as 1/Δt and long ones
self-average; the state Gaussian is an update failure and must be worst at short interval, where the
filter does the most tracking. The decomposition of §3.3 shows exactly that. This is a prediction of
SHAPE and not of position, which is why it escapes the circularity of §6: it needs no absolute scale.

**(5) Rank deficiency at high noise.** As the gating share vanishes the GIDM goes rank deficient, the
unitary-current and channel-count rows become identical, only the product N·i is identified, and the
emitted spectral scalars are computed over five eigenvalues instead of six (p_implied = 5.00 in those
rows). Predicted from the mean depending only on N·i, confirmed in the emitted spectrum. It is the
same mechanism as the paper's opening blow, appearing measured.

**And the negative prediction, which is the strongest thing here.** The breakdown was predicted to sit
in the per-sample cumulants, which is where any reader would look, and it measurably does not: 1.06
per sample where the total is 1.60. It sits in the cross-time term. A falsified hypothesis reported as
falsified, with a mechanism-level consequence, is stronger than a fitted curve that agrees, because a
clean refutation cannot be accused of parameter choice. It is also the fourth instance of the paper's
signature, now applied to the paper's own method: the available check passes and the failure lives one
level up.

## 8. What the figure should carry

Not a boundary. **The decomposition**, over the same plane: total, per-sample part, cross-time part.
All three are already emitted (GIDM, GSDM, CDM), so this costs no run. The figure then makes a claim
no boundary could: the part theory predicts is small, and the part that carries the boundary is the
part theory does not predict.

Plus the **estimation floor** drawn explicitly: anisotropy 1.043 ± 0.010 (flat, 1.040 to 1.049 across
intervals) and lambda_max 1.063. Anything below those is measuring the bootstrap, and the reader has
to be able to see it. Note that IR's ~1.04 readings elsewhere in the set sit at this floor.

The two surviving external anchors, (3) and (5) of §7, go beside it as evidence that the machinery is
tied to something outside itself. Neither is a boundary; both are non-circular.

## 9. Roads deliberately not taken, and why they were available

Recorded so the decision is not re-litigated, and so the cost is known if it is ever revisited.

**micro_IR as the exact reference. NOT TAKEN (author's decision, 2026-08-03).** It would have worked
and it already exists: `micro_full.h:1229` shows avg=2 doing plain Bayes on the joint (s_start, s_end)
space, with the state space being the OCCUPANCY (M = N+1 for two states, per the underflow note), so
the state side is exact and the only remaining approximation is the emission Gaussian. Running it
beside macro_IR on the same cells would have separated the two Gaussians experimentally: whatever
distortion remains in micro_IR is the emission Gaussian alone, and the difference to macro_IR is the
state Gaussian. Cost is (N+1)² per step, which is why it has runs at N_ch 5, 10 and 20 only, and that
range is exactly where the breakdown lives. It would also have supplied the exactness anchor the
code-verification claim otherwise lacks.

**Higher-cumulant emission. NOT TAKEN (author's decision, 2026-08-03).** The third and fourth cumulants
of the bridge-averaged conductance come from the same Van Loan block-exponential construction that
already produces `gmean_ij` and `gsqr_ij`: a block-bidiagonal matrix with Q on the diagonal and
diag(g) on the superdiagonal, whose (1, m+1) block gives the m-th moment integral over m!, so ONE
exponential of the five-block matrix yields moments one to four. For two states that is a 10x10
exponential per distinct Δt, negligible. With them, the per-sample distortion is a closed-form
function of γ3 and γ4 and the derivative vectors, with no free parameters, so GSDM becomes predicted
rather than fitted.

Two cautions were established and should be carried if this is revisited. The family matters more than
the moments: Gram-Charlier / Edgeworth goes NEGATIVE for moderate γ3 and γ4, precisely in the regime
where it is needed, and a likelihood that can go negative is unusable; Johnson S_B (bounded, and the
bridge average IS bounded) or the Pearson system are the valid four-moment alternatives. Better than
either, because the instrumental noise is convolved on afterwards, is to expand by TRANSITION COUNT
under the uniformization already in the code: the m = 0 term is an exact atom at g_i with weight
exp(-q_i Δt)/P_ii, the m = 1 term is uniform on the segment between g_i and g_j, both closed form
after convolution, and the m = 0 term is exactly what the Gaussian gets most wrong at short interval,
which is where the total distortion peaks.

**What the two together would have proved.** If a better emission removed micro_IR's residual
distortion, the responsible Gaussian would be identified by ELIMINATION rather than by decomposition:
an approximation removed and the number gone, which is causal where §3.3 is only correlational.

## 10. Still open

1. What produces the cross-time factor. It is the actual boundary and no candidate models it.
2. Why the N_ch exponent is -0.25 and why it drifts from -0.20 to -0.42 across the interval axis.
3. Which functional of the spectrum the paper is claiming about: anisotropy and lambda_max disagree by
   13x on the boundary in the same cell.
4. Whether the controlling group is N_ch or N_ch·p(1-p). Untestable as run: P_open is pinned at 0.5 in
   all 294 cells. **One cheap column at P_open = 0.1 would settle it**, and the mixture-kurtosis
   analysis predicts the SIGN of the magnitude distortion flips there because rare openings are
   leptokurtic. A sign-level prediction from an untested lens is the cheapest external validation
   available anywhere in this set, and it also settles whether "low P_open ≈ fewer channels" is true.
5. Reproducibility: seed = 0 is the `random_device` sentinel and the resolved seed is not logged, so no
   crossing is reproducible to better than the bootstrap width. The CI half-width on the far-field
   anisotropy is 0.011, which is 8% of the distance from the floor to 1.15.

## 11. Corrections owed to other documents

**`../_program/axes.md`, band boundaries.** The instrumental variance per averaged sample is
e = Current_Noise·f_s/n_samp = 0.1·z/d, and the single-channel gating variance of the interval average
is i²·p(1-p)·φ(2d) = 0.25·φ(2d). Equating them puts A/B at label = 2.5·d·φ(2d) and B/C at
label = 2.5·N_ch·d·φ(2d). The document carries 10·d and 10·N_ch·d: a factor 4 too high at short
interval (the missing 1/(p(1-p)) = 4) and 7 at d = 1 (the averaging bend). More of the existing grid
is already in band C than the coverage table says. TO VERIFY before editing.

**`../_program/axes.md` §1, the definition of τ.** It defines τ as 1/max|Re λ| of Q at the agonist
condition (1/200 s), while the dispatched `interval_in_tau` is Δt·k_off with τ = 1/k_off = 0.01 s. A
factor 2. Any criterion or caption written in "τ" must say which.

**The long-interval premise, stated in several planning documents and in the audio record.** On the
information distortion MacroIR does NOT degrade monotonically with interval; over most of the grid it
improves. At N_ch = 10, label 0.05 the anisotropy runs 1.358, 1.382, 1.342, 1.330, 1.270, 1.159, 1.176:
a shallow U with its floor near 0.5 τ and a small rebound at 1 τ. Over 42 cells the maximum sits at the
shortest interval in 13 and at the longest in 14; the minimum sits at 0.5 τ in 19. Median Spearman of
anisotropy against interval is -0.20 (sim) and -0.30 (pool), i.e. slightly IMPROVING. The genuine
monotone long-interval effect is r̄²_std and only at N_ch = 10 (§5). Ruled out as artefacts: segment
straddling (n_step_1 = 2/d and n_step_2 = n_step_3 = 4/d are exact integers at all seven grid
intervals, so no segment edge falls inside a sample) and estimator spread as the main driver. Real but
second-hand: Wald coverage degrades at long interval at every N_ch, tracking p/n_obs = 6d/10, which is
the textbook finite-sample rate and not misspecification. TO CONFIRM before this replaces the current
statement anywhere.

**Positioning against Münch 2022, and it is mandatory.** Their noise axis is state-dependent
open-channel noise, which the classical filter mismodels, so their error GROWS with noise. This paper's
axis is additive instrumental noise, which every member models correctly, so error FALLS. Without a
sentence saying so the two validity maps read as contradictory.

**The one control run that would separate the interval confound**, if it is ever wanted: one column
rerun with the observation count held fixed by lengthening the two plateau segments only, leaving the
transient physically identical. One column, not a grid.
