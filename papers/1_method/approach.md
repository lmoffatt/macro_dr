# The approach: how do you trust the tool?

> Updated: 2026-08-11 (§9 and §10 only; §1 to §8 and §11 are still the 2026-08-03 second pass and
> have not been re-read against the newer audios). Owns **the framing**: the reader's problem, the
> argument that answers it, the order it is told in, what is claimed, what is deliberately handed off,
> and what is still in doubt. It does not own the roster or the cells (`00_plan.md`), the diagnostic
> definitions (`../_program/machinery.md`), the section-by-section claim spine
> (`docs/manuscript-drafts/sections/README.md`), or what Figure 5 may say about the breakdown
> (`analysis_figure_5_breakdown_criterion.md`).
>
> SOURCE. The seven audios of 2026-08-03 (09.56 to 10.15), 2026-08-02 17.40.24 / 17.40.59 / 17.42.32
> / 17.52.43, and 2026-08-01 19.01.27 / 19.04.19, in `program/source-notes/audios/audios/Chat de
> WhatsApp con MacroIR 13/`; plus a working session the same day. Where the session CORRECTED an
> audio or an earlier draft of this file, both the old and the new statement are recorded, because a
> silently replaced claim is one that comes back.
>
> Two multi-agent runs stand behind §5 and §6: `wf_cd2f528f-2c3` (the breakdown criterion, five
> candidates, zero survivors) and `wf_c7d9da53-c4b` (the single-record variability of the diagnostic).

## 1. The spine

**How do you trust the tool?** The tool is not the algorithm. It is the algorithm **and the executable
code**, and a biophysicist has no way to check either one. That single question contains both halves
of what the paper does, and it is why the "algorithm validity" framing and the "code verification"
framing are one paper and not two.

The protagonist is the READER. The paper is the mentor, not the hero: it hands an instrument to
someone who has a problem. That fixes the voice of every section, and it decides where the climax
goes (§3, beat 8).

The villain is OPACITY, not least squares. With opacity as the villain the paper can recommend least
squares where it is adequate without contradicting itself, and least squares becomes a fallen ally
rather than a target.

## 2. The signature: local passes, non-local fails

One structure recurs, and naming it is what turns a series of measurements into an argument. The
first draft of this file called it "the available check passes", which was the weaker version: it
conflated checks a practitioner runs today with quantities this paper introduces, and a referee would
say so. The defensible statement is about LOCALITY.

**Everything local is right and everything non-local is wrong.** A point in parameter space is right
while the covariance around it is not; the fit in time is right while the curvature is not; a single
sample is right while the accumulation over samples is not. The mean is easy, the covariance is hard,
and the covariance across time is hardest.

| # | Where | The local object, and it is right | The non-local one, and it is wrong | Data |
|---|---|---|---|---|
| 1 | Figure 2 | the point estimates recover for every member | the reported covariance is wrong, and wrong in SHAPE, so no rescaling repairs it | measured |
| 2 | Figure 3 | per sample INR is calibrated (r² 0.9995, better than R) | accumulated 20.8, score autocorrelation 0.868 | measured |
| 3 | the pearl | the least-squares fit is visibly good in time | its numerical and Gaussian Fisher disagree | **no LSE column yet**, §10 item 1 |
| 4 | Figure 5 | MacroIR's own per-sample non-Gaussianity is small (1.06) | the total is 1.60, carried by the cross-time term | measured |

Instance 4 was found on 2026-08-03 after five candidate breakdown criteria were refuted, all of them
per-sample. That it arrived last matters and the paper must say so: a pattern presented as four
confirmations, with the fourth obtained by elimination, invites the reader to discount all four. The
refutation therefore goes in the paper visibly (§5), not only in the planning layer.

## 3. The beats

**Beat 0. The ordinary world.** Macroscopic currents are fitted by least squares. It works and you can
watch it work. *Figure 1.*

**Beat 1. First blow: information is lost.** Least squares fits the mean, and in the mean N_ch and i
enter only as the product N·i. What separates them lives in the fluctuations, which the method
discards. There are questions it cannot answer at all. *The empty LSE panels of Figure 2 are the loss,
drawn.* Stronger than opening on calibration, because it is a hard identifiability fact that needs
none of the paper's machinery and is the historical reason the field does fluctuation analysis.

**Beat 2. Second blow: what you can ask, you cannot trust.** Two failures at once, before any
machinery. The centre drifts off truth for the non-interval members, and the ellipse mis-sizes AND
mis-shapes the cloud: NR is off by 15x, 10x and 8.8x depending on which pair you ask about, and that
the three numbers DIFFER is the anisotropy. *Figure 2: the size and shape numbers per panel, and the
Fisher ellipse rescaled to the empirical area, which still misses.*

Carry one forward reference here, at the cost of a sentence: **there is a regime where none of this
matters, and beat 6 says which**. Without it the least-squares user leaves at beat 2 and never reaches
their pardon.

Two things this beat does beyond the blow. It introduces the interval axis without naming it, since
the members that model the interval average have no drift, and Figure 3 supplies the other half of the
dissociation. And the open marker is truth-plus-PREDICTED-bias, so the reader sees the machinery
anticipate the drift before they know they need it.

**Beat 3. The call and the refusal. END OF ACT I.** The fix is recursion. But a recursive filter's
prediction tracks the data because it knows the previous point, so it looks right whether or not it
is. In the author's words, it is like it cheats. Underneath: how do I trust the code, when it is not
easy to tell whether it is telling the truth? *Prose, plus Figure 3 row A, where predicted against
observed looks fine in all seven columns.*

**Beat 4. The instrument.** Two classical identities evaluated at the truth over simulated ensembles:
the score has mean zero, the Fisher equals the covariance of the score. Neither is new and the paper
says so.

On why one instrument serves, and the first draft of this file argued it BACKWARDS. It said the
identities are "indifferent to the cause, which is exactly why one instrument answers both halves".
That is wrong and a statistical referee will say so in a line: indifferent to the cause means it
answers NEITHER separately, it answers a joint question about (code, approximation, specification).
The defensible statement is about the READER'S DECISION, not about what the instrument verifies:
**the reader does not need the distinction, because their decision is the same either way.** Whether
the number in front of them is untrustworthy because of a coding error or because the closure has
failed, they act identically. Separating the causes is a different task and needs an exactness
anchor, which is what micro_IR and the closed-form simulator checks would supply (§10, §11).

**Beat 5. MIDPOINT. Transparency returns, and with it the second deception.** The identities are
visible at the trace level, so the reader gets eyes back, and the first thing they see is that the
easy check lies again: per sample INR is perfect and its accumulated information ratio is 20.8. The
double dissociation: interval averaging repairs the data level and nothing else, recursion repairs the
accumulation and nothing else, only both repair both. *Figure 3, rows in the order of the revelation,
now marked in three blocks with the E|F boundary carrying the reveal.*

**The double dissociation is undervalued in the current draft and should be promoted.** It is the most
TRANSFERABLE result in the paper: it holds for any state model with time-averaged observations, not
only for ion channels. It can be the title claim and still be the mechanism of the trust story.

**Beat 6. The reprieve, and which world you live in.** Noise equalises everyone. Outside that region
the interval members always win on bias and only MacroIR always wins on distortion. *Figure 4, LED BY
the noise result and not by the generalisation of the ladder.*

Close the false exit in the same paragraph: noise equalises because it BURIES the gating fluctuations,
so everyone is equal because everyone is blind to the same thing. The reprieve and the first blow are
one mechanism. Confirmed measured: at high noise the information matrix goes rank deficient and only
N·i is identified.

**Beat 7. MacroIR has measurable limits, and the paper maps them.** Priced carefully, because this
beat has been mis-set twice in this file: once demoted to "a complication about the paper", which
undersells it, and once called a first-order result, which oversells it.

It does two jobs, and both are real:

  * **It bounds the map.** Without a boundary there is no map, there is an endorsement.
  * **It is a control on the instrument.** If IR read 1.00 everywhere with no boundary anywhere, a
    reader would rightly suspect the diagnostic had been built to flatter it. That the identities
    catch the author's own filter departing is evidence they are not tuned to the answer. Modest but
    real: it shows the diagnostic is not rigged, which is not the same as showing it is right.

And what it is NOT worth, stated so nobody inflates it later:

  * **The departures are small.** Over 294 cells IR's maxima are anisotropy 1.42, magnitude 1.077 and
    lambda_max 1.78, against an estimation floor of 1.043. Dramatising that into a fatal flaw to look
    humble is dishonest in the other direction. The honest sentence is that the limits exist, are
    small, and are confined, and that this is most of the reason to trust the filter.
  * **It is measured and not explained.** The paper cannot say why, cannot predict where, and can
    claim the position on only one of the three axes, and there only as a bound (§5).
  * **Credibility is table stakes**, not a contribution. Mapping your own method's failure is what an
    honest paper does; it earns trust, it does not earn a result.

Keep this SEPARATE from the methodological note that accompanies it, which is that the paper looked
for that boundary where anyone would look, in the per-sample cumulants, and it measurably is not there
(§5). Merging the two makes the first read as an apology for the second.

**Beat 8. CLIMAX. The gift, and beat 7 is part of it.** The moment the instrument changes hands. In a
mentor structure this is the peak; the first draft of this file put the climax at beat 7, which
switches protagonist mid-story, and the fix is not to demote beat 7 but to recognise that **the
instrument that changes hands includes its limits**. "Here is the tool, here is where it works, here
is where it does not, and here is how you tell which you are in" is one deliverable. Strip the second
half and the gift becomes an endorsement; strip the first and there is nothing to give. *Figure 6, and
Figure 6B.* See §7 for what the gift is and §6 for what it deliberately is not.

**Beat 9. The lesson, as a mechanism and not a proverb.** "What is derived passes and what is patched
fails" is a proverb with n = 2 whose two instances the author built, which invites the reading that MR
and VR are straw men. The mechanism is in hand and is far stronger. MR uses the TOTAL conductance
variance per start state while IR uses the residual variance conditioned on both interval ends, and
because gvar_total = Σ_j gvar_ij + Var_j[gmean_ij], MR carries the end-state spread as plain
observation noise where IR carries the same magnitude in the state-coupled term, so only IR's routing
lets that spread reach the GAIN. The two therefore predict the SAME total variance from the same
prior; what MR cannot do is resolve it. <!-- CORRECTED 2026-08-06: this beat used to say MR "dumps the
end-state spread into the observation variance instead of resolving it", which reads as MR predicting
more variance than IR. It does not, at equal state: the boundary term enters gSg with a plus and ms
with a minus and cancels (verified term by term against legacy/qmodel.h:4568-4619, numerically on
random fields, and on the figure-1 dumps, where MR = IR = 1.048475625791748 while they still share a
prior). The misallocation is real and it is what the beat needs; it just shows up in the gain, and in
the variance only after the gains have driven the priors apart. Canonical: figures_build_plan.md
§195-245. --> **A partial
interval correction cannot interpolate because what is missing is not in the variance, it is in the
gain.** That is why MR is WORSE than R rather than sitting between R and IR: 1.30 to 1.56 in size
against R's 1.09 to 1.24, and 1.36 to 1.40 in shape against R's flat 1.18.

**Subplot: the two pearls.** The information ledger, and deceptive simplicity. The second is instance
3 of §2, so it sits LATE, near the lesson, where it reads as a rhyme.

**Open door.** Where MacroIR's breakdown comes from, stated as not understood, with the location of
the ignorance now identified (the cross-time term). And the limits look like a balance that has to
come out just right, which suggests a conserved quantity waiting to be named. As a conjecture, in
those words.

## 4. The declared scope

Three axes for the algorithms: recursion; interval; and whether the member is DERIVED or PATCHED,
which is a property of how it was arrived at. Three for the conditions: channel count, sampling
interval, instrumental noise.

**Outside, and named:** episodic against stationary recording, the number of states, the microscopic
regime, and evidence or model comparison.

**Also outside: the prior.** The study is at the level of the likelihood alone. That decision has a
second consequence and the first draft of this file OVERSTATED it. It said simulation-based
calibration is "unavailable by construction". False: the prior in SBC is an instrument of the check,
not an inferential commitment, so any prior may be placed on theta purely to calibrate and the paper
stays a likelihood paper. What is missing is the sampler running hundreds of replicates, and this repo
has parallel tempering. **It is a cost, not an impossibility**, and written as impossibility it invites
exactly the referee who knows better. The honest sentence is short: SBC is the standard, it is
expensive here, and the score identities are the likelihood-level analogue that needs only derivatives.

**Three restrictions that belong HERE and not buried in the doubts**, because §4 is what a referee
reads to judge scope:

- **P_open is pinned at 0.5** in all 294 cells, so the paper cannot distinguish whether the controlling
  group is N_ch or N_ch·p(1−p).
- **Two states.** This is the BEST case for any Gaussian approximation, since the occupancy is exactly
  binomial and its propagation is exact (qmodel.h:4389). Everything reported is therefore a LOWER
  bound on the problem, which turns the restriction into a strengthening if it is stated.
- **The interval axis is confounded with the observation count.** The recording is 10 tau long at every
  interval, so it holds 1000 observations at Δt·k_off = 0.01 and 10 at 1, for six parameters. Nothing
  in the design separates coarser resolution of the transient from fewer observations.

## 5. Beat 7 as it now stands: no predicted boundary

The original design put a predicted breakdown boundary on Figure 5, so the instrument would be
anchored to a calculation that did not come from the code under test. **Five candidate criteria were
derived, measured against the observed boundary and adversarially tested on 2026-08-03; none
survived.** Full reasoning in `analysis_figure_5_breakdown_criterion.md`. The one-line summary: every
criterion was a criterion on per-sample cumulants, and the per-sample factor is not what breaks.

Consequence for §9 of the first draft, which is now RETRACTED: it claimed that "the failure boundary
in channel count is predicted from the binomial without reference to the code, so its position is an
external check". That is refuted. The measured N_ch exponent is −0.25, neither the −0.5 of a
skewness mechanism nor the −1.0 of a kurtosis one, and it is not even a clean power law. What survives
is a one-sided EXCLUSION, which is a bound and not a position: occupancy non-normality alone cannot
produce a factor-2 distortion anywhere on this grid.

Figure 5 therefore carries **the decomposition** (total, per-sample, cross-time, all three already
emitted) and the **estimation floor** drawn explicitly, not a boundary. Two external anchors survive
beside it: the rank deficiency at high noise, and the square-root-of-noise relation for the worst
per-sample interval (+0.44 ± 0.04 against a derived 1/2).

## 6. The demarcation: what this paper does not settle

Decided 2026-08-03. The broad question, "when should I use MacroIR", **cannot be closed with the
elements of this paper**, and the honest move is to demarcate it rather than answer it badly. What the
paper cannot supply depends on facts about the reader's experiment:

1. **The cost comparison in general.** It depends on the number of states, the implementation and the
   optimizer. The paper has one model, one implementation, one machine.
2. **Whether a sweep-level bootstrap suffices.** It depends on how many sweeps the preparation
   tolerates and whether they are exchangeable, which is a rundown question.
3. **How much precision is lost with least squares, for YOUR question.** Measured here on this model;
   whether it matters is the reader's.
4. **How this composes with non-stationary fluctuation analysis** as a complete workflow.
5. Everything already named as outside in §4.

This demarcation is not a retreat, it is what makes the deliverable usable, and it changes what the
paper hands over: **a diagnostic, not a verdict.** The reader supplies their sweep count, their
rundown, their model size and their tolerance, and decides. The paper supplies the part they cannot
measure for themselves. A short subsection in the discussion, "what this does not settle", naming
these five, disarms the most likely objection, which comes from the reader who has plenty of sweeps.

## 7. The gift, in the reader's terms

Beat 8 must NOT say "and which method to use", which is the half §6 hands off. It says: **whether the
number your fit gave you is believable, and the one measurement that tells you.**

**Every coordinate the reader needs is obtainable BEFORE any kinetic fit**, and this is a positive
property that should be stated rather than a circularity to be defended. Peak current gives N_ch to a
factor of a few, the observed relaxation gives the time scale (they chose their sampling rate for it),
the baseline variance gives the instrumental noise, and the unitary conductance comes from the
literature or a variance-mean plot. On axes spanning decades, order-of-magnitude prior knowledge is
enough. Expressing the map in dimensionless groups is what makes it transferable to other channels;
indexed in raw pA² it would serve one preparation.

**The procedure, per sweep:** take the residual of the least-squares fit; divide by the constant SD
the fit assumes, the same one least squares uses (NOT a time-varying variance, or the number stops
being comparable to the axis); compute τ_int = 1 + 2·Σ_k ρ_k, equivalently N/N_effective, null 1 for a
white residual; repeat per sweep and **average the τ_int estimates, not the traces** (averaging traces
gains nothing, since the gating noise is independent across sweeps too, so the autocorrelation is
unchanged). On a stationary segment the fit is not needed at all: the residual is the trace minus its
own mean, which is classical fluctuation analysis.

**It is a K-sweep procedure and K must be stated.** Measured 2026-08-03 (`wf_c7d9da53-c4b`): a single
recording is NOT enough. At Δt·k_off = 0.1 the single-record spread makes an IR-versus-MR
discrimination fail 67 per cent of the time, worse than a coin flip, and a reader cannot distinguish
τ_int = 1.0 from 1.5 on one record at any interval coarser than 0.02. But K is small enough to
survive: for a 95 per cent half-width of 0.25 it needs K = 3, 6, 13, 27, 42, 93, 160 across the seven
intervals, and at 0.1 a measured K = 25 drops the error rate from 67 per cent to 3. A patch clamper
routinely has 5 to 50 sweeps. The self-centring objection is closed by measurement: a reader may run
the identical estimator centred on their own K-sweep mean, and the bias is 9 per cent at K = 2 and
1 per cent at K = 5.

**Two caveats that go in the caption.** The window is FIXED at max_lag = 10
(`legacy/moment_statistics.h:62-70`, and every production config sets it), so τ_int has a ceiling near
18.9 and the high end of the axis is compressed against it. And the two error directions run opposite
ways: a badly fitted mean INFLATES τ_int, which errs safe, while a short record DEFLATES it, which
does not, and that asymmetry is the reason to average sweeps.

**The reading, which is the point:** near 1, least squares is reporting honestly. At 3, the error bars
are wrong by about a factor 2. At 10 or more, by 4 or 5, and a nominal 95 per cent interval is not 95
per cent. In error-bar units take the square root of the y axis, since the distortion is a ratio of
informations: 1.15 is a 7 per cent error and nobody should care; 4 is a factor two.

## 8. What the alternatives cost, and what they answer

From the audios of 09.56.22 and 09.58.14, worked out in session. The recommendation is incomplete
without this, because a reader can answer "fine, I will bootstrap".

**Measured, and the only number worth quoting:** MacroIR costs 1.4x to 5x a least-squares likelihood
evaluation on this two-state model (1.012 s against 0.200 at Δt·k_off = 0.01, 0.0103 against 0.0072 at
1), and that ratio is **independent of the channel count**, which is the axis on which a reader will
assume a filter gets expensive. The ratio grows with the number of states as any filter's does; the
paper should NOT commit to where it crosses, since that needs a benchmark it has not run and invites a
fight about implementations.

| route | cost | what it answers |
|---|---|---|
| SEM across independent fits | free, you already have the fits | the uncertainty of the POPULATION MEAN across patches, mixing statistical and biological variability; does not remove bias, which stays while the SEM shrinks as √n |
| sweep-level bootstrap (resample whole sweeps, refit) | B refits, needs K ≥ 20-30 | **valid**, and it correctly preserves the within-sweep correlation. Gives honest uncertainty for whatever estimator you have; says nothing about whether that estimator is the right one |
| residual bootstrap | cheap | **invalid**: iid resampling destroys the correlation that is the problem |
| parametric bootstrap | B refits | assumes the model you were trying to check |
| block bootstrap | B refits | needs a block length much larger than the correlation time, i.e. it needs τ_int before it can be configured |
| simulation + sandwich | M simulations, and M must be ~10⁴ | the within-recording uncertainty, correctly |

**Why M must be about 10⁴ and not a few hundred**, measured on the same cell at nsim 200 / 256 / 1000
/ 1024 / 10000: the magnitude is unbiased and converges by nsim ~1000 (±3.7 per cent), but the
ANISOTROPY is biased UP and the bias decays as 1/√M. It reads 1.270, 1.265, 1.149, 1.127, 1.064 while
the floor is 1.043. **At nsim 1000 the bias alone equals the strict criterion being tested against**,
so anyone reproducing this measurement with fewer simulations manufactures an anisotropy that is not
there. This belongs in methods, and it is why 10 000 was run.

**The honest conclusion is not "MacroIR is cheaper"**, which reads as special pleading. It is: the
free route answers a different question and carries the bias forward; the valid cheap route (the sweep
bootstrap) gives an honest error bar around whatever the estimator converges to and is silent on
whether the estimator is right; and the route that answers the same question as the filter costs a
multiplier of order 10⁴ set by the precision needed on a covariance rather than by the model, applied
to the SIMULATION and not to the cheap evaluation, and it must simulate from the fitted model, so it
assumes what one wanted to verify.

**And where sweeps are plentiful, say so plainly.** With enough exchangeable sweeps, non-stationary
fluctuation analysis gives N and i with no likelihood at all, and a sweep bootstrap gives honest
intervals with no model. That is a real regime and it covers much of routine electrophysiology. The
case for the filter is strongest exactly where the classical route is weakest: few sweeps, rundown,
non-exchangeable conditions, a single precious recording, and wanting kinetics and amplitudes jointly
from the same record. Note also that Figure 6's own measurement puts the EFFICIENCY gap at a factor
1.9 on k_on at N_ch 10 and small elsewhere, so the recommendation cannot be "always use the filter".

## 9. How each idea from the audios landed

| Idea | Source | State |
|---|---|---|
| The tool is algorithm + code; the problem is trust | 10.11.53 + session | **the spine**, §1 |
| Emotional structure: disquiet, scheme, resolve, lesson, pretty birds | 09.59.31, 10.04.51 | **§3** |
| The paper as dimensions | 09.56.22 | **§4** |
| The identities are visible at trace level | session | **§3 beat 5**; already Figure 3 rows B–G |
| Per sample INR looks fine, you need the autocorrelation | session | **§2 row 2**, measured |
| Non-interval members also fail in BIAS | session | **§3 beat 2**, but see §11: the bias is smaller than claimed |
| Anisotropy: least squares' distortion is not isometric | 17.42.32, 17.52.43 | **implemented**, §10 |
| Figure 2 observes it, a Figure 4 supplement maps it | session | **implemented**, §10 |
| Noise equalises everyone | session | **§3 beat 6**, with the false exit closed |
| The practical rule: the residual autocorrelation | 17.40.59 | **§7 and Figure 6B**, now measured as a K-sweep procedure |
| Derived passes, patched fails | 10.11.53 | **§3 beat 9**, now as a mechanism |
| Deceptive simplicity | 10.04.51 | **§3 subplot**; the figure exists but has no LSE column (§10) |
| The Fisher as an information ledger | 10.01.07 | **§3 subplot** |
| MacroINR has no winning scenario | 17.40.24 | to be written into the recommendation |
| Cost of CIs by simulation, and LSE + bootstrap, against MacroIR | 09.56.22, 09.58.14 | **§8**, worked out |
| Least squares with all variability as instrumental noise | 10.04.51 | **not started** (§10) |
| Limits are conjecture; a conservation may be waiting | 10.15.18 | **§3 open door**, §11 |
| Control numbers, not figures | 19.04.19 | partly measured |

Added 2026-08-11, from the audios of 08-09, 08-10 and 08-11. Same table, later batch. **Read the
08-10 rows against the clock:** those audios are from 14:18 to 14:32, and the Introduction was
rewritten and committed at 23:25 that night, so several of them were acted on within hours and the
audio is the request, not the state.

| Idea | Source | State |
|---|---|---|
| No intrinsic signature any more: the 2007 flip state showed in a model-independent delay, the mechanisms at stake now depend on model validity | 08-10 14.32.13 | **already written**, and twice: `docs/manuscript-drafts/sections/01_introduction.tex:254-263`, the opening P2X2 paragraph, and `05_discussion.tex:229` |
| Reorient the paper: MacroIR and its validation only, the rest an addendum | 08-10 14.32.13 | **withdrawn by the author 2026-08-11**, in his terms a temporary weakness, back to the merged version. Tombstone at `../_program/program.md` §9, because the audio carries no retraction |
| A valid likelihood is a sine qua non for computing evidence | 08-10 14.25.33, 14.27.02 | **carried as necessity**, `01_introduction.tex:365-368`; the sufficiency version of it was removed as an overclaim on 08-10, see below |
| The rationale of the three ranges belongs in the Introduction | 08-10 14.25.33 | **written into its owner**, `../_program/axes.md` §2ab; still to be condensed into the paper, and one of the three needs a citation |
| Validation means simulator plus likelihood, checking the score and the Fisher; an implementation error shows as a departure | 08-10 14.18.17 | **already written**, `01_introduction.tex:358-380`, closing sentence |
| Apply the validation to the 2025 evidence results | 08-10 14.28.24 | not started, and it is the next paper rather than this one |
| Usage recommendations by regime: MacroR at low noise, MacroIR almost everywhere except few channels and low noise, non-recursive interval members unbiased but with no error bar | 08-10 14.25.33 | **in tension with §6**, which demarcates rather than recommends. The measured content is Figure 7's; the recommending voice is not the paper's |

| Idea | Source | State |
|---|---|---|
| Derive it twice, from MacroR and from the Kalman side, and show the two agree | 08-09 08.32.54 | **§10 item 9**, undecided |
| Say why MR fails: it marginalises too early | 08-09 08.36.28 | **§10 item 10**, and the wording trap is in it |
| A figure for the numerical Fisher being indefinite | 08-09 08.36.28 | **§10 item 11**, undecided |
| Automatic differentiation is an advance over the P2X2 build and must be said | 08-11 14.13.27 | **written**, `docs/manuscript-drafts/sections/06_methods.tex`, Reproducibility |
| The deliverable must carry the validation, not only the algorithm | 08-11 14.15.10 | **logged**, `../_program/carve_plan.md`, last section |
| bioRxiv first, then submit from the preprint | 08-11 14.05.44 | **logged**, `../_program/decisions.md` §1 |
| The pass over the figure roster and its supplements | 08-11 14.08.31, 14.13.27 | **closed 2026-08-11 by Luciano, nothing moves.** The audio is on the numbering from before the revision of 08-10, which is why its Figure 5 and Figure 6 read as one figure off against `decisions.md` "The figure set". Read it as a list of what should exist, not as a renumbering |
| The magnitude and anisotropy decomposition as a supplement | 08-11 14.08.31 | **already built, and it is a Figure 4 supplement** (Luciano, 2026-08-11). `figures/paper_both/figure_4_magnitude_anisotropy.Rmd`, built 2026-08-02, drawing `Figure_4_supplement_magnitude_anisotropy.pdf`; `04_results.tex:598,610` already quotes its medians. The audio's "me falta" is stale by nine days |
| The sample and correlation decomposition for the other members | 08-11 14.08.31 | open, but not a move: it exists for R against IR as `Figure_4_supplement_sample_corr.pdf` and the ask is to extend the roster, not to renumber |

## 10. What remains to implement

1. **[DESCOPED 2026-08-04, by Luciano.]** ~~The closed-form simulator checks.~~ Two judgements, both
   his. First, that the shared-specification argument is close to a technicality that nobody in the
   channel community will raise, and that for a two-state scheme the rate matrix, the conductance
   vector and the protocol are established objects that are nearly verifiable by inspection. That is
   right, and the remedy is now one sentence in Methods plus a bounding paragraph in the Discussion's
   demarcation, both written. Second, that the real way to settle it is property-based: check each
   object against the theorems it must satisfy, the semigroup and Chapman-Kolmogorov identities, the
   stationarity condition, the normalisation of the occupancy, and the moment identities of the
   emission. That is a piece of work with its own design and it is out of scope here. The Discussion
   now says so in those terms.

   **What survives as a live candidate, and why it is not ceremony.** The risk is not evenly spread.
   The kinetic specification is inspectable; the EMISSION specification is not, and it is where this
   repository has actually had errors: the `N*ms` interval-variance term went missing from the
   non-recursive path in a December refactor, and the Comm Biol validation was run with a defective
   `gvar_i`. Of the four checks, the one that covers that part is the non-stationary current
   covariance against Sigworth 1981 and Anderson & Stevens 1973, which is a published formula written
   by other people for exactly this quantity. One theorem is not the programme, but it is the one
   theorem that covers the part inspection cannot reach. **Undecided.**

   Original text, kept because the other three checks may still be wanted if the scheme ever grows:
   The
   code-verification half of the spine currently has no exactness anchor at all, and §11 says why that
   matters; these are the cheapest thing that closes it, they need no runs, and leaving them in the
   doubts while a figure column for a subplot sat at number one was a priority inversion. The mean
   relaxation against N·gᵀ·p0·exp(Qt); the non-stationary current covariance against Sigworth 1981 and
   Anderson & Stevens 1973 (both in `../../docs/bibliography/`); the occupancy against piQ = 0; the
   dwell times against their exponentials.
2. **Least squares at av = 0, and INR, with the non-`_G` emission**, so instance 3 of §2 has data.
3. **Figure 5 rebuilt as the decomposition**, with the estimation floor. No run needed.
4. **Figure 6B folded in as a panel**, with the K table and the two caption caveats of §7.
5. **The captions**: Figure 2 (the two numbers are different kinds; a 2D slice is a LOWER bound on the
   full anisotropy), Figure 3 (retold as the midpoint, naming the deception).
6. **Least squares with all variability assigned to instrumental noise** (10.04.51), one experiment, or
   left for a referee.
7. **The Münch 2022 positioning sentence**, mandatory: their noise is state-dependent open-channel
   noise so their error grows with noise; this paper's is additive instrumental so it falls. Without
   it the two validity maps read as contradictory.
8. **Governance:** `../_program/00_index.md` is still stamped 2026-07-20 and routes a three-paper world
   with paper 2 as a stub, when the merge of 2026-07-23 made this folder the fused paper; and the
   folder is still called `1_method` for something that is no longer the method paper.
9. **Derive it twice and show the two routes meet** (audio 2026-08-09 08.32.54). Present the member as
   a continuation of MacroR: a short synthesis of how MacroR is reached from the Bayesian side, then
   MR and IR derived from it taking MacroR as given; then the same destination reached by following
   the Kalman-filter logic; then the statement that the two formulations are equivalent, so a reader
   fluent in either one can enter through it. The precedent is his own: the 2007 paper presented one
   algorithm from the traditional route and from the Bayesian one. **Where it collides with what
   exists:** Appendix 1 already derives from MacroR, and the Kalman correspondence is currently a
   Discussion statement, quantified but not derived, which is the concession recorded in
   `../_program/decisions.md`. So this is not a new derivation, it is a second route through the
   existing one plus an explicit equivalence. **Undecided:** a second appendix, or an expansion of the
   Discussion correspondence. It is also the item most likely to cost length in a manuscript already
   under a length plan.
10. **Say why MR fails, as a mechanism and not as a measurement** (audio 2026-08-09 08.36.28). His
    words: it throws information away, because it marginalises before the interval's own measurement
    has been used. The paper currently characterises MR structurally (Table 1: conditioned on the
    state at the start of the window) and reports how badly it is calibrated, and it says what
    separates MR from IR (the gain). It does not say **why** in one sentence a reader can carry.
    **Checked against the canonical derivation before writing this**,
    `../../theory/macroir/notes/mr_vs_ir_from_macror.md`: MR marginalises `i_t`, the index of the
    state at the END of the window, and does it BEFORE running MacroR (:97, :112-113); IR conditions
    on the pair and marginalises `i_0` afterwards (:90); and the loss is named there, MR "conserva D²,
    pierde D", it marginalises precisely over the index whose profile carried the signal (:161, :166).
    **The wording trap, and the audio walks into it:** the audio says MR marginalises "al inicio",
    which reads as marginalising the initial state, and the index it marginalises is the final one.
    What is early is the timing, not the endpoint. Any sentence written for the paper has to say the
    timing. His follow-on question in the same audio, whether an MR that marginalised at the end would
    just be IR, is answered yes by the same note (:90), which is worth one clause because a reader
    will ask it.
11. **A display item for the numerical Fisher being indefinite** (audio 2026-08-09 08.36.28), as the
    thing that justifies the safeguards. The prose exists and is precise:
    `docs/manuscript-drafts/sections/06_methods.tex`, the second Fisher construction, says the
    differenced Fisher is widely indefinite replicate by replicate and that only pooling recovers a
    definite matrix. No figure shows it. **Undecided:** body or supplement, and on which data, since
    the battery that formed the differenced Fisher (`433ed13`) is the superseded lane and no number in
    the paper is taken from it, so a figure drawn on it needs a sentence saying why a superseded lane
    is the right place to show a defect of that lane.

    **The nearest thing that exists is not it, and checking that is what this entry is for.**
    `figures/paper_both/figure_4_gaussian_vs_numeric_fisher.Rmd` draws the analytic Gaussian Fisher
    against the differenced one in size and in shape, which is the agreement question, not the
    definiteness one; its own header calls it provisional and frames it as "how far apart are two
    positive definite matrices". Its three limits are the same ones any definiteness figure will hit:
    the differenced Fisher exists only in the non-`_G` runs, so least squares and `INR` are absent
    entirely, and the `NMR` column predates the restored `N*ms` term and is superseded. So item 11 is
    a new figure on a lane the paper otherwise does not use, and that, rather than the drawing, is
    what has to be decided first.

## 11. The doubts

**The specification is not tested, and this is now the largest one.** The identities compare the
likelihood code against the simulator code, and both are built from ONE model specification: Q from
theta, the mapping to rates, the conductance vector, the protocol, the sampling times. Everything
downstream of that split is tested; the specification is not. It is not a coincidence argument, it is
a shared dependency, and the failure mode that survives everything else is a specification error that
acts like a REPARAMETRISATION: records still look reasonable, the result structure is preserved, and
theta_sim is defined by the same code.

**One of its two supports has been withdrawn.** The first draft claimed the channel-count boundary was
an external check on its position; §5 retracts that. What is left is the structure-of-results argument
(a gross error would not produce a ladder that orders as theory predicts across five members and four
decades), which is soft, plus the one-sided exclusion, which is a bound and not a position. So §10
item 1 is not optional.

**The bias claim is weaker than this file previously stated**, measured 2026-08-03 over the source
data: median |bias| in log10 units is 0.0019 (ILSE, k_on), 0.0119 (NR), 0.0119 (R), 0.0009 (IR), and
the CI-aware `Bconf` median is zero for every member. The ORDERING follows the ladder, NR and R being
about 13x worse than IR on k_on, but the magnitudes are a few per cent in the median with large values
only in isolated cells. "The non-interval members are biased" is therefore true as an ordering and
false as a general property, and the bias cannot carry an argument against a sweep bootstrap.

**What produces the cross-time factor is not known.** It is the actual boundary and nothing in the
derivation set models it. This is the honest content of "we do not understand where MacroIR fails",
and it is sharper than before because the location of the ignorance is identified.

**The long-interval premise appears to be wrong** and has not been re-confirmed: the distortion is not
monotone in interval, and the genuine monotone effect is a 6.4 per cent over-prediction of the
predictive variance confined to ten channels. Numbers in
`analysis_figure_5_breakdown_criterion.md` §11, all marked TO CONFIRM.

**Which functional of the spectrum the paper claims about is a DECISION, not a doubt.** Anisotropy and
lambda_max disagree by 13x on the boundary in the same cell. The principled resolution is a mapping
from claim to functional rather than one choice: lambda_max answers "how wrong is my worst direction",
which is what matters if you report per-parameter intervals; the magnitude/anisotropy pair answers
"how wrong, and can one factor repair it", which is beat 2; log-det answers volume if a joint region is
ever reported. What is not allowed is one claim using one and its neighbour using another silently.

**Two roads deliberately not taken**, recorded with their costs in
`analysis_figure_5_breakdown_criterion.md` §9 so the decision is not re-litigated. micro_IR as the
exact reference: it already exists, does exact Bayes on the occupancy, and has runs at N_ch 5, 10 and
20, which is where the breakdown lives, so **cost was not the reason** and that file records the real
one, which is scope. And a higher-cumulant emission, which would have made the per-sample part
predicted rather than fitted.

**Reproducibility.** seed = 0 is the `random_device` sentinel and the resolved seed is not logged. The
two time configs additionally use different seeds despite a comment claiming otherwise, so per-record
pairing across algorithms through them is invalid.

**And the one the author names himself:** the limits look like a balance that has to come out just
right, which suggests a conserved quantity or a comparable conceptual tool that has not been
identified. Stated as a conjecture, in those words.
