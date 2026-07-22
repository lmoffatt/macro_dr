# Paper 1 (method): plan

> Updated: 2026-07-20. This is **paper 1**, the method paper of the three-paper program
> (`../_program/program.md`). It owns the interval-likelihood closure: given that you will compute a
> likelihood, what must it condition on? Roster `R`, `MR`, `VR`, `IR`.
>
> **Program-level framing is not owned here.** The three-paper map, the two-level structure (root
> question + ladder), the noise bands and the usage-map deliverable across all bands live in
> `../_program/`. This file cites them; it does not restate them. Where a section below was
> program-wide before the 2026-07-20 split, it now points out.
>
> Open items for the author are marked **[Q]**. Supersedes `00_master_plan.md` (retired stub).

## The argument

If we had an exact method, we could just demonstrate its validity mathematically.

On the other hand, if we have an approximation, a method we know is not exact, what do we do?

That's the question I was wondering myself.

The response I arrived at: find the situations where the method is very precise and the situations where it is not, and quantify both.

So the idea is to prove that the method works where it is supposed to work and breaks where it is supposed to break.

So, what is the method trying to solve?

The problem consists in approximating the likelihood function of time-averaged observations of a given Markovian process. We can simulate the unfolding of this Markovian process and its observations with very good accuracy; the inverse, going from an observed experiment to a likelihood, is the hard direction, and that asymmetry is exactly what lets us use the simulation as ground truth to judge any proposed approximation.

We have alternative approximations of the likelihood function and we want to see if they are good. After some introspection and reading the literature, we conclude that we can test whether a likelihood function faithfully represents a simulation process by a three-step analysis:

1. the normalized residuals should have mean zero, variance one, and no temporal autocorrelation;
2. the expected score should be zero, evaluated at the true parameters;
3. the variance of the score should equal the negative Hessian of the log-likelihood, evaluated at the optimum where the score is zero.

For a misspecified approximation, tests 2 and 3 land at different points, and that gap is the whole story of the bias.

Then we accept that we are dealing not with true likelihoods but with approximations, and we ask whether we can find a procedure to improve the approximation. We found such a procedure: an expression for the corrected bias of the MLE parameters and one for their corrected variance. Near the maximum, the same distortion propagates to the Bayesian evidence through two channels, a volume correction ½ log det C and an effective-sample rescaling α⋆ = p / tr C, so the correction is what keeps model comparison valid under an approximate likelihood. (Stated here as motivation; the derivation is deferred to a later program component.)

Now, what is the nature of the main algorithm we are presenting? The algorithm represents the predictive (prior) state of the channel population, the density that assigns a probability to each combination of channel states, by a multivariate Gaussian. That density is exactly multinomial when channels are independent (equivalently, the maximum-entropy closure given only the mean occupancies), and the Gaussian is its large-N limit. The higher moments the closure discards do not vanish; they reappear later as the correlation distortion the diagnostics measure.

## 0) Goal in one sentence

Show, for the recursive likelihood family in the **gating-dominated regime**, **what the likelihood
must condition on** to report its own uncertainty honestly: the interval-averaged conductance
conditioned on both interval endpoints (`IR`), and not less. The result is a within-family validity
map over `R`, `MR`, `VR`, `IR`, with the mechanism of each failure attributed, in a minimal two-state
channel model, in simulation.

Paper 1's place in the wider deliverable — the full usage map across least squares and the whole noise
axis — is `../_program/program.md`. Paper 1 is the middle band of that map: it establishes which
likelihood to use *when a likelihood is needed*, and hands the "is a likelihood needed at all"
question to paper 2 and the few-channel boundary to paper 3.

The frame that survives from the earlier single-paper plan: this is not "here is MacroIR", it is a
validation machinery plus the within-family map it produces, MacroIR's own failures included.

## 1) Core conceptual frame

**The two-level structure (root question + endpoint ladder) is owned by `../_program/program.md` §1.**
Paper 1 lives entirely on the lower level, the ladder: it takes as given that gating fluctuations are
being modelled and asks only what the Gaussian is conditioned on. `VR` splits the MR→IR step of that
ladder into two, isolating variance from gain (`../_program/decisions.md` §2).

Every method paper 1 studies makes **two Gaussian approximations**:

1. **Macro approximation.** The multinomial distribution of channel occupancies is replaced by a
   multivariate Gaussian. Valid for large N_ch; degrades for few channels (multinomial regime). **This
   is paper 3's subject, not paper 1's.** Paper 1 meets it only at its N_ch = 10 floor, where the one
   micro anchor cell attributes it (§2, `decisions.md`).
2. **Interval-likelihood approximation.** The conductance distribution over one measurement interval
   is replaced by a Gaussian, when the exact object is a stochastic-telegraph / Poisson-like
   distribution. Degrades for intervals much shorter than the relaxation time (telegraphic regime).
   **This is paper 1's subject.**

Three regimes result: multinomial (few channels, paper 3), telegraphic (very short intervals,
paper 1's failure edge), and Gaussian (many channels, moderate intervals, where IR is ideal). Because
IR is an approximation it must fail somewhere; paper 1 locates that failure and shows it is the
predicted degradation of the interval closure, not a coding artifact.

### 1a) The regime paper 1 occupies

**The noise bands and their two crossovers are owned by `../_program/axes.md` §2a.** Paper 1 occupies
**bands A and B**, the gating-dominated regime, where the answer to the root question is yes and the
open question is which rung of the ladder suffices. It does **not** run band C (instrumental noise
above the total gating noise); that is where a likelihood stops being needed, which is paper 2.

Paper 1 must **state that scope in band terms** and name the companion papers for the rest
(`../_program/program.md` §7). That stated scope, not the presence of a least-squares arm, is what
keeps the paper from reading as though it lived only where IR wins by construction.

## 2) Scope

### In scope
- Minimal **two-state** channel model (closed-open, `scheme_CO`), single K_on / K_off; **macroscopic**
  interval-averaged currents; **non-stationary** protocol. (All cross-paper; `../_program/decisions.md`
  §2.)
- **Roster `R`, `MR`, `VR`, `IR`** — the recursive-family ladder, in the gating-dominated regime.
- Control variables: N_ch (10 … 10⁴), interval Δ/τ, and instrumental noise **within bands A–B**.
- **Likelihood-only.** MLE / Gauss-Newton local maximum kept only to obtain the empirical parameter
  covariance.
- **Gaussian Fisher** as the distortion anchor (`../_program/machinery.md` §3).
- One **micro attribution anchor** (`micro_IR`, N_ch 10, noise 0.1) to assign IR's floor degradation
  to the occupancy closure; one or two annotated cells, not a column (`decisions.md`).

### Out of scope (owned by another paper or a later component)
- `LSE` and `NR` → **paper 2**. `NMR` → dropped from the program. Micro as a subject, and the
  multinomial boundary → **paper 3**.
- Band C (instrumental-noise-dominated) → paper 2.
- More than two states, the stationary regime, experimental data → later components.
- The posterior information-distortion framework and full model comparison → later; the likelihood-side
  evidence correction stays as motivation only (`../_program/machinery.md` §10).
- Taylor variance-correction variants (IRT, MRT) — cut. **Note the name collision:** `VR` is *not* a
  Taylor variant; the engine flag is `taylor_variance_correction` and Methods must say the two are
  unrelated (`../_program/program.md` §9).

## 3) Diagnostics

**Owned by `../_program/machinery.md`** (definitions, sign conventions, thresholds, the decomposition
identity, the hazards). Paper 1's Diagnostics *section* — how the machinery is presented and argued —
is `diagnostics.md`. Do not restate the definitions here.

## 4) The within-family validity map (the former ranking)

**The deliverable is a recommendation per regime, not a ranking.** For paper 1 the regime is bands A–B
and the roster is the ladder. The table below is the **band-A column**, the only band fully measured;
bands B and the boundary come from the runs.

| Algo | Bias | Variance distortion (band A) | Standing |
|---|---|---|---|
| **IR** | none | to ~1 in Gaussian regime; conservative (~0.5) in the few-channel corner | Calibrated; the default when gating noise dominates |
| **VR** | none apparent | **over-confident 1.8–2.2, worse than MR** (one cell, Gaussian anchor, 2026-07-22) | The control, and it fired: removing the variance without the boundary gain degrades calibration, so the gain is what recovers it |
| **MR** | none | over-confident 1.5–2.1 (drops the boundary cross-covariance N·γᵀΣγ that IR keeps) | Cautionary intermediate; the window axis is non-monotone, MR worse than R |
| **R** | none | residual distortion, size disputed (see `decisions/D-4_ranking_verdict.md`) | Usable; correlation distortion at short Δ |

Mechanism of MR's over-confidence: `theory/macroir/docs/Macro_MRT/macromrt_macromrt_paper_section.md`
and `theory/macroir/notes/gvar_i_overcount_audit.md` (use the mechanism, not the inverted May-2026
verdict). The audited finding: MR's error is genuine misspecification, not a code defect — it dumps
the end-state spread into the observation variance instead of resolving it through the gain. **`VR` is
the experiment that measures this**, and it is why the roster includes it.

**Open, and only the runs close them.** ~~VR's actual sign~~ **answered 2026-07-22: over-confident,
and more than MR** (`decisions.md`), so the mechanism claim holds at the headline cell and now needs
the grid to be stated at regime scope. R's corner behaviour against the two contested cells of
`decisions/D-4_ranking_verdict.md`; where the measured A/B boundary falls against the prediction.
**Do not pre-write these.** The previous table asserted verdicts the data then contradicted in those
two cells, which is how the four-copy drift started. NR, NMR and
LSE rows are deliberately gone: they are other papers.

## 5) Figures

**The arc is owned by `results.md`** (was `results_plan.md`; `../_program/00_index.md` routes it).
Paper 1's intended set, minus the program-wide map that belongs to paper 2:

- **Fig 1: Mechanism.** One filter step across R, MR, VR, IR (prior, observation, update, logL),
  showing where each approximation enters. Anchor-independent.
- **Fig 2: Recovery clouds.** MLE clouds with three ellipses (empirical, Gaussian-Fisher,
  sandwich-corrected), per method.
- **Fig 3: Time-resolved calibration.** logL gap, R² residual, score bias, per-interval and
  accumulated J_t/F_t, score autocorrelation.
- **Fig 4: Fisher profiles.** Per-step Fisher; the information about N_ch, k_on and i falls to zero
  once the open-channel count relaxes, while k_off stays informative through the decay.
- **Fig 5: The within-family validity map.** Over N_ch × noise **within bands A–B**, with interval
  carried inside each cell and the predicted A/B boundary drawn on it. (The full A–C map with LSE is
  paper 2's figure.) The old **N_ch × K_off** axes are dead: no K_off axis exists in the scripts
  (`figure_3_mle_G.macroir` hard-fixes `off = 100`) and building one was costed at ~35,000 CPU-hours.
- **Fig 6: Distortion decomposition.** Correlation vs sample distortion, and — the new job for this
  figure — the two steps MR → VR → IR made concrete: which step is variance, which is gain.
- **Supplement.** Two-path reconstruction (`IDM = K·CDM·Kᵀ`, `../_program/machinery.md` §5); the
  micro attribution anchor; band schematic; S1 score-mean; S3 corner; S4 bias and autocorrelation.

**The figure count is open** (`01_writing_plan.md` gates at six; carried to `../_program/program.md`
§9). Dropping NR/NMR/LSE columns and adding VR does not change the count; it changes what fills each
panel.

## 6) Data and reproducibility anchor

- Definitive figures anchor on the **Gaussian Fisher**; the all-algorithm Gaussian rerun becomes the
  reference basis. Provenance, hashes and the seed sentinel: `../_program/provenance.md`.
- Grid: **N_ch × noise within bands A–B**, interval inside every cell; `scheme_CO` only; Taylor
  variants dropped. K_off axis dead (§5).
- **VR must be run.** If the residual variance form is computable from existing Qdtm fields as the
  (inverted-verdict) note claims, it is a small change, not a new algorithm — **verify against current
  code before planning on it** (`../_program/decisions.md` references; `feedback_verify_dont_assume`).
- **The n_sims-uniformity hazard** is program-level and owned by `../_program/machinery.md` §6.1. For
  paper 1 specifically: do not pool the band-A cells (n_sims 10⁴) with anything at 1000, and do not
  pair the micro anchor (n_sims varies) except at its 10⁴ cell.

## 7) Manuscript target

Head manuscript `docs/manuscript-drafts/elife_paper.tex`. Rewrite the abstract and Results spine to the
within-family where-it-fails frame; keep the two-state theory; restage the MacroIR ≈ time-augmented
integrated Kalman link as a Discussion hypothesis (verified to ~1e-8 in the prior-art map), not a
Theory fact. Theory sources:
- Algorithm: `theory/macroir/docs/Macro_IR/macroir_macroir_paper_section.md`, `macroir_derivation.tex`.
- Distortion machinery: the supplements under `theory/macroir/docs/Likelihood_Information_Distortion/`,
  `Gaussian_Fisher_Distortion_Family.md`, `sample_correlation_distortion_analysis.md`.
- Results seed: `analysis_figure_S1_score_mean.md`.

## 8) Open decisions

Paper 1's own. Cross-paper decisions (venue, the N_ch partition, VR's name, one-repo-or-three) live in
`../_program/program.md` §9 and `../_program/decisions.md` §5.

**These are `Q-n`, not `D-n` (relabelled 2026-07-21).** The `D-n` labels belong to the
manuscript-production register: one brief per decision in `decisions/`, listed in `01_writing_plan.md`
§3. Both lists used `D-n` until now, so `D-0`, `D-1`, `D-3`, `D-4` and `D-5` each named two unrelated
things — inside this file, `D-4` in §4 above is the ranking brief while the old `D-4` here was the
threshold question. The old label is kept in parentheses so earlier references still resolve.

- **Q-1 [Q]** (was D-0) Confirm **NR** is out of paper 1 (moved to paper 2). Default: out.
- **Q-2 [Q]** (was D-1) MR and VR in main text or supplement? "Strawman" is retired; re-decide on the
  map footing — a method here is one with a domain, and VR's domain is the mechanism proof.
- **Q-3 [Q]** (was D-3) Placement of the Fisher-to-zero result: main Fig 4 as drawn, or supplement.
  Author: "me voló la cabeza" but "does not go in the abstract".
- **Q-4 [Q]** (was D-4, the threshold one) "Valid" thresholds on the map (candidate distortion < 1.1);
  the three-way threshold conflict is in `../_program/machinery.md` §8.
- **Q-5 [Q]** (was D-5) Anchor point for bias/distortion. Consensus: optimum / θ_pool; θ_sim exposes
  the bias.
- **Q-6 [Q]** (was D-7) Language of the planning docs: English (current) or Spanish. (Program-wide,
  but unrouted.)

## 9) Next actions

A pointer, not a list (a frozen plan cannot host a ledger that mutates weekly).

- **Writing order, section gates, definition of done:** `01_writing_plan.md`.
- **Engine work that must land in the frozen commit:** `../_program/carve_plan.md`, Freeze preconditions.
- **Who owns which topic across the four layers:** `../_program/00_index.md`.
