# MacroIR prior-art map

**Two parts.** Part I (this half) defends the **algorithm** (MacroIR, the integrated-measurement
filter). Part II (appended 2026-07-14, at the bottom) defends the **diagnostic** — the information
distortion matrix `C = H^{-1/2} J H^{-1/2}`, its decomposition, and the MLE/evidence corrections —
which is what the eLife paper now actually claims. **If you are positioning the eLife paper, read
Part II.** Part I is for the Comm Biol method.

Reference for positioning the MacroIR paper (exact integrated-measurement filter for
macroscopic ion-channel currents). Built from two adversarial literature workflows
(2026-06-26) plus a numerical equivalence check. Each entry gives the citation, the
local PDF (or paywall status), what the work actually does, its observation model, and
its **role** for our paper: `baseline-to-beat`, `concede` (device we did not invent),
`distinguish` (near-miss to cite and separate from), or `machinery`.

## Bottom line (what the searches settled)

1. **The channel domain cell is empty.** Every macroscopic / ensemble channel inference
   method models the observation as an *instantaneous* point sample of the open-channel
   count, `y_t = H n_t + noise`, never as the within-interval time-average
   `ȳ = (1/Δ) ∫ I dt`. Finite integration time is treated as a misspecification to be
   robust to, not modeled. MacroIR is genuinely first-in-domain.

2. **The device is not new.** The integrated-/averaged-measurement Kalman filter (augment
   the state with the running integral, propagate by a van Loan block matrix exponential,
   condition the endpoint on the observed integral) is textbook in control/econometrics
   and already used in *ensemble biophysics* (Folia & Rattray 2018, gene expression) and,
   with van Loan explicitly, in astronomy (Rubenzahl et al. 2026). It is a named current
   subfield (Integrated Measurement Kalman Filter, Yaghoobi & Särkkä 2024).

3. **The "exact CTMC vs their LNA" edge is weak and must not be oversold.** For the
   *linear* channel dynamics the linear-noise approximation (LNA) is itself exact in the
   first two moments. We verified numerically that the augmented filter built on the
   chemical-Langevin (LNA) diffusion equals MacroIR's exact-CTMC result to 1e-8
   (`claude/scripts/verify_IR_vs_augmentation.py`, against Monte Carlo). So a Folia &
   Rattray style LNA augmentation, specialized to channels, would coincide with MacroIR.
   The exact-CTMC realization is a **derivation / efficiency / channel-native** point, not
   a numerical-accuracy claim over LNA.

4. **No killer prior art.** The closest CTMC-with-integrated-observation filter (Bäuerle
   et al. 2014/2015) is a *single* telegraph chain (not an ensemble), exact only for two
   states, with no augmentation and no van Loan, applied to finance and queues. It is a
   cite-and-distinguish, not a refutation.

5. **The k² meta-state ⇄ (k+1) augmented-state equivalence is apparently unstated** as
   such, though both halves are standard and a matrix-analytic reviewer may call the
   reduction folklore. Frame it as a unification, openly conceding the ingredients.

**Defensible novelty:** (i) first integrated-measurement (time-averaged-observation)
treatment for macroscopic ion-channel currents; (ii) the exact, efficient, channel-native
realization for the many-channel finite-state continuous-time Markov chain (CTMC), where
the only prior CTMC integrated-observation filter, Bäuerle et al., was exact only for two
states; (iii) the empirical demonstration that the field's instantaneous methods are
misspecified (miscalibrated uncertainty, biased model selection). This is the same
contribution shape as Münch et al. 2022 (who brought the textbook Kalman filter to
channels): bring the established tool, show it is needed.

---

## A. In-domain channel baselines — the gap (all instantaneous; `baseline-to-beat`)

### Münch, Paul, Schmauder & Benndorf 2022 — the direct competitor
- eLife 11:e62714. doi:10.7554/eLife.62714. PDF: `Munch_2022_Bayesian_Kalman_IonChannels_eLife.pdf` (peer-review letter in `reviews/`).
- Generalized Bayesian / Kalman filter for macroscopic current and confocal patch-clamp
  fluorometry. Observation is **instantaneous** (Eq. 4): `y_t = H n_t + ν_t ~ N(·|H n_t, Σ_t)`,
  state-dependent variance `σ_m² + n_open σ_op²`, state propagated by a plain `exp(KΔt)`.
- State augmentation is mentioned **once**, only as an option for **colored** noise (citing
  Chang 2014); not adopted. Finite analog filtering / integration time is treated as a
  misspecification to be robust to ("higher analyzing frequency on minimally filtered data";
  finite integration time "induce[s] an intrinsic time scale ... not represented in the
  algorithms"). Notably it *cites* the ensemble integrated-Kalman/LNA literature (Komorowski;
  Finkenstädt; Folia & Rattray 2018) yet declines the integrated-observation construction.

### Moffatt 2007 — the baseline MacroIR extends, and the open problem it closes
- Biophys. J. 93(1):74-91. doi:10.1529/biophysj.106.101212. PDF: `Moffatt_2007_Macroscopic_Current_Fluctuations_BiophysJ.pdf`.
- Recursive Gaussian filter for macroscopic current capturing local time-correlation of
  intrinsic channel noise; instantaneous Gaussian readout of the current given occupancies;
  `exp(QΔt)` by eigendecomposition. No augmentation, no van Loan.
- **Names the gap explicitly:** "the measurements are considered to be instantaneous, where
  electrophysiological measurements always result from some time-averaging ... a formulation
  of the algorithms that can deal with time-averaged signals is therefore fundamental."
  MacroIR closes the author's own stated open problem — the strongest, least-attackable framing.

### Milescu, Akk & Sachs 2005
- Biophys. J. 88(4):2494-2515. doi:10.1529/biophysj.104.053256. PDF: `Milescu_2005_MLE_Macroscopic_Currents_BiophysJ.pdf` (already in repo).
- Macroscopic-current maximum likelihood. Gaussian in the **point** current `I_t` with mean
  and variance at the instantaneous occupancy; `A_t = exp(QΔt)` by spectral decomposition.
  Its van Loan citation is the Moler & Van Loan matrix-exponential *computation* survey, not
  the integrals method. Confirms instantaneous is the field default.

### Celentano & Hawkes 2004
- Biophys. J. 87(1):276-294. doi:10.1529/biophysj.103.036632. PDF: `CelentanoHawkes_2004_Covariance_Matrix_Kinetic_Fitting_BiophysJ.pdf`.
- Direct fit of kinetic parameters to macroscopic currents using the full inter-time-point
  covariance matrix of intrinsic noise (O(n³)). The covariance is **between time points**
  (autocorrelation of the process), not a model of the measurement's within-interval
  time-average. Observation still instantaneous. No augmentation/van Loan.

### Stepanyuk et al. 2011 / 2014 — faster successors to Celentano-Hawkes
- **CORRECTED 2026-07-14.** The old entry listed "Stepanyuk 2012, PLoS ONE 7(4):e35208". That DOI is
  **not** Stepanyuk: it is **Wang, Xiao, Zeng, Yao, Yuchi & Ding (2012)**, and it is a **least-squares**
  paper, not a likelihood successor (see the separate entry below). The local PDF has been renamed
  from `Stepanyuk_2012_...` to `Wang_2012_...`. The real Stepanyuk lineage is 2011 + 2014.
- 2011: PLoS ONE 6(12):e29731. doi:10.1371/journal.pone.0029731. **PDF not in repo.**
- 2014: Front. Cell. Neurosci. 8:303. doi:10.3389/fncel.2014.00303. PDF: `Stepanyuk_2014_Synaptic_Receptors_Macroscopic_FrontCellNeurosci.pdf`.
- Macroscopic current as a non-stationary Gaussian process (ML NSFA); exploit **semiseparability** of
  the autocovariance for near-linear cost. Mathematically a semiseparable-covariance Gaussian
  likelihood is filter-equivalent, but they never frame it as filtering and never mention Kalman; the
  covariance is the model's **prior** covariance, never conditioned on the observed trace. Sampling
  instants are point spacings, not averaging windows. No within-interval integration, no augmentation,
  no van Loan.
- Useful quotes on how unused this whole family is: *"PS NSFA [is] the most frequently used noise
  analysis approach"*; *"the unitary current is virtually the only parameter that can be reliably
  obtained from this type of analysis"*; **"To the best of our knowledge, kinetic rates have never
  been estimated for any synaptic receptors in their intrinsic environment."**

### Wang et al. 2012 — a macroscopic-fitting paper that declines the likelihood (ADDED 2026-07-14)
- PLoS ONE 7(4):e35208. doi:10.1371/journal.pone.0035208. PDF: `Wang_2012_Optimal_Estimation_Macroscopic_Currents_PLoSONE.pdf`.
- **Least squares on the deterministic mean current** (particle-swarm optimizer). No likelihood, no
  variance term, no gating fluctuation anywhere. Cites Milescu 2005 and **explicitly opts out of ML on
  computational grounds**.
- Role: direct evidence for the Introduction that even a dedicated macroscopic-fitting paper, seven
  years after the likelihood method existed, chose least squares on the mean.

### Del Core & Mirams 2025 — the newest competitor, still instantaneous
- Phil. Trans. R. Soc. A 383(2292):20240224. doi:10.1098/rsta.2024.0224. PDF: `DelCore_Mirams_2025_Stochastic_IonChannel_WholeCell_PhilTransA.pdf`.
- Parameter inference for stochastic ion-channel gating from whole-cell voltage-clamp data;
  LNA / moment state-space filter with a **point** Ohmic measurement `y_t = G_t x_t + ε`.
- **CORRECTED 2026-07-14 (verified against the PDF).** An earlier version of this map said they do
  not cite Moffatt 2007 or Münch 2022. **They do.** Refs 30-32 are Moffatt 2007, Milescu 2005,
  Münch 2022, and they name MacroR in the text. They do **not** cite Folia & Rattray or van Loan
  (that part of the old note was right). Do not repeat the old claim in a rebuttal.
- Three quotable gifts from their §on filtering:
  - **Third-party endorsement of the recursive point:** on Moffatt 2007, *"The authors proved that
    their method provides better estimates than approaches that do not consider statistical
    dependence of successive measurements [31]"* (i.e. beats Milescu 2005 precisely on time
    correlation).
  - **They argue against filtering on cost:** *"approaches that are based on filtering techniques
    are computationally expensive because an integration step of the differential moment equations
    (DMEs), and the corresponding updates in the correction step, are computed between every
    consecutive time points where the measurements are collected."* Our paper answers what you lose
    by skipping the filter.
  - **They declare the channel-number problem still open in 2025:** *"how to best infer the
    relationship between the noise observed in the ionic measurements and the total number of
    channels in the cell membrane contributing to gating is still unclear."* Directly relevant to
    our Fisher-information result on N_ch identifiability.
- Confirms the integrated-observation niche is still open in 2025 (their observation is instantaneous).

---

## B. The augmented / integrated-measurement Kalman filter — full bibliography (`concede`)

This is the column you must cite openly so the honesty is bulletproof. The filter MacroIR
re-derives is the **integrated-measurement Kalman filter via state augmentation**: append the
running integral `s = ∫ I dt` to the state, propagate `[x; s]` with a van Loan block matrix
exponential, condition the endpoint on the observed integral via `Cov(x_end, s)`. None of this
is ours. Cite the foundations (B1), the integrated/flow-observation pattern specifically (B2),
and the colored-noise augmentation that the channel field already uses (B3).

### B1. Foundations of Kalman filtering and state augmentation

- **Kalman, R. E. (1960).** A new approach to linear filtering and prediction problems.
  *Trans. ASME J. Basic Eng.* 82(1):35-45. doi:10.1115/1.3662552. PDF:
  `Kalman_1960_New_Approach_Linear_Filtering_JBasicEng.pdf`. The original filter; Münch 2022
  cites it as the basis. Cite as the root.
- **Jazwinski, A. H. (1970).** *Stochastic Processes and Filtering Theory.* Academic Press
  (Dover reprint 2007), ISBN 978-0486462745. **Book — not downloaded.** The textbook source for
  **state augmentation**, including correlated/colored noise and integrated/averaged measurements
  (append a state, condition the endpoint on the observed integral). The canonical augmentation cite.
- **Anderson, B. D. O. & Moore, J. B. (1979).** *Optimal Filtering.* Prentice-Hall (Dover reprint
  2005), ISBN 978-0486439389. PDF: `Anderson_Moore_1979_Optimal_Filtering.pdf`. The minimum-variance
  filter reference Münch 2022 cites for KF optimality; also covers augmentation and the
  continuous-discrete filter.
- **Särkkä, S. (2013).** *Bayesian Filtering and Smoothing.* Cambridge University Press, ISBN
  978-1107619289. PDF: `Sarkka_2013_Bayesian_Filtering_and_Smoothing_CUP.pdf` (free from the author).
  Modern textbook: continuous-discrete filtering, matrix-fraction decomposition, the
  Gaussian-process/SDE state-space view. 2nd ed.: **Särkkä, S. & Svensson, L. (2023)**, CUP, ISBN
  978-1108926645. The modern reference for the augmented/continuous-discrete machinery.
- **Van Loan, C. F. (1978).** Computing integrals involving the matrix exponential. *IEEE Trans.
  Automatic Control* 23(3):395-404. doi:10.1109/TAC.1978.1101743. **Paywalled (IEEE) — not
  downloaded.** The block-triangular matrix-exponential device for the integral of `exp(At)` and the
  process-noise covariance integral. The numerical primitive the augmentation rests on; cite by name.

### B2. Integrated / averaged / flow observation handled by augmentation (the exact pattern)

- **Zadrozny, P. (1988).** Gaussian likelihood of continuous-time ARMAX models when data are stocks
  and flows at different frequencies. *Econometric Theory* 4(1):108-124. doi:10.1017/S0266466600011890.
  **Paywalled (CUP) — not downloaded.** The cleanest econometric statement: **flow** variables
  (integrals over the period) handled by augmenting the state with the integral coordinate and running
  a Kalman filter on the exact-discretized SDE. Route 2, linear-Gaussian, decades old.
- **Harvey, A. C. (1989).** *Forecasting, Structural Time Series Models and the Kalman Filter.* CUP,
  §6.3 ("cumulator variable"), ISBN 978-0521405737. **Book — not downloaded.** The cumulator (running
  integral) variable for flow data. With **Phillips (1959/1978); Bergstrom (1984); Chambers** — the
  flow-vs-stock temporal-aggregation lineage. Cite to credit the general method's age.
- **Yaghoobi, F. & Särkkä, S. (2024).** Parallel state estimation for systems with integrated
  measurements. arXiv:2410.00627. PDF: `Yaghoobi_Sarkka_2024_Parallel_Integrated_Measurement_IMKF.pdf`.
  Explicit **"Integrated Measurement Kalman Filter (IMKF)"** for Slow-Rate inTegrated Measurement
  (SRTM) linear-Gaussian models. Cite to show the integrated-measurement filter is a current, *named*
  subfield, not exotic engineering.
- **Movassagh, Sh., Vahidi, A., Moshiri, B. et al.** Integrated-measurement Kalman filter for a process
  with slow-rate integrated (sample-and-accumulate) sensors (chemical-process control; IEEE).
  *Verify exact title/year/venue before citing.* A second explicit "integrated measurement Kalman
  filter" usage, linear-Gaussian, non-biophysics — corroborates the device is standard.
- **Rubenzahl, R. A., Hattori, S., Särkkä, S., Farr, B. F., Luhn, J. & Bedell, M. (2026).** Scalable
  Gaussian processes for integrated and overlapping measurements via augmented state-space models.
  *Astron. J.* doi:10.3847/1538-3881/ae4d0b (arXiv:2601.02527). PDF:
  `Rubenzahl_2026_Integrated_Overlapping_Measurements_Augmented_SSM_AJ.pdf`. Closest by **exact
  construction**: augments the SSM with an integral state that resets at exposure start and is observed
  as a time-average at exposure end, **using the van Loan (1978) block exponential explicitly**.
  Linear-Gaussian GP/SDE state, astronomy. Cite so you do **not** over-claim the
  van-Loan-on-augmented-integral pattern.
- **Folia, M. M. & Rattray, M. (2018).** Trajectory inference and parameter estimation in stochastic
  models with temporally aggregated data. *Stat. Comput.* 28:1053-1072. doi:10.1007/s11222-017-9779-x.
  PDF: `Folia_Rattray_2018_Trajectory_Inference_Aggregated_StatComput.pdf`. **The closest precedent in
  ensemble biophysics:** integrated/temporally-aggregated observation `y = P ∫X du + noise` handled by
  augmenting the LNA state with the running integral inside a continuous-discrete Kalman filter
  (single-cell luciferase integrated over minutes). LNA via moment-ODEs, no van Loan, gene expression.
  Note (Bottom line #3): for the *linear* channel dynamics this LNA approach coincides numerically with
  MacroIR. **Cite as the methodological precedent for integrated observation in biophysics.**
- **Calderazzo, S., Brancaccio, M. & Finkenstädt, B. (2019).** Filtering and inference for stochastic
  oscillators with distributed delays. *Bioinformatics* 35(8):1380-1388. doi:10.1093/bioinformatics/bty782.
  PDF: `Calderazzo_2019_Filtering_Inference_Oscillatory_LNA_Bioinformatics.pdf`. Integrated
  camera-exposure observation on the same LNA framework, built on Folia & Rattray.

### B3. Colored-noise augmentation — the only augmentation the channel field uses (context)

- **Chang, G. (2014).** On Kalman filter for linear system with colored measurement noise. *J. Geodesy*
  88(12):1163-1170. doi:10.1007/s00190-014-0751-7. **Paywalled (Springer) — not downloaded.** State
  augmentation / measurement-time-differencing to whiten **colored** noise. The one augmentation Münch
  2022 cites — and it is for noise coloring, **not** the time-average. Cite to show the channel field
  invokes augmentation only for colored noise, never for the integrated observation.

---

## C. The CTMC integrated-observation near-misses — cite and distinguish (`distinguish`)

### Blackwell 2018/2019 — integrated continuous-time HMMs (ADDED 2026-07-14)
- **Blackwell, P. G.** *Integrated continuous-time hidden Markov models.* arXiv:1807.11907 (stat.ME).
  **Not downloaded.** Found by the 2026-07-14 literature search; was missing from this map.
- Makes **the same conceptual move as MacroIR**, and says so in its own abstract: *"a new class of
  integrated continuous-time hidden Markov models in which each observation depends on the underlying
  state of the process over the whole interval since the previous observation, not only on its current
  state ... under appropriate conditioning, a model in this class can be regarded as a conventional
  hidden Markov model, enabling use of the Forward Algorithm for efficient evaluation of its
  likelihood."* Condition on the interval, recover a filter recursion.
- Why it does **not** sink the claim: it is **movement ecology**, a **single** chain (a switching
  diffusion for an animal's behavioural state), **not** a many-channel ensemble. No channel count, no
  open-channel variance, no Kronecker-sum collapse, no van Loan augmentation, no biophysics.
- **Why it must be cited anyway:** a statistician reviewer may well know it, and it is the closest
  statement of the general device for CTMC-with-integrated-observation outside our domain. Cite and
  distinguish, alongside Bäuerle.

### Bäuerle, Gilitschenski & Hanebeck 2014 / 2015 — the strongest near-miss
- arXiv:1411.0849; Statistics & Risk Modeling 32(3-4):159-176. doi:10.1515/strm-2015-0004. PDF:
  `Bauerle_Gilitschenski_Hanebeck_2015_Exact_Approximate_HMM_Filters_StochModels.pdf`.
- Hits the dangerous triple: a **finite-state CTMC** (not linear-Gaussian/LNA) + an **integrated**
  observation `Z_t = ∫ α(ε_s) ds + σ W_t` observed at discrete times + a closed-form treatment via
  the **endpoint-conditioned density** `g_ij(z,h)` of the integral given start `i` and end `j` — the
  same two-point object MacroIR's augmentation computes.
- Why it does **not** sink the claim: (1) it is a **single** CTMC (telegraph), **not** a many-channel
  ensemble — no channel count, no open-channel variance; (2) applied to **finance and fluid queues**,
  never biophysics/channels (Münch does not cite it); (3) **exact only for two states** (asymmetric
  telegraph / Bessel), with three approximations for k>2; (4) **no Kronecker-sum doubled generator, no
  `exp((Q⊕Q)t)=exp(Qt)⊗exp(Qt)` collapse, no (k+1) van Loan augmentation, no equivalence statement.**
  The retreat to approximations for k>2 is positive evidence the exact reduction was unavailable to
  expert practitioners. **The single most important cite-and-distinguish.**

---

## D. Mathematical machinery — concede as standard (`machinery`)

### CTMC reward / occupation-time moments (the moment engine)
- Hobolth & Jensen 2011, J. Appl. Probab. 48(4):911-924. PDF: `hobolth2011.pdf` (in repo).
  Endpoint-conditioned CTMC summary statistics via a van Loan auxiliary block — the (k+1) side,
  used in EM, conditioned on discrete endpoint *states* (not a time-averaged current).
- Tataru & Hobolth 2011, BMC Bioinformatics 12:465. PDF: `Tataru_Hobolth_2011_Endpoint_Sufficient_Statistics.pdf` (in repo).
- Pollett 2003, Math. Biosci. 182:213-225 (path integrals of a CTMC, start-conditioned). **Paywalled.**
- Bladt, Meini, Neuts & Sericola 2002, *Distributions of reward functions on CTMCs*, MAM4 proceedings;
  Sericola; Castella & Dujardin 2008. Reward distributions/moments via uniformization and block
  exponentials. **Proceedings/paywalled.** Marginal reward moments, never coupled to the end state as a
  filter.

### Kronecker-sum / doubled generator (the k² collapse identity)
- Zhou & Lange 2009, Adv. Appl. Probab. 41:no.1 (*Composition of Markov chains*). PDF: `Zhou_Lange_2009_Composition_Markov_Chains.pdf` (in repo). Kronecker composition of chains.
- Higham 2008, *Functions of Matrices*, SIAM: `exp(A⊕B)=exp(A)⊗exp(B)` as classical linear algebra. **Book.**
- Plateau & Stewart; Buchholz; Neuts (stochastic automata networks / Kronecker-structured CTMCs).
  **False friend:** the same `Q⊕Q` identity is used there to compose **independent** subsystems
  (product space, build-up), the opposite of MacroIR's correlated single-chain two-time (start,end)
  joint. Worth a footnote: even the experts who own the algebra were not pointing it at this object.

### Continuous-discrete local-linearization filtering (context)
- Carbonell / Jiménez, local linearization filters for continuous-discrete state-space models. PDF:
  `carbonell2008.pdf` (in repo). Discrete observation is an instantaneous point sample `z_{t_k}=Cx(t_k)+e`.
  Establishes the continuous-discrete filtering frame; not integrated-observation.

---

## E. Colored-noise augmentation — the only augmentation in the channel lineage (context)

### Chang 2014 (paywalled, Springer)
- J. Geodesy 88(12):1163-1170. doi:10.1007/s00190-014-0751-7. **Paywalled — not downloaded.**
- State augmentation (and measurement-time-differencing) to whiten **colored** measurement noise.
  This is the one augmentation reference Münch 2022 cites, and it is for noise coloring, **not** for a
  time-averaged/integrated observation — confirming the integrated-observation use is absent in the
  channel-Kalman lineage.

---

## Claim / concede / cite — one-line summary

| What | Status | Key references |
|------|--------|----------------|
| Integrated-measurement (time-averaged) treatment **for macroscopic ion channels** | **CLAIM (first in domain)** | gap vs Münch 2022, Moffatt 2007, Milescu 2005, Celentano-Hawkes 2004, Stepanyuk 2012/2014, Del Core-Mirams 2025 |
| Exact, efficient, all-k realization for the **many-channel finite-state CTMC** | **CLAIM** (vs Bäuerle's two-state-only) | Bäuerle 2014/2015 (distinguish) |
| Empirical demonstration that instantaneous methods are **misspecified** | **CLAIM** | figures 2/3/4 |
| Integrated-/averaged-measurement Kalman filter as a **device** | CONCEDE | Kalman 1960; Jazwinski 1970; Anderson-Moore 1979; Van Loan 1978; Särkkä 2013/2023; Zadrozny 1988; Harvey 1989; Yaghoobi-Särkkä 2024; Folia & Rattray 2018; Rubenzahl 2026 |
| van-Loan-on-augmented-integral pattern | CONCEDE | Rubenzahl 2026 |
| "Exact CTMC" as accuracy over LNA (for channels) | DO NOT CLAIM (they coincide; verified 1e-8) | `claude/scripts/verify_IR_vs_augmentation.py` |
| k² meta-state = (k+1) augmented-state equivalence | unification, openly conceded ingredients | Higham 2008; SAN literature (false friend) |
| CTMC reward-moment / Kronecker machinery | CONCEDE as standard | Hobolth-Jensen 2011; Bladt-Meini-Neuts-Sericola 2002; Zhou-Lange 2009 |

*Sources: two adversarial literature workflows (2026-06-26), task IDs wg0zu622u (k²⇄augmentation
equivalence) and wxa9rb7v8 (van Loan / augmented Kalman in macrocurrents), plus the numerical
equivalence check in `claude/scripts/verify_IR_vs_augmentation.py`. Confidence caveat: the surveys
are web + standard-literature, not line-by-line reads of every matrix-analytic monograph (Bladt &
Nielsen; Bladt-Meini-Neuts-Sericola were not fully retrieved). No killer prior art was found by either
adversarial critic.*

---
---

# PART II — prior art for the DIAGNOSTIC (the eLife paper's actual claim)

Added 2026-07-14 from a four-cluster adversarial search. Part I defends the algorithm; the eLife
paper claims a **diagnostic**. Almost none of the diagnostic is new. This part says exactly which
parts are conceded, to whom, and what is left.

**Notation bridge (memorize this; the reviewer has).** In the misspecification literature:
`H` = **sensitivity matrix** (the model's own information / expected negative Hessian);
`J` = **variability matrix** (covariance of the score); `G = H J⁻¹ H` = **Godambe information**;
`G⁻¹ = H⁻¹ J H⁻¹` = the **sandwich**. Our `C = H^{-1/2} J H^{-1/2}` has exactly the eigenvalues
`λᵢ` of `H⁻¹J`, and `tr C = tr(H⁻¹J)`. Symmetrizing adds nothing mathematically. **`C = I` is
Lindsay's "information unbiasedness"; the failure of `C = I` has been the definition of likelihood
misspecification since 1977.**

## II.0 Bottom line

1. **Every one of the three tests is a named classical test.** (i) residual whiteness = innovation
   whiteness / filter consistency (Mehra 1970; Bar-Shalom 2001; Ljung 1999). (ii) zero-mean score =
   **first Bartlett identity**. (iii) score covariance = reported Fisher = **second Bartlett identity**
   = White's **information matrix equality**. Chesher, Dhaene, Gouriéroux & Scaillet (1999) is the
   single cite for "test the Bartlett identities".
2. **`C` is not a new object.** Its spectrum is the weighted-χ² weight set of the misspecified LRT
   (Foutz & Srivastava 1977; Kent 1982). Its determinant and trace are **published test statistics**:
   the Log-Det and Log-Trace **Generalized Information Matrix Tests** (Golden, Henley, White & Kashner
   2016). Its trace is the **effective number of parameters** (TIC; Varin-Reid-Firth 2011).
3. **α⋆ = p / tr C is published verbatim, twice, as a Bayesian power on the likelihood.**
   Pauli, Racugno & Ventura (2011) eq. 2.3; Ribatet, Cooley & Davison (2012) §2.1 ("magnitude
   adjustment", `k = p / Σλᵢ`). Frequentist roots: Rotnitzky & Jewell (1990), Geys, Molenberghs &
   Ryan (1999). **Claiming this is novel is a first-paragraph kill.**
4. **½ log det C as an evidence correction is published**, as Lv & Liu (2014, JRSS-B) Theorem 1 and
   their **GBIC**, where `H_n = A_n⁻¹ B_n` is called the *covariance contrast matrix*. Their **GBIC_p**
   carries `tr(Ĥ_n) − log det(Ĥ_n)`, i.e. **both** invariants of `C` in one evidence criterion. It is
   also **contested** (Bhattacharya & Pati 2020) on the grounds that the sandwich does not appear in the
   posterior limit (Kleijn & van der Vaart 2012).
5. **The correlation/sample decomposition is the HAC long-run-variance decomposition**, and its scalar
   `κ_det = det(R)^{1/d}` is *literally* the Vats-Flegal-Jones multivariate-ESS ratio.
6. **In-domain: Münch et al. 2022 (eLife) already published diagnostic (i) plus an N_ch validity rule
   of thumb for a Kalman ion-channel likelihood** (verified from the paper: residual ACF is white;
   *"one should be careful with both algorithms for time traces with N_ch ∈ [10¹, 10²]"*;
   `error(N_ch) ∝ a/N_ch`). They do **not** compute the score, the Fisher information, the score
   covariance, or any sandwich. That is the delta.

## II.1 The concession table — what we call it vs what statistics calls it

| Our name | Standard name | Owner (cite these) |
|---|---|---|
| The three-part test | **Bartlett identities test** (IM test = 1st-order case) | Chesher, Dhaene, Gouriéroux & Scaillet 1999, LIDAM/IRES DP 1999019 |
| (i) standardized residuals white, unit variance | **innovation whiteness / filter consistency**; NIS/NEES | Mehra 1970 *IEEE TAC* 15(2):175-184; Bar-Shalom, Li & Kirubarajan 2001; Ljung 1999; Desroziers et al. 2005 *QJRMS* 131:3385. **In our own domain: Münch et al. 2022 eLife 11:e62714, Box 1-fig 2** |
| (ii) E[score] = 0 at truth | **first Bartlett identity**; QMLE pseudo-true value θ⋆ ≠ θ₀ | Huber 1967; White 1982 *Econometrica* 50(1):1-25 |
| (iii) Cov(score) = reported Fisher | **information matrix equality** / second Bartlett identity | White 1982; White 1994 ch. 11 |
| Information distortion matrix `C = H^{-1/2} J H^{-1/2}` | **the sandwich, symmetrically normalized.** Spectrum = eig(H⁻¹J) = the weights of the misspecified LRT (Kent; Vuong's `W`); `C = I` = **information unbiasedness** | Godambe 1960; Foutz & Srivastava 1977 *AoS* 5:1183; Kent 1982 *Biometrika* 69:19-27; Lindsay 1982; Vuong 1989 *Econometrica* 57:307. **No standard name for the symmetric form** — closest are the "information ratio" (Zhou, Song & Thompson 2012 *JASA* 107:205) and the "covariance contrast matrix" (Lv & Liu 2014) |
| eigenvalues `λᵢ(C)` | **generalized design effects** (survey statistics!) — `Δ̂ = Cov_model⁻¹ Ĉov_true`, "the eigenvalues of Δ̂ are termed generalized design effects"; Rao-Scott 1st-order correction divides χ² by the **mean** eigenvalue (= our α⋆), 2nd-order is a Satterthwaite match on Σλ, Σλ² | **Rao & Scott 1981** *JASA* 76:221-230; **Rao & Scott 1984** *AoS* 12:46-60 |
| `ĉ` variance-inflation factor (QAIC) | the **isotropic collapse** of the same spectrum: QAIC = TIC under the restriction `H⁻¹J = ĉ·I` | Burnham & Anderson 2002 |
| `log det C`, `tr C` as scalar summaries | **Log-Det GIMT and Log-Trace GIMT** (published test statistics; GIMT_Det = (1/q) log det(A⁻¹B)) | Golden, Henley, White & Kashner 2016 *Econometrics* 4(4):46; Golden et al. 2013 (Springer, White festschrift) |
| DCC = `H⁻¹ J H⁻¹` | **Godambe / sandwich covariance** | Godambe 1960 *AoMS* 31:1208; Huber 1967; White 1982 |
| DIB = `H⁻¹ E[score]` | first-order Newton step to the **pseudo-true value** θ⋆ (KL minimizer) | White 1982 |
| `tr C` as "effective sample" | **effective number of parameters**; TIC penalty; CLIC/CLAIC/CLBIC | Takeuchi 1976; Varin & Vidoni 2005 *Biometrika* 92:519; Varin, Reid & Firth 2011 *Stat. Sinica* 21:5-42 (p. 11: `dim(θ) = tr{H G⁻¹}`); Ng & Joe 2014 *Bernoulli* 20:1738 |
| Peak correction `α⋆ = p / tr C` | **magnitude adjustment / calibrated composite likelihood** — same constant | **Pauli, Racugno & Ventura 2011** *Stat. Sinica* 21:149-164 eq. (2.3); **Ribatet, Cooley & Davison 2012** *Stat. Sinica* 22:813-845 §2.1; Rotnitzky & Jewell 1990 *Biometrika* 77:485; Geys, Molenberghs & Ryan 1999 *JASA* 94:734 |
| Volume correction `½ log det C` | **GBIC**: the misspecification term in the Laplace expansion of the log marginal likelihood | **Lv & Liu 2014** *JRSS-B* 76(1):141-167, Thm 1 + GBIC/GBIC_p |
| Correlation distortion `R = J_sample^{-1/2} J_total J_sample^{-1/2}` | **HAC / long-run variance** `Σ = Γ₀ + Σ(Γ_k + Γ_kᵀ)` | Newey & West 1987; Andrews 1991 |
| `κ_det = det(R)^{1/d}` | **multivariate ESS ratio** `n/mESS`, `mESS = n(|Λ|/|Σ|)^{1/p}` | **Vats, Flegal & Jones 2019** *Biometrika* 106(2):321-337 |
| Sample distortion `C_sample = H^{-1/2} J_sample H^{-1/2}` | **information ratio** (variability vs sensitivity), OPG-vs-Hessian form of the IM test | Zhou, Song & Thompson 2012 *JASA* 107(497):205-213 |
| Power-likelihood framing | **generalized Bayes** | Bissiri, Holmes & Walker 2016 *JRSS-B* 78:1103; Grünwald & van Ommen 2017 *Bayesian Anal.* 12:1069 |
| "Simulation as ground truth for an approximate likelihood" | closest genres: SBC (algorithm-vs-model, **cannot** see this), ABC/SBI misspecification diagnostics (summary-level) | Talts et al. 2018; Frazier, Robert & Rousseau 2020 *JRSS-B* 82:421; Hermans et al. 2022 *TMLR* |

## II.2 What is actually left (defensible novelty, ranked)

1. **First-in-domain application.** Nobody has run the information-matrix-equality machinery on an
   ensemble-CTMC / macroscopic-current likelihood. Münch 2022 got as far as residual whiteness and
   stopped. This is the same contribution shape as Part I: bring the established tool, show it is needed.
2. **The population quantities are computed exactly, not estimated from one dataset.** The IM test's
   famous defect is its finite-sample behaviour (true vs nominal level off by 10× — Taylor 1987; Kennan
   & Neumann 1988; Orme 1990 *J. Econometrics* 46:309; fixed by bootstrap, Horowitz 1994 *J. Econometrics*
   61:395). Because the forward process is exactly simulable, we get `H` analytically and `J` by Monte
   Carlo over replicates. We are not *testing* misspecification against a null; we are **measuring** a
   known misspecification. Different problem, and it dodges the defect. **This is the strongest framing
   and should be stated in the Introduction.**
3. **The multiplicative decomposition `C = C_sample^{1/2} R C_sample^{1/2}` as a diagnostic pair.**
   Both factors are standard separately (IM test; HAC). Attributing an algorithm's distortion to
   *within-interval non-Gaussianity* vs *cross-interval score correlation* — and reading the second as
   the fingerprint of the discarded higher moments — appears to be unstated. Small, but real. Scope it
   tightly. Bonus: Milescu 2005 names *"the local time correlation of the current"* as his dominant
   error source, which is exactly the `R` factor.
4. **The full spectrum, not a scalar collapse.** Every consumer of `H⁻¹J` in the literature collapses
   it: TIC takes the trace, QAIC's `ĉ` assumes `C ∝ I`, Rao-Scott takes the mean eigenvalue,
   Satterthwaite takes `(Σλ)²/Σλ²`. Reporting the anisotropy (which *parameter directions* are
   inflated) is a legitimate, if modest, methodological point — and the data support it (the
   inflate/deflate directionality claim in the abstract *is* the anisotropy).
5. **The physical identification of `C_sample` with the third and fourth cumulants of the interval
   current.** White (1982, p. 7) proves that for a Gaussian model the IM equality holds iff skewness = 0
   and kurtosis = 3; our own cumulant expansion says the same. So `C_sample` is not an abstract
   discrepancy, it is a *measurement of the non-Gaussianity of the channel-ensemble current*, and it
   maps onto the multinomial and telegraphic corners for a reason. That reading is ours.
6. **The regime map.** A 3-axis (N_ch, Δ/τ, noise) quantitative map for five likelihood approximations.
   Genre exists (Schnoerr, Sanguinetti & Grima 2014 for moment closure; Grima et al. 2011 for
   CFPE/LNA; Dalmasso et al. 2024 for parameter-indexed coverage), and Münch 2022 published a 1-axis
   verbal version for this exact model class. An increment, not a discovery. **Do not lead with it.**

## II.3 Must-cite (omitting any of these is an unforced error)

Huber 1967 (Corollary to Thm 3, p. 231 — the sandwich); **White 1982** (Thm 3.3 = information matrix
equivalence; §4 = the IM test; p. 7 = the Gaussian skew/kurtosis case); White 1994 (ch. 11);
**Godambe 1960**; Foutz & Srivastava 1977; **Kent 1982**; Vuong 1989; **Rao & Scott 1981/1984**
(generalized design effects); Takeuchi 1976; Watanabe 2010 (ν = tr(H⁻¹J));
Chesher, Dhaene, Gouriéroux & Scaillet 1999; Orme 1990 + Horowitz 1994 (IM-test finite
sample); Freedman 2006; **Golden, Henley, White & Kashner 2016** (log-det / trace GIMT); Zhou, Song & Thompson 2012;
**Varin, Reid & Firth 2011**; Rotnitzky & Jewell 1990; Geys, Molenberghs & Ryan 1999;
**Pauli, Racugno & Ventura 2011**; **Ribatet, Cooley & Davison 2012**; Chandler & Bate 2007;
Pace, Salvan & Sartori 2011; **Lv & Liu 2014**; Bissiri, Holmes & Walker 2016; Grünwald & van Ommen 2017;
Lyddon, Holmes & Walker 2019; **Kleijn & van der Vaart 2012**; **Müller 2013** (*Econometrica* 81:1805);
Bhattacharya & Pati 2020; Newey & West 1987; **Vats, Flegal & Jones 2019**; Bollerslev & Wooldridge 1992;
Douc & Moulines 2012; Mehra 1970; Bar-Shalom, Li & Kirubarajan 2001; Ljung 1999; Talts et al. 2018;
Frazier, Robert & Rousseau 2020; **Münch et al. 2022**; Milescu, Akk & Sachs 2005; Moffatt 2007.

## II.4 The five objections a statistician reviewer will raise

**O1 (fatal if unanswered). `H` is the wrong matrix for the sandwich.** White's sandwich is
`A⁻¹ B A⁻¹` with `A = −E_true[∇²ℓ]`, the expected Hessian **under the true law**, at the pseudo-true
`θ⋆`. Our `H` is the model's **own** Gaussian Fisher `Σ_t [(∂μ)ᵀ(∂μ)/σ² + (∂σ²)ᵀ(∂σ²)/2σ⁴]`, i.e.
`−E_model[∇²ℓ]`. Under misspecification these differ, and the difference is *itself* part of what is
being measured. So `C = H^{-1/2} J H^{-1/2}` is a legitimate **IM-equality test statistic** (under the
null all three matrices coincide) but `H⁻¹ J H⁻¹` is **not** the asymptotic covariance of θ̂ unless
`A ≈ H`. **Answerable, and cheaply:** the codebase already computes both the numerical FD Hessian `F_b`
and the bridge `GFD = G_b^{-1/2} F_b G_b^{-1/2}`. Report GFD; if it is ≈ I, the substitution is
harmless and say so; where it is not, the DCC must be anchored on `F_b`. **Do not ship the paper
without this check.**

**O2. The anchor point.** White's theory holds at `θ⋆` (the KL minimizer), not at `θ_sim`. That
`E[score] ≠ 0` at `θ_sim` is precisely the statement `θ_sim ≠ θ⋆`. So the DIB is an estimate of
`θ⋆ − θ_sim` (the misspecification displacement), and a sandwich anchored at `θ_sim` is not the
asymptotic covariance of the estimator. Two different objects, both wanted, and the paper must say
which is which. (`battery_sim` vs `battery_pool` in the data is exactly this split.)

**O3. The evidence correction is an imposition, not an approximation.** The Laplace expansion of
`∫L(θ)π(θ)dθ` is already *correct* with the observed Hessian, misspecified or not; under
Kleijn & van der Vaart (2012) the misspecified posterior concentrates with covariance `H⁻¹`, **not**
the sandwich, so credible sets are not confidence sets. Bhattacharya & Pati (2020) make this objection
in print against Lv & Liu's GBIC. Adding `½ log det C` does not compute the same integral better; it
computes a **different, generalized-Bayes object**. That is defensible, but the warrant must be
decision-theoretic (Müller 2013: the sandwich posterior has lower asymptotic frequentist risk;
Bissiri et al. 2016), not "we improved the Laplace approximation". Also: McAlinn & Takanashi (2026,
arXiv:2602.01573) argue a tempered normalizer is not evidence and Bayes factors from it are not
well defined. State the warrant, or cut the claim.

**O4. Double counting / the scalar power is known to be the weak member.** Applying `α⋆` already
rescales the mode curvature to `α⋆H`, shifting the Laplace volume by `−(p/2) log α⋆`; the peak and
volume corrections coincide only when all `λᵢ(C)` are equal (Pace, Salvan & Sartori 2011 make exactly
this point: the moment-matched, Satterthwaite, and Chandler-Bate adjustments coincide only when
`d = 1`). And a **scalar** power provably cannot deliver the sandwich covariance: the magnitude-adjusted
posterior is `N(θ₀, (np)⁻¹ tr(H⁻¹J) H⁻¹)`, still proportional to `H⁻¹` (Ribatet et al. 2012 eq. 13;
Miller 2021 *JMLR* 22(168)). Only an **affine/curvature** adjustment does (Chandler & Bate 2007:
`CᵀHC = HJ⁻¹H`). If we keep the scalar, we must say why (cost, PSD, identifiability of the affine map).
Note also that `α⋆ ≠ ` Lyddon-Holmes-Walker's `w⋆ = tr(HJ⁻¹Hᵀ)/tr(H)`: they are different rules, equal
only when `J ∝ H`. Do not cite LHW as the source of `α⋆`; the sources are Pauli et al. / Ribatet et al.

**O5. The decomposition is HAC.** `J_total = J_sample + cross-lag terms` is the textbook long-run
variance. Frame it as: *the per-step score of a correctly specified filter is a martingale difference
sequence, so any non-zero cross-lag term is itself a misspecification signal, and we use the HAC
decomposition to attribute it.* That framing survives; "we introduce a decomposition" does not.

**O6 (the one nobody saw coming). For a Gaussian model, the IM test IS a skewness/kurtosis test.**
White 1982 works the normal case explicitly (p. 7): for `N(μ,σ²)`, `A(θ₀) = −B(θ₀)` **iff** skewness
`√β₁ = 0` and kurtosis `β₂ = 3`. Our emission density is Gaussian by construction, so a referee can
say *"your `C_sample = I` check is a Jarque-Bera normality test in a costume."* **This is true, and it
is good news, not bad** — it is exactly what the cumulant expansion in
`figures/paper/sample_correlation_distortion_analysis.md` already found (`J_sample − G_b` = skew +
kurtosis terms, exactly). Say it first: the sample distortion **is** the third and fourth cumulants of
the interval current, and that is why it tracks the multinomial/telegraphic corners. Owning this
converts an attack into the physical interpretation of the diagnostic.

**O7. Freedman's objection, which we can turn.** Freedman 2006 (*Am. Statist.* 60:299-302): *"It
remains unclear why applied workers should care about the variance of an estimator for the wrong
parameter."* His target is people who patch the variance with a sandwich and carry on. **Cite him
approvingly**: our point is that `C ≠ I` signals `θ⋆ ≠ θ_true`, i.e. bias, which is exactly the thing
he says everyone ignores — and the DIB measures it. A diagnostic beats a sandwich patch.

**Two numerical land mines, both already in print.** (1) White 1982 §4: the IM test's covariance
`V(θ⋆)` includes the estimation-error correction `∇D A⁻¹ ∇log f`; a naive covariance of the indicators
is wrong, and White himself warns the simplified estimator "is neither consistent nor necessarily
positive semi-definite when the null fails." (2) Ward 2023 (*Entropy* 25:512): TIC is not used in
practice because `tr(H⁻¹J)` is numerically hopeless above ~20 parameters (98-100% singular Hessians).
Report `p` and the conditioning; with `p = 6` we are safe, and saying so loudly is the defence.

**Bonus name collision.** "Information distortion" is already an established term in *neural coding*
(Dimitrov & Miller, rate-distortion clustering of stimulus-response classes) — a different meaning, in
an audience eLife shares. Not fatal; be aware.

*Sources: four adversarial search clusters (2026-07-14) covering composite-likelihood adjustment,
White/Huber/IM-test lineage, power-likelihood + evidence under misspecification, and filter
diagnostics / HAC / simulation-based-calibration; plus direct fetches of Münch et al. 2022, the GIMT
definitions, and the Chesher et al. Bartlett-identities paper.*
