# MacroIR prior-art map

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

### Stepanyuk et al. 2012 / 2014 — faster successors to Celentano-Hawkes
- 2012: PLoS ONE 7(4):e35208. doi:10.1371/journal.pone.0035208. PDF: `Stepanyuk_2012_Optimal_Estimation_Macroscopic_Currents_PLoSONE.pdf`.
- 2014: Front. Cell. Neurosci. 8:303. doi:10.3389/fncel.2014.00303. PDF: `Stepanyuk_2014_Synaptic_Receptors_Macroscopic_FrontCellNeurosci.pdf`.
- Macroscopic current as a non-stationary Gaussian process; exploit (semi)separability of the
  autocovariance for near-linear cost. Sampling instants are point spacings, not averaging
  windows. No within-interval integration, no augmentation, no van Loan.

### Del Core & Mirams 2025 — the newest competitor, still instantaneous
- Phil. Trans. R. Soc. A 383(2292):20240224. doi:10.1098/rsta.2024.0224. PDF: `DelCore_Mirams_2025_Stochastic_IonChannel_WholeCell_PhilTransA.pdf`.
- Parameter inference for stochastic ion-channel gating from whole-cell voltage-clamp data;
  LNA / moment state-space filter with a **point** Ohmic measurement `y_t = G_t x_t + ε`.
  Does not cite Moffatt 2007, Münch 2022, Folia & Rattray, or van Loan. Confirms the niche is
  still open in 2025.

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

## C. The CTMC integrated-observation near-miss — cite and distinguish (`distinguish`)

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
