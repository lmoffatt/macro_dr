# MacroIR prior-art map

**Rewritten 2026-07-14** after an exhaustive four-agent bibliographic review that checked the previous
version's claims against the actual PDFs. The previous version is preserved verbatim at
`archive/MacroIR_prior_art_map_v1_2026-06-26.md`. **It contains roughly twenty factual errors and one
false headline claim.** Do not cite from it.

**Non-destructive rule for this file.** Nothing is deleted. When a claim turns out to be wrong, it moves
to Part 0 (Errata) with the old wording quoted, the corrected fact, and the consequence for the paper.
The errata section is the most valuable part of this document: it is the list of things a referee could
have caught and we would have had no answer for.

**Three parts.**

- **Part 0. Errata.** What the old map asserted, and what is actually true.
- **Part I. Prior art for the ALGORITHM** (MacroIR, the integrated-measurement filter). This is the
  Comm Biol method. Positioning here matters for the method paper and for the eLife Introduction.
- **Part II. Prior art for the DIAGNOSTIC** (the information distortion matrix `C = H^{-1/2} J H^{-1/2}`,
  its decomposition, the MLE and evidence corrections). **This is what the eLife paper actually claims.
  If you are positioning the eLife paper, read Part II.**
- **Part III. What to claim, in what order**, after both parts.

Conventions. Each entry gives citation, local PDF (or paywall status), what the work actually does, its
observation model, and its **role**: `baseline-to-beat`, `concede` (a device we did not invent),
`distinguish` (near-miss to cite and separate from), `machinery`, or `gift` (a quotable that helps the
Introduction). `[VERIFY]` marks a claim that has **not** been checked against the source and must be
before it enters a manuscript.

---
---

# PART 0 — ERRATA

Fourteen corrections. E1 through E5 are citation hygiene (embarrassing, cheap to fix). E6 through E9
are mischaracterizations that would have made us look like we had not read the papers we were
distinguishing ourselves from. **E10 through E14 change what the paper can claim.**

## Citation hygiene

**E1. The Stepanyuk PDF was the wrong paper.**
Old map: *"Stepanyuk 2012, PLoS ONE 7(4):e35208"*, filed as a likelihood successor to Celentano-Hawkes.
That DOI is not Stepanyuk. It is **Wang, Xiao, Zeng, Yao, Yuchi & Ding (2012)**, and it is a
**least-squares** paper. The local PDF has been renamed (`git mv`) from `Stepanyuk_2012_...` to
`Wang_2012_...`. The real Stepanyuk lineage is **2011 + 2014**. Both now have their own entries, and
Wang 2012 turns out to be a *gift* (see A.6).

**E2. Zhou & Lange 2009: wrong title, and cited for an identity it never states.**
Old map: *"Adv. Appl. Probab. 41:no.1 (Composition of Markov chains)"*, cited for the Kronecker-sum
collapse. The paper is about composition Markov chains **of multinomial type**, and it **never states**
`exp(A ⊕ B) = exp(A) ⊗ exp(B)`. Correct home for that identity: **Horn & Johnson (1991), *Topics in
Matrix Analysis*, §4.4** (Kronecker sums and the exponential). Higham 2008 (*Functions of Matrices*) was
also cited for it; it is a fine book but not the source of this identity either. **Cite Horn & Johnson.**

**E3. "Movassagh, Vahidi, Moshiri et al." does not exist.**
Old map listed it as a second explicit "integrated measurement Kalman filter" usage, with a
*"verify exact title/year/venue"* note attached. It was searched for and **not found**. **Delete it. Do
not cite it.** The real second usage, and the actual root of the IMKF name, is **Fatehi & Huang (2017)**
(see E4).

**E4. Jazwinski 1970 was credited with the integrated/averaged measurement. Unsupported.**
Old map: *"the textbook source for state augmentation, including correlated/colored noise and
integrated/averaged measurements"*. The colored-noise half is right. The
**integrated/averaged-measurement half was not found in the text.** Cite Jazwinski for augmentation and
colored noise only. For the integrated-measurement filter as a *named* construction, the root is
**Fatehi & Huang (2017)** `[VERIFY venue/pages]`, with **Yaghoobi & Särkkä (2024)** as the current
statement (but see E9).

**E5. `carbonell2008.pdf` is not the paper the entry describes.**
Old map filed it under *"Carbonell / Jiménez, local linearization filters for continuous-discrete
state-space models"*. The actual PDF is **Carbonell, Jiménez & Pedroso**, a **generalization of Van Loan's
method for computing multiple integrals involving matrix exponentials**. This is *better* for us than
what the entry claimed: it belongs in the **machinery/concede** column (Part I, section D), directly
alongside Van Loan 1978, not in a "local linearization filtering context" section. Read it before
citing: it may already contain the exact block-exponential generalization MacroIR uses.

**E6. Rubenzahl et al. 2026: incomplete author list.** **Foreman-Mackey** is omitted from the old
entry's list. Fix before citing.

## Mischaracterizations

**E7. Del Core & Mirams 2025 is NOT a filter, and their argument against filtering is not only about cost.**
The old map called it an *"LNA / moment state-space filter"* and quoted them arguing against filtering on
**computational expense**. Two errors. (1) **Their method is not a filter.** Avoiding the filter is their
*selling point*: they integrate differential moment equations and fit, without a per-observation
correction step. (2) **They also argue against filtering on accuracy**, comparing their DME method
favourably against an EKF, and reporting roughly **0.5 h vs 36 h**. So they are a *sharper* adversary
than the old map made them, not a softer one. Our paper's answer must be a like-for-like accuracy
comparison, not "they only cared about speed". (The old map's other Del Core correction stands and is
worth repeating: **they do cite Moffatt 2007, Milescu 2005 and Münch 2022** (refs 30-32) and name
"MacroR" in the text. The pre-correction version of this map said they did not. Never repeat that.)

**E8. At least one "Münch 2022" quote is from the peer-review letter, not the paper.**
The old entry quotes Münch on finite integration time *"induc[ing] an intrinsic time scale ... not
represented in the algorithms"* and on *"higher analyzing frequency on minimally filtered data"*. At
least one of these lives in the **public peer-review response** (`reviews/`), not the paper body. Both are
citable (eLife publishes its reviews and responses, and a quote from the authors' own response is fair
game) but they must be attributed to the right document. **Action: re-verify each quote against the PDF
with a page number in hand before it enters the manuscript.** `[VERIFY]`

**E9. Yaghoobi & Särkkä 2024 explicitly *decline* state augmentation.**
The old map used them as evidence that "the integrated-measurement filter by augmentation is a standard,
named subfield". Half right. They do name the **Integrated Measurement Kalman Filter (IMKF)**, which is
the useful part. But their contribution is a **parallel-in-time** formulation, and they **avoid**
augmenting the state (that is what their construction is for). Cite them for the *name and the
subfield*, not as an instance of the augmentation device. The augmentation instances we actually have are
**Zadrozny 1988**, **Harvey 1989** (cumulator), **Folia & Rattray 2018**, **Calderazzo 2019**, and
**Rubenzahl 2026** (the only one with van Loan explicit).

**E10. Section D's blanket "never coupled to the end state" is false.**
The old map wrote of the CTMC reward-moment literature: *"Marginal reward moments, never coupled to the
end state as a filter."* **Hobolth & Jensen 2011 is endpoint-conditioned**, and the old map says so two
lines earlier. The blanket claim contradicts our own entry. What is actually true, and what we should say:
that literature conditions on **observed discrete endpoint states** (for EM in phylogenetics); we condition
on an **observed continuous integral** and marginalize the endpoint state. The distinction is the
observation model, not the conditioning.

## Claim-changing

**E11 (the big one). Bottom line #1 was FALSE. The anti-alias filter has been inside the likelihood
since 1992 — in the single-channel branch.**
Old map, bottom line #1: *"The channel domain cell is empty. Every macroscopic / ensemble channel
inference method models the observation as an instantaneous point sample ... Finite integration time is
treated as a misspecification to be robust to, not modeled."* The first sentence ("the channel domain
cell is empty") is **false as written**, and a single-channel referee will know it. The
**single-channel HMM lineage has modeled the recording filter inside the likelihood for thirty years**:

- **Fredkin & Rice (1992)** — filtered single-channel likelihood.
- **Michalek, Lerche, Wagner & Timmer (1999)** — literally *"moving-average filtered hidden Markov
  models"*.
- **Venkataramanan & Sigworth (1998-2002)** — the metastate/vector-HMM formulation for filtered data.
- **Qin, Auerbach & Sachs (2000)** — *"The filtering is modeled using a finite impulse response (FIR)
  filter"*, **and**: *"Extension of the algorithm to data containing multiple channels is described."*

The last one is the dangerous sentence. Someone can hold it up and say the multi-channel filtered
likelihood was published in 2000.

**What survives, and it survives cleanly, is the scaling argument.** The FIR/metastate construction
builds a hidden state over the filter's memory: cost `k^L` for `L` filter taps, and for `N` channels the
aggregated state space grows combinatorially on top of that. It is usable for one channel and a handful
of taps. It is not usable for `N = 10³` or `10⁴` channels. **MacroIR's cost is `O(k³)` per interval,
independent of `N`.** That is the claim: not "nobody modeled the filter", but **"the exact treatment
existed only at a cost that scales out of reach of the macroscopic regime, and we give one that does
not."** Rewrite claim (i) accordingly (see Part III).

Adjacent, and worth having: **Schröder & Hansen (2005)**, **Almanjahie et al. (2015/2019)**, and the
beta-distribution family for filtered single-channel amplitudes. Same message.

**E12. Kilic, Sgouralis & Pressé et al. (2021) is the strongest near-miss anywhere in the literature, and
it was missing.** See Part I, section C.2. It writes the observation as an **exact within-window integral
of the state-dependent rate, for an arbitrary number of states**. Claim (ii) as previously worded ("the
only prior CTMC integrated-observation filter, Bäuerle, was exact only for two states") **is no longer
true and must be rewritten.**

**E13. The Bäuerle "category error".** The old map's most attackable sentence: *"The retreat to
approximations for k>2 is positive evidence the exact reduction was unavailable to expert
practitioners."* **This is wrong, and it is wrong in a way that reveals we did not understand their
problem.** Bäuerle et al. need the **full density of the integral, and that density has an atom** (the
integral of a single telegraph process concentrates mass at the extremes). Moments cannot produce an
atom. That is why they retreat to approximations for `k > 2`: not because the "exact reduction" was
beyond them, but because **the object they need is not a moment problem.** MacroIR gets away with two
moments **only because the ensemble is Gaussian by the central limit theorem** — the very assumption
whose breakdown is this paper's subject. Two further corrections in the same entry: their PDE route for
`k > 2` is **exact up to numerics** (so "exact only for two states" is also wrong), and their noise
structure (`σ W_t`, a Brownian increment with variance ∝ `h`) **agrees** with MacroIR's `1/Δ` scaling
rather than contradicting it. **Delete the "positive evidence" sentence. Replace it with the atom
argument, which is a genuinely good distinction and makes us look like we read them.**

**E14. Kronecker structure is already in the channel field.**
**Albertsen & Hansen (1994)** use the product/Kronecker construction for multi-channel recordings. The
`Q ⊕ Q` doubled-generator move should be conceded as an available idea in this domain, not presented as
an import from matrix-analytic queueing theory. (The "false friend" observation from the old map still
stands and is still worth a footnote: the stochastic-automata-networks literature uses the same `⊕` to
compose **independent** subsystems, the opposite of our correlated two-time joint of a **single** chain.)

---
---

# PART I — prior art for the ALGORITHM

## I.0 Bottom line (REVISED 2026-07-14)

1. **The *macroscopic* cell is empty; the *channel* cell is not.** Every macroscopic/ensemble channel
   inference method models the observation as an instantaneous point sample of the open-channel count,
   `y_t = H n_t + noise`, never as the within-interval time-average `ȳ = (1/Δ) ∫ I dt`. That part
   survives every check. But **the single-channel HMM branch has modeled the recording filter inside the
   likelihood since 1992** (E11), and Qin et al. 2000 explicitly extend to multiple channels. The claim
   must be stated as **scaling**, not as absence: the exact treatment existed at `k^L` cost and does not
   reach the macroscopic regime; MacroIR is `O(k³)` per interval, independent of `N`.

2. **The device is not new.** The integrated/averaged-measurement Kalman filter (augment the state with
   the running integral, propagate by a van Loan block matrix exponential, condition the endpoint on the
   observed integral) is old in control and econometrics (Zadrozny 1988; Harvey 1989), already used in
   *ensemble biophysics* (Folia & Rattray 2018), and, with van Loan explicit, in astronomy (Rubenzahl et
   al. 2026). It is a named subfield (IMKF; Fatehi & Huang 2017; Yaghoobi & Särkkä 2024). **Concede it
   loudly.**

3. **The "exact CTMC beats their LNA" edge is weak and must not be sold.** For the *linear* channel
   dynamics the linear-noise approximation is exact in the first two moments. We verified numerically
   that the augmented filter built on the chemical-Langevin diffusion equals MacroIR's exact-CTMC result
   to 1e-8 (`claude/scripts/verify_IR_vs_augmentation.py`). A Folia & Rattray style LNA augmentation,
   specialized to channels, **would coincide with MacroIR.** The exact-CTMC realization is a
   derivation/efficiency/channel-native point, not a numerical-accuracy claim.

4. **The strongest near-miss is Kilic et al. 2021, not Bäuerle.** Kilic writes the exact within-window
   integral for arbitrary `K`. It is single-molecule, the observation is Poisson photon counts, and they
   sample trajectories by MCMC instead of marginalizing analytically. That is the distinction, and it is
   a real one, but it is much narrower than the old map believed.

5. **The k² meta-state ⇄ (k+1) augmented-state equivalence is apparently unstated as such**, though both
   halves are standard and a matrix-analytic reviewer may call it folklore. Frame as unification,
   concede the ingredients (and see E14: Kronecker is already in this field).

**Defensible novelty, revised:**
(i) The first **macroscopic-regime-feasible** exact treatment of the time-averaged observation: prior
exact treatments exist for single channels at `k^L` cost; MacroIR is `O(k³)` per interval, independent of
channel number.
(ii) The exact, efficient, channel-native realization for the **many-channel ensemble** CTMC, where the
existing exact integrated-observation constructions are for **single** chains (Bäuerle; Blackwell) or for
**single molecules under Poisson emission with MCMC trajectory sampling** (Kilic).
(iii) The empirical demonstration that the field's instantaneous methods are misspecified (miscalibrated
uncertainty, biased model selection). **This is the eLife paper, and its prior art is Part II.**

Same contribution shape as Münch et al. 2022 (who brought the textbook Kalman filter to channels): bring
the established tool, show it is needed.

---

## A. In-domain macroscopic baselines: the gap (all instantaneous; `baseline-to-beat`)

### A.1 Münch, Paul, Schmauder & Benndorf 2022 — the direct competitor
- eLife 11:e62714. doi:10.7554/eLife.62714. PDF: `Munch_2022_Bayesian_Kalman_IonChannels_eLife.pdf`
  (peer-review letter in `reviews/`).
- Generalized Bayesian / Kalman filter for macroscopic current and confocal patch-clamp fluorometry.
  Observation is **instantaneous** (Eq. 4): `y_t = H n_t + ν_t ~ N(·|H n_t, Σ_t)`, state-dependent
  variance `σ_m² + n_open σ_op²`, state propagated by a plain `exp(KΔt)`.
- State augmentation appears **once**, only as an option for **colored** noise (citing Chang 2014); not
  adopted. Finite analog filtering / integration time is treated as a misspecification to be robust to.
  They *cite* the ensemble integrated-Kalman/LNA literature (Komorowski; Finkenstädt; Folia & Rattray
  2018) and still decline the integrated-observation construction.
- **Quote hygiene: see E8.** Verify page numbers, and separate paper from peer-review response.
- **They are also Part II prior art.** Münch already published diagnostic (i), residual whiteness, plus
  an `N_ch` validity rule of thumb. See II.0.6. That is the closest thing to this paper that exists.

### A.2 Moffatt 2007 — the baseline MacroIR extends, and the open problem it closes
- Biophys. J. 93(1):74-91. doi:10.1529/biophysj.106.101212. PDF:
  `Moffatt_2007_Macroscopic_Current_Fluctuations_BiophysJ.pdf`.
- Recursive Gaussian filter for macroscopic current capturing local time-correlation of intrinsic channel
  noise; instantaneous Gaussian readout of the current given occupancies; `exp(QΔt)` by
  eigendecomposition. No augmentation, no van Loan.
- **Names the gap explicitly:** *"the measurements are considered to be instantaneous, where
  electrophysiological measurements always result from some time-averaging ... a formulation of the
  algorithms that can deal with time-averaged signals is therefore fundamental."* MacroIR closes the
  author's own stated open problem. Strongest, least attackable framing.

### A.3 Milescu, Akk & Sachs 2005
- Biophys. J. 88(4):2494-2515. doi:10.1529/biophysj.104.053256. PDF:
  `Milescu_2005_MLE_Macroscopic_Currents_BiophysJ.pdf`.
- Macroscopic-current maximum likelihood. Gaussian in the **point** current `I_t`, mean and variance at
  the instantaneous occupancy; `A_t = exp(QΔt)` by spectral decomposition. Its van Loan citation is the
  Moler & Van Loan matrix-exponential *computation* survey, not the integrals method.
- **Gift for Part II:** Milescu names *"the local time correlation of the current"* as his dominant error
  source. That is precisely the `R` (correlation-distortion) factor. Use it.

### A.4 Celentano & Hawkes 2004
- Biophys. J. 87(1):276-294. doi:10.1529/biophysj.103.036632. PDF:
  `CelentanoHawkes_2004_Covariance_Matrix_Kinetic_Fitting_BiophysJ.pdf`.
- Direct fit of kinetic parameters to macroscopic currents using the full inter-time-point covariance
  matrix of intrinsic noise (O(n³)). The covariance is **between time points** (autocorrelation of the
  process), not a model of the measurement's within-interval time-average. Observation still
  instantaneous. No augmentation, no van Loan.

### A.5 Stepanyuk et al. 2011 / 2014 — faster successors to Celentano-Hawkes
- **See E1: the old "Stepanyuk 2012 / PLoS ONE e35208" entry was a different paper.** The real lineage:
- 2011: PLoS ONE 6(12):e29731. doi:10.1371/journal.pone.0029731. **PDF not in repo.**
- 2014: Front. Cell. Neurosci. 8:303. doi:10.3389/fncel.2014.00303. PDF:
  `Stepanyuk_2014_Synaptic_Receptors_Macroscopic_FrontCellNeurosci.pdf`.
- Macroscopic current as a non-stationary Gaussian process (ML NSFA); exploits **semiseparability** of the
  autocovariance for near-linear cost. A semiseparable-covariance Gaussian likelihood is
  filter-equivalent mathematically, but they never frame it as filtering and never mention Kalman; the
  covariance is the model's **prior** covariance, never conditioned on the observed trace. Sampling
  instants are point spacings, not averaging windows.
- **Gifts, on how little this whole family is used:** *"PS NSFA [is] the most frequently used noise
  analysis approach"*; *"the unitary current is virtually the only parameter that can be reliably obtained
  from this type of analysis"*; **"To the best of our knowledge, kinetic rates have never been estimated
  for any synaptic receptors in their intrinsic environment."**

### A.6 Wang, Xiao, Zeng, Yao, Yuchi & Ding 2012 — a macroscopic-fitting paper that declines the likelihood (`gift`)
- PLoS ONE 7(4):e35208. doi:10.1371/journal.pone.0035208. PDF:
  `Wang_2012_Optimal_Estimation_Macroscopic_Currents_PLoSONE.pdf`.
- **Least squares on the deterministic mean current** (particle-swarm optimizer). No likelihood, no
  variance term, no gating fluctuation. Cites Milescu 2005 and **explicitly opts out of maximum
  likelihood on computational grounds.**
- Role: direct evidence for the Introduction that a dedicated macroscopic-fitting paper, seven years
  after the likelihood method existed, chose least squares on the mean.

### A.7 Del Core & Mirams 2025 — the newest competitor, and a sharper one than we thought
- Phil. Trans. R. Soc. A 383(2292):20240224. doi:10.1098/rsta.2024.0224. PDF:
  `DelCore_Mirams_2025_Stochastic_IonChannel_WholeCell_PhilTransA.pdf`.
- Parameter inference for stochastic ion-channel gating from whole-cell voltage-clamp data. Differential
  moment equations with a **point** Ohmic measurement `y_t = G_t x_t + ε`.
- **See E7. They are not a filter, and not filtering is the point.** They compare against an EKF and beat
  it on accuracy as well as cost (roughly 0.5 h vs 36 h). Our reply has to be a like-for-like accuracy
  comparison.
- **They do cite Moffatt 2007, Milescu 2005 and Münch 2022** (refs 30-32) and name "MacroR" in the text.
  They do not cite Folia & Rattray or van Loan.
- Three quotable gifts:
  - **Third-party endorsement of the recursive point**, on Moffatt 2007: *"The authors proved that their
    method provides better estimates than approaches that do not consider statistical dependence of
    successive measurements [31]"* (i.e. beats Milescu 2005 precisely on time correlation).
  - **Their argument against filtering:** *"approaches that are based on filtering techniques are
    computationally expensive because an integration step of the differential moment equations (DMEs),
    and the corresponding updates in the correction step, are computed between every consecutive time
    points where the measurements are collected."* Our paper answers what you lose by skipping the filter.
  - **They declare the channel-number problem open in 2025:** *"how to best infer the relationship between
    the noise observed in the ionic measurements and the total number of channels in the cell membrane
    contributing to gating is still unclear."* Directly relevant to our Fisher-information result on
    `N_ch` identifiability.

### A.8 The macroscopic field mostly does not fit a likelihood at all (NEW, and it is the Introduction's premise)
This is the evidence that killed the earlier draft sentence *"several likelihood approximations are in
routine use"*. It is false. The field fits the deterministic mean by least squares.

- **Clerx, Beattie, Gavaghan & Mirams (2019), "Four Ways to Fit an Ion Channel Model"**, *Biophys. J.*
  117(12):2420-2437. **All four ways are deterministic.** No stochastic likelihood among them.
- **Owen & Mirams (2025), IonBench** `[VERIFY exact venue/year]`: the benchmark states plainly, *"We do
  not include optimisation approaches for stochastic models."*
- Non-stationary fluctuation analysis (NSFA) is widely used but does a **different job**: it extracts
  unitary conductance and channel number from the variance-mean relation, not kinetic rates (see
  Stepanyuk's own quote in A.5).
- **Consequence.** Likelihood-based fitting of macroscopic currents is a niche of four or five groups
  (Milescu; Celentano-Hawkes; Stepanyuk; Moffatt; Münch). The Introduction must say **that**, and the
  paper's motivation is not "fix the approximations everyone uses" but **"the approximations that the few
  who do this rely on are distorted, and here is by how much, and in which regime."** Honest, verifiable,
  and still worth publishing.

### A.9 Requadt et al. 2025 `[VERIFY — get the PDF, read it, extract the quote]`
Flagged by the 2026-07-14 search as an Introduction-grade quotable on the current state of
macroscopic-current inference. **Not yet read. Do not cite until it is.** Action item: retrieve, read,
and either promote to a full entry or delete this stub.

---

## B. The augmented / integrated-measurement Kalman filter — full bibliography (`concede`)

Cite this column openly so the honesty is bulletproof. The filter MacroIR re-derives is the
**integrated-measurement Kalman filter via state augmentation**: append the running integral `s = ∫ I dt`
to the state, propagate `[x; s]` with a van Loan block matrix exponential, condition the endpoint on the
observed integral via `Cov(x_end, s)`. **None of this is ours.**

### B1. Foundations

- **Kalman, R. E. (1960).** A new approach to linear filtering and prediction problems. *Trans. ASME J.
  Basic Eng.* 82(1):35-45. doi:10.1115/1.3662552. PDF:
  `Kalman_1960_New_Approach_Linear_Filtering_JBasicEng.pdf`. The root; Münch 2022 cites it as the basis.
- **Jazwinski, A. H. (1970).** *Stochastic Processes and Filtering Theory.* Academic Press (Dover 2007).
  **Book, not downloaded.** The textbook source for **state augmentation**, including correlated/colored
  noise. **See E4: do not credit it with the integrated/averaged measurement.**
- **Anderson, B. D. O. & Moore, J. B. (1979).** *Optimal Filtering.* Prentice-Hall (Dover 2005). PDF:
  `Anderson_Moore_1979_Optimal_Filtering.pdf`. The minimum-variance filter reference Münch 2022 cites for
  KF optimality. **It is deliberately discrete-time.** Do not cite it for the continuous-discrete filter.
- **Särkkä, S. (2013).** *Bayesian Filtering and Smoothing.* CUP. PDF:
  `Sarkka_2013_Bayesian_Filtering_and_Smoothing_CUP.pdf`. Cite for the general Bayesian filtering frame.
  The **matrix-fraction decomposition** attribution in the old map was not found in this text as
  described; that material is in **Särkkä & Solin (2019), *Applied Stochastic Differential Equations***.
  Cite the right one. `[VERIFY which of the two you need]`
- **Van Loan, C. F. (1978).** Computing integrals involving the matrix exponential. *IEEE Trans. Automatic
  Control* 23(3):395-404. doi:10.1109/TAC.1978.1101743. **Paywalled.** The block-triangular device for the
  integral of `exp(At)` and the process-noise covariance integral. The numerical primitive the
  augmentation rests on. **Cite by name.**
- **Horn, R. A. & Johnson, C. R. (1991).** *Topics in Matrix Analysis*, §4.4. The actual source of
  `exp(A ⊕ B) = exp(A) ⊗ exp(B)`. **See E2.**

### B2. Integrated / averaged / flow observation handled by augmentation (the exact pattern)

- **Zadrozny, P. (1988).** Gaussian likelihood of continuous-time ARMAX models when data are stocks and
  flows at different frequencies. *Econometric Theory* 4(1):108-124. **Paywalled.** The cleanest
  econometric statement: **flow** variables (integrals over the period) handled by augmenting the state
  with the integral coordinate and running a Kalman filter on the exact-discretized SDE. Linear-Gaussian,
  decades old.
- **Harvey, A. C. (1989).** *Forecasting, Structural Time Series Models and the Kalman Filter.* CUP, §6.3
  (**"cumulator variable"**). **Book.** With Phillips (1959/1978), Bergstrom (1984), Chambers: the
  flow-vs-stock temporal-aggregation lineage. Cite to credit the method's age.
- **Fatehi, A. & Huang, B. (2017).** The root of the **"integrated measurement Kalman filter"** name and
  the slow-rate integrated (sample-and-accumulate) sensor construction. `[VERIFY venue, volume, pages]`
  **This is the citation that E3's phantom was standing in for.**
- **Yaghoobi, F. & Särkkä, S. (2024).** Parallel state estimation for systems with integrated
  measurements. arXiv:2410.00627. PDF: `Yaghoobi_Sarkka_2024_Parallel_Integrated_Measurement_IMKF.pdf`.
  Names the **Integrated Measurement Kalman Filter (IMKF)** for Slow-Rate inTegrated Measurement models.
  **See E9: they explicitly avoid augmentation.** Cite for the name and the subfield's currency, not as an
  augmentation instance.
- **Rubenzahl, R. A., Hattori, S., Särkkä, S., Farr, B. F., Luhn, J., Foreman-Mackey, D. & Bedell, M.
  (2026).** Scalable Gaussian processes for integrated and overlapping measurements via augmented
  state-space models. *Astron. J.* doi:10.3847/1538-3881/ae4d0b (arXiv:2601.02527). PDF:
  `Rubenzahl_2026_Integrated_Overlapping_Measurements_Augmented_SSM_AJ.pdf`. **Closest by exact
  construction:** augments the SSM with an integral state that resets at exposure start and is observed
  as a time-average at exposure end, **using van Loan (1978) explicitly**. Linear-Gaussian GP/SDE state,
  astronomy. **Cite so we do not over-claim the van-Loan-on-augmented-integral pattern.** (Author list
  corrected, E6.)
- **Folia, M. M. & Rattray, M. (2018).** Trajectory inference and parameter estimation in stochastic models
  with temporally aggregated data. *Stat. Comput.* 28:1053-1072. PDF:
  `Folia_Rattray_2018_Trajectory_Inference_Aggregated_StatComput.pdf`. **The closest precedent in ensemble
  biophysics:** integrated observation `y = P ∫X du + noise` handled by augmenting the **LNA** state with
  the running integral inside a continuous-discrete Kalman filter (luciferase integrated over minutes).
  Moment-ODE LNA, no van Loan, gene expression. **Per bottom line #3, for linear channel dynamics this
  coincides numerically with MacroIR.** Cite as the methodological precedent.
- **Calderazzo, S., Brancaccio, M. & Finkenstädt, B. (2019).** Filtering and inference for stochastic
  oscillators with distributed delays. *Bioinformatics* 35(8):1380-1388. PDF:
  `Calderazzo_2019_Filtering_Inference_Oscillatory_LNA_Bioinformatics.pdf`. Integrated camera-exposure
  observation on the same LNA framework, built on Folia & Rattray.

### B3. Colored-noise augmentation — the only augmentation the channel field uses (context)

- **Chang, G. (2014).** On Kalman filter for linear system with colored measurement noise. *J. Geodesy*
  88(12):1163-1170. **Paywalled.** State augmentation / measurement-time-differencing to whiten
  **colored** noise. The one augmentation Münch 2022 cites, and it is for noise coloring, **not** the
  time-average. Cite to show the channel field invokes augmentation only for colored noise.

---

## C. Integrated-observation near-misses — cite and distinguish (`distinguish`)

Ordered by how much they hurt. **C.2 (Kilic) is the one to worry about.**

### C.1 Blackwell 2018/2019 — integrated continuous-time HMMs
- **Blackwell, P. G.** *Integrated continuous-time hidden Markov models.* arXiv:1807.11907 (stat.ME).
  **Not downloaded.**
- Makes **the same conceptual move as MacroIR**, and says so in its abstract: *"a new class of integrated
  continuous-time hidden Markov models in which each observation depends on the underlying state of the
  process over the whole interval since the previous observation, not only on its current state ... under
  appropriate conditioning, a model in this class can be regarded as a conventional hidden Markov model,
  enabling use of the Forward Algorithm for efficient evaluation of its likelihood."* Condition on the
  interval, recover a filter recursion.
- Why it does not sink the claim: **movement ecology**, a **single** chain (a switching diffusion for an
  animal's behavioural state). No channel count, no open-channel variance, no Kronecker-sum collapse, no
  van Loan augmentation.
- **Cite anyway.** A statistician referee may know it, and it is the cleanest general statement of the
  device for CTMC-with-integrated-observation outside our domain.

### C.2 Kilic, Sgouralis, Pressé et al. 2021 — the strongest near-miss (ADDED 2026-07-14; see E12)
- **Cell Rep. Phys. Sci.** 2:100409, with a companion in **Biophys. J.** 120:409. `[VERIFY exact author
  list and both DOIs; PDFs not in repo]`
- **The observation is an exact within-window integral of the state-dependent rate**, for an **arbitrary
  number of states**:
  `w_n ~ Poisson( μ_back · τ + ∫ μ_{T(t)} dt )`,
  where the integral runs over the exposure window and `T(t)` is the hidden CTMC trajectory. That is
  structurally our observation model, and it is not restricted to two states.
- **This is why the old claim (ii) is dead.** "The only prior CTMC integrated-observation filter was exact
  only for two states" is false: Kilic is exact for any `K`.
- **What still distinguishes MacroIR, and it is real:**
  1. **Single molecule, not an ensemble.** No channel count `N`, no `N`-dependent gating variance, no CLT
     regime, and therefore none of this paper's actual subject.
  2. **Poisson photon emission, not Gaussian current with a noise floor.** Different emission family; the
     conjugacy and the failure modes are different.
  3. **They sample the trajectory by MCMC; they do not marginalize it analytically.** MacroIR's whole point
     is the closed-form marginalization: `O(k³)` per interval, deterministic, differentiable, embeddable in
     a likelihood a Newton optimizer can use. Kilic's cost is an MCMC over paths.
  4. **No filter recursion, no van Loan, no Kronecker collapse.** The construction is Bayesian nonparametric
     sampling, not a Gaussian filter.
- **Role:** cite prominently and distinguish on (1) and (3). Getting caught not knowing this paper would be
  worse than any of its content. **Retrieve both PDFs.**

### C.3 Bäuerle, Gilitschenski & Hanebeck 2014 / 2015 (REVISED; see E13)
- arXiv:1411.0849; *Statistics & Risk Modeling* 32(3-4):159-176. doi:10.1515/strm-2015-0004. PDF:
  `Bauerle_Gilitschenski_Hanebeck_2015_Exact_Approximate_HMM_Filters_StochModels.pdf`.
- Hits the dangerous triple: a **finite-state CTMC** (not linear-Gaussian/LNA) + an **integrated**
  observation `Z_t = ∫ α(ε_s) ds + σ W_t` at discrete times + a closed-form treatment via the
  **endpoint-conditioned density** `g_ij(z,h)` of the integral given start `i` and end `j`, the same
  two-point object MacroIR's augmentation computes.
- **The correct distinction (the atom argument).** Bäuerle et al. need the **full conditional density** of
  the integral, and for a single telegraph process **that density has an atom** (mass concentrated at the
  extremes, corresponding to no-jump paths). **Moments cannot represent an atom.** That is why they work
  with special functions for `k = 2` and go to a PDE for `k > 2`. **MacroIR needs only the first two
  moments, and it is entitled to them only because the ensemble is Gaussian by the central limit theorem.**
  In other words: their problem is harder than ours, and our shortcut is legitimate exactly in the regime
  where `N` is large. **The regime where our shortcut fails is the subject of this paper.** This is a much
  better paragraph than the old map's.
- Three further corrections to the old entry:
  - **"Exact only for two states" is wrong.** Their PDE route for `k > 2` is exact up to numerics. What is
    special about `k = 2` is the closed form (asymmetric telegraph / Bessel), not the exactness.
  - **Delete "the retreat to approximations for k>2 is positive evidence the exact reduction was
    unavailable to expert practitioners."** It is a category error (E13) and it is the single most
    attackable sentence in the old map.
  - **Their noise structure agrees with ours.** `σ W_t` is a Brownian increment over the window: variance ∝
    `h`, so the *average* has variance ∝ `1/h`. That is MacroIR's `1/Δ` scaling. Do not present it as a
    contrast.
- Still true, and still the reason it does not sink us: it is a **single** telegraph chain, **not** a
  many-channel ensemble; applied to **finance and fluid queues**, never biophysics (Münch does not cite
  it); **no Kronecker-sum doubled generator, no `exp((Q⊕Q)t) = exp(Qt) ⊗ exp(Qt)` collapse, no `(k+1)` van
  Loan augmentation, no equivalence statement.**

### C.4 The single-channel filtered-likelihood lineage (NEW; this is E11, and it is the one that bites)

**The old map's headline claim ("the channel domain cell is empty") dies here.** The single-channel HMM
branch has put the recording filter **inside the likelihood** for thirty years.

- **Fredkin, D. R. & Rice, J. A. (1992).** Filtered single-channel likelihood. `[VERIFY exact title/venue]`
- **Michalek, S., Lerche, H., Wagner, M. & Timmer, J. (1999).** *"Moving-average filtered hidden Markov
  models"* for single-channel currents. `[VERIFY venue]` The name says it.
- **Venkataramanan, L. & Sigworth, F. J. (1998-2002).** The **metastate / vector-HMM** formulation:
  expand the hidden state over the filter's memory so that the filtered observation is again Markov.
- **Qin, F., Auerbach, A. & Sachs, F. (2000).** *"The filtering is modeled using a finite impulse response
  (FIR) filter"* **and** *"Extension of the algorithm to data containing multiple channels is described."*
  **The most dangerous sentence in the prior art.** Read this paper. `[VERIFY the multi-channel section:
  how many channels, at what cost]`
- Adjacent: **Schröder & Hansen (2005)**; **Almanjahie et al. (2015/2019)**; the beta-distribution family
  for filtered single-channel amplitude distributions.

**How to survive it: make the claim about scaling, not absence.**
The metastate/FIR construction is exact and it is *published*. Its cost is `k^L` for `L` filter taps, on
top of an aggregated state space that grows combinatorially in the number of channels. It works for one
channel and a short filter. It does not reach `N = 10³` or `10⁴`. **MacroIR costs `O(k³)` per interval,
independent of `N`.** State it that way and the claim is unassailable; state it as "nobody did this" and
a single-channel referee kills the paper in one sentence.

**Bonus:** this is not a weakness in the story, it is a *better* story. The exact treatment of the
recording filter is a known necessity in the single-channel world (nobody there would dream of ignoring
it), and the macroscopic world simply dropped it because the exact route did not scale. That is a clean,
true, and interesting thing to say in an Introduction.

---

## D. Mathematical machinery — concede as standard (`machinery`)

### CTMC reward / occupation-time moments (the moment engine)
- **Hobolth & Jensen 2011**, *J. Appl. Probab.* 48(4):911-924. PDF: `hobolth2011.pdf`.
  **Endpoint-conditioned** CTMC summary statistics via a van Loan auxiliary block: the `(k+1)` side, used
  in EM, conditioned on discrete endpoint **states**. **See E10:** the old map's blanket "never coupled to
  the end state" is false; this *is* endpoint-conditioned. The real distinction is the **observation
  model**: they condition on an observed endpoint state, we condition on an observed continuous integral
  and marginalize the endpoint.
- **Tataru & Hobolth 2011**, *BMC Bioinformatics* 12:465. PDF: `Tataru_Hobolth_2011_Endpoint_Sufficient_Statistics.pdf`.
- **Pollett 2003**, *Math. Biosci.* 182:213-225 (path integrals of a CTMC, start-conditioned). **Paywalled.**
- **Bladt, Meini, Neuts & Sericola 2002**, *Distributions of reward functions on CTMCs*, MAM4; Sericola;
  Castella & Dujardin 2008. Reward distributions/moments via uniformization and block exponentials.
  **Proceedings/paywalled; not fully retrieved.**

### Kronecker-sum / doubled generator (the k² collapse identity)
- **Horn & Johnson 1991**, *Topics in Matrix Analysis*, §4.4. **The identity's actual source** (E2).
- **Albertsen & Hansen 1994.** `[VERIFY exact citation]` **Kronecker/product-state construction for
  multi-channel recordings, in this field** (E14). Concede: the doubled-generator idea is not an import
  from queueing theory, it is already at home in channel work.
- **Zhou & Lange 2009**, *Adv. Appl. Probab.* 41(1), *composition Markov chains of multinomial type*. PDF:
  `Zhou_Lange_2009_Composition_Markov_Chains.pdf`. **Title corrected; it does not state the identity**
  (E2). Keep only if you actually use composition chains.
- **Plateau & Stewart; Buchholz; Neuts** (stochastic automata networks / Kronecker-structured CTMCs).
  **False friend, worth a footnote:** there the same `Q ⊕ Q` composes **independent** subsystems (product
  space, build-up), the opposite of MacroIR's correlated single-chain two-time (start, end) joint. Even the
  experts who own the algebra were not pointing it at this object.

### Matrix-exponential integrals
- **Van Loan 1978** (see B1). The primitive.
- **Carbonell, Jiménez & Pedroso.** PDF: `carbonell2008.pdf`. **See E5: this is a generalization of Van
  Loan's method for computing multiple integrals involving matrix exponentials**, not the
  local-linearization filtering paper the old entry described. **Read it before citing:** it may already
  contain the exact block-exponential generalization MacroIR uses, in which case it becomes a `concede`
  and a convenience, not a competitor. `[ACTION: read]`

---

## E. Claim / concede / cite — Part I summary (REVISED)

| What | Status | Key references |
|------|--------|----------------|
| Integrated-measurement (time-averaged) treatment **at macroscopic scale, at O(k³) independent of N** | **CLAIM** (scaling, not absence: see E11) | gap vs Münch 2022, Moffatt 2007, Milescu 2005, Celentano-Hawkes 2004, Stepanyuk 2011/2014, Del Core-Mirams 2025; **vs Qin 2000 / Venkataramanan-Sigworth on cost** |
| Exact, efficient, all-k realization for the **many-channel ensemble** CTMC | **CLAIM**, narrowed | distinguish from Bäuerle (single chain, needs a density with an atom), Blackwell (single chain, ecology), **Kilic 2021 (single molecule, Poisson, MCMC)** |
| Empirical demonstration that instantaneous methods are **misspecified** | **CLAIM** (this is the eLife paper; see Part II for its prior art) | figures 2/3/4 |
| The recording filter modeled inside the likelihood | **CONCEDE — since 1992** | Fredkin & Rice 1992; Michalek 1999; Venkataramanan & Sigworth; **Qin, Auerbach & Sachs 2000** |
| Integrated-/averaged-measurement Kalman filter as a **device** | CONCEDE | Kalman 1960; Van Loan 1978; Zadrozny 1988; Harvey 1989; **Fatehi & Huang 2017**; Yaghoobi & Särkkä 2024; **Folia & Rattray 2018**; Calderazzo 2019; Rubenzahl 2026 |
| van-Loan-on-augmented-integral pattern | CONCEDE | Rubenzahl 2026 |
| Kronecker/product structure in channel modeling | **CONCEDE — already in the field** | **Albertsen & Hansen 1994** |
| `exp(A⊕B) = exp(A)⊗exp(B)` | CONCEDE as classical | **Horn & Johnson 1991 §4.4** (not Zhou-Lange, not Higham) |
| "Exact CTMC" as accuracy over LNA (for channels) | **DO NOT CLAIM** (they coincide; verified 1e-8) | `claude/scripts/verify_IR_vs_augmentation.py` |
| k² meta-state = (k+1) augmented-state equivalence | unification, ingredients conceded | Horn & Johnson; SAN literature (false friend) |
| CTMC reward-moment machinery | CONCEDE as standard | Hobolth-Jensen 2011 (endpoint-conditioned!); Bladt-Meini-Neuts-Sericola 2002 |
| "Likelihood approximations are in routine use" | **FALSE, DO NOT WRITE** | Clerx 2019 (four ways, all deterministic); Owen & Mirams IonBench ("we do not include ... stochastic models") |

*Sources: two adversarial literature workflows (2026-06-26), task IDs wg0zu622u and wxa9rb7v8, plus the
numerical equivalence check in `claude/scripts/verify_IR_vs_augmentation.py`; corrected and extended by a
four-agent exhaustive review (2026-07-14) that read the PDFs.*

---
---

# PART II — prior art for the DIAGNOSTIC (the eLife paper's actual claim)

Written 2026-07-14 from a four-cluster adversarial search. Part I defends the algorithm; the eLife paper
claims a **diagnostic**. **Almost none of the diagnostic is new.** This part says exactly which parts are
conceded, to whom, and what is left.

**Notation bridge (memorize this; the reviewer has).** In the misspecification literature:
`H` = **sensitivity matrix** (the model's own information / expected negative Hessian);
`J` = **variability matrix** (covariance of the score); `G = H J⁻¹ H` = **Godambe information**;
`G⁻¹ = H⁻¹ J H⁻¹` = the **sandwich**. Our `C = H^{-1/2} J H^{-1/2}` has exactly the eigenvalues `λᵢ` of
`H⁻¹J`, and `tr C = tr(H⁻¹J)`. Symmetrizing adds nothing mathematically. **`C = I` is Lindsay's
"information unbiasedness"; the failure of `C = I` has been the definition of likelihood misspecification
since 1977.**

## II.0 Bottom line

1. **Every one of the three tests is a named classical test.** (i) residual whiteness = innovation
   whiteness / filter consistency (Mehra 1970; Bar-Shalom 2001; Ljung 1999). (ii) zero-mean score =
   **first Bartlett identity**. (iii) score covariance = reported Fisher = **second Bartlett identity** =
   White's **information matrix equality**. Chesher, Dhaene, Gouriéroux & Scaillet (1999) is the single
   cite for "test the Bartlett identities".
2. **`C` is not a new object.** Its spectrum is the weighted-χ² weight set of the misspecified LRT (Foutz
   & Srivastava 1977; Kent 1982). Its determinant and trace are **published test statistics**: the Log-Det
   and Log-Trace **Generalized Information Matrix Tests** (Golden, Henley, White & Kashner 2016). Its trace
   is the **effective number of parameters** (TIC; Varin-Reid-Firth 2011).
3. **α⋆ = p / tr C is published verbatim, twice, as a Bayesian power on the likelihood.** Pauli, Racugno &
   Ventura (2011) eq. 2.3; Ribatet, Cooley & Davison (2012) §2.1 ("magnitude adjustment", `k = p / Σλᵢ`).
   Frequentist roots: Rotnitzky & Jewell (1990), Geys, Molenberghs & Ryan (1999). **Claiming this is novel
   is a first-paragraph kill.**
4. **½ log det C as an evidence correction is published**, as Lv & Liu (2014, JRSS-B) Theorem 1 and their
   **GBIC**, where `H_n = A_n⁻¹ B_n` is called the *covariance contrast matrix*. Their **GBIC_p** carries
   `tr(Ĥ_n) − log det(Ĥ_n)`, i.e. **both** invariants of `C` in one evidence criterion. It is also
   **contested** (Bhattacharya & Pati 2020) on the grounds that the sandwich does not appear in the
   posterior limit (Kleijn & van der Vaart 2012).
5. **The correlation/sample decomposition is the HAC long-run-variance decomposition**, and its scalar
   `κ_det = det(R)^{1/d}` is *literally* the Vats-Flegal-Jones multivariate-ESS ratio.
6. **In-domain: Münch et al. 2022 (eLife) already published diagnostic (i) plus an N_ch validity rule of
   thumb for a Kalman ion-channel likelihood** (verified from the paper: residual ACF is white; *"one
   should be careful with both algorithms for time traces with N_ch ∈ [10¹, 10²]"*; `error(N_ch) ∝
   a/N_ch`). They do **not** compute the score, the Fisher information, the score covariance, or any
   sandwich. **That is the delta.**

## II.1 The concession table — what we call it vs what statistics calls it

| Our name | Standard name | Owner (cite these) |
|---|---|---|
| The three-part test | **Bartlett identities test** (IM test = 1st-order case) | Chesher, Dhaene, Gouriéroux & Scaillet 1999, LIDAM/IRES DP 1999019 |
| (i) standardized residuals white, unit variance | **innovation whiteness / filter consistency**; NIS/NEES | Mehra 1970 *IEEE TAC* 15(2):175-184; Bar-Shalom, Li & Kirubarajan 2001; Ljung 1999; Desroziers et al. 2005 *QJRMS* 131:3385. **In our own domain: Münch et al. 2022 eLife 11:e62714, Box 1-fig 2** |
| (ii) E[score] = 0 at truth | **first Bartlett identity**; QMLE pseudo-true value θ⋆ ≠ θ₀ | Huber 1967; White 1982 *Econometrica* 50(1):1-25 |
| (iii) Cov(score) = reported Fisher | **information matrix equality** / second Bartlett identity | White 1982; White 1994 ch. 11 |
| Information distortion matrix `C = H^{-1/2} J H^{-1/2}` | **the sandwich, symmetrically normalized.** Spectrum = eig(H⁻¹J) = the weights of the misspecified LRT (Kent; Vuong's `W`); `C = I` = **information unbiasedness** | Godambe 1960; Foutz & Srivastava 1977 *AoS* 5:1183; Kent 1982 *Biometrika* 69:19-27; Lindsay 1982; Vuong 1989 *Econometrica* 57:307. **No standard name for the symmetric form**; closest are the "information ratio" (Zhou, Song & Thompson 2012 *JASA* 107:205) and the "covariance contrast matrix" (Lv & Liu 2014) |
| eigenvalues `λᵢ(C)` | **generalized design effects** (survey statistics) — `Δ̂ = Cov_model⁻¹ Ĉov_true`, "the eigenvalues of Δ̂ are termed generalized design effects"; Rao-Scott 1st-order correction divides χ² by the **mean** eigenvalue (= our α⋆), 2nd-order is a Satterthwaite match on Σλ, Σλ² | **Rao & Scott 1981** *JASA* 76:221-230; **Rao & Scott 1984** *AoS* 12:46-60 |
| `ĉ` variance-inflation factor (QAIC) | the **isotropic collapse** of the same spectrum: QAIC = TIC under the restriction `H⁻¹J = ĉ·I` | Burnham & Anderson 2002 |
| `log det C`, `tr C` as scalar summaries | **Log-Det GIMT and Log-Trace GIMT** (published test statistics; GIMT_Det = (1/q) log det(A⁻¹B)) | Golden, Henley, White & Kashner 2016 *Econometrics* 4(4):46; Golden et al. 2013 (Springer, White festschrift) |
| DCC = `H⁻¹ J H⁻¹` | **Godambe / sandwich covariance** | Godambe 1960 *AoMS* 31:1208; Huber 1967; White 1982 |
| DIB = `H⁻¹ E[score]` | first-order Newton step to the **pseudo-true value** θ⋆ (KL minimizer) | White 1982 |
| `tr C` as "effective sample" | **effective number of parameters**; TIC penalty; CLIC/CLAIC/CLBIC | Takeuchi 1976; Varin & Vidoni 2005 *Biometrika* 92:519; Varin, Reid & Firth 2011 *Stat. Sinica* 21:5-42 (p. 11: `dim(θ) = tr{H G⁻¹}`); Ng & Joe 2014 *Bernoulli* 20:1738 |
| Peak correction `α⋆ = p / tr C` | **magnitude adjustment / calibrated composite likelihood**, the same constant | **Pauli, Racugno & Ventura 2011** *Stat. Sinica* 21:149-164 eq. (2.3); **Ribatet, Cooley & Davison 2012** *Stat. Sinica* 22:813-845 §2.1; Rotnitzky & Jewell 1990 *Biometrika* 77:485; Geys, Molenberghs & Ryan 1999 *JASA* 94:734 |
| Volume correction `½ log det C` | **GBIC**: the misspecification term in the Laplace expansion of the log marginal likelihood | **Lv & Liu 2014** *JRSS-B* 76(1):141-167, Thm 1 + GBIC/GBIC_p |
| Correlation distortion `R = J_sample^{-1/2} J_total J_sample^{-1/2}` | **HAC / long-run variance** `Σ = Γ₀ + Σ(Γ_k + Γ_kᵀ)` | Newey & West 1987; Andrews 1991 |
| `κ_det = det(R)^{1/d}` | **multivariate ESS ratio** `n/mESS`, `mESS = n(|Λ|/|Σ|)^{1/p}` | **Vats, Flegal & Jones 2019** *Biometrika* 106(2):321-337 |
| Sample distortion `C_sample = H^{-1/2} J_sample H^{-1/2}` | **information ratio** (variability vs sensitivity), OPG-vs-Hessian form of the IM test | Zhou, Song & Thompson 2012 *JASA* 107(497):205-213 |
| Power-likelihood framing | **generalized Bayes** | Bissiri, Holmes & Walker 2016 *JRSS-B* 78:1103; Grünwald & van Ommen 2017 *Bayesian Anal.* 12:1069 |
| "Simulation as ground truth for an approximate likelihood" | closest genres: SBC (algorithm-vs-model, **cannot** see this), ABC/SBI misspecification diagnostics (summary-level) | Talts et al. 2018; Frazier, Robert & Rousseau 2020 *JRSS-B* 82:421; Hermans et al. 2022 *TMLR* |

## II.2 What is actually left (defensible novelty, ranked)

1. **First-in-domain application.** Nobody has run the information-matrix-equality machinery on an
   ensemble-CTMC / macroscopic-current likelihood. Münch 2022 got as far as residual whiteness and
   stopped. Same contribution shape as Part I: bring the established tool, show it is needed.
2. **The population quantities are computed exactly, not estimated from one dataset.** The IM test's famous
   defect is its finite-sample behaviour (true vs nominal level off by 10×: Taylor 1987; Kennan & Neumann
   1988; Orme 1990 *J. Econometrics* 46:309; fixed by bootstrap, Horowitz 1994 *J. Econometrics* 61:395).
   Because the forward process is exactly simulable, we get `H` analytically and `J` by Monte Carlo over
   replicates. **We are not *testing* misspecification against a null; we are *measuring* a known
   misspecification.** Different problem, and it dodges the defect. **This is the strongest framing and
   should be stated in the Introduction.**
3. **The multiplicative decomposition `C = C_sample^{1/2} R C_sample^{1/2}` as a diagnostic pair.** Both
   factors are standard separately (IM test; HAC). Attributing an algorithm's distortion to
   *within-interval non-Gaussianity* vs *cross-interval score correlation*, and reading the second as the
   fingerprint of the discarded higher moments, appears to be unstated. Small, but real. Scope it tightly.
   **Bonus: Milescu 2005 names "the local time correlation of the current" as his dominant error source,
   which is exactly the `R` factor.**
4. **The full spectrum, not a scalar collapse.** Every consumer of `H⁻¹J` in the literature collapses it:
   TIC takes the trace, QAIC's `ĉ` assumes `C ∝ I`, Rao-Scott takes the mean eigenvalue, Satterthwaite
   takes `(Σλ)²/Σλ²`. Reporting the anisotropy (which *parameter directions* are inflated) is a legitimate,
   if modest, methodological point, and the data support it (the inflate/deflate directionality claim in
   the abstract *is* the anisotropy).
5. **The physical identification of `C_sample` with the third and fourth cumulants of the interval
   current.** White (1982, p. 7) proves that for a Gaussian model the IM equality holds iff skewness = 0
   and kurtosis = 3; our own cumulant expansion says the same. So `C_sample` is not an abstract
   discrepancy, it is a *measurement of the non-Gaussianity of the channel-ensemble current*, and it maps
   onto the multinomial and telegraphic corners for a reason. **That reading is ours.**
6. **The regime map.** A 3-axis (N_ch, Δ/τ, noise) quantitative map for five likelihood approximations.
   The genre exists (Schnoerr, Sanguinetti & Grima 2014 for moment closure; Grima et al. 2011 for
   CFPE/LNA; **Dalmasso et al. 2020** for parameter-indexed coverage — *year corrected, see II.5*), and
   **Münch 2022 published a 1-axis verbal version for this exact model class.** An increment, not a
   discovery. **Do not lead with it.**

## II.5 ADDENDUM (second pass, 2026-07-14) — the martingale argument, and the one that hurts

Added after a targeted check of three questions Part II did not cover: who owns the martingale-difference
property of the prequential score; how much of the non-recursive story the composite-likelihood literature
already owns; and whether the regime-map genre has reached time-integrated observations. Derivation:
`theory/macroir/docs/Likelihood_Information_Distortion/score_martingale_argument.md`.

### II.5.1 The MDS proposition is textbook, and it has an owner in our exact model class

O5's framing recommendation ("the per-step score of a correctly specified filter is a martingale difference
sequence") had **no citation attached**. It has one now, and it is squarely in hidden-Markov territory:

- **Cappé, Moulines & Rydén (2005), *Inference in Hidden Markov Models*, Springer, ch. 12** (statistical
  properties of the MLE). Verbatim: *"we want to prove a central limit theorem (CLT) for the score function
  evaluated at the true parameter. A quite general way to do that is to recognize that the corresponding
  score increments form, under reasonable assumptions, a martingale increment sequence with respect to the
  filtration generated by the observations."* Their limiting Fisher information is a **single per-step
  second moment with no cross-lag terms**, which is exactly `J_total = J_sample`, i.e. `R = I`.
- Primary sources: **Bickel, Ritov & Rydén (1998)** *AoS* 26(4):1614; **Douc, Moulines & Rydén (2004)** *AoS*
  32(5):2254. General (non-HMM): **Hall & Heyde (1980)**, *Martingale Limit Theory*, ch. 6;
  **Barndorff-Nielsen & Sørensen (1994)** *Int. Statist. Rev.* 62:133; **Godambe (1985)** *Biometrika*
  72:419. Prediction-error decomposition: **Schweppe (1965)**; **Harvey (1989)**. Prequential:
  **Dawid (1984)** *JRSS-A* 147:278.
- **The converse** (non-zero cross-lag score covariance signals a wrong one-step predictive) is the null of
  **White (1987), "Specification testing in dynamic models"**, in Bewley (ed.), *Advances in Econometrics:
  Fifth World Congress* Vol. 1, CUP, pp. 1-58 (the **dynamic information matrix test**) `[VERIFY: full text
  not retrieved]`, and of the m-testing / conditional-moment framework (**Newey 1985** *Econometrica*
  53:1047; **Tauchen 1985**; **Wooldridge 1990** *Econometric Theory* 6:17). In the prequential frame:
  **Seillier-Moiseiwitsch & Dawid (1993)** *JASA* 88:355.

**Claim status:** state the proposition **with citation, as a proposition, not as a result.** What is ours is
using the cross-lag score covariance as a *quantitative per-algorithm measure* instead of a test statistic,
and the exact identity `Cov(s_i, s_j) = Cov(s_i, b_j)` with `b_j = E[s_j | F_{j-1}]` (the conditional score
bias), which says **all correlation distortion is the error of the one-step predictive**.

### II.5.2 O8 (NEW, and it is the one that hurts). Milescu, Akk & Sachs 2005 already owns the phenomenon

Part II's entry for Milescu (II.2, item 3) says only that he *"names the local time correlation of the
current as his dominant error source"*. **That badly understates it, and the paper is in the repo**
(`docs/bibliography/Milescu_2005_MLE_Macroscopic_Currents_BiophysJ.pdf`, *Biophys. J.* 88(4):2494-2515).
He does four things:

1. **States the independence assumption explicitly** (p. 2497): *"We assume that the state probability
   vector P_t is a function of the initial state probabilities P_0, but does not depend on the state
   probabilities at any other time … which is equivalent to saying that the data have no autocorrelation
   beyond that contained in the A matrix. Thus, the conditional probability of the entire current trace is
   equal to the product of the individual conditional probabilities."* That is our `NR`, and he wrote down
   its defect.
2. **Diagnoses the failure as autocorrelation.** Fig. 10 caption, verbatim: **"The bias of the estimates is
   caused by autocorrelation."**
3. **Proves it by a decorrelation experiment** (p. 2504): *"we hypothesized it might be a result of the
   algorithm ignoring the local time correlation of the data. To test this, we generated data without
   autocorrelation (Fig. 10). … The parameter estimates and the mean currents obtained from the
   noncorrelated traces showed no bias, and higher precision."*
4. **Predicts the fix** (p. 2504): *"We speculate that only a Bayesian filtering algorithm … could obtain
   unbiased estimates from macroscopic data, but this would entail a high computational cost."*

**Consequence for the manuscript. Do NOT write "we show that ignoring temporal correlation degrades
inference" as a discovery sentence.** Milescu owns the phenomenon, in this model class, twenty years ago.
Concede it in the Introduction, in his own words, and then say what is actually new:

> Milescu et al. (2005) identified the temporal correlation of the current as the dominant error source of
> the non-recursive likelihood, demonstrated it by decorrelating the data, and predicted that a filtering
> algorithm would remove it. What was not available then, and is what we supply, is a *measurement*: how
> much information such a likelihood over-reports (10 to 16 fold), the fact that the per-interval
> information identity Cov(s_t) = F_t holds exactly for every member of the family so that the entire
> failure is carried by cross-time score correlation, and how the inflation depends on channel number,
> sampling interval and noise.

**And we can correct him, which is better than agreeing with him.** He attributes the *bias* to
autocorrelation. Our decomposition says bias and information inflation are governed by **different axes**:
the per-interval predictive moments (the averaging axis) govern the score bias, and the conditioning
(recursion) axis governs the cross-time correlation and hence the information inflation. `NMR` is the
counterexample that separates them: non-recursive, fitted to autocorrelated data, and **unbiased** — because
its per-interval density is normalized with the correct moments, so E[s_t] = 0 exactly, autocorrelation or
not. So the honest answer to Milescu's speculation is **no**: a non-filtering method need not be biased. It
will, however, be over-confident. `[VERIFY: confirm NMR's unbiasedness on the final Gaussian run before
this sentence is written; it is load-bearing.]`

### II.5.3 The non-recursive members are independence composite likelihoods. Concede it.

Missing from the II.1 concession table: Chandler & Bate and Varin-Reid-Firth appear there only as sources
for the *adjustment*, not as the owners of the **framing**.

- **Chandler & Bate (2007)** *Biometrika* 94(1):167-183 coined the **"independence loglikelihood"** and the
  sandwich-based adjustment.
- **Varin, Reid & Firth (2011)** *Statistica Sinica* 21:5-42, §2.2, p. 8, verbatim and decisive:
  *"Composite likelihoods may be seen as misspecified likelihoods, where misspecification occurs because of
  the working independence assumption among the likelihood terms forming the pseudolikelihood. Consequently,
  the second Bartlett identity does not hold, and we need to distinguish between the sensitivity matrix H(θ)
  … and the variability matrix J(θ) … and the Fisher information needs to be substituted by the Godambe
  information matrix."*
- **Ribatet, Cooley & Davison (2012)**, on the overconfidence, verbatim: *"they mention that this
  substitution may lead to overly precise inferences, they do not describe how to correct this"* — treated
  as common knowledge.

**What the objection eats:** the qualitative overconfidence, its sign, and the sandwich fix.
**What it does not eat:** the magnitude; the exact per-interval calibration result (Cov(s_t) = F_t for all
five, so the failure is *entirely* cross-time); the regime dependence; and the identification of the
published ion-channel likelihoods as members of this class. **Not found anywhere:** the overconfidence of an
independence likelihood as a *quantified law* with a scaling. The composite-likelihood literature states
H ≠ J and stops. The scaling κ = 1 + 2Σρ_k is the HAC long-run-variance identity (Newey & West 1987; the
multivariate ESS form is Vats, Flegal & Jones 2019), not a composite-likelihood theorem. **Our prediction
κ ≈ 2τ/Δ, and its test against the interval sweep, is ours** (see the derivation note, §4).

### II.5.4 The regime-map claim survives, narrowly

- **Schnoerr, Sanguinetti & Grima (2014)** *JCP* 141:084103 and **(2015)** *JCP* 143:185101 map moment-closure
  validity over **mean molecule number only**: no sampling-interval axis, no noise axis, and the metric is
  *accuracy of moments*, not calibration of a likelihood.
- **Gillespie & Golightly (2014)**, arXiv:1409.1096, *SAGMB*, is the closest genuine hit: standardized
  prediction errors and **confidence-interval coverage** of the approximation against exact simulation over a
  space-filling design. Still predictive coverage, not the information-matrix equality, and no (N, Δ, noise)
  grid.
- **Dalmasso, Lee, Izbicki, Pospisil, Kim & Lin (2020)**, AISTATS, PMLR 108:3349-3361 (**not 2024**; fix the
  citation in II.2 item 6): parameter-indexed validation of approximate likelihoods, local two-sample tests.
- **No study in this genre uses a time-integrated observation.** That axis is unoccupied.

**Paste-able claim, stated as narrowly as it must be:**
> Validity maps for approximate likelihoods of partially observed Markov jump processes exist, but they map
> the accuracy of moments over system size (Schnoerr et al. 2014, 2015) or the coverage of an approximate
> predictive over parameter space (Gillespie & Golightly 2014; Dalmasso et al. 2020). None maps the
> information-matrix equality, none compares a family of likelihood approximations on a common axis set, and
> none uses a time-integrated observation.

**The claim dies if it is stated as "nobody has produced validity maps".** State it narrowly, and cite
Gillespie & Golightly and Dalmasso as the nearest genre.

## II.3 Must-cite (omitting any of these is an unforced error)

Huber 1967 (Corollary to Thm 3, p. 231, the sandwich); **White 1982** (Thm 3.3 = information matrix
equivalence; §4 = the IM test; p. 7 = the Gaussian skew/kurtosis case); White 1994 (ch. 11);
**Godambe 1960**; Foutz & Srivastava 1977; **Kent 1982**; Vuong 1989; **Rao & Scott 1981/1984**
(generalized design effects); Takeuchi 1976; Watanabe 2010 (ν = tr(H⁻¹J)); Chesher, Dhaene, Gouriéroux &
Scaillet 1999; Orme 1990 + Horowitz 1994 (IM-test finite sample); Freedman 2006; **Golden, Henley, White &
Kashner 2016** (log-det / trace GIMT); Zhou, Song & Thompson 2012; **Varin, Reid & Firth 2011**; Rotnitzky
& Jewell 1990; Geys, Molenberghs & Ryan 1999; **Pauli, Racugno & Ventura 2011**; **Ribatet, Cooley &
Davison 2012**; Chandler & Bate 2007; Pace, Salvan & Sartori 2011; **Lv & Liu 2014**; Bissiri, Holmes &
Walker 2016; Grünwald & van Ommen 2017; Lyddon, Holmes & Walker 2019; **Kleijn & van der Vaart 2012**;
**Müller 2013** (*Econometrica* 81:1805); Bhattacharya & Pati 2020; Newey & West 1987; **Vats, Flegal &
Jones 2019**; Bollerslev & Wooldridge 1992; Douc & Moulines 2012; Mehra 1970; Bar-Shalom, Li &
Kirubarajan 2001; Ljung 1999; Talts et al. 2018; Frazier, Robert & Rousseau 2020; **Münch et al. 2022**;
Milescu, Akk & Sachs 2005; Moffatt 2007.

## II.4 The objections a statistician referee will raise

**O1 (fatal if unanswered). `H` is the wrong matrix for the sandwich.** White's sandwich is `A⁻¹ B A⁻¹`
with `A = −E_true[∇²ℓ]`, the expected Hessian **under the true law**, at the pseudo-true `θ⋆`. Our `H` is
the model's **own** Gaussian Fisher `Σ_t [(∂μ)ᵀ(∂μ)/σ² + (∂σ²)ᵀ(∂σ²)/2σ⁴]`, i.e. `−E_model[∇²ℓ]`. Under
misspecification these differ, and the difference is *itself* part of what is being measured. So
`C = H^{-1/2} J H^{-1/2}` is a legitimate **IM-equality test statistic** (under the null all three matrices
coincide) but `H⁻¹ J H⁻¹` is **not** the asymptotic covariance of θ̂ unless `A ≈ H`. **Answerable, and
cheaply:** the codebase already computes both the numerical FD Hessian `F_b` and the bridge
`GFD = G_b^{-1/2} F_b G_b^{-1/2}`. Report GFD; if it is ≈ I, the substitution is harmless and say so; where
it is not, the DCC must be anchored on `F_b`. **Do not ship the paper without this check.**

**O2. The anchor point.** White's theory holds at `θ⋆` (the KL minimizer), not at `θ_sim`. That
`E[score] ≠ 0` at `θ_sim` is precisely the statement `θ_sim ≠ θ⋆`. So the DIB is an estimate of
`θ⋆ − θ_sim` (the misspecification displacement), and a sandwich anchored at `θ_sim` is not the asymptotic
covariance of the estimator. Two different objects, both wanted, and the paper must say which is which.
(`battery_sim` vs `battery_pool` in the data is exactly this split.)

**O3. The evidence correction is an imposition, not an approximation.** The Laplace expansion of
`∫L(θ)π(θ)dθ` is already *correct* with the observed Hessian, misspecified or not; under Kleijn & van der
Vaart (2012) the misspecified posterior concentrates with covariance `H⁻¹`, **not** the sandwich, so
credible sets are not confidence sets. Bhattacharya & Pati (2020) make this objection in print against Lv
& Liu's GBIC. Adding `½ log det C` does not compute the same integral better; it computes a **different,
generalized-Bayes object.** That is defensible, but the warrant must be decision-theoretic (Müller 2013:
the sandwich posterior has lower asymptotic frequentist risk; Bissiri et al. 2016), not "we improved the
Laplace approximation". Also: McAlinn & Takanashi (2026, arXiv:2602.01573) argue a tempered normalizer is
not evidence and Bayes factors from it are not well defined. **State the warrant, or cut the claim.**

**O4. Double counting, and the scalar power is known to be the weak member.** Applying `α⋆` already
rescales the mode curvature to `α⋆H`, shifting the Laplace volume by `−(p/2) log α⋆`; the peak and volume
corrections coincide only when all `λᵢ(C)` are equal (Pace, Salvan & Sartori 2011 make exactly this point:
the moment-matched, Satterthwaite and Chandler-Bate adjustments coincide only when `d = 1`). And a
**scalar** power provably cannot deliver the sandwich covariance: the magnitude-adjusted posterior is
`N(θ₀, (np)⁻¹ tr(H⁻¹J) H⁻¹)`, still proportional to `H⁻¹` (Ribatet et al. 2012 eq. 13; Miller 2021 *JMLR*
22(168)). Only an **affine/curvature** adjustment does (Chandler & Bate 2007: `CᵀHC = HJ⁻¹H`). If we keep
the scalar, we must say why (cost, PSD, identifiability of the affine map). **Note also that
`α⋆ ≠` Lyddon-Holmes-Walker's `w⋆ = tr(HJ⁻¹Hᵀ)/tr(H)`:** different rules, equal only when `J ∝ H`. **Do
not cite LHW as the source of `α⋆`; the sources are Pauli et al. / Ribatet et al.**
**ACTION: `theory/macroir/docs/Likelihood_Information_Distortion/supplement_evidence_correction.tex`
currently mis-cites LHW as the source. Fix it.**

**O5. The decomposition is HAC.** `J_total = J_sample + cross-lag terms` is the textbook long-run variance.
Frame it as: *the per-step score of a correctly specified filter is a martingale difference sequence, so
any non-zero cross-lag term is itself a misspecification signal, and we use the HAC decomposition to
attribute it.* That framing survives; "we introduce a decomposition" does not.

**O6 (the one nobody saw coming). For a Gaussian model, the IM test IS a skewness/kurtosis test.** White
1982 works the normal case explicitly (p. 7): for `N(μ,σ²)`, `A(θ₀) = −B(θ₀)` **iff** skewness `√β₁ = 0`
and kurtosis `β₂ = 3`. Our emission density is Gaussian by construction, so a referee can say *"your
`C_sample = I` check is a Jarque-Bera normality test in a costume."* **This is true, and it is good news,
not bad.** It is exactly what the cumulant expansion in
`figures/paper/sample_correlation_distortion_analysis.md` already found (`J_sample − G_b` = skew +
kurtosis terms, exactly). **Say it first:** the sample distortion **is** the third and fourth cumulants of
the interval current, and that is why it tracks the multinomial and telegraphic corners. Owning this turns
an attack into the physical interpretation of the diagnostic.

**O7. Freedman's objection, which we can turn.** Freedman 2006 (*Am. Statist.* 60:299-302): *"It remains
unclear why applied workers should care about the variance of an estimator for the wrong parameter."* His
target is people who patch the variance with a sandwich and carry on. **Cite him approvingly:** our point
is that `C ≠ I` signals `θ⋆ ≠ θ_true`, i.e. bias, which is exactly the thing he says everyone ignores, and
the DIB measures it. A diagnostic beats a sandwich patch.

**Two numerical land mines, both already in print.** (1) White 1982 §4: the IM test's covariance `V(θ⋆)`
includes the estimation-error correction `∇D A⁻¹ ∇log f`; a naive covariance of the indicators is wrong,
and White himself warns the simplified estimator "is neither consistent nor necessarily positive
semi-definite when the null fails." (2) Ward 2023 (*Entropy* 25:512): TIC is not used in practice because
`tr(H⁻¹J)` is numerically hopeless above ~20 parameters (98-100% singular Hessians). Report `p` and the
conditioning; with `p = 6` we are safe, and saying so loudly is the defence.

**Name collision (affects the title).** "Information distortion" is already an established term in *neural
coding* (Dimitrov & Miller: rate-distortion clustering of stimulus-response classes), in an audience eLife
shares. Not fatal, but the chosen title
(*"Information distortion in likelihood approximations for macroscopic ion-channel currents"*) inherits it.
See `papers/macroir-elife-2025/title_options.md`. **Open decision.**

*Sources: four adversarial search clusters (2026-07-14) covering composite-likelihood adjustment, the
White/Huber/IM-test lineage, power-likelihood + evidence under misspecification, and filter diagnostics /
HAC / simulation-based calibration; plus direct fetches of Münch et al. 2022, the GIMT definitions, and the
Chesher et al. Bartlett-identities paper.*

---
---

# PART III — what to claim, in what order

The two parts above pull in opposite directions, and the paper has to resolve them.

**What got weaker (Part I).** The algorithm's novelty is narrower than we thought: the recording filter has
been inside single-channel likelihoods since 1992 (E11), the exact integrated observation for arbitrary `K`
exists in single-molecule work (E12, Kilic 2021), and the augmented-integral device is standard (B2). What
survives is **scaling** (`O(k³)`, independent of `N`, versus `k^L` and combinatorial) and **the ensemble**.

**What got weaker (Part II).** The diagnostic is a rediscovery of the White/Godambe information-matrix
machinery, almost item for item, including both corrections (`α⋆` and `½ log det C`) as published objects.

**What survives, and it is enough for one good paper.** Two things, and they should be said in this order:

1. **We are measuring a known misspecification, not testing an unknown one.** The forward process is
   exactly simulable, so `H` is analytic and `J` is Monte Carlo over replicates. The literature's version of
   this object is a *test statistic* with a notorious finite-sample defect; ours is a *measurement* with no
   null hypothesis in sight. Nobody has been able to do this before because nobody had a likelihood whose
   data-generating process they could sample exactly. **This is the sentence the Introduction is built on.**
2. **The distortion has a physical identity.** `C_sample` is the third and fourth cumulants of the interval
   current (O6); `R` is the cross-interval score correlation, which is exactly the error source Milescu
   named in 2005. So the diagnostic does not merely say "the likelihood is wrong", it says **which physical
   feature of the channel ensemble the approximation threw away, and in which parameter direction.**

**And what to concede, in the Introduction, before anyone asks:** the sandwich is White's; the corrections
are Pauli's and Lv & Liu's; the whiteness test is Münch's; the integrated-measurement filter is a device
from control theory; the exact filtered likelihood is Qin's and Sigworth's for one channel.
**Conceding all of that in one honest paragraph costs a hundred words and buys the referee's trust for the
whole paper.**

**Open items before submission.**
- Report GFD (`G_b^{-1/2} F_b G_b^{-1/2}`) to answer O1. **Blocking.**
- Fix the LHW mis-citation in `supplement_evidence_correction.tex` (O4).
- Decide the evidence-correction warrant (O3) or cut the claim.
- Verify the **MR sign**: `00_master_plan_v2.md` says MacroMR *overestimates* variance; the Figure_S3
  caption reports MR as *over*-confident. One of them has the sign backwards, and the abstract's
  bidirectionality sentence depends on it.
- Retrieve and read: **Kilic 2021** (both papers), **Qin, Auerbach & Sachs 2000**, **Requadt 2025**,
  **`carbonell2008.pdf`** (which is not what we thought it was).
- Revisit the title given the "information distortion" name collision.
