# Use of the Covariance Matrix in Directly Fitting Kinetic Parameters: Application to GABA_A Receptors

**Citation.** Celentano JJ, Hawkes AG (2004). *Use of the Covariance Matrix in Directly Fitting Kinetic Parameters: Application to GABA_A Receptors.* Biophysical Journal 87(1):276–294. DOI: 10.1529/biophysj.103.036632. Submitted 29 Oct 2003, accepted 30 Mar 2004. Repo PDF: `docs/bibliography/CelentanoHawkes_2004_Covariance_Matrix_Kinetic_Fitting_BiophysJ.pdf`.

**Verification status.** Verified. Title, authors, journal, volume, pages, and DOI confirmed against both the repo PDF title page and external indexers (ScienceDirect, PMC1304350). All quotations below are verbatim from the PDF with page locators. The one identifiability limit that is most load-bearing for the caveat (linear vs branched topology) is shown empirically, not proven; see Caveats.

## Core result

The paper introduces **covariance fitting (CVF)**: fit macroscopic (many-channel) desensitization currents to a kinetic model by maximum likelihood using not only the mean current trajectory but the full covariance between pairs of time points. The mean and covariance come from the Colquhoun-Hawkes (1977) relaxation-and-fluctuation theory for a population of N independent channels:

- Mean (Eq. 3, p. 292): `Current(t) = N · p₀ · e^{Qt} · Γ · u`, where `p₀` is the initial state distribution, `Q` the rate matrix, `Γ` a diagonal conductance matrix (**0 for closed and desensitized states**), `u` a column of ones.
- Covariance (Eq. 4, p. 293), for `0 ≤ t₁ ≤ t₂`: `Cov(t₁,t₂) = N·[ (p₀ e^{Qt₁} Γ e^{Q(t₂−t₁)} Γ u) − (p₀ e^{Qt₁} Γ u)(p₀ e^{Qt₂} Γ u) ]`.

The variance (diagonal) is not constant: it grows as the open fraction nears 50%, so the mean-only sum-of-squares error model is misspecified. Coupled with a likelihood-ratio test, CVF discriminates models with up to 11 free parameters. Applied to α1β2γ2S GABA_A traces it favors a model with two open and three desensitized states.

The result relevant here is the split between **what the covariance adds to the mean** and **what stays unidentifiable even with it**.

## What is / is not identifiable

**What the mean already gives (sum-of-squares / variance-weighted SS = wSS).** Magnitude of the current and the number of exponential decay phases, hence roughly the number of desensitized states. wSS "only takes into account the magnitude of the current at each point in time" and "fails to discriminate between differences in the arrangement of states or in the pattern of desensitization" (p. 290).

**What the covariance adds (identifiable by CVF but not by mean-only fitting):**
- The **arrangement (topology) of states** at equal state count and equal parameter count. "Comparing models with the same number of states and free parameters is made possible by the use of the full covariance matrix" (p. 290).
- The **pattern of desensitization**, i.e. which desensitized state produces which kinetic phase (proximal-vs-distal assignment). CVF "heavily favors the correct pattern"; "When applied to 3D-I,pdd, wSS fails to discriminate between the two patterns" (p. 284).
- The **number of desensitized states on each side of the gating isomerization** (proximal count vs distal count). "CVF is very sensitive to changes in the number of desensitized states on each side of the gating isomerization" (p. 290). Also the number of states overall, and (with the correct error model) approximately-normal MLEs that make the LR test valid.

**What remains unidentifiable from macroscopic data even with the covariance:**
- **Linear vs branched arrangement at equal state count per side.** "It has been impossible to discriminate between linear and branched arrangements that have the same number of desensitized states on each side of the gating isomerization ... There is no difference in the total L[LH] and for each data set, L[LR] is close to zero" (p. 290). This includes topologies differing in the **number of gateway states** (states directly connected to an open state): "In fitting to either 3D-I or 3D-II, it is impossible to discriminate between branched and linear forms" (p. 283).
- **Which nonconducting state is which (label-swapping of indistinguishable desensitized states).** Data from a symmetric 2D-I model (D1==C==O==D2) fits "the correct model with at least two distinct local maxima depending on which of the two desensitized states is responsible for each of the two phases of desensitization" (p. 281). The likelihood is multimodal because the two nonconducting states are exchangeable in the observable.

**How macroscopic vs single-channel observation changes the picture.** The number of gateway states "dramatically affects the correlations observed when fitting adjacent events derived from single-channel data (Ball et al., 1989; Colquhoun and Hawkes, 1987)" (p. 283), yet it is invisible to the macroscopic covariance. So single-channel adjacent-interval correlations carry a topological discriminant (gateway-state count, linear vs branched) that the macroscopic ensemble covariance does not. The trade the other way: a single macroscopic trace of several hundred channels carries the statistical power that single-channel analysis would need "several hundred" repetitions to match, and CVF "does not require correction for missed events" (p. 291). Macroscopic covariance thus sits strictly between mean-only fitting (below it) and full single-channel dwell-time analysis (above it, for topology).

## Relevance to the theory-model-observable picture (A mechanistic vs B aggregation)

The paper is almost entirely about source **(B), aggregation (rate matrix -> observable)**, and gives a clean, quantitative instance of it.

- The observable current sees the state space only through the diagonal conductance matrix Γ, which is **zero for every closed and desensitized state**. All nonconducting states are aggregated into "invisible," and their kinetic footprint reaches the data only through their effect on open-state occupancy. Distinct Q matrices that yield the same mean (Eq. 3) and the same covariance (Eq. 4) are macroscopically indistinguishable. Empirically, linear and branched topologies with the same per-side state counts give identical total log-likelihood, and symmetric-desensitized-state models give a multimodal likelihood. This is the classic aggregated-Markov equivalence, here for the many-channel/interval observable rather than single-channel dwell times: the covariance narrows the equivalence class (it resolves state number, per-side counts, pattern) but does not collapse it to one topology.
- The label-swapping multimodality (p. 281) is structurally the same phenomenon as the left-right coupling bimodal posterior in the P2X2 work (Moffatt & Pierdominici-Sottile, Comm Biol 2025): a macroscopic symmetry that leaves an equivalence class the data cannot break, resolvable only by importing an independent modality. Celentano-Hawkes even flags the analogous escape hatch: multi-perturbation protocols (deactivation, recovery) "may prove more powerful in discriminating between such models" (p. 290).

Source **(A), mechanistic (theory -> rate matrix)**, is largely out of scope. The paper works at the rate-matrix level and does not discuss molecular stories mapping to a matrix. It touches (A) only in that the state labels are treated as abstract: "The distinction between the C and D states is otherwise arbitrary" (p. 277). A rate matrix encodes topology and rates, not which physical conformation each node is; the observable adds no information about node identity beyond conducting-vs-not, so the mechanistic labeling ambiguity is untouched by, and compounds on top of, the aggregation ambiguity this paper measures.

Bottom line for the program: adding the covariance to the mean is a real identifiability gain (topology arrangement, per-side counts, phase-to-state pattern), but it does not make macroscopic data topologically complete. Linear-vs-branched and gateway-state count stay in the equivalence class, and exchangeable nonconducting states stay multimodal, which is exactly where an external modality has to enter.

## Verbatim (sourced quotes)

- (Abstract, p. 276) "Unlike conventional sum-of-squares minimization, CVF fits both the magnitude of the recorded current and the strength of the correlations between different time points."
- (Abstract, p. 276) "Coupled with the likelihood ratio test, it accurately discriminates between models with different numbers of states, discriminates between most models with the same number but a different arrangement of states, and extracts meaningful information on the relationship between the desensitized states and the phases of macroscopic desensitization."
- (p. 281) "Because the two desensitized states are indistinguishable, data generated using model 2D-I (D1==C==O==D2) fit to the correct model with at least two distinct local maxima depending on which of the two desensitized states is responsible for each of the two phases of desensitization."
- (p. 283) "In fitting to either 3D-I or 3D-II, it is impossible to discriminate between branched and linear forms. The value of the L[LR]s for these comparisons are rarely >0.001 or below 0. This is surprising given that, in the case of 3D-I, the comparison involves models with differing numbers of gateway states (states connected directly to an open state). The number of gateway states dramatically affects the correlations observed when fitting adjacent events derived from single-channel data (Ball et al., 1989; Colquhoun and Hawkes, 1987)."
- (p. 290) "Comparing models with the same number of states and free parameters is made possible by the use of the full covariance matrix. Fitting by wSS, which only takes into account the magnitude of the current at each point in time, fails to discriminate between differences in the arrangement of states or in the pattern of desensitization."
- (p. 290) "CVF is very sensitive to changes in the number of desensitized states on each side of the gating isomerization. However, it has been impossible to discriminate between linear and branched arrangements that have the same number of desensitized states on each side of the gating isomerization ... There is no difference in the total L[LH] and for each data set, L[LR] is close to zero."
- (Appendix I, p. 292) "Each kinetic model generates a unique Q matrix. A change in the number of states changes the dimensions of the Q matrix, and a change in the arrangement of states changes the arrangement of the elements within the Q matrix."

## Caveats

- **Empirical, not a theorem.** The linear-vs-branched and gateway-state non-identifiability is a demonstrated failure to discriminate (L[LR]≈0, identical total log-likelihood) on simulated **single-perturbation** data, not a proved manifest-equivalence result. The authors explicitly say "It remains to be seen if different application paradigms with multiple perturbations ... prove more powerful in discriminating between such models" (p. 290). Treat it as strong evidence of a macroscopic blind spot under one protocol, not as a universal impossibility.
- **Observable differs from MacroIR.** Here the observable is point-sampled multichannel current with pairwise-time covariance (100-200 selected points per trace); MacroIR uses interval-integrated currents. Same relaxation+fluctuation family, but the exact sufficient statistics differ, so the identifiability boundary need not transfer verbatim.
- **Scope.** Two-state-open-and-desensitized gating models fit to GABA_A traces; no explicit agonist-binding step in most models; single conductance class. Conclusions about "arrangement" concern desensitized-state topology around a gating isomerization, not arbitrary Markov graphs.
- **(A) mechanistic non-identifiability is not addressed** beyond the remark that state labels are arbitrary; do not cite this paper for theory-to-matrix multiplicity.

## Sources (URLs)

- Repo PDF: `docs/bibliography/CelentanoHawkes_2004_Covariance_Matrix_Kinetic_Fitting_BiophysJ.pdf` (pp. 276-294 read directly).
- ScienceDirect record: https://www.sciencedirect.com/science/article/pii/S000634950473516X
- PubMed Central full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC1304350/
- DOI: https://doi.org/10.1529/biophysj.103.036632