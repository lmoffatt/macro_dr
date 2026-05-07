# Literature review: differentiable binary-doubling construction of the symmetric-subspace propagator with boundary-conditioned moments

**Date drafted:** 2026-05-01
**Status:** working memo, supersedes the three exploratory reports (`Novel CTMC Ensemble Inference Algorithm.md` [Gemini], `compass_artifact_…md` [Claude.ai], `deep-research-report.md` [GPT]). Quotes below are verified directly from the PDFs in this folder.

---

## Bottom line up front

The construction is best positioned as **algorithmic synthesis**, not foundational invention. Each ingredient has prior art; the integration into a single forward-mode-AD-stable pipeline targeting exact ion-channel micro-likelihood at experimental N is novel.

The strongest novelty claim defensible from the verified readings:

> The first algorithm to compute the exact boundary-conditioned moment bundle (transition probability *m*, conditional mean *μ*, conditional variance *σ²* of an integrated observable) for an ensemble of N independent identical CTMCs by binary doubling on N — exploiting the Kronecker-sum factorization of the symmetric-subspace generator established in the composition-Markov-chain framework of **Zhou & Lange (2009)** — in a form that is trivially differentiable under forward-mode automatic differentiation. This makes derivative-based exact-microscopic Bayesian inference (gradient MLE, Fisher information, HMC) feasible at the experimental ion-channel N where Gaussian/Kalman moment closure (the standard route since **Moffatt 2007** and most recently **Münch et al. 2022**) is the documented practical bottleneck.

---

## 1. The closest mathematical home: Zhou & Lange (2009)

**Citation:** Hua Zhou & Kenneth Lange, "Composition Markov chains of multinomial type," *Adv. Appl. Probab.* 41(1):270–291, 2009. DOI 10.1239/aap/1240319585. PDF: `Zhou_Lange_2009_Composition_Markov_Chains_AdvApplProb.pdf`.

### What they prove

Given n identical particles each evolving under the same single-particle Markov chain on d states, the *composition chain* recording counts (n₁, …, n_d) on the multinomial state space C(n+d−1, d−1):

- **Multinomial equilibrium** (Proposition 1): the composition chain has multinomial equilibrium with cell probabilities determined by the single-particle equilibrium π.
- **Eigenstructure transfer** (Propositions 1–2, Facts 1–9): all eigenvalues of the composition chain are obtained by projecting eigenvectors of the product chain to the symmetric subspace; column eigenvectors are the **multivariate Krawtchouk polynomials of Griffiths**.
- **Continuous-time carry-over** (Section 5, page 280): "Propositions 1 and 2 continue to hold for continuous-time Markov chains provided that the n identical particles move independently." For the continuous-time generator, the composition-chain eigenvalue is β_ñ = Σ nᵢ γᵢ where γᵢ are the eigenvalues of the single-particle generator Q.
- **Kronecker structure** is central and explicit. Page 274: the product kernel for s = n is M ⊗ ⋯ ⊗ M; the composition chain inherits Kronecker structure through the lumping map.

### What they do *not* do

- **No endpoint-conditioned conditional moments** of integrated observables — the paper is purely about the chain itself (equilibrium, convergence rates via χ² bounds, spectral structure). No (μ, σ²) bundle for any path-integrated quantity.
- **No algorithm for finite-time inference.** The continuous-time finite-time transition matrix is given as a spectral expansion (equation 3): exp(tΩ*) = Σ exp(γᵢt) wⁱ (vⁱ)ᵀ diag(π). Closed-form via single-particle spectral decomposition — *exactly the spectral path the user's algorithm avoids* due to derivative pathology near degenerate eigenvalues.
- **No autodiff angle, no Bayesian inference, no HMM filtering.** Pure mathematical analysis: ergodicity, reversibility, χ² convergence rates.
- **No binary doubling on N.** They treat the composition chain at fixed N via lumping + spectral decomposition.

### Verdict

The named mathematical home for the count-space construction. **The user's algorithm operates on exactly the object Zhou-Lange characterize**, but provides what they don't: a forward-arithmetic recursive *algorithm* for the propagator with derivatives, a *bundle* of boundary-conditioned moments, and a path to *inference*. Cite as the foundational mathematical reference; do not claim novelty for the count-space framework itself.

---

## 2. The closest competitor on the differentiable-propagator axis: Rupp et al. (2024)

**Citation:** Kevin Rupp, Rudolf Schill, Jonas Süskind, Peter Georg, Maren Klever, Andreas Lösch, Lars Grasedyck, Tilo Wettig, Rainer Spang, "Differentiated uniformization: a new method for inferring Markov chains on combinatorial state spaces including stochastic epidemic models," *Computational Statistics* 39:3643–3663, 2024. DOI 10.1007/s00180-024-01454-9. arXiv 2112.10971. PDF: `Rupp_2024_Differentiated_Uniformization_ComputStat_published.pdf`.

### What they do

Extend Grassmann's (1977) uniformization with derivative tracking. Algorithm 2 ("Differentiated Uniformization") computes both p(t) = exp(tQ)p(0) **and** p′(t) = ∂p(t)/∂θ in a single Poisson-weighted series:

- Uniformization: P := Q/γ + I with γ ≥ max_x |Q_xx|. Then exp(tQ) = e^{−γt} Σ (γt)ⁿ/n! · Pⁿ.
- Derivative recursion (equation 10): ∂Pⁿ/∂θ = (∂P/∂θ)Pⁿ⁻¹ + P(∂Pⁿ⁻¹/∂θ).
- All operations are matrix–vector products and additions on non-negative objects. **No eigendecomposition**, no spectral pathology.
- Operators P, ∂P/∂θ, γ, ∂γ/∂θ only need to support matrix–vector products; the M×M matrix never has to be stored if Q has a tensor representation.

### What they apply it to

Stochastic SIR for the COVID-19 first wave in Austria. State space (S, I) ∈ {0,…,N}² with N = 9 million, so |X| = 81 trillion. Q is represented as a sum of 4 tensor products of (N+1)×(N+1) banded matrices (equation 21), reducing storage from O(N⁴) to O(N) and matrix-vector cost from O(N⁴) to O(N²). HMC for full Bayesian inference of infection rate β and recovery rate α.

### Direct quotes establishing scope

Page 3645:

> "In this paper, we provide an algorithm that directly computes exp(tQ) and, crucially, ∂exp(tQ)/∂θ at the same time. For the SIR model it scales cubically in the population size but is practical."

Page 3645:

> "However, Ho et al. (2018) have recently provided an algorithm that solves the Kolmogorov equation in the Laplace domain… for which their algorithm scales quadratically in the population size."

### What Rupp does *not* do

- **Output is m_ij and ∂m_ij/∂θ only.** No conditional mean μ_ij, no conditional variance σ²_ij of any integrated observable. The full bundle the user's algorithm carries is absent.
- **Does not exploit identical-particle / multinomial-symmetric-subspace structure.** Their tensor representation of Q for SIR is for *distinguishable species* (S vs I are different roles, not exchangeable particles). The state space is a Cartesian product (S, I) ∈ ℕ², not a multinomial composition.
- **Cost on a fixed-N evaluation is O(γN²)**, with γ ~ N². Total ~O(N⁴) for SIR. Same order as the user's binary-doubling final stage.

### Verdict

The closest existing prior on the "differentiable propagator without eigendecomposition" axis. Same general problem class (HMC on combinatorial-state-space CTMCs), genuinely complementary mechanism (uniformization series vs. binary-doubling convolution), strictly smaller output (m only). **The user's algorithm distinguishes itself by (a) exploiting identical-particle multinomial-subspace structure, (b) carrying the full conditional moment bundle.** Both are substantive distinguishing features, neither is gilding.

The right framing in the paper: complementary methodologies serving different scope. Rupp handles arbitrary tensor-sum CTMCs (including coupled species); the user's method specializes to non-interacting ensembles where the additional symmetry permits a richer bundle and a different algorithmic route.

---

## 3. The closest ion-channel conceptual neighbor: Münch et al. (2022)

**Citation:** Jan L Münch, Fabian Paul, Ralf Schmauder, Klaus Benndorf, "Bayesian inference of kinetic schemes for ion channels by Kalman filtering," *eLife* 11:e62714, 2022. DOI 10.7554/eLife.62714. PDF: `Munch_2022_Bayesian_Kalman_IonChannels_eLife.pdf`.

### What they do

Bayesian (HMC) inference of ion-channel kinetic schemes from macroscopic patch-clamp + cPCF fluorescence data. Hidden state is the count vector **n**(t) = (n₁(t), …, n_K(t))ᵀ, the same object the user's algorithm operates on. Dynamics modeled as a multivariate-Gaussian Markov process (equation 3):

- n_{t+1} = T n_t + ω_t, ω ~ N(0, Q_t)
- T = exp(KΔt), Q_t = Q(T, n_t) (state-dependent covariance per Ball 2016)
- Observations y_t = H n_t + ν_t, multivariate Gaussian.

Apply Kalman filtering on (mean, covariance) of **n**. Generalizes the user's Moffatt 2007 paper by adding state-dependent fluctuations (open-channel noise) and Poisson fluorescence noise, and switching from ML to Bayesian.

### Direct quotes — these are the paper's load-bearing gap statements

Page 4:

> **"On the one hand, a complete HMM analysis (forward algorithm) would deliver the most exact likelihood of macroscopic data. On the other hand, the computational complexity of the forward algorithm limits this type of analysis in ensemble patches to no more than a few hundred channels per time trace (Jahnke and Huisinga, 2007). To tame this computational complexity, we approximate the solution of the CME with a Kalman filter (KF), thereby remaining in a stochastic framework Kalman, 1960."**

Page 4:

> "due to the state-dependent noise the central Bayesian update equation loses its analytical solution. We derived an approximation which is correct for the first two moments of the probability distributions."

Page 2:

> "Previous rigorous attempts to incorporate the autocorrelation of the intrinsic noise of current data into the estimation (Celentano and Hawkes, 2004) suffer from cubic computational complexity (Stepanyuk et al., 2011), rendering the algorithm non-optimal or even impractical for a Bayesian analysis of larger data sets."

Page 4:

> "Our approach generalizes the work of **Moffatt, 2007** by including state-dependent fluctuations such as open-channel noise and Poisson noise in additional fluorescence data."

### What Münch does *not* do

- The exact HMM forward algorithm at experimental N. They acknowledge it as the gold standard but cite it as computationally infeasible past "a few hundred channels per time trace."
- Carry the (m, μ, σ²) bundle — they propagate (mean, covariance) only, with the explicit "correct for the first two moments" approximation caveat.
- Differentiable construction without spectral decomposition. T = exp(KΔt) is computed conventionally; derivatives go through standard expm machinery.

### Verdict

This is the gap statement the user's paper directly fills. Münch et al. (2022):

1. Explicitly identify "complete HMM forward algorithm" as the exact-but-infeasible alternative to their KF.
2. Cite Jahnke & Huisinga (2007) for the "few hundred channels" computational ceiling.
3. Generalize Moffatt (2007) — the user's own paper — into the Bayesian Kalman framework.
4. Acknowledge their core approximation is "correct for the first two moments" (i.e., it is moment closure, not exact).

**The natural opening for the follow-up paper writes itself**: cite Münch's gap statement directly, then state that the binary-doubling construction makes the exact HMM forward algorithm tractable past the "few hundred channels" ceiling Münch identifies, by specializing to the non-interacting-channel regime that already underlies their Kalman setup.

This is the most defensible, highest-impact framing identified in the entire literature search.

---

## 3a. The immediate predecessor: Moffatt & Pierdominici-Sottile (2025) — MacroIR

**Citation:** L. Moffatt and G. Pierdominici-Sottile, "Bayesian inference of functional asymmetry in a ligand-gated ion channel," *Communications Biology* 2025 (DOI 10.1038/s42003-025-09056-x — supplementary material in `42003_2025_9056_MOESM2_ESM.pdf`, peer review file in `42003_2025_9056_MOESM1_ESM.pdf`). **Author of this codebase.**

This is the *immediate published predecessor* of the binary-doubling Qdt construction. It is also the paper whose limitations the new algorithm is designed to address, by the same author. The follow-up paper writes itself as a methodological extension of this one.

### What MacroIR (the published method) does

For a Markov model with N_states, transition rate matrix Q, state-wise current vector γ:

- **Single-channel boundary-conditioned conditional mean** (eq. S1 of the SI):

  γ̄_{i→j}(t) = (1/P_{i→j}(t)) Σ_{k, n_1, n_2} V_{i n_1} V⁻¹_{n_1 k} γ_k V_{k n_2} V⁻¹_{n_2 j} E_2(λ_{n_1} t, λ_{n_2} t)

  where V is the eigenvector matrix of Q, λ_n are its eigenvalues, and E_2(x, y) = (eˣ − eʸ)/(x − y) for x ≠ y, E_2(x, x) = eˣ.

  **This is the eigendecomposition route to μ_ij at single-CTMC level — the same construction Hobolth–Jensen Section 2 / Tataru–Hobolth give for endpoint-conditioned mean time/transition counts, applied here to integrated current.**

- **Macroscopic prediction (eqs. S2–S3):**

  ȳ^pred_{0→t} = N_ch · (μ_0^prior · γ̄_0)
  σ²_{ȳ_{0→t}^pred} = ε²/t + ν² + N_ch [γ̄_0ᵀ (Σ_0^prior − diag(μ_0^prior)) γ̄_0 + μ_0^prior · 𝐄]

  Macro state is (μ^prior, Σ^prior) — mean and covariance of the count vector. **This is the Gaussian moment closure** on the multinomial subspace.

- **Update (eqs. S4–S6):** standard Kalman correction with a Kalman-gain-like term.

### What MacroIR (2025) **does not** do — i.e., what the new algorithm changes

1. **Eigendecomposition is required.** Equation S1 explicitly uses V, V⁻¹, and the eigenvalues λ_n. This inherits the spectral pathology near degenerate eigenvalues (Schranz et al. 2008) and is hostile to forward-mode AD because of the V⁻¹ term and the Sylvester-equation structure of ∂V/∂θ.

2. **Only the single-channel μ_{i→j}, not σ²_{i→j}.** MacroIR computes only the conditional mean of integrated current at single-channel level. The conditional variance is collapsed into the macroscopic ensemble Gaussian variance via moment closure; it is not propagated as a per-(i, j)-pair single-channel object.

3. **Gaussian moment closure on the count vector.** The lift from single-channel to ensemble uses (μ, Σ) of the count vector with Gaussian update — i.e., it is a Kalman filter, not an exact micro-state filter. The "correct for the first two moments" approximation flagged explicitly by Münch et al. (2022, page 4) is the same approximation MacroIR makes, by construction.

### How the new algorithm extends MacroIR

The new binary-doubling Qdt construction extends MacroIR (2025) along three orthogonal axes:

1. **Replaces the eigendecomposition route to γ̄_{i→j}** (eq. S1) with the augmented-block-triangular Van Loan / Carbonell route (see §5.3 below), eliminating V, V⁻¹, λ_n, E_2 and the spectral derivative pathology entirely. The single-channel primitives are computed by the same arithmetic operations used in the rest of the pipeline.

2. **Adds the σ²_{i→j} bundle at single-channel level**, via the Carbonell augmented-matrix-exponential block-triangular trick (eq. 13 of Carbonell et al. 2008). This is the missing piece that lets variance propagate exactly through the multinomial convolution under the variance-additivity-given-split-endpoints observation.

3. **Replaces Gaussian moment closure with exact independence-factorized binary doubling** on the multinomial symmetric subspace. The macro-state (μ, Σ) of the count vector is replaced by the full conditional moment bundle on the C(N+S−1, S−1)-dimensional symmetric subspace, propagated by linear-additive convolution rather than Kalman update.

### Implications for the follow-up paper

The natural opening sentence of the follow-up paper is now:

> "Our previous MacroIR method (Moffatt & Pierdominici-Sottile 2025) computed the single-channel boundary-conditioned conditional mean via eigendecomposition and lifted to the macroscopic ensemble via Gaussian moment closure on the count vector — the same approximation strategy independently identified by Münch et al. (2022) as bounded to 'no more than a few hundred channels per time trace' for exact analysis. Here we extend MacroIR to the full conditional moment bundle (m, μ, σ²), eliminate the eigendecomposition by using the augmented-matrix-exponential representation of the integrated-observable moments at single-channel level (Carbonell et al. 2008; Hobolth & Jensen 2011 route 3), and replace the Gaussian moment closure with exact independence-factorized binary doubling on the multinomial symmetric subspace (Zhou & Lange 2009). The resulting construction is forward-mode-AD-stable end to end and makes derivative-based exact-microscopic Bayesian inference tractable at experimental N for the first time."

This is a much cleaner and more defensible novelty story than positioning against Rupp 2024 or Münch 2022 alone. The lineage:

- **Moffatt 2007** (Biophys J) — KF for ML estimation from macroscopic currents.
- **Münch et al. 2022** (eLife) — Bayesian/HMC generalization, state-dependent fluctuations, fluorescence; explicitly cites the exact-HMM ceiling.
- **Moffatt & Pierdominici-Sottile 2025** (Comms Biol) — MacroIR, recursive Bayesian likelihood for time-averaged currents, conformational models for P2X2; eigendecomposition + Kalman.
- **New algorithm (this work)** — exact micro likelihood, full (m, μ, σ²) bundle, eigen-free, binary-doubling on N.

Each step extends the previous in exactly one specifiable direction, and the whole arc is by the same author or close collaborator group.

---

## 4. Verified novelty claim, calibrated against all three readings

### What is genuinely novel (defensible from this search)

1. **The integration**: independence-factorized binary-doubling-on-N construction of the symmetric-subspace propagator, packaged as the full (m, μ, σ²) bundle, with forward-mode AD as a built-in property — applied to ion-channel HMM inference.
2. **The variance-tracking choice**: explicitly motivating the propagation of g_var (centered second moment) rather than g_total (raw second moment) by the variance-additivity-under-conditional-independence-given-split-endpoints observation. None of Zhou-Lange, Rupp, or Münch carry conditional variance of an integrated observable.
3. **The "few hundred channels" ceiling**: Münch's documented practical bottleneck for exact HMM inference of patches is the headline use case the algorithm relieves.

### What is not novel (folklore)

- Kronecker-sum identity exp(Q⊕Q) = P⊗P and the symmetric-subspace factorization of non-interacting many-body propagators. Standard physics (Negele–Orland, Doi–Peliti) and explicit in Zhou-Lange Section 5.
- The composition-Markov-chain / multinomial-state framework as such. Established in Zhou-Lange 2009 and earlier (Darvey & Staff 1966; Karlin–McGregor; Griffiths Krawtchouk).
- Multinomial-Poisson convolution structure for monomolecular CMEs. Jahnke & Huisinga 2007.
- Endpoint-conditioned (μ, σ²) of integrated functionals at single-CTMC level. Hobolth–Jensen 2011, Tataru–Hobolth 2011, Minin–Suchard 2008, Carbonell–Jiménez–Biscay 2008.
- Differentiable matrix-exponential construction without eigendecomposition for HMC on combinatorial state spaces. Rupp et al. 2024 (uniformization route); Al-Mohy & Higham 2009 (Fréchet machinery).
- Kalman-filter / Gaussian moment-closure approach to N-channel macroscopic inference. Moffatt 2007, Münch 2022, the entire RE/KF lineage.

### Honest framing recommendation

> "We bring together four established lines: (i) the composition-Markov-chain framework for multinomial-state lumping of N independent identical particles (Zhou & Lange 2009, Griffiths Krawtchouk, Jahnke & Huisinga 2007); (ii) endpoint-conditioned moment computation for integrated observables at the single-CTMC level (Hobolth & Jensen 2011, Minin & Suchard 2008); (iii) differentiable propagator construction that bypasses eigendecomposition (Rupp et al. 2024 differentiated uniformization, Al-Mohy & Higham 2009); and (iv) ion-channel ensemble inference where Gaussian/Kalman moment closure is the documented practical alternative (Moffatt 2007, Milescu 2005, Münch et al. 2022). Our contribution is the integration: a forward-mode-AD-stable construction of the full (m, μ, σ²) bundle by binary doubling on N, exploiting the Kronecker-sum factorization on the multinomial symmetric subspace. The construction extends the exact HMM forward algorithm, which Münch et al. (2022) identified as currently restricted to 'no more than a few hundred channels per time trace' (citing Jahnke & Huisinga 2007), into the experimental-N regime relevant for typical patch-clamp recordings."

---

## 5. Single-CTMC primitives the construction lifts to the population level

The user's algorithm gets the single-CTMC building blocks (m_ij, μ_ij, σ²_ij) from established methods and lifts them to the N-particle multinomial subspace via independence-factorized binary doubling. The three references below provide those primitives. All three operate on a *single* chain — endpoint conditioning is x(0) = a, x(T) = b, not a population/symmetric-subspace conditioning.

### 5.1 Hobolth & Jensen (2011): three routes for endpoint-conditioned summary statistics

**Citation:** A. Hobolth and J. L. Jensen, "Summary statistics for endpoint-conditioned continuous-time Markov chains," *J. Appl. Probab.* 48:911–924, 2011. doi:10.1239/jap/1324046009. PDF: `hobolth2011.pdf`.

Surveys and extends three computational routes for the endpoint-conditioned moments needed in CTMC EM:

- E[D_α 1(x(T)=b) | x(0)=a] — mean time in state α (eq. 2a)
- E[N_αβ 1(x(T)=b) | x(0)=a] — mean number of α→β transitions (eq. 2b)
- E[D_α D_β 1(x(T)=b) | x(0)=a], E[N_αβ N_γδ 1(x(T)=b) | x(0)=a], E[N_αβ D_γ 1(x(T)=b) | x(0)=a] — second-order moments (eqs. 2c–2e)

The integral representations (eqs. 3–4):

  I(a,b,α,β) = ∫₀^T P_aα(t) P_βb(T−t) dt
  I(a,b,α,β,γ,δ) = ∫₀^T ∫₀^t P_aα(u) P_βγ(t−u) P_δb(T−t) du dt

are the single-CTMC objects the user's algorithm lifts.

Three computational routes:
1. Eigenvalue decomposition (Section 2) — assumes Λ is diagonalizable.
2. **Uniformization** (Section 3, Theorem 1) — their main contribution: a unified framework for general statistics H = ψ(z₀) + Σ φ(z_{i−1}, z_i) f(T_{i+1}−T_i), expanded as a Poisson-weighted series in the discrete-time chain R = Q/μ + I.
3. Van Loan integrals (Section 4) — augmented block-triangular matrix exponential trick.

Cites Ball & Milne (2005) for ion-channel applications. **Single chain only, no population structure.** This is the most directly cited reference for the user's μ_ij and σ²_ij at single-CTMC level.

### 5.2 Minin & Suchard (2008): the augmented-generator generating-function trick

**Citation:** V. N. Minin and M. A. Suchard, "Counting labeled transitions in continuous-time Markov models of evolution," *J. Math. Biol.* 56:391–412, 2008. doi:10.1007/s00285-007-0120-8. PDF: `minin2007.pdf`.

For a single CTMC X_t and a labeled set R of transitions, the counting process N_t = number of R-transitions in (0, t]. The matrix probability generating function is (eq. 30):

  G(r, t) = exp((Λ_R̄ + r Λ_R) t)

where Λ_R is Λ restricted to labeled transitions and Λ_R̄ is everything else. Restricted factorial moments fall out by differentiation at r = 1 (eq. 31): M^[k](t) = ∂^k G(r,t)/∂r^k |_{r=1}, with the recursive integral formula (eq. 33):

  M^[k](t) = k ∫₀^t M^[k−1](θ) Λ_R exp(Λ(t−θ)) dθ.

This is the augmented-block-triangular trick (same structure as Van Loan / Carbonell). Page 12 explicitly compares to Hobolth-Jensen 2005 / Narayana-Neuts uniformization — Hobolth-Jensen 2011 (above) generalizes the framework.

The construction assumes diagonalizability for the spectral route (page 11 quote: "almost all evolutionary models derive reversible CTMCs. Reversibility implies that the infinitesimal generator is similar to a symmetric matrix and hence is diagonalizable"). Ion-channel CTMCs with absorbing states (e.g. inactivation) are not in this regime, so the spectral route inherits the same eigendecomposition pathology the user's algorithm avoids.

**Single chain only, no population, no autodiff.** Provides the parametric-matrix-exponential generating-function machinery for second-moment objects at single-CTMC level.

### 5.3 Carbonell, Jiménez & Pedroso (2008): block-triangular trick generalized to arbitrary multiplicity

**Citation:** F. Carbonell, J. C. Jiménez, L. M. Pedroso, "Computing multiple integrals involving matrix exponentials," *J. Comput. Appl. Math.* 213:300–305, 2008. doi:10.1016/j.cam.2007.01.007. PDF: `carbonell2008.pdf`.

Extends Van Loan's 1978 block-triangular trick from multiplicity k ≤ 4 to arbitrary k. Theorem 1 gives, for a block triangular A = [(A_lj)]_{l,j=1:n}, that exp(At) recovers all the iterated integrals B_{lj}(t) = ∫⋯∫ exp(A_ll·) A · ⋯ exp(A_jj·) in the upper-triangular blocks of a single matrix exponential.

The application that matters for this work (eq. 13):

  ∫₀^t exp(F(t−s)) G(s) G^T(s) exp(F^T(t−s)) ds = B_{1, 2p+2}(t) B_{11}^T(t)

— the *covariance* matrix of an integrated observable, computed as the upper-right block of an augmented matrix exponential. Cited for the system noise covariance of the extended Kalman filter for continuous-discrete state-space models with additive noise.

**This is the σ²_ij primitive at single-CTMC level.** No algorithm for ensemble, no autodiff, no identical-particle structure. The user's algorithm uses this (or its equivalent) once at N = 1, then composes via binary doubling. The Carbonell route avoids eigendecomposition entirely at the single-CTMC level if the user wants — paired with the binary doubling, the entire pipeline can be eigen-free end to end.

### Direct connection to MacroIR (Moffatt & Pierdominici-Sottile 2025)

MacroIR's equation S1 computes γ̄_{i→j}(t) using V, V⁻¹, λ_n, and E_2(λ_a t, λ_b t) = (e^{λ_a t} − e^{λ_b t})/(λ_a − λ_b). This is the *spectral evaluation* of exactly the kind of integral Carbonell handles: ∫₀^t e^{Q(t−s)} γ e^{Qs} ds and similar. The E_2 function is the closed-form result when Q is diagonalizable; near eigenvalue degeneracy E_2 → te^{λt} and the formula becomes ill-conditioned in the V⁻¹ multiplication.

Carbonell's block-triangular trick replaces V Λ V⁻¹ with the augmented matrix exponential exp(A·t) for a block A built from Q, γ, and identity blocks. The same integral, no eigendecomposition, no E_2 special-case branch. **For the new algorithm, this is the recommended replacement for equation S1 of MacroIR.** Then σ²_{i→j} via eq. 13 of Carbonell, and the entire single-channel (m, μ, σ²) bundle is constructed by two augmented matrix exponentials at single-CTMC level — no eigendecomposition, fully AD-stable.

### 5.4 Jahnke & Huisinga (2007): the algebraic ancestor for marginal distributions

**Citation:** T. Jahnke and W. Huisinga, "Solving the chemical master equation for monomolecular reaction systems analytically," *J. Math. Biol.* 54:1–26, 2007. doi:10.1007/s00285-006-0034-x. PDF: `jahnke2006.pdf`.

For any monomolecular reaction system (only first-order reactions: S_j → S_k conversion, * → S_k inflow, S_j → * degradation) with deterministic initial condition X(0) = ξ, **Theorem 1**:

  P(t, ·) = P(·, λ(t)) ⋆ M(·, ξ₁, p^(1)(t)) ⋆ ⋯ ⋆ M(·, ξ_n, p^(n)(t))

where 𝒫 is the product Poisson distribution (for inflow-generated molecules), 𝓜 is the multinomial distribution, and the parameter vectors λ(t) and p^(k)(t) evolve under the deterministic reaction-rate equation ṗ = A(t)p, λ̇ = A(t)λ + b(t).

**Closed system + multinomial initial → multinomial stays multinomial** (Proposition 1).
**Open system + Poisson initial → Poisson stays Poisson** (Proposition 2).

The proof (page 13) is exactly the independence-factorization argument: split all molecules into n+1 disjoint subsets indexed by species-at-t=0 (plus a 0-th subset for inflow-generated molecules); each subset evolves independently; the total joint is the convolution of independents. **This is the marginal-distribution version of the same independence argument the user's binary doubling exploits at the propagator level.**

**What Jahnke-Huisinga do NOT do:**
- They give P(t, x) — the *marginal* distribution at time t. No endpoint conditioning, no boundary moments of an integrated observable.
- The multinomial parameter ODE is ṗ = A(t)p — a *deterministic* rate equation, not a propagator construction. They never construct a transition propagator m_ij; the multinomial structure lives at the marginal level.
- No HMM filtering, no inference, no autodiff. Pure forward analysis.

**On the "few hundred channels" claim:** Münch et al. 2022 attribute this specific number to Jahnke & Huisinga 2007. The Jahnke-Huisinga paper does *not* contain that exact number — instead, the introduction (page 2) says "a rather small system of only three species with molecule numbers varying between, say, 0 and 99, contains 100³ different states, and 1,000,000 coupled differential equations have to be solved." Same complexity worry, different specific quantification. **The "few hundred channels" framing is Münch's**, consistent with Jahnke-Huisinga's complexity discussion but not a verbatim quote from them.

### What this collection of single-CTMC primitives establishes for the user's algorithm

The four references above provide:
- (m_ij at single CTMC) standard expm — done a thousand times.
- (μ_ij at single CTMC) Hobolth-Jensen route 1/2/3, Minin-Suchard generating-function trick, or Carbonell augmented-matrix-exponential.
- (σ²_ij at single CTMC) Carbonell augmented-block-triangular Van Loan trick, equivalently second derivative of Minin-Suchard's G(r,t).
- (multinomial structure at marginal level) Jahnke-Huisinga 2007 — the algebraic ancestor of the count-space convolution structure, but only for marginals.

The user's contribution is the **lift**: take any of these single-CTMC routes for (m, μ, σ²) at N=1, and compose to (m, μ, σ²) at arbitrary N by binary doubling on the multinomial symmetric subspace, in a forward-arithmetic AD-trivial manner. Each ingredient is established at N=1; the ensemble lift via binary-doubling-with-conditional-moments is what the literature does not contain.

This sharpens the framing: the user's paper is positioned as the *bridge* between a well-developed single-CTMC summary-statistics literature (Hobolth, Jensen, Minin, Suchard, Carbonell) and a well-developed marginal-distribution count-space literature (Jahnke-Huisinga, Vastola, Zhou-Lange) and a well-developed ion-channel ensemble-inference literature (Moffatt 2007, Münch 2022, Milescu 2005). Each literature on its own is missing one of the three pieces; the construction integrates all three.

---

## 6. Status and verification

### All eight principal references verified directly from PDFs in this folder
- **Zhou & Lange 2009** (sections 1–6, complete spectral framework, continuous-time carry-over confirmed page 280, Kronecker structure confirmed page 274).
- **Rupp et al. 2024** (sections 1–3, full algorithm, tensor SIR application; m_ij and ∂m_ij/∂θ only — confirmed no μ or σ²).
- **Münch et al. 2022** (introduction + Kalman derivation; "few hundred channels" gap statement confirmed page 4; "correct for the first two moments" approximation caveat confirmed page 4; cite Jahnke-Huisinga 2007 as the complexity reference, generalize Moffatt 2007).
- **Jahnke & Huisinga 2007** (full multinomial-Poisson convolution Theorem 1 confirmed; marginal-distribution-only scope confirmed; the specific "few hundred channels" quantification is Münch's, not theirs).
- **Hobolth & Jensen 2011** (three routes confirmed; single-CTMC scope confirmed; Theorem 1 uniformization framework documented).
- **Minin & Suchard 2008** (G(r,t) = exp((Λ_R̄ + rΛ_R)t) generating-function trick confirmed; single-CTMC scope confirmed; reversibility/diagonalizability assumption noted).
- **Carbonell, Jiménez & Pedroso 2008** (block-triangular trick generalized to arbitrary multiplicity; eq. 13 covariance integral formula confirmed; single-shot matrix exponential, no recursion). Connection to MacroIR (2025) eq. S1 made explicit: Carbonell's augmented matrix exponential is the eigen-free replacement for V Λ V⁻¹ + E_2 in the user's published method.
- **Moffatt & Pierdominici-Sottile 2025** (`42003_2025_9056_MOESM2_ESM.pdf`, SI fully read): MacroIR algorithm equations S1–S6 verified; eigendecomposition route to single-channel γ̄_{i→j} confirmed in eq. S1; Gaussian moment closure on count vector confirmed in eqs. S2–S5. **The immediate predecessor of the new construction.**

### Other PDFs in folder, not yet read in detail (do not change headline claims)
- Ho, Crawford & Suchard 2018 — listed by Claude.ai report as next-closest competitor for HMC-on-stochastic-SIR. Continued-fraction route, no Kronecker. Worth verifying for completeness.
- Vastola 2021 — Doi-Peliti rederivation of Jahnke-Huisinga. Worth a skim.
- Milescu 2005 — earlier reports flagged "microscopic-recursive likelihood." Worth confirming the user's prior framing matches this paper's actual content.
- Schranz et al. 2008, Tataru & Hobolth 2011 — already cited, in folder, not load-bearing for the headline claim.

---

## 7. Citations

### Verified directly from PDFs in this folder
- Zhou H., Lange K. (2009). Composition Markov chains of multinomial type. *Adv. Appl. Probab.* 41(1):270–291. doi:10.1239/aap/1240319585
- Rupp K., Schill R., Süskind J., Georg P., Klever M., Lösch A., Grasedyck L., Wettig T., Spang R. (2024). Differentiated uniformization: a new method for inferring Markov chains on combinatorial state spaces including stochastic epidemic models. *Comput. Stat.* 39:3643–3663. doi:10.1007/s00180-024-01454-9
- Münch J.L., Paul F., Schmauder R., Benndorf K. (2022). Bayesian inference of kinetic schemes for ion channels by Kalman filtering. *eLife* 11:e62714. doi:10.7554/eLife.62714
- Jahnke T., Huisinga W. (2007). Solving the chemical master equation for monomolecular reaction systems analytically. *J. Math. Biol.* 54:1–26. doi:10.1007/s00285-006-0034-x
- Hobolth A., Jensen J.L. (2011). Summary statistics for endpoint-conditioned continuous-time Markov chains. *J. Appl. Probab.* 48:911–924. doi:10.1239/jap/1324046009
- Minin V.N., Suchard M.A. (2008). Counting labeled transitions in continuous-time Markov models of evolution. *J. Math. Biol.* 56:391–412. doi:10.1007/s00285-007-0120-8
- Carbonell F., Jiménez J.C., Pedroso L.M. (2008). Computing multiple integrals involving matrix exponentials. *J. Comput. Appl. Math.* 213:300–305. doi:10.1016/j.cam.2007.01.007

### The author's own published lineage (in folder)
- Moffatt L. (2007). Estimation of Ion Channel Kinetics from Fluctuations of Macroscopic Currents. *Biophys. J.* 93:74–91. doi:10.1529/biophysj.106.101212 — Kalman filter for ML estimation; foundation of the lineage.
- Moffatt L., Pierdominici-Sottile G. (2025). Bayesian inference of functional asymmetry in a ligand-gated ion channel. *Communications Biology* (DOI 10.1038/s42003-025-09056-x; SI in `42003_2025_9056_MOESM2_ESM.pdf`) — **MacroIR algorithm**: recursive Bayesian likelihood for time-averaged currents, conformational models for P2X2. **The immediate predecessor of the new binary-doubling Qdt construction**, with eigendecomposition + Gaussian moment closure as the two specific limitations the new algorithm addresses.

### Other PDFs in folder, cited but not yet read in detail
- Milescu L.S., Akk G., Sachs F. (2005). Maximum Likelihood Estimation of Ion Channel Kinetics from Macroscopic Currents. *Biophys. J.* 88:2494–2515. doi:10.1529/biophysj.104.053256
- Vastola J.J. (2021). Solving the chemical master equation for monomolecular reaction systems and beyond: a Doi-Peliti path-integral view. *J. Math. Biol.* 83:48. doi:10.1007/s00285-021-01670-7
- Ho L.S.T., Crawford F.W., Suchard M.A. (2018). Direct likelihood-based inference for discretely observed stochastic compartmental models of infectious disease. *Ann. Appl. Stat.* 12(3):1993–2021. doi:10.1214/18-AOAS1141
- Schranz H.W., Yap V.B., Easteal S., Knight R., Huttley G.A. (2008). Pathological rate matrices: from primates to pathogens. *BMC Bioinformatics* 9:550. doi:10.1186/1471-2105-9-550
- Tataru P., Hobolth A. (2011). Comparison of methods for calculating conditional expectations of sufficient statistics for continuous time Markov chains. *BMC Bioinformatics* 12:465. doi:10.1186/1471-2105-12-465
- Chowdhury D. (2013). Modeling stochastic kinetics of molecular machines at multiple levels: from molecules to modules. *Biophys. J.* 104:2331–2341. doi:10.1016/j.bpj.2013.04.042 — peripheral background on stochastic kinetics of molecular machines (motors, polymerases); not directly competing prior art for ion-channel Qdt construction, but useful framing context.

### Cited but not in folder (peripheral, easy institutional access)
- Al-Mohy A.H., Higham N.J. (2009). Computing the Fréchet derivative of the matrix exponential. *SIAM J. Matrix Anal. Appl.* 30:1639–1657. doi:10.1137/080716426
- Van Loan C.F. (1978). Computing integrals involving the matrix exponential. *IEEE Trans. Automat. Control* 23:395–404. doi:10.1109/TAC.1978.1101743
- Ball F., Milne R.K. (2005). Simple derivations of properties of counting processes associated with Markov renewal processes. *J. Appl. Probab.* 42:1031–1043.
- Grassmann W.K. (1977). Transient solutions in Markovian queueing systems. *Comput. Oper. Res.* 4:47–53.
