# Sample vs correlation distortion: definition, cumulant analysis, and numerical closure

macro_IR Gaussian-likelihood diagnostics (figure_7 family). This note derives what the
**sample** and **correlation** distortions measure, why the sample-distortion peak sits at an
intermediate sampling interval Δ and moves the way it does with noise and channel count, and
closes the argument against the simulation data. Claims are tagged by how well they survived
adversarial + numerical verification:

- **[VERIFIED]** independently re-derived and/or confirmed against source / data.
- **[PARTIAL]** the structure holds but a quantitative sub-claim was corrected.
- **[REFUTED]** the original claim was wrong; the corrected statement is given.
- **[OPEN]** not settled by the current data.

---

## 1. What the two distortions are (from the source)

**[VERIFIED — `src/core/likelihood.cpp`]** The code builds three information matrices per cell
and chains them:

```
G_b   = mean( Σ_t Gaussian_Fisher_Information )     # analytic Gaussian FIM the model self-assigns
J_sample = covariance< dlogL >                      # empirical per-sample score covariance
J     = covariance< Σ_t dlogL >                     # empirical FULL score covariance (with cross-sample terms)
```

and defines (with `apply_normalized_congruence(W_A, B) = A^{-1/2} · B · A^{-1/2}`):

```
Gaussian_Sample_Distortion  (gsdm) = G_b^{-1/2} · J_sample · G_b^{-1/2}     # likelihood.cpp:3553
Likelihood_Correlation_Distortion (cdm) = J_sample^{-1/2} · J · J_sample^{-1/2}   # likelihood.cpp:3438
```

So the built-in decomposition is the chain

> **G_b  →[ sample distortion ]→  J_sample  →[ correlation distortion ]→  J**

- **Sample distortion** = does the analytic **Gaussian** Fisher information match the **empirical
  per-sample** score covariance? It is 1 (identity matrix) iff the per-sample emission is Gaussian
  for scoring purposes. It is a **per-sample non-Gaussianity** measure.
- **Correlation distortion** = does the per-sample information match the full information? It is 1
  iff consecutive samples' scores are uncorrelated. It is a **temporal-correlation** measure.

Note: a sibling `Likelihood_Sample_Distortion` (line 3424) whitens `J_sample` by the *numerical*
Fisher instead of `G_b`; it is empty in the gaussian-only run and is **not** the quantity plotted
in figure_7 (which is `Gaussian_Sample_Distortion`).

The per-sample Gaussian FIM has the textbook form **[VERIFIED — likelihood.cpp:3170-3173, 3214-3217,
3719-3722]**, with σ² = `y_var` the moment-matched predictive variance:

```
G_b = (∂μ/∂θ)² / σ²  +  (∂σ²/∂θ)² / (2σ⁴)
```

A consequence that kills one alternative explanation up front: because σ² **is** the model's
per-sample emission variance by construction, the standardized residual has E[u²]=1 exactly, so
there is **no first/second-moment mis-calibration term** in the sample distortion. Whatever
`gsdm` measures is genuinely a **higher-cumulant** effect, not a variance-calibration bump. (This
is also why the variance-calibration diagnostic r̄²_std peaks at *long* Δ while the sample
distortion peaks at *intermediate* Δ — they are different objects.)

---

## 2. The right expansion: cumulants of the per-sample emission

**[VERIFIED — independent re-derivation]** Write the standardized Gaussian score for one parameter
with a = ∂μ/∂θ, b = ∂σ²/∂θ, u = (y−μ)/σ:

```
s = a·u/σ + ½·b·(u² − 1)/σ²
```

`G_b = E_gauss[s²]` uses u∼N(0,1) (E[u³]=0, E[u⁴]=3), giving `G_b = a²/σ² + b²/(2σ⁴)` — matching the
code. `J_sample = E_true[s²]` uses the **true** per-sample law, which (with the mean/variance
matched, E[u]=0, E[u²]=1) has E[u³]=S (standardized skewness) and E[u⁴]=3+K (standardized excess
kurtosis). Using E[(u²−1)²]=2+K and E[u(u²−1)]=S:

> **J_sample − G_b  =  (∂μ/∂θ)(∂σ²/∂θ)·S/σ³  +  (∂σ²/∂θ)²·K/(4σ⁴)**

This is the whole content of "which Taylor analysis": it is an **Edgeworth/cumulant expansion of
the score's second moment around the Gaussian**. The leading correction is the emission's
standardized **skewness S** (a mean×variance cross term) plus its **excess kurtosis K** (a
pure-variance term). Δ, N and noise enter only through S and K.

Dominance depends on how θ enters:
- a **mean-driven** parameter (a ≫ b) → the **skewness** cross-term dominates the distortion ratio
  (∝ S).
- a **variance-driven** parameter (b ≫ a) → the **kurtosis** term dominates, ratio → K/2.

---

## 3. How gSg and noise enter — dilution by the Gaussian pedestal

**[VERIFIED structurally]** The emission is y = μ + δ_gating + δ_noise, with δ_noise the additive
**Gaussian** instrumental noise (`Current_Noise`) and δ_gating the **non-Gaussian** channel-gating
current. Cumulants add across independent components, and a Gaussian contributes nothing above
order 2:

```
σ²   = v_g + v_n           v_g = gating variance (∝ N·gSg),  v_n = noise variance
κ₃(y) = κ₃(gating)          κ₄(y) = κ₄(gating)
```

so the **standardized** cumulants dilute by the gating variance fraction f = v_g/(v_g+v_n):

```
S = S_g · f^{3/2}          K = K_g · f²          f = v_g / (v_g + v_n)
```

The Gaussian noise is a **pure-variance pedestal**: it adds to σ² but not to κ₃, κ₄, so it shrinks
the *standardized* non-Gaussianity. This is where `gSg` and `noise` enter — **only through the
ratio v_g/v_n that sets f**. More noise lowers f (bigger denominator); more channels/interval
raises v_g and raises f.

**[VERIFIED — exact emission variance, `legacy/qmodel.h:3472, 3489, 3742`]** The predictive
variance the model uses is

```
y_var = e + N·gSg
  e    = Current_Noise · fs / number_of_samples          # noise floor; ∝ 1/Δ (n_samples ∝ Δ), N-independent
  gSg  = g' · Σ · g  (+ within-state shot term)           # gSg = g'(P_Cov − diag(P_mean))g + P_mean·gvar
```

with `Σ` the occupancy covariance `P_Cov`. So `v_n = e` and `v_g = N·gSg`. At the **prior** level
`Σ` is the **per-channel** state-indicator covariance (`qmodel.h:858`, "Cov(X) = diag(P_mean) −
P_mean·P_meanᵀ"), N-independent — which is why the explicit `×N` in `N·gSg` gives the total gating
variance of N independent channels.

**[VERIFIED — recursive N-dependence, `qmodel.h:~4218`]** BUT in the recursive (IR) filter the
measurement update shrinks the occupancy covariance by an amount **∝ N**:

```
Σ_new = Σ_pre − (α · N / y_var) · gSᵀ gS
```

so the **steady-state** `Σ` — and therefore `gSg` — **does depend on N** (more channels → more
informative current → tighter occupancy estimate → smaller `Σ`). The Riccati balance (process
noise `Q` from occupancy diffusion vs the ∝N measurement reduction), `Q = (N/y_var)·(gΣ)²`, has two
regimes:

- **gating-dominated** (N·gSg ≫ e, f→1): y_var ≈ N·gSg ⟹ `Σ ≈ Q` (N-independent) ⟹ `v_g = N·gSg ∝ N`.
- **noise-influenced** (e ≳ N·gSg): y_var ≈ e ⟹ `Σ² = Qe/(Ng²)` ⟹ **`Σ ∝ 1/√N`** ⟹ **`v_g = N·gSg ∝ √N`**.

The sample distortion lives precisely in the **noise-influenced** regime (it is what the noise
dilutes), so there `Σ_ss ∝ 1/√N` and `v_g ∝ √N`. This is a **recursive-filter** effect: a
non-recursive filter leaves Σ at the (N-independent) prior and would give `v_g ∝ N`.

---

## 4. The interval peak and its shift

**[PARTIAL]** Combining §2–§3, the sample distortion for one parameter is a product of two factors
with opposite Δ-trends:

```
sampDist − 1  ≈  (bare gating non-Gaussianity, ↓ with N·Δ by CLT) × (f-factor, ↑ with Δ and N, ↓ with noise)
```

- bare non-Gaussianity S_g, K_g decreases with the number of channel transitions per sample (CLT),
- the f-factor increases toward 1 as v_g grows (long Δ, many channels).

Their product peaks at an **intermediate** Δ — the "mixture" regime where the gating current is
comparable to the noise floor yet still non-Gaussian. This matches the data (peak at Δ≈0.1 at
N=10, noise 0.1). At short Δ the gating is drowned by the noise pedestal (f→0, per-sample ≈
Gaussian); at long Δ the CLT has Gaussianized the gating (S_g, K_g → 0).

**Peak shift — direction [VERIFIED], magnitude [PARTIAL]:**
- **More noise → peak moves to LONGER Δ.** Raising v_n lowers f everywhere; since f grows with Δ,
  the product's maximum moves right. Confirmed: centroid peakΔ 0.072 → 0.70 as noise 0.05 → 10
  (N=10), log-log slope **+0.46** (theory +0.5 → good).
- **More N → peak moves to SHORTER Δ.** Confirmed in direction: centroid peakΔ 0.104 → 0.022 as N
  10 → 10000 (noise 0.1), slope **−0.20** — correct sign but only ~40% of the −0.5 magnitude. Part
  of this is grid truncation (Δ floor 0.01), but the numerical closure below shows there is a
  **genuine √N anomaly**, so the clean symmetric law Δ* ∝ √(noise/N) is **not** claimed; only the
  noise axis is quantitatively √.

---

## 5. Numerical closure

At fixed (N, Δ), varying only the noise, the theory predicts the sample distortion depends on noise
**only through f**, with a single amplitude C = C(N,Δ):

```
sampDist(N,Δ,noise) − 1  =  C / (1 + noise/b)^m ,     b ≡ N·gSg (v_g),   m the f-exponent
```

Fitting C, b (and optionally m) per (N,Δ) to the six noise levels {0.05,0.1,0.2,0.5,1,10} for
k_off (independently reproduced by two extractors):

**[VERIFIED] the noise-dilution shape fits.** R² on (sampDist−1) is **0.98–0.9996** across the
focus grid. Sanity cell N=10, Δ=0.1, k_off: data (sampDist−1) = {.05:0.240, .1:0.146, .2:0.070,
.5:0.023, 1:0.007}, fit **C=0.458, b=0.130, R²=0.9996**. So "noise acts through a single saturating
fraction f" is solid.

**[RESOLVED — was mislabelled] b ∝ √N, C ∝ 1/√N — and this is genuine recursive-filter physics.**
The fitted coefficients each scale as **N^{±0.5}**, not N^{±1}:

| Δ | fitted b vs N | fitted C vs N |
|---|---|---|
| 0.10 | b ∝ N^{+0.51} (b/N falls 8×, N:10→1000) | C ∝ N^{−0.52} |
| 0.05 | b ∝ N^{+0.39…0.56} | C ∝ N^{−0.50} |
| 0.02 | b ∝ N^{+0.43} | C ∝ N^{−0.44} |

The √N is **not** a fitting artifact and **not** the "naive b∝N": it is the **recursive Kalman
steady-state** of §3. In the noise-influenced regime where the sample distortion lives, the IR
filter's occupancy covariance is `Σ_ss ∝ 1/√N` (Riccati balance `Q = (N/y_var)·(gΣ)²` with
y_var≈e), so `v_g = N·gSg ∝ √N` ⟹ crossover `b ∝ √N`, and the bare non-Gaussianity amplitude
`C ∝ 1/√N`. Both match. (The user identified the missing N-dependence of `Σ`/`gSg`; my earlier
"code forces b∝N" reading looked only at the per-step per-channel formula and missed that the
recursive measurement update makes the steady-state `Σ` N-dependent.) A **cross-check that would
seal it**: a non-recursive algorithm (no measurement update) should show `b ∝ N` (Σ frozen at the
N-independent prior); only recursive filters give the √N. **[RESOLVED for IR; NR cross-check OPEN]**

**[REFUTED] "kurtosis-driven / ½K_g·f²".** The f-exponent is **not identified** by a noise sweep of
5–6 points (m and b are degenerate; R²≈1 for a range of m). A free-exponent fit on the sanity cell
gives **m ≈ 1.7 ± 0.14**, and m=1.5 fits marginally *better* than m=2. This is exactly what a
**skewness (f^{3/2}) + kurtosis (f²) mixture** produces — consistent with §2's algebra, which keeps
BOTH terms. Moreover k_off was never shown to be variance-driven (at equilibrium p_open≈0.09 both
∂μ/∂k_off and ∂σ²/∂k_off are nonzero), so there is no basis for dropping the skewness term. The
amplitude C ∝ 1/√N independently points to skewness (standardized skewness of a sum of N channels
∝ 1/√N; kurtosis ∝ 1/N).

> **Corrected leading statement:** the sample distortion is a **per-sample higher-cumulant
> (skewness+kurtosis) non-Gaussianity**, diluted by the Gaussian noise fraction f with an effective
> exponent **between 3/2 and 2 (empirically ≈1.7)**, and — for k_off — dominated by the **skewness**
> cross-term (amplitude ∝ 1/√N), not pure kurtosis.

---

## 6. Sample vs correlation: the physical picture (holds)

**[VERIFIED — data]** The two distortions live at opposite ends of the Δ axis with different peak
mobility:

- **Correlation distortion** peaks at **short Δ** (state persistence: e^{−Δ·k}, samples correlated).
  Its centroid is nearly **immobile** in Δ across noise/N (0.047→0.122 over noise 0.05→10, ×2.6);
  noise/N only scale its magnitude. Signature of a kinetics/persistence effect (set by Δ·k).
- **Sample distortion** peaks at **intermediate Δ** and is **mobile**: peak → longer Δ with noise
  (~√), → shorter Δ with N. Signature of the gating/noise mixture (set by v_g/v_n).

In figure_7 this reads directly: the correlation band is ~**vertical** across the noise/N facets;
the sample band **tilts**. And the earlier r_std diagnostics (r̄, r̄², τ_int) are subsumed — the
sample/correlation pair is the interpretable decomposition of the same information distortion.

One correctness point this also settles: at short Δ the residuals are white and per-sample
calibration is fine, yet the correlation and information (F≠J) distortion are large. That is **not**
a model failure — it is the (correct) information-redundancy cost of oversampling; the macro_IR
likelihood is well specified there.

---

## 7. Status summary

| Claim | Status |
|---|---|
| Definitions gsdm = G_b^{−1/2} J_sample G_b^{−1/2}, cdm = J_sample^{−1/2} J J_sample^{−1/2}; chain G_b→J_sample→J | **VERIFIED** (code) |
| G_b = (∂μ)²/σ² + (∂σ²)²/(2σ⁴); σ² moment-matched ⇒ no variance-miscalibration term | **VERIFIED** (code) |
| J_sample − G_b = (∂μ)(∂σ²)S/σ³ + (∂σ²)²K/(4σ⁴) (cumulant expansion) | **VERIFIED** (re-derived) |
| Noise dilutes standardized non-Gaussianity via f = v_g/(v_g+v_n); S=S_g f^{3/2}, K=K_g f² | **VERIFIED** (structure) |
| sampDist−1 = C/(1+noise/b)^m fits (noise enters only through f) | **VERIFIED** (R² 0.98–0.9996) |
| Dominant cumulant = kurtosis; exponent = 2 | **REFUTED** → skewness+kurtosis mixture; f-exponent m unidentified by noise sweep (1↔2 all fit R²≈0.99); C∝1/√N points to skewness for k_off |
| b ∝ N (gSg const), C ∝ 1/N | **REFUTED** → b ∝ √N, C ∝ 1/√N |
| The √N is a fitting artifact | **REFUTED** → genuine; recursive Kalman steady-state Σ_ss∝1/√N in the noise-influenced regime |
| v_g = N·gSg with gSg N-independent | **REFUTED** → gSg N-dependent via the recursive covariance update (Σ_new = Σ_pre − (αN/y_var)gSᵀgS) |
| Peak → longer Δ with noise | **VERIFIED** (exponent +0.46 ≈ +0.5) |
| Peak → shorter Δ with N | **VERIFIED in direction**; magnitude weak (−0.20) |
| Clean symmetric Δ* ∝ √(noise/N) | **PARTIAL** — √ on noise, weaker on N |
| Correlation immobile (Δ·k), sample mobile (v_g/v_n) | **VERIFIED** (data) |
| Origin of the √N | **RESOLVED** — IR steady-state occupancy covariance (Riccati), not a bug |

### Open items
1. **NR cross-check.** The √N is predicted to be recursive-specific: a non-recursive algorithm
   should show `b ∝ N` (Σ frozen at the N-independent prior). Confirm against NR/NMR sample
   distortion if a noise sweep exists at enough channels.
2. **Fix the f-exponent m (skewness f^{3/2} vs kurtosis f²).** The per-(N,Δ) noise sweep is
   degenerate (m trades with the crossover; free-b→m=2, b∝N→m≈1, all R²≈0.99). Break it either by
   (a) a global (N,Δ,noise) fit that imposes the Riccati `Σ_ss(N)` form and the CLT scalings
   S_g∝1/√(N·n_trans), K_g∝1/(N·n_trans); or (b) instrumenting the C++ to emit the two terms of
   `J_sample − G_b` (skewness cross vs kurtosis) separately, or the raw per-sample emission
   skewness/kurtosis. Neither is in the current battery CSVs.

*Verification: independent algebra re-derivation, source read (likelihood.cpp), two independent
numerical extractors, and an adversarial critic. The definitions and the cumulant identity are
clean; the corrections above are entirely in the interpretation of the coefficients.*
