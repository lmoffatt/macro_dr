# The score as a martingale difference: why the correlation distortion is exactly the failure of the one-step predictive

> Written 2026-07-14. Short derivation, kept at hand. It is the theoretical content behind Figure 4 and it explains the whole five-algorithm ranking from one exact identity.
> Status of novelty: the martingale-difference property of a correctly specified prequential score is textbook (prediction-error decomposition; see the prior-art map, Part II, objection O5, which already recommends this framing). **The exact identity in §3, and the reading of the two algorithm axes it produces (§5), are what we add.** A literature check is running; do not claim novelty until it lands.

## 1. Setup

Data y_1, …, y_T from the true law P (the exactly simulated channel ensemble), true parameters θ₀. Let F_t = σ(y_1, …, y_t) be the data filtration.

A candidate likelihood is a family of densities. The five algorithms fall into two structural classes:

- **Prequential (recursive):** L(θ) = ∏_t q_t(y_t | y_{1:t−1}, θ). The hidden occupancy distribution is conditioned on the data as it arrives. This is `R`, `MR`, `IR`.
- **Open loop (non-recursive):** L(θ) = ∏_t m_t(y_t | θ), with m_t the marginal predictive propagated from the initial condition, never conditioned on the observed trace. This is `NR`, `NMR`.

In both cases write the per-interval score s_t = ∇_θ log (the t-th factor), evaluated at θ₀, and S = Σ_t s_t.

The two information matrices the diagnostics use:

  J_sample = Σ_t Cov(s_t),  J_total = Cov(S) = J_sample + Σ_{i≠j} Cov(s_i, s_j),

and the correlation distortion R_c = J_sample^(−1/2) J_total J_sample^(−1/2). **R_c = I if and only if the cross-interval score covariances vanish.**

The only regularity needed below: each factor is a normalized density **in y_t for every θ**, and differentiation under the integral sign is licit. Both hold for every member of the family (all factors are Gaussian in y_t).

## 2. Proposition 1 (the classical part)

**If the prequential factors are the true one-step conditionals at θ₀**, that is q_t(· | y_{1:t−1}, θ₀) = p(· | y_{1:t−1}) for every t, **then the score is a martingale difference sequence with respect to F_t**, and therefore J_total = J_sample, R_c = I.

*Proof.* Conditional on F_{t−1},

  E[s_t | F_{t−1}] = ∫ ∇_θ log q_t(y | y_{1:t−1}, θ)|_{θ₀} · p(y | y_{1:t−1}) dy
                   = ∫ (∇_θ q_t / q_t) · q_t dy      (using p = q_t at θ₀)
                   = ∇_θ ∫ q_t(y | y_{1:t−1}, θ) dy |_{θ₀}
                   = ∇_θ 1 = 0.

So {s_t, F_t} is a martingale difference sequence. For i < j, by the tower property,

  Cov(s_i, s_j) = E[ s_i · E[s_j | F_{j−1}]ᵀ ] = 0.

Hence J_total = J_sample. ∎

Two corollaries fall out at once: E[S] = 0 (no score bias), and applying the per-step second Bartlett identity, J_sample = H, so the full distortion C = I. This is the correctly specified case, and it is the benchmark everything is measured against.

## 3. Proposition 2 (the exact identity, and the useful one)

Drop the assumption that the predictive is correct. Define the **conditional score bias**

  **b_t(F_{t−1}) := E[s_t | F_{t−1}]**.

Then, without any approximation,

  **b_t = ∫ ∇_θ log q_t(y | y_{1:t−1}, θ₀) · [ p(y | y_{1:t−1}) − q_t(y | y_{1:t−1}, θ₀) ] dy**,

because the q_t-weighted integral of the score vanishes identically (the step in §2). So b_t is the pairing of the score with the **error of the one-step predictive density**, and b_t ≡ 0 whenever the predictive is right.

And for i < j, again by the tower property,

  **Cov(s_i, s_j) = Cov(s_i, b_j).**

Therefore

  **J_total − J_sample = Σ_{i≠j} Cov(s_i, b_j)   (symmetrized),**

which says:

> **All correlation distortion is the covariance of the past score with the conditional score bias, and the conditional score bias is exactly the error of the one-step predictive. A likelihood has R_c = I if and only if its one-step predictive is the true conditional.**

This is the statement worth having. It converts the correlation distortion from a measured quantity into a **test of the one-step predictive**, and it unifies two of the paper's three diagnostics: the standardized-residual test and the correlation distortion probe the *same* failure, the first through the first two moments of y, the second through the score. They are not independent checks; the residual test is the cheap, parameter-free shadow of the score test.

## 4. Proposition 3 (the open-loop case: unbiased but overconfident)

For the non-recursive members, each factor m_t(y_t | θ) is a normalized density in y_t, so **the first two Bartlett identities still hold per interval**:

  E[s_t] = 0  and  Cov(s_t) = F_t   (provided m_t has the correct marginal moments at θ₀).

That is exactly what Figure 4 shows: **the per-step identity J_t = F_t holds for all five algorithms.** Nothing is wrong per interval.

But there is no filtration argument, because m_t does not condition on y_{1:t−1}. The cross terms do not vanish: y_i and y_j are dependent under the true law (the ensemble current is autocorrelated over the relaxation time), s_i and s_j are functions of them, and so Cov(s_i, s_j) ≠ 0. In the open-loop case the conditional bias is b_t = E[s_t | F_{t−1}], which is generally large, because m_t ignores everything the data have already said about the hidden state.

**Scaling, and a testable prediction.** If the per-interval score has an autocorrelation ρ_k (lag k in intervals) with a finite integral time, then to leading order

  J_total ≈ κ · J_sample,  **κ = 1 + 2 Σ_{k≥1} ρ_k**,

so R_c ≈ κ I, C ≈ κ I, and the reported information is inflated by the factor κ, which is exactly the ratio of nominal to effective sample size. Because the score decorrelates on the process's relaxation time τ, oversampling the interval means

  **κ ≈ 2 τ_eff / Δ**, i.e. the overconfidence of the non-recursive likelihoods should grow linearly with the number of intervals per relaxation time, and should be removed by sampling at Δ ≈ τ.

**Numerical sanity check against what we already have.** The figure-3/4 cell is Δ = 0.1 τ, so this predicts κ of order 10 to 20. Measured (Figure_S3): empirical-to-Fisher variance ratio **10 to 16** for NR and NMR. The order and the magnitude match.

**This is a prediction the existing data can test and nobody has run yet.** The score autocorrelation ρ_k is already plotted (Figure 3, row F), so κ can be computed directly from it and compared, cell by cell, against the measured distortion. And the interval sweep already exists, so the predicted 1/Δ scaling of the non-recursive inflation is checkable across the whole map. If it holds, the paper stops reporting a number and starts explaining it.

## 5. What this says about the five algorithms (the payoff)

The two axes of the family map onto **two different identities**, and that is the cleanest statement of the paper's structure:

| Axis | What it controls | Which identity | Symptom when it fails |
|---|---|---|---|
| **Averaging** (`av` = 0, 1, 2: how much of the interval the conductance is conditioned on) | whether the per-interval predictive has the right moments | the **per-step** Bartlett identities | score bias, and per-step Cov(s_t) ≠ F_t |
| **Recursion** (conditioning the occupancy on the data) | whether b_t = E[s_t | F_{t−1}] vanishes | the **martingale-difference** property | cross-time score correlation, information inflation |

Read the ranking through it:

- **`NR`** gets the averaging wrong (instantaneous conductance, so the per-interval predictive moments are wrong) **and** does not condition. Biased score **and** inflated information. Both identities broken.
- **`NMR`** averages correctly, so its per-interval predictive has the right moments and its score is unbiased; but it does not condition, so b_t ≠ 0 and the information is inflated. **Unbiased and overconfident**, which is exactly the observed combination, and which no purely empirical account explained.
- **`R`** conditions but with a wrong predictive (instantaneous conductance), so b_t ≠ 0 and a residual correlation distortion survives (measured ≈ factor 2).
- **`MR`** conditions with a start-conditioned predictive: better than `R`'s, still not the true conditional, and it drops the boundary cross-covariance term.
- **`IR`** conditions on both interval endpoints, which is the closest available approximation to the true one-step conditional; b_t ≈ 0 and R_c ≈ I, up to the Gaussian closure. Its residual distortion is the residual error of the closure, which is precisely the multinomial and telegraphic degradation the paper maps.

**So Milescu's 2005 question gets a precise answer.** He asked whether estimates must be biased if the method is not a Bayesian filtering algorithm. The answer: **recursion is not what removes the bias, it is what makes the score a martingale difference.** Bias is governed by the other axis (whether the per-interval predictive has the right moments), which is why `NMR` is unbiased without filtering. What filtering buys is the *credibility of the reported uncertainty*. And Milescu's own diagnosis of his dominant error source, "the local time correlation of the current", is b_t, named twenty years before it was computed.

## 6. It also resolves the score-mean puzzle in the supplement

`papers/macroir-elife-2025/analysis_figure_S1_score_mean.md` reports a pattern that looked wrong and was flagged as a possible power artifact: counting the intervals whose mean score is significantly non-zero (out of 40 in-pulse, chance level 2), **`R` and `MR` show 24 to 38, while `NR` and `NMR` show only 5 to 8, and `IR` sits at chance (2 to 5)**. The recursive algorithms looked *more* biased than the non-recursive ones, which contradicts every other result in the paper.

Proposition 2 predicts exactly this, and it is not an artifact.

- For an **open-loop** member, the per-interval factor is a normalized density in y_t with (for `NMR`) the correct moments, so **E[s_t] = 0 exactly**, interval by interval, whatever the autocorrelation of the data. The per-step score mean is zero by construction. Hence 5 to 8 out of 40, at chance.
- For a **recursive** member, E[s_t] = E[b_t], and b_t = 0 only if the one-step predictive is the true conditional. `R` and `MR` condition on the data with a *wrong* predictive, so b_t ≠ 0 and the per-step score mean is systematically non-zero. Hence 24 to 38 out of 40.
- `IR`'s predictive is right up to the Gaussian closure, so b_t ≈ 0 and it returns to chance. Hence 2 to 5.

So the counts are not a bug and not a power effect: **they are a direct measurement of b_t, the conditional score bias**, and they rank the algorithms by how wrong their one-step predictive is. That supplement panel is more valuable than it looked, and it should be re-captioned as a measurement of b_t. **[Verify: rule out the power alternative by reporting the standardized bias alongside the counts. If both readings agree, the panel is a headline.]**

## 7. Honest scope, with the citations (literature check, 2026-07-14)

- **Proposition 1 is textbook, and it has an owner in our exact model class.** **Cappé, Moulines & Rydén (2005), *Inference in Hidden Markov Models*, Springer, ch. 12**: *"the corresponding score increments form, under reasonable assumptions, a martingale increment sequence with respect to the filtration generated by the observations"*, with the limiting Fisher information a single per-step second moment and no cross-lag terms. Also **Hall & Heyde (1980)** ch. 6; **Bickel, Ritov & Rydén (1998)**; **Douc, Moulines & Rydén (2004)**. Prediction-error decomposition: Schweppe 1965; Harvey 1989. **State it as a proposition with citation. Do not claim it.**
- **The converse is also owned**, as a test: it is the null of **White (1987)**'s dynamic information matrix test, and of the m-testing / conditional-moment framework (Newey 1985; Tauchen 1985; Wooldridge 1990); in the prequential frame, Seillier-Moiseiwitsch & Dawid (1993). What is ours is using the cross-lag score covariance as a **quantitative per-algorithm measure** rather than a test statistic, and the identity Cov(s_i, s_j) = Cov(s_i, b_j) read as *all correlation distortion is a wrong one-step predictive*.
- **The non-recursive members are independence composite likelihoods**, and that literature owns the qualitative result. **Chandler & Bate (2007)** *Biometrika* 94:167 coined the term; **Varin, Reid & Firth (2011)** §2.2: *"Composite likelihoods may be seen as misspecified likelihoods, where misspecification occurs because of the working independence assumption … Consequently, the second Bartlett identity does not hold"*. **Concede it in the Introduction.** What is *not* in that literature: the magnitude, the exact per-interval calibration result, and any scaling law. The composite-likelihood papers state H ≠ J and stop.
- **κ = 1 + 2Σρ_k** is the HAC long-run-variance identity (Newey-West 1987; the multivariate ESS form is Vats, Flegal & Jones 2019). **The prediction κ ≈ 2τ/Δ, and its test against the interval sweep, is ours.**
- **Milescu, Akk & Sachs (2005) own the phenomenon in-domain.** See the prior-art map, II.5.2. They state the independence assumption, diagnose the failure as autocorrelation, prove it by decorrelating the data, and predict that only a Bayesian filter would fix it. **Concede this loudly.** And then correct it: they conflate bias with information inflation, and `NMR` separates them (non-recursive, autocorrelated data, unbiased score, over-confident). The answer to their speculation is **no**, a non-filtering method need not be biased; it will be over-confident.

## 8. Related correction

## 7. Related correction

The reconstruction identity in `supplement_information_distortion_main.tex` §10 is printed wrong (C = C_s^(1/2) R C_s^(1/2)); the exact form is C = K R Kᵀ with K = H^(−1/2) J_sample^(1/2). Independent of everything above, but it lives in the same supplement. See `papers/macroir-elife-2025/correction_idm_reconstruction.md`.
</content>
