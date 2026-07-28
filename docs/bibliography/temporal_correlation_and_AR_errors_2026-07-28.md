# Correlated residuals: the AR/ARMA route, and who used temporal correlation before MacroR

**Written 2026-07-28**, answering the two questions left open in the audio notes of the same day
(`program/source-notes/audios/audios/Chat de WhatsApp con MacroIR 13/`, 15.12 and 15.17):

1. What are the other methods for handling non-independent residuals (autoregressions), and why do
   they not apply here?
2. Is MacroR really *the* classical method for temporal correlation in macroscopic currents, or is
   there another one?

**Short answers.** (1) There is a full classical ladder, all of it whitening; one paper has applied
ARMA to ion-channel calibration and it is worth citing. (2) **No.** The correlation route in this
field is from 1973, and the direct-fit covariance likelihood (Celentano & Hawkes) predates MacroR by
three years. The Introduction must not claim absence. It should claim that nobody characterized when
these methods are valid, which is true and is what the paper does.

All PDFs added to `docs/bibliography/`; all bib entries added to
`papers/1_method/docs/manuscript-drafts/biblio_full.bib` (113 -> 124 entries).

---

## Part 1 — the AR/ARMA route

### 1.1 The classical ladder

| Job | Method | Reference |
|---|---|---|
| Detect correlation | Durbin-Watson (AR(1)); Ljung-Box portmanteau; Breusch-Godfrey LM | `durbin1950testing`, `ljung1978measure`; Breusch-Godfrey `[VERIFY, not in bib]` |
| Fix the estimator by whitening | GLS. AR(1) quasi-differencing (Cochrane-Orcutt; Prais-Winsten keeps observation 1); general ARMA(p,q) errors, Box-Jenkins | `cochrane1949application`, `prais1954trend` |
| Fix only the standard errors, keep the OLS point estimate | HAC / long-run variance | `newey1987simple`, Andrews 1991 (already in bib) |
| Same thing in state-space form | any ARMA error model has a state-space representation, estimated by a Kalman filter | `harvey1989forecasting` (already in bib) |

### 1.2 The structural point (this is the Introduction in one sentence)

**Every one of these whitens. What separates them is where the whitening covariance comes from.**

- AR/ARMA/HAC: fitted as **free nuisance parameters**, estimated from the residuals themselves.
- MacroR / MacroIR: **derived from the same rate constants that generate the mean**.

That is exactly why the mechanistic route gains information instead of merely widening error bars.
Under an ARMA error model the correlation is a nuisance to be neutralized. Under the Markov model the
correlation is a second, independent measurement of the same theta. Least squares throws that
measurement away; an ARMA correction throws it away too, and additionally spends parameters doing so.

### 1.3 Why ARIMA specifically fails on macroscopic currents

The audio gives two reasons (finite channel population, concentration/voltage jumps). Both are right;
here they are in the form a referee cannot argue with, all four derivable from machinery already in
the paper.

1. **Nonstationarity.** The gating covariance is
   `C(t1,t2) = N i^2 p(t1) [ p11(t2|t1) - p(t2) ]` (Sigworth 1981, Eq. 4). Its diagonal
   `N i^2 p(1-p)` tracks the mean current and vanishes when the current does. An ARMA error model is
   **stationary by construction**: one covariance function for the whole trace. After a concentration
   or voltage jump the residual process restarts; an ARIMA fit has to average across the restart.
2. **No `N`.** Signal scales as `N`, gating noise as `sqrt(N)`. That ratio is what identifies
   `N_ch` and the unitary conductance. An ARMA error model contains no `N` and no unitary current, so
   it cannot deliver either. Those are precisely the parameters MacroIR gets for free, and they are
   the ones Del Core & Mirams 2025 declare still unclear.
3. **Untied timescales.** ARMA timescales are free parameters; in the Markov model they are the
   eigenvalues of `Q`, already constrained by the mean current. Free timescales absorb exactly the
   information the mechanism would have contributed. **This is the information-budget figure stated in
   words.**
4. **Whiteness alone is not falsification.** Adding enough ARMA terms will whiten anything. So residual
   whiteness cannot certify that the kinetics are right, which is the argument for the score and
   information-equality tests on top of test (i).

### 1.4 The one in-domain precedent: Lei et al. 2020

`lei2020considering` — *Considering discrepancy when calibrating a mechanistic electrophysiology
model*, Phil. Trans. R. Soc. A 378:20190349. PDF:
`Lei_Mirams_2020_Discrepancy_ARMA_GP_Electrophysiology_PhilTransA.pdf`.

ARMA(2,2) plus two Gaussian-process variants (`GP(t)` and `GP(O,V)`) on the residuals of a hERG
kinetic model, synthetic data, i.i.d. Gaussian observation noise sigma = 25 pA. Verbatim, on the ARMA:
*"we do not condition on the observed discrepancy sequence ... but only use it to correlate the
discrepancy structure in time."*

Findings, in their words:

- ignoring discrepancy gives *"model parameters that are wrong, and yet we are certain about this
  wrong value"* -> *"spuriously certain parameter inference and overly-confident and wrong
  predictions"*;
- *"the choice of the discrepancy model can shift the posterior distribution significantly, both in
  terms of its location and spread"*;
- the ARMA(2,2) *"increases the width of the posterior (compared to i.i.d. noise), but its posterior
  mean prediction does not follow the data as closely as the two GP models"*;
- and the killer: *"in some cases, the best predictions were still made by ignoring discrepancy."*

**How to use it.** It is the empirical demonstration that the phenomenological route has been tried in
this exact domain and buys uncertainty inflation without a reliable gain. That sets up the mechanistic
alternative better than any argument from first principles.

**Caveat, state it or a referee will.** Their discrepancy is *model* discrepancy (wrong kinetic
scheme, structural error), not finite-channel gating noise. Different target. Their ARMA is a patch
over a misspecified mean; ours is a correctly specified second moment. The honest sentence is that
their result shows what a correlation model without a mechanism can and cannot buy.

---

## Part 2 — who used the temporal correlation before MacroR

**The audio's claim ("the classical method for temporal correlation is MacroR, I don't know if there
is another") is not correct, and it is the kind of sentence a reviewer of this field kills in one
line.** Three tiers, oldest first.

### Tier 1 — spectral, stationary conditions (1970s)

- `katz1970membrane`, `katz1972statistical` — the origin: the fluctuations of the macroscopic current
  carry single-channel information.
- `anderson1973voltage` — **Anderson & Stevens 1973**, J. Physiol. 235:655-691. The channel closing
  rate `alpha` read off the **Lorentzian power spectrum** of the end-plate current, with its voltage
  and temperature dependence (`alpha = B e^{AV}`, `B = 0.17 +- 0.04 ms^-1`, `A = 0.0058 +- 0.0009
  mV^-1` at 8 C). PDF in repo.
- `neher1977conductance` — the review that consolidated it.

The power spectrum **is** the Fourier transform of the autocovariance. So the field has been reading
kinetics out of the temporal correlation of macroscopic currents since 1973. What it required:
stationarity, and a model simple enough that one corner frequency identifies one rate.

### Tier 2 — the nonstationary two-time covariance (1980-81)

- `conti1980conductance` — the two-conductance-level covariance expression.
- `sigworth1981covariance` — **Sigworth 1981**, Biophys. J. 34:111-133. PDF in repo.
  **This is the closest in-domain ancestor and it was missing from the prior-art map.**
  Derives the exact two-time covariance for `N` identical independent Markov channels,
  `C(t1,t2) = N sum_k sum_j p_j(t1) p_jk(t2|t1) [i_j i_k - i_j mu(t2)]` (Eq. 3), reducing for two
  conductance levels to `C(t1,t2) = N i^2 p(t1)[p11(t2|t1) - p(t2)]` (Eq. 4), whose diagonal is the
  NSFA variance `N i^2 p (1-p)` (Eq. 5). Estimated from ensembles of 256-504 repeated sweeps and used
  to discriminate kinetic hypotheses (one open state vs several vs several channel populations).

  **This is the same object MacroR propagates recursively.** The differences, and they are the whole
  gap: it is estimated **empirically from replicate ensembles** rather than predicted from the model
  inside a likelihood; it is used as a **discriminator**, not to fit rate constants; and there is no
  uncertainty statement on anything.

### Tier 3 — likelihoods that use the correlation (2004 onwards)

- `celentano2004use` — **Celentano & Hawkes 2004**, "covariance fitting" (CVF). Verbatim from the
  abstract: *"Unlike conventional sum-of-squares minimization, CVF fits both the magnitude of the
  recorded current and the strength of the correlations between different time points."*
  Q-matrix covariance + maximum likelihood, `O(n^3)`. **Predates MacroR by three years and is the
  direct-fit statement of the idea.**
- `stepanyuk2011efficient` / `stepanyuk2014maximum` — same likelihood, semiseparable covariance,
  near-linear cost.
- `milescu2005maximum` — Gaussian in the point current: has the gating **variance** but not the
  cross-time correlation. (And names *"the local time correlation of the current"* as his dominant
  error source.)
- `moffatt2007estimation` — MacroR, the recursive form.
- `munch2022bayesian` — Kalman, Bayesian.
- MacroIR — adds the within-interval integration.

### 2.1 Consequence for the Introduction

Do **not** write "MacroR is the classical method for temporal correlation" and do **not** write that
the correlation was unused. Both are false and both are cheap to refute.

**Write this instead** (it is true, verifiable, and a better story):

> The kinetic information in the temporal correlation of macroscopic currents has been recognized
> since the 1970s, first through the power spectrum under stationary conditions (Katz & Miledi 1970;
> Anderson & Stevens 1973), then through the nonstationary two-time covariance (Conti et al. 1980;
> Sigworth 1981), and finally inside a likelihood (Celentano & Hawkes 2004; Milescu et al. 2005;
> Moffatt 2007; Stepanyuk et al. 2011; Münch et al. 2022). Yet the field still fits macroscopic
> currents by least squares on the deterministic mean (Clerx et al. 2019; Wang et al. 2012; Owen &
> Mirams 2025). The obstacle is not that the correlation was unknown: it is that no criterion exists
> for deciding when a method that exploits it is working.

The gap moves from *"nobody used the correlation"* (false, killable) to *"nobody said when using it
is valid"* (true, and it is exactly this paper). It also gives the region figure its meaning: the
regions are the answer to a question the field has never been able to ask.

---

## Loose ends

- **Sigworth 1980**, J. Physiol. 307:97-129 (*The variance of sodium currents at the node of
  Ranvier*, PMID 6259340, PMC1283036) — **not retrieved**, PMC blocked the download. It is the
  variance/NSFA companion to the 1981 covariance paper. Get it if the Introduction needs the NSFA
  citation from the source rather than from Stepanyuk's summary.
- **Breusch 1978 / Godfrey 1978** (LM test for autocorrelation, valid with lagged dependent
  variables) — citation not verified, deliberately kept out of the bib. Verify before use.
- **Lei et al. 2020 author list** — taken from the arXiv version. Verify against the published Phil.
  Trans. version before it enters a manuscript.
- `MacroIR_prior_art_map.md` should get Sigworth 1981 and Anderson & Stevens 1973 into Part I section
  A. Their absence was the map's one real in-domain hole: it had the *likelihood* users of the
  correlation but none of the *pre-likelihood* ones.
