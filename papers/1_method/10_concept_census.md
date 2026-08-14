# Concept census: where each concept is defined, and what it rests on

Written 2026-08-12, after the carve of `09_carve_plan.md`. One row per concept the manuscript uses:
the symbol if it has one, the place that **defines** it, and the concepts it depends on. Sorted into
layers, so the table reads as a dependency order and not as a glossary.

**What this does not duplicate.** `_program/nomenclature.md` owns the letters and what they mean;
`_program/machinery.md` owns the diagnostic definitions and sign conventions;
`_program/notation_map.md` owns symbol↔code↔CSV. None of them says **where in the manuscript** a
concept is defined, which is what makes a forward reference findable. That is this file's only job.

**How the first-use column was measured, so it can be re-run.** A script over
`sections/` in reading order (abstract, introduction, framework, results, discussion, methods,
appendices), comments stripped, whitespace collapsed so that a definition broken across two source
lines still matches. First use is the first occurrence in that order. Anything whose first use
precedes its definition is a forward reference, and §6 lists them.

Reading order after the merge: **00** abstract, **01** introduction, **02** framework
(= the merged Theory + Diagnostics), **04** results, **05** discussion, **06** methods,
**08** Appendix 1 derivation, **09** Appendix 2 members, **10** Appendix 3 diagnostics,
**11** Appendix 4 repairs.

---

## Layer 0. The physical model

| Concept | Symbol | Defined in | Depends on |
|---|---|---|---|
| macroscopic current | — | 02 opening | — |
| channel count | $N_{\mathrm{ch}}$ | 02 opening | — | <!-- pat: N_\{\\mathrm\{ch\}\} | def: 02 -->
| kinetic states, state count | $K$ | 02 opening | — |
| generator / rate matrix | $\mathbf{Q}$ | 02 opening | $K$ | <!-- pat: \\mathbf\{Q\} | def: 02 -->
| transition matrix | $\mathbf{P}(t)$, $P_{i\to j}$ | 02 opening | $\mathbf{Q}$ |
| per-state conductance | $\bm{\gamma}$ | 02 opening | $K$ | <!-- pat: \\bm\{\\gamma\} | def: 02 -->
| conductance diagonal | $\bm{\Gamma}=\mathrm{diag}(\bm{\gamma})$ | 08 §boundary-conditioned moments | $\bm{\gamma}$ | <!-- pat: \\bm\{\\Gamma\} | def: 08 -->
| occupancy mean, covariance (per channel) | $\bm{\mu},\bm{\Sigma}$ | 02 §interval update; algebra in 08 | $K$ |
| ensemble counts | $\mathbf{N}$ | 08 §interval average | $\bm{\mu},\bm{\Sigma},N_{\mathrm{ch}}$ |
| initial occupancy | $\bm{\pi}$, Eq. init | 06 §the initial condition | $\bm{\mu},\bm{\Sigma}$ |
| two-state scheme C/O, the six parameters | $k_{\mathrm{on}},k_{\mathrm{off}},i,b$ | 06 §minimal model | $\mathbf{Q},\bm{\gamma}$ |
| open probability of the sweep | $p=0.5$ | 06 §minimal model | the six parameters | <!-- pat: open probability | def: 06 | ok: 01 -->
| instrumental noise, white | $\epsilon^2$, `Current_Noise` | 06 §minimal model | — |
| closing time constant | $\tau=1/k_{\mathrm{off}}$ | 06 §minimal model | $k_{\mathrm{off}}$ |
| exact simulation by uniformization | — | 06 §minimal model | the model |

## Layer 1. The observable

| Concept | Symbol | Defined in | Depends on |
|---|---|---|---|
| acquisition window | $\Delta$ | 02 §the observable | — | <!-- pat: \\Delta\b | def: 02 -->
| interval average (the observable) | $\bar y_t$, Eq. obs-avg | 02 §the observable | $\Delta$, the current, $\eta_t$ | <!-- pat: interval average | def: 02 -->
| instrumental term of one sample | $\eta_t$ | 02 §the observable | $\epsilon^2$ |
| uniform window (boxcar) | — | 02 §the observable | $\bar y_t$ |
| amplifier kernel, and what the window misses | $f_c\Delta$ | **08 §the acquisition kernel** (bound, table, four properties); scope sentence in 02 | uniform window | <!-- pat: f_c\\Delta | def: 08 -->
| block-averaging a record | — | 08 §the acquisition kernel | $f_c\Delta$ |

## Layer 2. The exact objects, and why they are unavailable

| Concept | Symbol | Defined in | Depends on |
|---|---|---|---|
| microscopic description / filter | — | 02 opening (cost), 08 §interval average | $\mathbf{N}$ |
| integrated conductance over one window | $A_{0\to\Delta}$, Eq. integrated-conductance | 08 §interval average | $\bm{\gamma}$, the path | <!-- pat: A_\{0\\to\\Delta\} | def: 08 -->
| telegraph-like law, its atom and series | — | 02 box (claim), 08 §interval average (law) | $A_{0\to\Delta}$, $K$ |
| exact likelihood, at the top of the ladder | — | 02 §two axes | the two above |

## Layer 3. The two closures, and the regimes

| Concept | Symbol | Defined in | Depends on |
|---|---|---|---|
| closure (the operation) | — | 02 §two Gaussian approximations | — |
| occupancy closure, and its two validity conditions | — | 02 box | $\bm{\mu},\bm{\Sigma}$, $N_{\mathrm{ch}}$, $\epsilon^2$ | <!-- pat: occupancy closure | def: 02 -->
| interval-signal closure, and its validity | — | 02 box | $A_{0\to\Delta}$, $\Delta$ | <!-- pat: interval-signal closure | def: 02 -->
| microscopic regime | — | 02 §two Gaussian approximations | occupancy closure |
| telegraphic regime | — | idem | interval-signal closure |
| Gaussian regime | — | idem | both closures |

## Layer 4. The family

| Concept | Symbol | Defined in | Depends on |
|---|---|---|---|
| window axis | — | 02 §least squares and the two axes | interval average | <!-- pat: window axis | def: 02 -->
| recursion axis (prequential / open loop) | — | idem | the likelihood factorization | <!-- pat: prequential | def: 02 -->
| variance form (total / residual) | — | 02 naming rule; algebra in 09 switch 3 | interval variance |
| naming rule (prefix, suffix, V) | — | 02 naming rule | the three axes above |
| the eight members | `LSE ILSE NR INR R MR VR IR` | 02 Table 1 (words), 09 `tab:appendix-members` (equations), 06 Table 2 (flags) | naming rule |
| least squares as the anchor, off the lattice | — | 02 §least squares; 09 §least squares | window axis |
| cost ladder | — | 02 §least squares | the members | <!-- pat: cost ladder | def: 02 | ok: 01 -->
| boundary state | $(i\to j)$ | 02 §interval update; device reading in 08 | $\mathbf{P}(\Delta)$ | <!-- pat: boundary state | def: 02 -->
| the update, gain, predicted mean and variance | Eq. macror, $\mathbf{g}$ | 02 §interval update; general form in 08 | $\bm{\mu},\bm{\Sigma},\bm{\gamma},\epsilon^2$ |
| sufficiency, and what the update gives up | — | 02 §interval update; concession in 05 | the update |
| posterior-as-prior closure | — | 02 §interval update | the update, recursion axis |
| boundary-conditioned single-channel moments | $\bar\gamma_{i\to j},\bar v_{i\to j}$, Eq. gammabar-vbar | 08 §boundary-conditioned moments | $A_{0\to\Delta}$, boundary state | <!-- pat: \\bar\\gamma_\{i\\to j\}|\\bar\{\\gamma\}_\{i\\to j\} | def: 08 | ok: 02 -->
| tilted-semigroup derivatives | $\mathbf{F}_1,\mathbf{F}_2$ | 08 idem | $\mathbf{Q},\bm{\Gamma},\Delta$ | <!-- pat: \\mathbf\{F\}_1 | def: 08 | ok: 06 -->
| spectral evaluation, divided differences | $E_2,E_3$ | 08 idem | $\mathbf{F}_1,\mathbf{F}_2$ |
| the tilde contraction on pairs | $\widetilde{\;\cdot\;}$, Eqs. tilde, vector-tilde | 08 §three steps; restated 09 §the tilde operator | boundary state, $\bm{\Sigma}$ |
| the three steps of the construction | Eqs. mu-bnd…sig-post | 08 §three steps | all of the above |
| MR against IR: same variance, different gain | Eq. mr-ir-identity | 08 §what separates them; 09 | the update, variance form |
| linear-filtering correspondence | Eqs. kf… | 08 §linear-filtering frame | the update |
| pseudo-count regularizer | $\varepsilon$, Eq. pseudocount | 06 §safeguards | $\bar\gamma_{i\to j}$, $P_{i\to j}$ | <!-- pat: pseudo-count | def: 06 -->
| trust coefficient | $\alpha_\mu$ | 06 §safeguards | the update | <!-- pat: trust coefficient | def: 06 -->
| binomial floor, minimum occupancy | $5$, $10^{-12}$ | 06 §safeguards | the update |

## Layer 5. Estimation and the two anchors

| Concept | Symbol | Defined in | Depends on |
|---|---|---|---|
| the log-likelihood each member scores with | Eq. loglik | 06 §the likelihood | predicted mean and variance |
| simulation truth | $\theta_{\mathrm{sim}}$ | 02 §what is measured; 06 §two anchors | the model | <!-- pat: \\theta_\{\\mathrm\{sim\}\} | def: 02 -->
| pooled optimum | $\theta_{\mathrm{pool}}$ | idem | the likelihood | <!-- pat: \\theta_\{\\mathrm\{pool\}\} | def: 02 -->
| misspecification bias | $\theta_{\mathrm{pool}}-\theta_{\mathrm{sim}}$ | idem | the two anchors |
| per-group MLE cloud, group size | — | 06 §parameter estimation | the likelihood |
| optimiser, warm start, Newton decrement | — | 06 §parameter estimation | — |
| percentile bootstrap over groups | — | 06 §diagnostics | the cloud |
| retained subspace of the anchor | $\lambda_i>\lambda_{\max}10^{-10}$ | 06 §diagnostics | $\mathbf{H}$ |

## Layer 6. The diagnostics

| Concept | Symbol | Defined in | Depends on |
|---|---|---|---|
| standardized residual, its three properties | $r_t$ | 02 §what is measured | predicted mean and variance | <!-- pat: standardized residual | def: 02 -->
| integrated autocorrelation of the residual | $\tau_{\mathrm{int}}$ | **04 prose and the Figure 5 caption only** | $r_t$ | <!-- pat: integrated autocorrelation | def: 02 | ok: 01 -->
| score, total and per interval | $S$, $s_t$ | 02 §what is measured | the log-likelihood |
| Gaussian Fisher information (the anchor) | $\mathbf{H}$, $\mathrm{GFI}_t$ | 02 (named), **06 §diagnostics (the display)** | predicted moments and their $\theta$-derivatives | <!-- pat: Gaussian Fisher information | def: 02 -->
| per-interval and accumulated information | $F_t,F_T$ | 02 §what is measured | $\mathbf{H}$ |
| score covariance, total and within-interval | $\mathbf{J},\mathbf{J}_s$ | 02 §what is measured | $S$, $s_t$ |
| information distortion | $\mathbf{C}$ | 02 §what is measured | $\mathbf{H},\mathbf{J}$ | <!-- pat: \\mathbf\{C\}=\\mathbf\{H\} | def: 02 -->
| sign convention | $C_{ii}>1$ over-confident | 02 §what is measured | $\mathbf{C}$ |
| sample distortion | $\mathbf{C}_s$ | 02 §what is measured | $\mathbf{H},\mathbf{J}_s$ | <!-- pat: sample distortion | def: 02 -->
| correlation distortion | $\mathbf{R}$ | 02 §what is measured | $\mathbf{J}_s,\mathbf{J}$ | <!-- pat: correlation distortion | def: 02 -->
| how the two compose | $\mathbf{C}=\mathbf{K}\mathbf{R}\mathbf{K}^\top$ | 10 §how the two parts compose | $\mathbf{C},\mathbf{C}_s,\mathbf{R}$ |
| effective sample size | $\kappa$, $T_{\mathrm{eff}}$ | 02 §what is measured | $\mathbf{R}$ |
| magnitude and anisotropy | $m$, $a$, Eq. magnitude-anisotropy | 02 §what is measured | eigenvalues of $\mathbf{C}$, retained subspace | <!-- pat: multiplicative spread | def: 02 -->
| affine-invariant distance | $D_{\mathrm{AI}}$ | 10 §the two scalars | $m,a$ | <!-- pat: D_\{\\mathrm\{AI\}\} | def: 10 -->
| finite-sample floor of the anisotropy | $1.081$ | 10 §the two scalars | $a$ |
| sandwich covariance | $\bm{\Sigma}=\mathbf{H}^{-1}\mathbf{J}\mathbf{H}^{-1}$ | 02 §what is measured; expansion in 10 | $\mathbf{H},\mathbf{J}$ | <!-- pat: sandwich | def: 02 -->
| distortion-induced bias | $b=\mathbf{H}^{-1}\mathbb{E}[S]$ | 02 §what is measured; expansion in 10 | $\mathbf{H}$, $\mathbb{E}[S]$ |
| coverage of the nominal region | — | 04 §recovery at one cell | sandwich |
| least-squares premise (constant variance) | — | 02 (clause), 04 (in full) | $\mathbf{H}$, LSE |

## Layer 7. The design space and the map

| Concept | Symbol | Defined in | Depends on |
|---|---|---|---|
| dimensionless acquisition interval | $\widetilde{\Delta}=\Delta k_{\mathrm{off}}$ | **06 §regime grid**; glossed at 04 opening | $\Delta,\tau$ | <!-- pat: \\widetilde\{\\Delta\} | def: 02 -->
| dimensionless instrumental noise | $\widetilde{S}=Sk_{\mathrm{off}}/i^2$ | **06 §regime grid**; glossed at 04 opening | noise, $k_{\mathrm{off}}$, $i$ | <!-- pat: \\widetilde\{S\} | def: 02 -->
| the noise ladder, the swept grid | Eq. noise-ladder | 06 §regime grid | $\widetilde{S}$ |
| design plane (channel count against noise) | — | 02 box (why separable), 04 | the two axes |
| the two criteria | $1.15$ and $2$ | 02 (error-bar units), 04 and figure captions | $\mathbf{C}$ |
| the four boundaries and five regions | — | 04 §a usage map | the criteria, the corrected standard error | <!-- pat: five regions | def: 04 -->

---

## 5b. State: all six repaired, 2026-08-12

Every defect in §6 was fixed the same day, on Luciano's authorisation. What changed, in the
manuscript and not only here:

1. Table 1 gained a fifth column, "Called, in the Results", and the naming paragraph gained one
   sentence binding each periphrasis to its code. The two mappings that were not obvious were read
   off their use sites first: MR is *the one-endpoint member* and VR *the variance-corrected rung*,
   from the paragraph that runs both against R over four channel counts, where "the member it
   corrects" is R. *Window-ignoring* was deliberately left unbound: it names a set, the three
   members that predict an instant, and the sentence that uses it names them.
2. $\tau_{\mathrm{int}}$ is now defined in "What is measured", with the form the producer computes:
   $1+2\sum_k\rho_k$, read from `legacy/moment_statistics.h:1680-1730`, whose own header calls it
   the Kish/Sokal integral autocorrelation time with positive and negative lags counted
   symmetrically. Methods was corrected in the same pass: its max-lag sentence said "for the score
   autocorrelation" when `max_lag = 10` is one battery parameter governing both.
3. $\widetilde{\Delta}$ and $\widetilde{S}$ moved from the Results opening into the framework. Net
   cost zero.
4. The display of $\mathbf{H}$ moved back into the body as Eq. gaussian-fisher and Methods now
   points at it. Costs no counted words, and it is a half-reversal of a cut made in
   `09_carve_plan.md`.
5. $\kappa$ and $T_{\mathrm{eff}}$ are gone as symbols; the idea stays in words.
6. $\mathbf{F}_1,\mathbf{F}_2$ gained a gloss in Methods, where they are used, ahead of the appendix
   that derives them. This was the one live hit of the new check.

**And the census now has a keeper.** `_program/concept_firstuse.py` re-measures first use for every
row carrying a `<!-- pat: … | def: … -->` marker and reports any concept used before its defining
section; `check.sh` runs it as item 10 and prints how many rows are covered, so partial coverage
never reads as full coverage. 33 rows carry a pattern today and the check is at zero. A row may
declare `ok: NN`, the earliest section where an undefined mention is deliberate: the abstract is
exempt by construction, and four rows are marked `ok: 01` because the Introduction names them
before the body defines them, which is what an introduction is for.

## 6. What the census found

Six defects, each of them a place where a reader meets a thing before the paper has told them what
it is. Ordered by how much they cost.

**(1) The Results' name for the paper's own member is never bound to its code.** Counted, not
guessed: *boundary-conditioned* occurs 12 times in the Results, 6 in Methods and 4 in the
Discussion, and **zero times in the framework section**. *Open-loop member* is 8 and 0,
*window-ignoring* 3 and 0, *the classical fit* 4 and 0. Only *recursive instantaneous* is bound in
the framework, twice, and by accident of the prior-art sentence rather than on purpose. So a reader
who has learned the naming rule and Table 1 meets, one page later, four names that appear in
neither. The repair is one sentence in the naming paragraph binding each periphrasis to its code, or
a column in Table 1.

**(2) $\tau_{\mathrm{int}}$ is defined in a caption and in passing.** The Introduction uses the
integrated autocorrelation as the reader-facing statistic before any definition exists; the symbol
first appears in the Figure 5 caption; the prose gloss is a subordinate clause in the Results. It is
the one quantity the paper offers an experimentalist to compute on their own recording, and it is
the least well defined thing in it. It belongs in "What is measured", beside $r_t$, in one sentence.

**(3) Two dimensionless design axes are defined after the section that uses them.**
$\widetilde{\Delta}$ and $\widetilde{S}$ are defined in Methods, which eLife prints after the
Results. The Results opening glosses both and points at Methods, which is the mitigation and is
probably enough; recorded here so the decision is deliberate.

**(4) The anchor $\mathbf{H}$ is named in the body and displayed only in Methods.** By decision of
the carve (the display was printed twice); the framework now says what it is and where the formula
is. Recorded because it is the one diagnostic whose formula a reader cannot see where it is used.

**(5) Concepts defined and then used only in words.** $\kappa$ and $T_{\mathrm{eff}}$ are defined
with symbols in the framework and used nowhere else in the manuscript: the Results say "no single
effective count serves both" and "how many intervals a likelihood effectively has", never the
symbols. Verified by grep, both symbols, all sections. Either they earn their place downstream or
the framework defines the idea in words and drops them. $s_t$ is the same case; $F_t$ is not, being
used in the Results and in three figure captions.

**(6) $\bm{\Gamma}$, $\bm{\pi}$, $E_2$, $E_3$ and $D_{\mathrm{AI}}$ live entirely outside the body.**
That is correct after the carve and is listed so that nobody "restores" them to the body on the
grounds that they are undefined there.

## 7. What the census confirms

The dependency order matches the section order. Every layer 0 to 4 concept the Results use is
defined in the framework or in Methods, no layer-6 concept depends on anything defined after it, and
the appendices depend on the body and never the reverse. The one genuine cycle the merge could have
created, the framework's "what is measured" depending on Methods for the anchor display, is a
pointer and not a definition.
