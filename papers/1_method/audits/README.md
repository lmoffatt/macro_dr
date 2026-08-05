# Audits

> Multi-agent audits of the manuscript, kept because their **negative results are as valuable as
> their findings**: the "what passed" list is a record of what has already been checked, so the same
> ground does not get re-audited, and the refuted claims are a record of what looks like a defect and
> is not. Raw agent output is beside this file as `raw_*.json`; the curated reading is below.
>
> Both runs are dated 2026-08-04. Neither replaces a human read.

---

## 1. Completeness audit (`raw_completeness_20260804.json`)

**The question asked:** is the information a competent reader needs in order to UNDERSTAND the results
actually present in the manuscript? Not readability, not order, not style.

**Design.** Six auditors, one lens each (the design and the data; the seven members; the instrument;
every number in Results and Discussion; the six figures; claim against evidence). Each read the nine
section files directly. Every reported gap was then attacked by an independent agent whose default was
to refute it, since the dominant failure mode of a completeness audit is reporting as missing from
Results what is defined in Methods. Then one synthesis over the survivors. Thirteen agents, 2.06M
tokens, 23 minutes, no errors.

**Verdict: complete with named holes.** 76 gaps reported, 11 survived refutation.

### The two structural ones, each found independently by four of the six lenses

**A. The information-budget subsection has no figure in the `.tex`.** `04_results.tex` says "A
calibrated member's reported covariance can be spent, and Figure~\ref{fig:design} spends it", but that
float pulls `Figure_5.pdf`, whose caption is about the departure corner and its per-sample/cross-time
decomposition. The five fan factors (1.24, 3.11, 288.19, 60.50, 3236.85) and the design advice drawn
from them have no display item.
**Cheap fix, and the audit did not know it:** `Figure_5_budget.pdf` exists in `figures/paper_both/`,
built, with `Figure_5_budget_caption.md` and a supplement, and its caption opens on exactly what the
subsection describes. The figure is not missing; the float is. Its own caption says "budget, demoted",
so it was demoted at some point and the Results prose was never told. **The `.tex` carries six body
floats and the canonical set (the figures with a `_caption.md` beside them) has seven pieces.**

**B. The coverage ensemble is internally contradictory.** `04_results.tex:54` measures coverage "over
a thousand fits in six parameters at group size 100", while `06_methods.tex:322-323` states that
10,000 recordings give "1,000 groups of 10 and 100 groups of 100". At group size 100 there are 100
fits, not 1,000. Since the Figure 2 caption scales reported covariances "by the reciprocal of the
group size", a reader cannot tell what the coverage numbers 0.914 and 0.947-0.949 are covariances of.
**This is the only survivor that needs someone to look at data rather than to write.**

### The rest, all fixable by stating a condition already known

3. **The Introduction's collapse claim has no measurement in the body** (two lenses). It says the
   channel count and the noise enter through the residual autocorrelation "and through nothing else".
   That is measured, in `figures/paper_both/figure_6B.Rmd` (header, points 1 and 2), and it is not in
   the manuscript. Move the measurement into Results with its `% src:`; no run needed.
4. **Which parameter the accumulated information ratios belong to** (two lenses). 14.96 / 13.63 /
   1.07 / 1.08 at `04_results.tex:91`, and the quantity is per-parameter everywhere else.
5. **The conditions of the residual autocorrelations** (two lenses). 1.008 / 1.30 / 3.93 / 9.26 with
   no channel count, no interval, and no statement that they are medians.
6. **"Pulse information"** at `04_results.tex:102` (a sentence added the same day): the quantity, its
   segment, its aggregation and its denominator are undefined, and "the two non-recursive members"
   does not identify a set when three are drawn.
7. **Which one-endpoint member** reaches 20% at `04_results.tex:177`; MR and VR both qualify. The
   Discussion instance was refuted, because `05_discussion.tex:32` pins MR by a property only MR has.
8. **The Kalman agreement** at `05_discussion.tex:58`, "to about $10^{-8}$", with no quantity, no norm
   and no case count.
9. **The numerical cutoff** constant at `03_diagnostics.tex:54` is referred to and never given.
10. **The affine-invariant distance** at `04_results.tex:213`: no formula, no null, no citation. No
    body number depends on it.
11. **Small and local:** `p` unbound in $\alpha^\star=p/\operatorname{tr}\mathbf{C}$; the cost ordering
    in "better calibrated interval by interval than any of the three recursive members that cost more"
    rests on no cost measurement; "a lower bound on the error for anything richer" has no stated link
    from state count to closure quality; the three preparation footprints lack ranges and sources;
    "the single-molecule integrated filters" is uncited; Figure 1's interval and noise level are
    absent; Figure 6's drawn axis is in $\tilde S$ while its caption says label units.

### Two arithmetic errors the audit did NOT surface, verified by hand afterwards

- **The τ_int gloss is inverted.** `04_results.tex:187` says it "counts how many intervals a likelihood
  effectively has for every interval it believes it has" and then reads 9.26 for least squares. τ_int
  is N/N_eff, so 9.26 is nine believed per effective, the other way round, and the next sentence draws
  the correct conclusion, contradicting the gloss. The same inverted phrasing is in `figure_6B.Rmd`
  and was copied into `01_introduction.tex` on 2026-08-04.
- **Area quoted as width.** "under-reports the width of the parameter distribution by a factor of ten
  to fifteen" follows two sentences that measure a ratio of ellipse AREAS. The width factor is the
  square root, about three and a half. Same square-root family as the τ_int error-bar scale, so treat
  it as systematic and sweep for it.

### What passed, so it is not re-audited

The kinetic scheme, generating parameters, protocol, observable, simulator exactness and the seed
caveat; the noise-label conversion and the interval axis; the eleven channel counts, seven intervals
and the honestly non-rectangular per-member coverage; all seven members defined three times
consistently, with the MR/VR/IR algebra written out and the MR-IR identity proved; both least-squares
configurations and the 4x4 versus 6x6 incommensurability; score, Bartlett identities, **J**, **H**,
**C** with its sign convention, the sandwich and its verification against the clouds, the sample and
correlation split with the log-det-only caveat, kappa and T_eff, both anchors and which figure uses
which; the five usage regions; the bootstrap and its half-widths.

**Refuted, and worth recording because they look like defects and are not:** the alleged
inconsistency of the noise axis between prose and captions; the alleged contradiction in the +-0.184
band; the two corrected bands; the region numbering; and the Discussion's one-endpoint member.

### Not checkable by this run

The declared Supplementary Information is outside the nine section files, so the tilde derivation and
the sandwich Taylor expansion could not be verified as present. Supplement graphics beyond Figures 4,
5 and 6 were not rendered. Arithmetic correctness was out of scope by design, which is why the two
errors above had to be found by hand.

---

## 2. Introduction thread bake-off (`raw_intro_bakeoff_20260804.json`)

Five Introductions written on five different narrative threads from one shared fact packet, at fixed
length, then three critique lenses: one structural over all five, one adversarial referee over three,
and five blind readers, one per draft.

**Use it for the structural and adversarial findings. Discount the readers.** The five "readers" were
agents role-playing an electrophysiologist; their testimony is generated, not measured, and the five
drafts shared the same fact packet, so their convergence partly measures the imposed material rather
than the five architectures. That caveat is the most important thing in the file.

**What survives on its own logic:** the events-per-sample ratio does not separate the single-channel
from the macroscopic regime (the channel count does), which killed one candidate thread; under
misspecification the score vanishes at the pseudo-true value and not at the simulator's parameters,
which no draft stated; Godambe 1960 is optimum estimating functions and not the score identity; a
channel count is dimensionless, so "three quantities, two of them dimensionless" is wrong; Hodgkin and
Huxley fitted by hand; and IonBench "does not cover optimisation approaches for stochastic models",
which is narrower than "excludes stochastic models". All six were applied to the manuscript the same
day.
