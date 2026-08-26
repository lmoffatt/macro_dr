# Shortening plan for the revised manuscript

> **RECONCILIATION, 2026-08-26, before execution.** The plan was written against the PDF of
> 2026-08-25 23:12 and does not see the commits of 2026-08-26 (`5491400` figure citation order,
> `fdb2b9a` caption overflow). Four rulings, decided with the author, override the items below:
>
> 1. **Locate every passage by its TEXT, never by the line numbers.** The extraction is stale by two
>    commits and by every edit this pass lands.
> 2. **Ruling (m) is now settled the other way, and the Results item is CANCELLED.** The plan rejects
>    cutting the Figure 2 caption's magnitude/anisotropy decode because the Results item cuts the
>    prose gloss naming the caption as survivor. The caption decode was already cut on 2026-08-26 to
>    bring the caption inside its page, which is physical and cannot be reverted. So the PROSE gloss
>    (Results, "its magnitude, the typical factor by which...") is now the sole survivor and must
>    NOT be cut. Item "COMPRESS ~47 words, lines 586-594" is void; ~47w come off the total.
> 3. **Figure 5—figure supplement 3 is DROPPED (author sign-off given).** Beyond the plan's own
>    reasoning, its page carries the document's worst overfull box (121.9 pt, recorded in
>    `changes.md`), so dropping it also removes that defect. The fallback is not taken.
> 4. **The abstract appositive is KEPT; that item is DECLINED.** "the smallest that isolates the
>    approximation's own error" entered on 2026-08-25 because five independent readers raised the
>    two-state objection; eight words do not pay for reopening it.
>
> Also stale, and harmless: the Figure 3 caption sentence the plan names as the surviving statement
> of N_ch/i inseparability was cut on 2026-08-26 (survivors: that figure's supplement 2 parenthesis
> and the Results identifiability sentence); "and the difference is not cosmetic" in Figure 4's
> caption is already gone, so that item counts zero.
>
> **Word counts:** the plan measures with `pdftotext` (main text ~24,300); `papers/_program/wordcount.py`
> measures ~13,100 by a different definition. The two arithmetics are not comparable. Pages are the
> comparable unit (69 -> ~55); re-measure with `wordcount.py` against `08_length_plan.md` after the pass.

Target document: `elife_paper.pdf`, revision of 2026-08-25 23:12 (69 pages, about 38,600 words, main text about 24,300 words).

Line numbers below refer to the plain-text extraction produced by `pdftotext -layout elife_paper.pdf elife_paper.txt` on that revision; page numbers refer to the PDF itself and are the stable reference if the text is re-extracted after edits.

## Summary

The plan holds 117 items, each proposed by a section editor and then checked against the whole text for orphaned references, double counting, and cases where two sections each cut their own copy of a shared passage. After that check, 91 items execute as proposed, 24 execute with the stated adjustment, and 2 are rejected and must not be executed. The verified saving is about 7,750 words, which is 10 to 11 pages of running text plus two to three figure and table pages. Expected outcome: 69 pages down to roughly 55, main text from about 24,300 words to about 17,000.

Savings by mechanism, largest first:

1. Methods material moved to Supplementary File 1 or the repository documentation, about 3,100 words. Methods keeps model, simulation, grid, noise conversion, estimation, anchors, bootstrap, and the seed caveat.
2. Removal of second and third statements of ideas the paper repeats, about 1,900 words across Introduction, the family section, Results, and Discussion. The canonical copy of each repeated passage is fixed in the rulings below.
3. Caption and supplement-legend compression, about 1,250 words plus two whole pages (one supplement dropped, two merged). Rule: the caption decodes the figure once, the Results text reads it once.
4. Appendix pedagogy, about 1,300 words, lowest priority since eLife appendices carry less reader cost.

## Ground rules

Do not cut, in any item: the two-closures box, Table 1, the three diagnostics definitions with Eqs. 3-4, the m and a summaries, the main figures, the numerical-safeguard disclosures (Eq. 7, the smooth trust coefficient with its constants, the count-floor fallback sentence), and the limitations subsection. Where an item touches one of these, the residue requirement in the item says what must survive.

Suggested execution order: Methods moves first (largest, least entangled), then captions and legends, then Discussion, then Results, then the family section and Introduction (their cuts depend on rulings about which copy of a passage survives), then appendices. Recompile and re-extract after each section pass; later line numbers shift as earlier cuts land.

Two kinds of decision stay with the author. Dropping Figure 5, figure supplement 3 outright is marked SIGN-OFF on the item itself, with a stated fallback (compress its legend instead). And every Methods MOVE names a destination (Supplementary File 1 or repository documentation); those destinations are proposals, so confirm each before executing the move.

## Canonical-copy rulings

These rulings resolve every passage that appears more than once. When an item below conflicts with a ruling, the ruling wins.

- (a) natural-units/tilde definition: the body p.5 copy (lines 236-243) survives IN FULL, interpretive tail included (family editor's trim of the tail is reversed, since the Methods copy sheds its re-derivation); the Results p.11 copy (576-579) becomes a pointer; the Methods copy (2083-2102) keeps only the Current_Noise-to-S-tilde engine conversion, with the S-tilde=1 gloss surviving once at 2113-2114. The scalar-vs-boundary-pair accent disambiguation survives ONLY in the Appendix 2 tilde box (2894-2896) - not relocated to p.5, and the Methods copy of it (2100-2102) is cut.
- (b) 'the average is not a function of the occupancy at any one time... on the K^2 pairs one does': body p.7 (348-352, compressed by the family item) and Appendix 1 'The three steps that build the interval update' (2456-2467) survive; the Appendix 1 preview copy (2302-2312) is the cut. Verified no text between 2312 and 2456 references 'the three steps'; the 'step two' uses at 2517 and 2567 come after the surviving subsection.
- (c) tau_int gloss 'counts how many intervals a fit effectively has for every interval it believes it has': the protected p.10 definition (497-499) is the SOLE surviving full gloss. Cut from the Introduction (124-126), the family echo (552-556), Results p.18 (1037-1038), the Figure 5 caption (which keeps only 'integrated autocorrelation of its standardized residual... reads 1 for a white residual' as axis decode), and the Fig 5-supp 3 legend (leaves with the moved supplement). The intro editor's claim that the Figure 5 caption copy survives is superseded.
- (d) Milescu local-time-correlation diagnosis: p.10 (541-543, inside the protected diagnostics definitions) and Results p.14 (861-863, with the measured -0.004/0.191/0.864 numbers, verified uncut) survive; the Discussion copy (1443-1451) is the cut, compressed to keep 'Conditioning on both ends leaves a white score' as the Mehra sentence's antecedent.
- (e) factor-2217 white-residual example: Results p.18 (1061-1064, verified uncut) survives as the full measured statement; the Fig 5-supp 1 legend keeps it as a census datum (3839-3840); the Discussion copy (1467-1469) is cut back to the principle.
- (f) 'Chaining intervals leaves a single state variable per junction': body p.7 (326-328, verified untouched by the family compress at 344-358) survives; the Appendix 1 copy (2471-2472) is the cut, with the static-condensation analogy retained there.
- (g) simulator-shares-specification: the protected limitations paragraph (1591-1597) is canonical; Methods (1681-1691) keeps the bare facts with its 'symbol by symbol below' pointer re-aimed at Supplementary File 1; the Discussion closing echo (1621-1622) is deleted; the family section's clause at 175-177 survives as one clause.
- (h) grouping/record-reduction result (both editors cut their copy): Results 1213-1221 SURVIVES UNTOUCHED - the results editor's compress is rejected; the Discussion copy (1475-1492) is the cut, compressed to the design implication plus the three details that exist only there (geometric growth schedule; proportionately smaller instrumental variance with no sample discarded; cost linear in the groups) and the '(Results)' pointer.
- (i) washout/conditioning mechanism (both editors cut their copy): the mechanism sentence and the 15-16%-vs-0.5-0.8% numbers SURVIVE in Results 866-876 (the results compress shrinks to its framing sentence only); the Discussion copy (1470-1474) is the cut, keeping only the callback and the unique design advice 'spend its samples on the jump and its rise'.
- (j) per-member per-sample/correlation sweep numbers (Results 1014-1019 vs Fig 4-supp 1 legend 'Third' reading): the LEGEND copy survives and becomes the numbers' home - the body paragraph MOVEs into it and the legend's 'which is the paragraph above' is deleted; the captions editor must NOT cut the Third reading (that item is reduced accordingly).
- (k) Figure 5 caption vs Results readings: Results text 1034-1074 survives as the readings' single home; the caption is trimmed only via the captions editor's three finer items (the results editor's whole-caption compress is superseded and counted at zero).
- (l) 'ribbons showing where each moves as the interval is swept...' fragment (1427-1428): one deletion executed once - counted under the results item, the identical captions item counted at zero; the Figure 7 caption copy (1381-1382) is the survivor.
- (m) m/a meaning gloss (three statements): the protected p.10/Eq. 4 definitions AND the Figure 2 caption decode (668-670) survive; the Results-text gloss (588-592) is the cut. The captions item cutting the Figure 2 caption decode is rejected - it named the Results-text copy as a survivor while the results editor cut that copy naming the caption as survivor.

## Items by section

Tiers: DELETE removes the passage outright; COMPRESS shrinks it and the item states what the residue must retain; MOVE relocates it and the item states the destination.

A section editor's note at the end of a section records that editor's original reasoning. Where a note disagrees with a canonical-copy ruling or with an item's binding adjustment, the ruling and the adjustment win; the notes were written before the cross-section check.

### Materials and methods, Data availability

Span: pp. 28-38, lines 1627-2204. Section size about 7176 words. Items: 26.

- [ ] **MOVE**, about 24 words, lines 1635-1638, p.28
  - Passage: "read from the production run scripts rather than from the parameter tables in the repository (the tables hold stale values the figures never used)"
  - Why: A deposit-consistency warning for repository users, not paper content; the sentence keeps its parameter values and drops the provenance aside.
  - Do: Repository documentation (README data-provenance note).
  - Cross-references: NONE; nothing else in the paper mentions the stale tables.

- [ ] **DELETE**, about 27 words, lines 1668-1669, p.29
  - Passage: "That window is the acquisition interval the whole study is about, and its length is set by the number of raw samples entering each recorded value."
  - Why: Both facts are stated in the sentences immediately before (mean of n_samp,t consecutive raw samples, fixed at acquisition) and after (uniform averaging over the window).
  - Cross-references: NONE.

- [ ] **MOVE**, about 36 words, lines 1679-1681, p.29
  - Passage: "The production scripts do pass a number_of_substeps argument, and on this branch it has no effect..."
  - Why: A defensive note for someone reading the scripts, which is exactly what repository documentation is for.
  - Do: Repository documentation.
  - Cross-references: NONE; the Figure 1 substep-sampler sentence (lines 1695-1697) is about a different run and stays.

- [ ] **COMPRESS**, about 100 words, lines 1681-1691, p.29
  - Passage: "This is what licenses treating the simulation as the reference ... stated with the others in the Discussion." [duplication (g)]
  - Why: Duplication (g): the Discussion copy (lines 1592-1594 and 1621, in the protected limitations material) survives; the Methods copy sheds the justificatory prose ("licenses", "by design and not by convenience") and keeps only the facts.
  - Do: Residue (~55 words) must retain: inference is likelihood-only (no prior, no evidence computation); likelihood and simulator are independent implementations sharing only the upstream specification (rate matrix, conductance vector, protocol, sampling times); the cost of that sharing is treated with the limitations in the Discussion.
  - Cross-references: The clause "writes that specification out symbol by symbol below" points at the symbol table being moved (lines 1821-1846); the residue's pointer must aim at Supplementary File 1 instead.

- [ ] **MOVE**, about 384 words, lines 1704-1735, pp.29-30
  - Passage: Dispatch-flag prose: "The members are the two axes ... the I prefix records the window setting" plus the grid-cell inventory "Both least-squares arms are dispatched ... below the range this paper reports" and the macroscopic flag definitions
  - Why: Flag values, data keys, the 9e-16 mean-agreement check, and superseded-cell bookkeeping are code-mapping material the body family section (pp.3-10, protected) already covers conceptually.
  - Do: Supplementary File 1, as the prose accompanying the moved Table 2; body retains ~70 words: least squares = family approximation 2, macroscopic = 0, the three flags (recursive, endpoint count av in {0,1,2}, variance form), the shared-mean fact, and a pointer to the supplementary table.
  - Cross-references: Residue must keep the fact that each least-squares arm shares its mean model exactly with the macroscopic member of the same window flag (LSE-NR, ILSE-INR) and differs only in predictive variance, which Results rely on when comparing arms; body line 483 and Appendix 2 lines 2986-2987 point readers to "Table 2" for these flags and must be re-pointed.
  - **Adjustment (binding)**: Residue must keep one sentence stating the LSE-NR / ILSE-INR shared-mean identity WITH the 9e-16 verification number: the family editor deletes the body's 10^-15 demonstration naming line 1715 as the survivor. Savings 384 to ~364w.

- [ ] **MOVE**, about 438 words, lines 1737-1765, p.30
  - Passage: Table 2 ("Member / Data key / Recursive / av / Variance / Conditioning") with its full caption "The members as flags. ..."
  - Why: The span guidance's heaviest MOVE candidate: an implementation dispatch table whose conceptual content (member definitions, ladder, INR = MacroINR identification) already lives in the protected body family section, Table 1, and Appendix 2.
  - Do: Supplementary File 1, as a supplementary table (renumber; re-point the two citations).
  - Cross-references: Cited at body lines 481-483 ("the dispatched flags and data keys are in Table 2") and Appendix 2 lines 2986-2987 - both citations must be re-pointed to the supplementary table; the caption's INR = MacroINR identification (also stated at lines 46, 481, 2986) and the MR/IR same-total-variance note (cross-referenced to Appendix 2) must survive in the moved caption.

- [ ] **COMPRESS**, about 135 words, lines 1769-1784, pp.30-31
  - Passage: "Table 2 lists the members. Reading the endpoint flag as a ladder ... the names the implementation carries."
  - Why: The ladder reading duplicates both the Table 2 caption and the body family section; the paragraph's unique content is three citations' worth of anchoring for IR.
  - Do: Residue (~60 words) must retain: IR conditions on the boundary state (the pair of channel states at the two window ends) and marginalises the interior analytically, specified by Appendix 1; first applied to macroscopic recordings of the P2X2 receptor (Moffatt and Pierdominici-Sottile, 2025); coincides with a time-integrated Kalman filter (Discussion); every member scores with Eq. 5 and the regularizer of Eq. 7.
  - Cross-references: Line 1783's "the regularizer of Eq. 7 below" must survive in the residue (Eq. 7 stays in the body); the Kalman correspondence is detailed in the Discussion (line 1507) and the P2X2 priority in the Introduction, so pointers suffice; the symbol-correspondence clause re-points to Supplementary File 1.

- [ ] **COMPRESS**, about 29 words, lines 1816-1819, p.31
  - Passage: "Eq. 6 is also the one moment at which Sigma genuinely is a single channel's covariance ... comes from that distinction."
  - Why: Three sentences of conceptual aside compress to one without losing the Cov(N)/N_ch distinction, which Appendix 1 derives anyway.
  - Do: Residue (one sentence) must retain: from the first update onward, conditioning on a shared recorded value correlates the channels, so the recursion carries Cov(N)/N_ch, no single channel's covariance.
  - Cross-references: NONE.

- [ ] **MOVE**, about 208 words, lines 1821-1846, p.31
  - Passage: "Every symbol against the name it carries in the implementation." - the full symbol-to-implementation table (P, gmean_ij, gtotal_ij, gsqr_ij, SmD, gSg, ms, gS, y_var) and its intro
  - Why: A code-mapping table for reimplementers; the equations of Appendix 1 stand alone in symbols, so the body needs only a pointer.
  - Do: Supplementary File 1 (keep it archival rather than repo-only, since it is presented as part of the specification); body retains one sentence: the one-to-one correspondence between Appendix 1 symbols and implementation names is given in Supplementary File 1.
  - Cross-references: Referenced by the ladder paragraph (lines 1783-1784) and by the "symbol by symbol" clause in lines 1689-1690 - both residues must re-point to Supplementary File 1; no implementation name in this table is used anywhere else in the paper (checked).

- [ ] **MOVE**, about 55 words, lines 1855-1860, p.32
  - Passage: "Throughout the reported runs the remaining build flags were held fixed: taylor variance correction off ... MacroMR and not MacroMRT."
  - Why: Build-flag inventory named by the span guidance as a MOVE; the V-in-VR versus Taylor-flag disambiguation and the MRT/IRT exclusion matter only to someone matching deposited filenames.
  - Do: Supplementary File 1 (build-flag inventory); body retains ~15 words: remaining build flags were held fixed in every reported run, settings in Supplementary File 1.
  - Cross-references: NONE outside the span (no MacroMRT/MRT/IRT reference elsewhere in the paper, checked); the supplement must state the MacroMR-not-MacroMRT identification for readers of the deposit.

- [ ] **COMPRESS**, about 680 words, lines 1861-1921, pp.32-33
  - Passage: The three numerical safeguards: "The filter carries three numerical safeguards ... a minimum occupancy probability of 10^-12."
  - Why: Per guidance: keep the disclosures and the formulas, move the design-rationale prose (why epsilon is keyed to the condition number, the truncate-then-renormalise history, the epsilon/2 bias arithmetic, covariance-step and untied-trusts detail) to Supplementary File 1.
  - Do: Residue (~200 words) must retain: (i) one framing sentence that the filter carries three safeguards and that without the latter two the recursive members return no value at small channel counts; (ii) safeguard 1: the boundary-conditioned moments are quotients by P_ij(Delta), identically 0/0 over the zero-agonist segments, regularized by the conjugate pseudo-count of Eq. 7 (display kept) with epsilon = machine epsilon times the Frobenius condition number of the eigenvector matrix and endpoint-average prior means; (iii) safeguard 2: the mean step is damped onto the simplex by the smooth coefficient alpha = (1/2)(1 + c*alpha_p - sqrt((1 - c*alpha_p)^2 + eps^2)) with c = 0.9, eps = 1e-4, smooth because the log-likelihood is differentiated in theta (a hard minimum reproduces the mean current and corrupts every quantity the paper measures); (iv) safeguard 3 verbatim: binomial-count floor 5 with fallback to the non-recursive form, occupancy floor 1e-12. Excised rationale goes to Supplementary File 1.
  - Cross-references: Eq. 7 is cross-referenced at line 1783, so Eq. 7 must stay in the body; "the trust coefficient above" at lines 2147-2148 requires the smooth-alpha formula to remain; the safeguard-3 fallback disclosure is on the protected list and is retained verbatim.
  - **Adjustment (binding)**: Protected disclosures: residue must keep Eq. 7 with the epsilon definition, the smooth-alpha formula with its constants (required by 'the trust coefficient above' at 2147-2148), safeguard 3's fallback disclosure verbatim, plus one sentence 'damping the step is what replaced correcting the state after it' and the reimplementation warning (the Discussion compress at 1515-1518 assumes Methods keeps the contrasts). Savings 680 to ~640w.

- [ ] **COMPRESS**, about 55 words, lines 1952-1961, p.33
  - Passage: Six-parameter configuration: "Both arms of it are produced by ops/local/figure_3_time.macroir ... The script's own header records why ... misalign them without any visible error."
  - Why: Script name, "the script's own header records why", and the parameter-index-alignment rationale are repository material; the statistical justification (Fisher accumulated, never inverted) and the two visible flat directions are load-bearing and stay.
  - Do: Residue must retain: the six-parameter configuration feeds Figure 3; both arms are produced by the same script that simulates the ensemble and scores the macroscopic members, so all eight members score one common set of recordings; all six parameters free in log10 coordinates; the time-resolved figures accumulate per-interval Fisher information and never invert it, so the flat directions are harmless and two are visible (Current_Noise zero score/Fisher row; unitary current and channel number rank-one degenerate along N_ch*i).
  - Cross-references: Figure 3 comparisons rely on "all eight members of that figure score one common set of recordings" - the residue keeps that fact; the moved index-alignment note goes to repository documentation.

- [ ] **MOVE**, about 60 words, lines 1969-1971 and 1982-1983, p.34
  - Passage: Four-parameter configuration script names "is ops/local/figure_4_LSE.macroir ... figure_3_mle_LSE.macroir" and the file index dictionary "Its parameter indices are 0 = kon ... index 2 is the baseline."
  - Why: Script filenames and output-file column indices document the deposit, not the method; the statistical content of the configuration (which parameters are fixed, why fixing is forced, fixed at truth so least squares is helped) stays in full at lines 1973-1981.
  - Do: Repository documentation (script inventory and output-file parameter-index dictionary); body retains ~15 words: the four-parameter configuration produces every grid cell and, where needed, the per-group estimate cloud.
  - Cross-references: NONE in the body; the index dictionary must land in repository documentation for anyone parsing the deposited LSE files, and the retained sentence "This is the configuration behind every ILSE cell of Figures 4 and 6" (line 1981) must survive.

- [ ] **COMPRESS**, about 95 words, lines 1997-2005, p.34
  - Passage: Optimiser convergence machinery: "Its damping schedule and its three convergence limits are in Supplementary File 1. What stops most fits ... The Newton-decrement constants are compiled in."
  - Why: The criterion hierarchy keeps its numbers but sheds the two rationale clauses (noise-floor chasing; relative-gradient-norm anti-inflation) and the compiled-in aside, which join the schedule already promised to Supplementary File 1.
  - Do: Residue (~55 words) must retain: Gauss-Newton with Levenberg damping, warm-started at the simulation truth; primary stop is the Newton decrement (1/2)g'H^-1 g below 1e-8 for the first ten iterations and 1e-4 after; the three script-set limits act as fallbacks; schedules and limits in Supplementary File 1.
  - Cross-references: Lines 2009-2013 ("decided by the criteria above") must still resolve - the residue keeps the Newton-decrement criterion as primary with the script-set limits as fallbacks.

- [ ] **COMPRESS**, about 43 words, lines 2028-2032, p.35
  - Passage: "Faithfulness of each approximate likelihood to the simulator is measured with likelihood-level diagnostics: standardized residuals, the score (the gradient ...) ... here we describe only how they were computed."
  - Why: Re-defines the three diagnostics that the protected definitions section (pp.3-10) already defines; a naming sentence suffices.
  - Do: Residue (~30 words) must retain: the three diagnostics as named items (standardized residuals; the score; the information-distortion matrix comparing the covariance J of the score across replicates with the model's Fisher information), definitions and thresholds as given where introduced, this section giving only the computation.
  - Cross-references: The symbol J is used at line 2059, so the residue must keep the phrase "the covariance J of the score across replicates".

- [ ] **COMPRESS**, about 55 words, lines 2043-2049, p.35
  - Passage: "A second Fisher construction exists in the codebase, a central difference of the analytic score ... An earlier battery anchored on it (data directory 433ed13) is retained in the deposit ..."
  - Why: The disclosure that a finite-difference cross-check exists and that no reported number uses it stays; the deposit-archaeology sentence about directory 433ed13 belongs with the data documentation.
  - Do: Residue (~45 words) must retain: a second, central-difference Fisher construction exists; it is widely indefinite replicate by replicate (one reason the body anchors on the analytic construction); averaged over recordings it is the check the analytic anchor is measured against, over the full roster, in Supplementary File 1; no number reported in the paper is taken from it.
  - Cross-references: NONE in the body; directory 433ed13 must be recorded in repository documentation so the deposit stays interpretable.

- [ ] **COMPRESS**, about 83 words, lines 2062-2072, p.35
  - Passage: Bootstrap paragraph internals: "Despite the internal column names, no probit transform is applied ... verified equal to the R quantile of type 6" and "The program does run a group bootstrap of it ... held in memory and never written ..."
  - Why: The protocol (nonparametric bootstrap over groups, the five quantiles, 100 replicates, lag 10, 200 for the capstone) is protected and stays; the column-name caveat and the in-memory-never-written archaeology are repository documentation.
  - Do: Residue must retain: full bootstrap protocol with all quoted numbers, plus one sentence: the per-group estimate cloud is deposited raw and every interval quoted from it is computed in the figure code from the raw estimates.
  - Cross-references: Supplements that quote intervals from the estimate cloud rely on the retained sentence that the cloud is deposited raw and intervals are computed downstream in the figure code; the probit/type-6 note moves to repository documentation.

- [ ] **COMPRESS**, about 150 words, lines 2083-2102, p.36
  - Passage: Noise-axis conversion and tilde re-definition: "The noise axis needs its conversion stated once ... the two uses never meet in one expression." [duplication (a)]
  - Why: Duplication (a): the surviving natural-units/tilde definition is the body p.5 copy (lines 236-245, protected section); the Methods copy keeps only the engine-to-units mapping and drops the re-derivation, the "only dimensionless combination" justification, and the S-tilde = 1 gloss (which survives once at lines 2113-2114).
  - Do: Residue (~90 words) must retain: the engine parameter Current_Noise is a white-noise power spectral density in units of g^2 s, reaching the likelihood only as e_t = Current_Noise/Delta_t with no correlation time; the sweep holds the density fixed so per-point variance rises as 1/Delta_t; in the natural units defined earlier (Delta-tilde = Delta/tau, S-tilde = S*koff/i^2), the reference cell is S-tilde = 0.01.
  - Cross-references: The accent-disambiguation sentence (tilde as boundary-pair contraction, Appendix 1 Eqs. A18-A19) should be relocated to the surviving definition at body p.5 - coherence pass must check p.5 does not already carry it; the Results p.11 copy of the definition is another agent's span.
  - **Adjustment (binding)**: Keep only the engine-to-units mapping as proposed, but do NOT relocate the accent-disambiguation to p.5: the Appendix 2 tilde-box copy (2894-2896) is retained instead (item 95 adjusted accordingly). S-tilde=1 gloss survives once at 2113-2114. Savings 150w.

- [ ] **DELETE**, about 25 words, lines 2118-2119, p.36
  - Passage: "The noise range is deliberately not rectangular: the highest levels were run only at the larger channel numbers, where the crossover the map measures sits."
  - Why: Verbatim-degree duplicate, within one paragraph, of line 2112 ("denser below one, where the crossover the map measures sits") and lines 2115-2117 ("Coverage is deliberately not a rectangle and it differs by member").
  - Cross-references: The unique fact (highest noise levels only at larger channel numbers) is preserved by the per-member cell manifest in Supplementary File 1 (cited at line 2117).

- [ ] **MOVE**, about 90 words, lines 2103-2108, p.36
  - Passage: "The results reported here are read from three data directories, each named by the git commit hash ... 1c2ae6f, 87889e6 and 0ffbda7, searched in that order ... by matching nsim_10000 in the filename ..."
  - Why: Directory hashes, search order, and filename matching are deposit navigation, not method.
  - Do: Repository documentation; body retains ~20 words: every cell entering a body figure is at n_sim = 10^4, enforced by the figure code.
  - Cross-references: NONE in the body; the hashes and search order must appear in repository documentation since Data availability gestures at per-file provenance.

- [ ] **MOVE**, about 34 words, lines 2119-2121, p.36
  - Passage: "The figure code auto-detects which cells exist on disk at render time and prints the detected grid into the render log ..."
  - Why: Render-log mechanics document the pipeline, not the science.
  - Do: Repository documentation.
  - Cross-references: NONE.

- [ ] **COMPRESS**, about 18 words, lines 2122-2124, p.36
  - Passage: "Each grid cell is reduced through five stages ... each writing its own family of comma-separated-value files; the stages and their outputs are listed in Supplementary File 1."
  - Why: The five-stage pipeline named by the guidance is already delegated to Supplementary File 1; the CSV-family clause is the only part left to shed.
  - Do: Residue (~20 words) must retain: the reduction stages from per-group estimate cloud to the empirical-against-theoretical distortion capstone, and their outputs, are listed in Supplementary File 1.
  - Cross-references: NONE; Supplementary File 1 already holds the stage list.

- [ ] **COMPRESS**, about 50 words, lines 2127-2135, pp.36-37
  - Passage: "The engine is macro_dr ... Provenance is self-documenting at three levels: on every invocation the binary snapshots the assembled script ... documented per file."
  - Why: Keep the engine identification and the commit-hash stamping that the rest of the paper leans on; the run-ledger snapshot and "multi-commit provenance per file" detail is repository documentation.
  - Do: Residue (~45 words) must retain: the engine is macro_dr, C++20, each analysis a script in a small typed domain-specific language (Moffatt, 2025); every output stamps the build's short git commit hash and data directories are named by it, so results from different code versions cannot overwrite one another.
  - Cross-references: NONE; Data availability (line 2166) repeats that each output records the engine commit hash.

- [ ] **COMPRESS**, about 115 words, lines 2136-2149, p.37
  - Passage: Automatic-differentiation paragraph: "The engine that produced these results is not the build behind the P2X2 analysis ... rather than carried over from that work."
  - Why: The clauses restating that the hand-wired Fisher is positive semidefinite and needs no second pass duplicate lines 2040-2042 of the Diagnostics subsection (which survive), and the MCMC history compresses to one clause.
  - Do: Residue (~100 words) must retain: not the build behind the P2X2 analysis; every score comes from automatic differentiation through the recursion itself, value and derivative in one templated C++ type, the spectral reconstruction differentiated in place with its divided differences and coincidence limits (Appendix 1); the Fisher information is wired by hand as the Gaussian form of Eq. 3 over the predictive moments and first derivatives; the earlier study's MCMC fits needed no gradient, so the score-based diagnostics are new here, and the dependence is what forces the trust coefficient to be infinitely differentiable.
  - Cross-references: The residue must keep the AD-forces-smooth-trust-coefficient link, which the safeguard-2 residue (lines 1894-1918 item) depends on; "positive semidefinite by construction" remains stated at line 2042.
  - **Adjustment (binding)**: Residue must keep the one-templated-type / derivative-through-the-recursion sentence (the intro drops its copy naming this as survivor) and the AD-forces-smooth-trust-coefficient link (safeguard-2 residue depends on it). Savings ~115w stand.

- [ ] **COMPRESS**, about 60 words, lines 2153-2158, p.37
  - Passage: "The one reproducibility limit worth stating twice is the simulation seed. Because the seed was the random sentinel ... is stated in the Data and code availability section."
  - Why: The seed caveat is protected but is stated in full twice within Methods; the first copy (lines 1698-1702, p.29) survives, and this second copy keeps only its one new fact.
  - Do: Residue (~22 words) must retain: all resampling that operates on a fixed dataset (the bootstraps, the group partitioning) is deterministic and reproducible given the data.
  - Cross-references: Guidance protects the seed caveat - it remains fully stated at lines 1698-1702 and briefly again in Data availability (lines 2167-2170); the resampling-determinism fact appears nowhere else, so the residue must keep it.

- [ ] **COMPRESS**, about 35 words, lines 2171-2177, p.37
  - Passage: macroir library paragraph: "A separate library, macroir, carries ... by tools/cross_language_check.py, which is the evidence for the statement that the three interfaces return the same numbers."
  - Why: The library's existence and scope stay in Data availability; the verification-tool filename and the bindings-not-reimplementations elaboration are repository documentation.
  - Do: Residue (~65 words) must retain: the macroir library carries the boundary-conditioned likelihood, its score, its Gaussian Fisher information and the distortion diagnostic in a self-contained C++ core with R and Python packages bound over the one core, verified against it on a reference cell.
  - Cross-references: NONE in the body; the tool name tools/cross_language_check.py must move to repository documentation so the same-numbers claim keeps its evidence trail.

Section editor's note: Span = Materials and methods + Data availability (lines 1627-2204, pp.28-38), 7,176 words including ~30 words of page-number artifacts. Proposed saving 3,084 words (~43% of the span), of which roughly 1,750 are MOVEs preserved in Supplementary File 1 or repository documentation rather than deletions - appropriate since this is the designated MOVE-heavy section. Duplication survivors, as instructed: (a) natural-units/tilde - the body p.5 copy (lines 236-245, protected section) survives; the Methods p.36 copy is reduced to the engine-to-units mapping (the Results p.11 copy is outside this span). (g) simulator-shares-specification - the Discussion copy (lines 1592-1594 and 1621, in/near the protected limitations) survives; the Methods copy compresses to one factual clause. Span-internal duplications found beyond the orchestrator's list: the seed caveat is stated in full twice in Methods (first copy at lines 1698-1702 survives; second copy compressed) and a third time briefly in Data availability (kept); the S-tilde = 1 gloss appears twice on p.36 fifteen lines apart (lines 2113-2114 copy survives); "coverage is deliberately not a rectangle / not rectangular ... where the crossover the map measures sits" appears twice in one paragraph (lines 2115-2117 copy survives, lines 2118-2119 deleted). Kept untouched per guidance: model + parameter values, emission model, stimulus, uniformization simulation, run sizes, seed caveat (first copy), Eq. 5 scoring and Eq. 6 initial condition, the two prior-work reductions (lines 1848-1854), Eqs. 7-10 (all remain in the body), the marginalized-scale least-squares passage (lines 1931-1951), the forced-fixing/helped-not-handicapped and 6x6-vs-4x4 scoping passages, the two-stage estimation protocol, group sizes, dropped-fit disclosure, warm-start rationale, the two anchors, the analytic Gaussian-Fisher anchor and eigenvalue-tolerance/grey-cell handling, the bootstrap protocol numbers, the grid sweep and noise ladder, SLURM note, and all end-matter boilerplate. Biggest coherence dependency for the synthesis pass: moving Table 2 requires re-pointing its two citations (body lines 481-483 and Appendix 2 lines 2986-2987) and preserving the INR = MacroINR identification in the moved caption.

### Figure-supplement legends and caption audit

Span: pp. 55-69, lines 3139-4026, plus main captions. Section size about 4296 words. Items: 19.

- [ ] **COMPRESS**, about 75 words, lines 3284-3296 (p.57)
  - Passage: Figure 2—figure supplement 1 legend, 'The reported covariance against the empirical one...'
  - Why: 209-word legend; the method-context sentences (why centred on the cloud mean, why this cell has 1000 fits and others 100, chi-squared validity note) restate rationale a supplement legend does not need.
  - Do: Residue (~134 words) retains: what-is-plotted (squared Mahalanobis distance from the cloud mean vs chi-squared-6 quantiles, under the reported Fisher and the sandwich; identity = calibrated, above = under-covers), the coverage numbers 0.914 vs 0.947-0.949, one clause that centring on the cloud mean tests the second moment alone, and the cell identity (IR alone, Nch=10, S=0.005, group size 100, the one cell with a thousand fits). Cut the 'every other cell being simulated at ten thousand recordings' justification and the closing chi-squared-validity sentence.
  - Cross-references: Results lines 700-702 cite this supplement for the 0.914-to-0.947 coverage lift; those numbers stay in the residue.

- [ ] **COMPRESS**, about 65 words, lines 3353-3366 (p.58)
  - Passage: Figure 3—figure supplement 1 legend, 'Per-step information against score variance, the opening rate and the closing rate...'
  - Why: 241-word legend; the six-intervals-to-floor fact is stated in Results 866-870 (which cites this supplement), the 'so the cheapest member...is better calibrated than any of the three recursive members' clause is interpretation the seven medians already show, and 0.183889 is spurious precision.
  - Do: Residue retains: layout (seven members, kon band A/B, koff band C/D; upper row Ft vs Jt magnitudes, lower row log10(Jt/Ft) with 95% bootstrap), 'curves stop at the information mask', the head-clause 'per-interval calibration is not monotone along the cost ladder' with all seven koff medians, and the least-squares ratio-row explanation compressed to one clause with median 0.184. DELETE the sentence 'In the two recursive columns the kon band reaches that floor within six intervals...' (survivor: Results 866-870).
  - Cross-references: Results 866-870 is the surviving copy of the six-intervals fact; keep it there. Fig 3—supp 3 legend references 'figure supplement 2', not this one.

- [ ] **COMPRESS**, about 35 words, lines 3425-3433 (p.59)
  - Passage: Figure 3—figure supplement 2 legend, 'The other half of the preceding supplement...'
  - Why: 132-word legend; the six-intervals fact duplicates Results 866-875, and the Nch-inseparable-from-i explanation duplicates Figure 3's main caption (lines 825-826).
  - Do: Residue: 'The other half of the preceding supplement, read the same way: the unitary current (A, B) and Nch (C, D)'; one clause 'the Nch band reaches the floor within six intervals of washout, the unitary-current band does not'; the least-squares absence reduced to a parenthesis '(no least-squares panel: Nch is not separable from i there)'. Cut 'which is the contrast the Discussion draws on how a recording is spent'.
  - Cross-references: Discussion draws on this supplement's contrast ('how a recording is spent'); the pointer direction Discussion-to-figure survives, the legend's pointer to the Discussion goes. Figure 3 caption lines 825-826 is the surviving inseparability statement.

- [ ] **COMPRESS**, about 55 words, lines 3475-3488 (p.60)
  - Passage: Figure 3—figure supplement 3 legend, 'The checks this figure has no room for...'
  - Why: 220-word legend; the noise-as-control-direction sentence (~50 words) duplicates Results 875-876, and the least-squares plug-in sentence (~40 words) duplicates Results 726 and 854-855.
  - Do: Residue retains: A/B/C what-is-plotted, 'the first residual moment separates nobody', the informative-interval counts (IR 0.067/0.050/0.050, INR 0.100/0.100/0.000, recursive members 0.409-0.975), one clause 'the noise level, measured before the agonist arrives, keeps a white score even where the kinetic directions do not', and 'Least squares carries no noise curve: its noise level is a plug-in, not a free direction.' Cut the two inference tails (memory-lives-in-gating-directions; same-fact-as-residual-pinned-at-one).
  - Cross-references: Results 875-876 explicitly cites this supplement for the white noise-score; the residue keeps that observation as a clause so the citation still lands on something.

- [ ] **COMPRESS**, about 70 words, lines 3543-3564 (p.61)
  - Passage: Figure 4—figure supplement 1 legend, 'The distortion factored...'
  - Why: 348-word legend, the longest; the 'Third' reading is a near-verbatim duplicate of Results 1014-1019 (the legend itself says 'which is the paragraph above'), and the framing sentence plus caveat can shrink.
  - Do: Residue retains: the two congruence formulas with short glosses (per-sample = single-interval departure from the assumed Gaussian; correlation = memory across intervals), 'same plane, same rows, same scale, same six members as the body figure', the full NR-vs-INR reading with its number strings and the conclusion that NR's divergence is a one-interval variance failure and not temporal (unique to this page), the INR-to-IR 'recursion buys the other factor and only the other factor' reading, and a ~40-word caveat (exact in the log-determinant, K C K-transpose reconstruction, 1.468 vs 1.481 example). DELETE 'Three readings, and the first is what this page is for.' and the whole 'Third, R and IR are equally faithful...the paragraph above.'
  - Cross-references: Results 1014-1019 is the surviving R/IR per-sample-vs-correlation statement (verify it stays); Figure 6's caption describes the same decomposition independently and is unaffected.
  - **Adjustment (binding)**: Do NOT cut the 'Third' reading: it becomes the home of the per-member numbers moved out of Results 1014-1019 (and its 'which is the paragraph above' is deleted). Cut the framing sentence and tighten the caveat only. Savings 70 to ~30w.

- [ ] **COMPRESS**, about 55 words, lines 3627-3641 (p.62)
  - Passage: Figure 4—figure supplement 2 legend, 'The magnitude of the information distortion against its anisotropy...'
  - Why: 226-word legend; the m/a meaning-glosses are the third statement of definitions given at p.10 and reprised in Results 588-592, and the affine-invariant-distance geometry duplicates the diagnostics section; the rescaling decision rule and the unique LSE/ILSE finding stay.
  - Do: Residue: roster; '(A) the magnitude m and (B) the anisotropy a of the distortion (Eq. 4), m signed about a null of one, a one-sided'; the rule 'a at one means a single effective-count rescaling is exact and m is its factor; a above one means none exists'; the coarsest-interval LSE 0.074 vs ILSE 1.00 finding with 'the mean model failing, not the variance model' and the conservative-cells census (11 of LSE's 210, 4 of R's, all at the two coarsest intervals).
  - Cross-references: m/a definitions (p.10) are protected survivors; Results 1067-1074 cites this supplement for the no-single-effective-count argument and remains the interpretive copy; the m/a numbers themselves are kept (protected m/a summaries).

- [ ] **COMPRESS**, about 95 words, lines 3706-3717 (p.63)
  - Passage: Figure 4—figure supplement 3 legend, 'The same plane for the recursive family...'
  - Why: 179-word legend; the rung-by-rung member definitions duplicate Table 1, the why-R-and-IR-are-included sentence is design rationale, and the 1.32mm-vs-1.00 cell-width comparison is layout trivia.
  - Do: Residue (~80 words): 'Figure 4 unchanged in every respect but the roster: the four rungs of the recursive family, R, MR, VR and IR (Table 1), the same two moments on the same shared factor scale with the same line families and nesting. The reading is in the text: neither partial correction improves on the member it corrects.'
  - Cross-references: Results 1157-1161 is the surviving full reading of MR/VR and cites this supplement; Table 1 (protected) carries the member definitions the legend drops. Note: this supplement is NOT droppable - it is the only plane-wide evidence for VR (absent from supp 2's roster).

- [ ] **COMPRESS**, about 100 words, lines 3763-3777 (p.64)
  - Passage: Figure 4—figure supplement 4 legend, 'What each member achieves, as opposed to what it reports...'
  - Why: 259-word legend; the two design-rationale 'because' sentences (why the nesting is transposed, why the scale is sequential) and the padded justification for drawing the corrected error are commentary, not what-is-plotted.
  - Do: Residue retains: the definition (the standard error each member delivers = square root of the diagonal of the distortion-corrected covariance, same plane and orientation as Figure 4), 'the reading is in the text', 'its level set at a factor of two draws two of the boundaries of Figure 7', bare structure statements ('blocks are per parameter, members side by side'; 'the scale is sequential: pale is a parameter pinned, dark one barely determined'), the two marks (an error bar of 15% and one of a factor two), one clause 'the corrected error is drawn rather than the reported one, which Figure 4 shows wrong by up to a factor of 74', and the least-squares-arms-on-rate-blocks note.
  - Cross-references: Load-bearing supplement: Figure 7's caption, Results 1396-1400 and the whole 'Sampling a hundred times faster' subsection (1182-1221) read off it - the definition and the level-set-draws-Figure-7 sentences must stay. The factor-74 number appears only here (and possibly Discussion); keep it.

- [ ] **COMPRESS**, about 75 words, lines 3829-3841 (p.65)
  - Passage: Figure 5—figure supplement 1 legend, 'Each member against the memory left in its own residual...'
  - Why: 208-word legend; the IR-0.98-to-1.07-vs-17.4 sentence duplicates Results 1044-1046 verbatim (which cites this supplement), and the one-directional moral duplicates Results 1058-1064.
  - Do: Residue retains: axis = each member's own-residual tau_int, 'computable only by running the member', the question the panel asks, the rank correlations (+0.92 ILSE, +0.96 INR, +0.93 MR, +0.91 VR; about +0.99 within one interval), the NR census as bare data ('NR carries 73 cells with tau_int < 1.05; an eighth are distorted, the worst by a factor of 2217'), and 'each panel is on its own scales'. DELETE 'What conditioning on both ends of the interval does...1.0 to 17.4 for the least-squares arms' (survivor Results 1044-1046) and 'The implication is one-directional...memory implies distortion, while a white residual certifies nothing' (survivor Results 1058-1064).
  - Cross-references: Duplication (e), the factor-2217 white-residual example: Results p.18 (lines 1061-1064) is the surviving measured copy; the Discussion p.25 copy is the Discussion agent's call; here it shrinks to the census datum. Results 1058-1061 cites this supplement for the mechanism claim - the rank correlations must stay as its evidence.

- [ ] **COMPRESS**, about 90 words, lines 3924-3941 (p.67)
  - Passage: Figure 5—figure supplement 2 legend, 'The instrumental noise and the channel count act only through their ratio...'
  - Why: 305-word legend; 'a reader who measures tau_int on their own record has measured r without knowing it' duplicates Results 1049-1051, and the bootstrap-width and ill-conditioned-cell discussions carry padding around numbers Results already quotes (1%, 3.5%, 13.6 at lines 1047-1055).
  - Do: Residue retains: same cells drawn twice; line construction (one channel count at one interval swept over noise, twenty-eight per panel); left axis S vs right axis r = S/Nch; four-curves-merge-to-seven; 'nothing is fitted, rescaled or normalized between the columns'; top row = the body figure's own axis and collapses too; 'bars are bootstrap quantiles, about 3% per point; the four channel counts land within 1%'; meaning of the right-panel number; 'the one large value, R's, is the single ill-conditioned cell at Nch=10^4, Delta=1 that Figure 4 greys out; without it R agrees to 23%'; 'IR's bottom-row departures are within its own error bar'.
  - Cross-references: Figure 5's main caption and Results 1054-1055 point here; Results is the surviving quantitative statement of the collapse. Keep 'nothing is fitted, rescaled or normalized' (honesty disclosure).

- [ ] **MOVE**, about 185 words, legend lines 4021-4026, extraction truncated - full page is PDF p.69; pointer at Results lines 1065-1066 (p.18); listing line at 1140-1141 (p.19)  
  **SIGN-OFF required.**
  - Passage: Figure 5—figure supplement 3 entire (figure page + legend 'Where a recording leaves memory in its own residual...'), plus its two pointers
  - Why: Cited exactly once, by a one-sentence pointer making no load-bearing claim; its content is carried elsewhere (the least-squares tau_int is Figure 5's own axis, each member's own-residual range is in supp 1, and the ratio collapse of the same statistic is supp 2's top row) - dropping it saves one full page.
  - Do: The MacroIR code/data repository (Zenodo archive) or a source-data file of Figure 5; delete Results lines 1065-1066 ('Where on the design plane a recording leaves that memory is itself a map...') and the supplement-listing line under Figure 5's caption (lines 1140-1141) with it.
  - Cross-references: Duplication (c), third copy: this legend's 'counts how many intervals a fit effectively has for every interval it believes it has' gloss and the computed-with-no-parameter sentence vanish with it (survivors: the protected p.10 definition and Results 1037-1041). If the author insists on keeping it, fallback is COMPRESS the ~140-word legend to ~75 by deleting those two duplicated sentences.
  - **Adjustment (binding)**: Count 161w here, not 185: the 24w Results pointer (1065-1066) is the results editor's separately counted deletion. Also edit the listing line at 1140-1141. Dropping a whole supplement needs the author's sign-off; the stated fallback (compress the legend by its two duplicated sentences) is sound if declined.

- [ ] **MOVE**, about 55 words, pp.58-59; legends lines 3353-3366 and 3425-3433; listing lines 827-829 (p.14)
  - Passage: Merge Figure 3—figure supplement 1 and Figure 3—figure supplement 2 into a single supplement
  - Why: Identical layout (seven members, information + log-ratio band per parameter) and identical reading instructions; supp 2's legend already calls itself 'the other half of the preceding supplement' - one figure page with four bands and one legend saves a page and the second title plus scaffolding.
  - Do: One figure page: four parameter bands (kon, koff, i, Nch) over the same seven columns; one merged legend built from the two compressed residues; one listing line under Figure 3's caption.
  - Cross-references: Results 869-870 cites 'Figure 3-figure supplements 1 and 2' - becomes one citation; Fig 3-supp 3's legend says 'the preceding supplements' and 'the ordering repeats the one in figure supplement 2' - renumber to the merged figure (supp 3 becomes supp 2). Words counted here are only the scaffolding beyond the two per-legend compressions above, to avoid double counting.
  - **Adjustment (binding)**: Renumber list to execute in full: Results 869-870 and 875-876, Discussion 1472, Methods 2052, the Figure 3 caption listing (827-830), and old supp 3's internal references ('the preceding supplements', 'the ordering repeats the one in figure supplement 2').

- [ ] **DELETE**, about 17 words, lines 1121-1122 (p.19)
  - Passage: Figure 5 caption sentence 'It is what a reader holds before choosing a method, the classical fit having been made anyway.'
  - Why: Interpretive sentence duplicated in Results: survivor is Results 1040-1041 ('a quantity an experimentalist can compute on a real recording without an ensemble and without a known truth') plus 1056 ('the one other number the reader already has').

- [ ] **COMPRESS**, about 94 words, lines 1119-1121, 1126-1132 (p.19)
  - Passage: Figure 5 caption, axis/colour/line-construction block ('tau_int, the integrated autocorrelation...which counts how many intervals...', 'Colour is the acquisition interval...because the first does not suffice...', 'Each line joins one channel count...reads that ratio off the record...')
  - Why: Three caption passages restate Results text: the tau_int counting gloss is duplication (c) (survivors: protected p.10 definition and Results 1037-1038); the colour justification with the 8.9-vs-1.6 point duplicates Results 1055-1058; the ratio mechanism duplicates Results 1047-1055 and supp 2's legend.
  - Do: Residue: 'tau_int, the integrated autocorrelation of its standardized residual (reads 1 for a white residual)'; 'Colour is the acquisition interval Delta.'; 'Each line joins one channel count at one interval swept over the instrumental noise; the four channel counts land on one another (Figure 5-figure supplement 2).'
  - Cross-references: Duplication (c) is the named triplicate - this is one of the two cut copies (the other is the supp-3 legend, item above). The residue must still let a reader parse both axes unaided; it does.
  - **Adjustment (binding)**: Residue must keep a self-sufficient axis decode ('tau_int, the integrated autocorrelation of its standardized residual in normalized form; reads 1 for a white residual'), colour = acquisition interval, and one clause of the lines-collapse construction (needed to parse 28-drawn-7-visible); cut the counting gloss (dup c), the 8.9-vs-1.6 justification and the ratio mechanism. Savings 94 to ~75w.

- [ ] **COMPRESS**, about 52 words, lines 1122-1126, 1132-1136 (p.19)
  - Passage: Figure 5 caption, bar/roster/least-squares block ('...its length is what no single correction factor can remove', 'The four members are the two axes of the family crossed, ILSE averaging...', 'The least-squares distortion matrix has four dimensions...nothing is pooled across the two.')
  - Why: The no-single-correction clause duplicates Results 1067-1074; the member glosses duplicate Table 1; the four-vs-six-dimensions sentence duplicates Results 1001-1002 (all survivors).
  - Do: Residue: bar definition kept through 'best- and worst-determined directions'; 'The four members are the two axes of the family crossed: ILSE, INR, R and IR (Table 1).'; 'Each panel reports what its own fit says about its own parameter set (Methods); nothing is pooled.'
  - Cross-references: Table 1 is protected and carries the member definitions being dropped; Results 1001-1002 must stay as the surviving dimensionality statement.
  - **Adjustment (binding)**: Keep the bar construction (point = magnitude, ends = best/worst directions), 'The four members are the two axes of the family crossed (Table 1)', and one nothing-is-pooled clause; cut the no-single-correction clause, the per-member glosses and the four-vs-six detail (Results 1001-1002 verified kept). Savings 52 to ~40w.

- [ ] **DELETE**, about 20 words, lines 1427-1428 (p.23)
  - Passage: Results fragment 'ribbons showing where each moves as the interval is swept from 0.01 to 1 in units of tau.'
  - Why: Dangling sentence fragment duplicating Figure 7's caption sentence 'Ribbons are the same boundary swept over Delta from 0.01 to 1...' - survivor: the caption (lines 1381-1382), where the how-to-read belongs; deleting also fixes the grammar.
  - Cross-references: Sits at the end of the review-demanded thinner-evidence paragraph (1422-1428): the four disclosures before it must stand untouched; the fragment discloses nothing not in the caption.
  - **Adjustment (binding)**: Identical to the results editor's deletion; execute once, savings counted there (0 here).

- [ ] **COMPRESS**, about 9 words, lines 1386-1388 (p.23)
  - Passage: Figure 7 caption closing pair 'The map is a concept map, not a phase diagram. Its boundaries are level sets of continuous diagnostics, so a looser criterion moves each of them...'
  - Why: The level-sets construction is stated with its pointer in Results 1396-1400 (survivor); the caption keeps the caveat itself.
  - Do: Residue: 'The map is a concept map, not a phase diagram: a looser criterion moves each boundary by up to a decade in noise without changing the layout.'

- [ ] **COMPRESS**, about 42 words, lines 948-950, 954-957 (p.16)
  - Passage: Figure 4 caption, least-squares clause and factor-meaning sentence ('Both least-squares arms are on the grid and each holds one column per half...', 'What a factor means differs between them and the difference is not cosmetic...')
  - Why: The four-parameter-configuration fact is stated as the scope disclosure at Results 1414-1416 (survivor) and again in supp 4's legend; the factor-meaning sentence keeps its point at half the length if the (A) example and the 'not cosmetic' flourish go.
  - Do: Residue: 'Both least-squares arms hold one column per half, their four-parameter configuration fixing i and the noise level (Methods).' and 'In (A) a factor is on the parameter; in (B) it is a ratio of variances, so 2 means the reported error bar is root-2, about 40%, too narrow.'
  - Cross-references: Results 1414-1416 is part of the scope statement demanded by review - it survives; supp 4's legend keeps its own one-line version (item above). Keep the root-2/40% example - it prevents the commonest misreading of panel B.
  - **Adjustment (binding)**: Merged with the results caption item (48): resolve this item's internal contradiction by KEEPING the root-2/40% example; cut 'and the difference is not cosmetic' once and halve the layout clause while keeping the four-parameter scope disclosure. Savings 42 to ~25w after removing the overlap.

- **REJECTED, do not execute**: COMPRESS, lines 668-670 (p.12)
  - Passage: Figure 2 caption orange-numbers gloss 'the orange numbers are the magnitude of the distortion, the typical factor by which the reported uncertainty is wrong, and, in parentheses, its anisotropy, how much that factor changes with direction'
  - Why: Third statement of what m and a mean: defined at p.10 (protected survivor) and glossed with the Eq. (4) anchor in Results 588-593 (surviving prose copy).
  - Do: Residue: 'the orange numbers are the distortion's magnitude and, in parentheses, its anisotropy (m and a of Eq. (4)); the blue number is the empirical ellipse over the corrected one. All three read 1.0 for a calibrated member.'
  - Cross-references: Protected 'm/a summaries': the numbers, their naming and the Eq. (4) pointer all stay; only the repeated meaning-gloss goes. Optional item - lowest priority in this list.
  - Reason: Circular survivorship: the results editor cuts the Results-text gloss (586-594) naming this caption decode as its survivor, and Figure 2 is the reader's first encounter with the orange numbers - cutting both would leave only the p.10 definition two pages back. The item was flagged optional/lowest priority; keep the caption copy.

Section editor's note: Span inventory: lines 3139-3256 are references (843 words, untouched); the eleven figure-supplement pages run pp.57-69 with blanks at pp.66 and 68 (verified in the PDF: 2-s1 p.57, 3-s1 p.58, 3-s2 p.59, 3-s3 p.60, 4-s1 p.61, 4-s2 p.62, 4-s3 p.63, 4-s4 p.64, 5-s1 p.65, 5-s2 p.67, 5-s3 p.69). The extraction truncates the 5-s3 legend at line 4026 mid-sentence; its full text is on PDF p.69. Legend text totals ~2425 words; the per-legend compressions (items 1-10, 12) remove ~770 of it (~32%), and dropping 5-s3 removes ~140 more. PAGES SAVED: 2 of the 13 supplement pages - merging 3-s1+3-s2 (one page) and moving 5-s3 to the repository (one page); the two blank pages vanish on reflow regardless. Legend compression alone saves words but no pages, since each legend shares a page with its full-page figure and none overflows. No other supplement is droppable: 4-s3 is the only plane-wide evidence for VR (absent from 4-s2's roster), 4-s4 feeds Figure 7's boundaries and the sampling-faster subsection, 2-s1 backs the sandwich-coverage claim (Results 700-702), 5-s1 and 5-s2 are cited with specific numbers. MAIN-TEXT CONTRIBUTION: the caption/Results cuts (items 13-19 plus the two caption blocks) total ~252 words of main text. CAPTION-AUDIT PAIRS (survivor first): (1) Results 1040-1041 > Fig 5 caption 1121-1122 (reader-holds-it sentence, DELETE); (2) p.10 definition + Results 1037-1038 > Fig 5 caption 1119-1121 tau_int gloss AND Fig 5-s3 legend gloss - this is the assigned duplication (c), both unprotected copies cut here; (3) Results 1055-1058 > Fig 5 caption 1126-1128 (colour justification); (4) Results 1047-1055 > Fig 5 caption 1128-1132 (ratio mechanism); (5) Results 1001-1002 > Fig 5 caption 1133-1136 (four-vs-six dimensions); (6) Results 1067-1074 > Fig 5 caption 'no single correction factor' clause; (7) Fig 7 caption 1381-1382 > Results 1427-1428 ribbons fragment (rare case where the caption survives - the Results copy is a dangling fragment); (8) Results 1396-1400 > Fig 7 caption level-sets clause; (9) Results 1414-1416 > Fig 4 caption 948-950 least-squares-configuration clause; (10) Results 1014-1019 > Fig 4-s1 legend 'Third' reading (DELETE); (11) Results 1044-1046 and 1058-1064 > Fig 5-s1 legend sentences (DELETE) - duplication (e)'s legend copy reduced to the bare NR census, with Results p.18 the survivor and the Discussion p.25 copy left to that span's agent. Protected items verified untouched: Table 1, the diagnostics definitions (p.10 survives as the unique definition after these cuts), the m/a summary numbers, all main figures, Figure 6's caption (already defers to text), and the thinner-evidence disclosures at 1422-1426 (only the trailing fragment 1427-1428 goes). Renumbering fallout for the coherence pass: merging 3-s1/3-s2 renumbers 3-s3; dropping 5-s3 shortens Figure 5's supplement listing; Results citations at 869-870, 1054-1055, 1065-1066 need the corresponding updates.

### Discussion

Span: pp. 24-27, lines 1431-1626. Section size about 2681 words. Items: 10.

- [ ] **COMPRESS**, about 47 words, lines 1443-1451, pp.24-25
  - Passage: "The diagnosis is Milescu's, now measured. Milescu et al. (2005) named the local time correlation ... appears in how the intervals accumulate."
  - Why: Third statement of duplication (d) plus a near-verbatim restatement of the Results finding at lines 847-849 (per-interval within tens of per cent, accumulated by an order of magnitude); the p.10 definitional copy (lines 541-543, inside the protected diagnostics definitions) and the Results p.14 measured copy (lines 861-863, with the -0.004/0.191/0.864 numbers) survive; this Discussion copy is the one cut.
  - Do: Residue (~12 words) must retain: attribution to Milescu et al. (2005), the phrase 'local time correlation', 'now measured', and a (Results) pointer, e.g. "The diagnosis is Milescu et al.'s (2005) local time correlation, now measured (Results)."
  - Cross-references: The following Mehra (1970) sentence (lines 1451-1452) depends on 'Conditioning on both ends leaves a white score', which is retained; line 1513 cites Milescu again for a different point (analytic derivatives) and must stay; Mehra appears nowhere else in the paper, so that sentence is untouchable.

- [ ] **COMPRESS**, about 18 words, lines 1454-1457, p.25
  - Passage: "When we proposed that P2X2 passes through a flip state before opening (Moffatt and Hume, 2007), the claim carried a signature ... (Jiang et al., 2012)."
  - Why: Third telling of the flip-state story, already told in the Abstract (line 15) and twice in the Introduction (lines 41, 55); the Discussion needs only the callback, not the narrative.
  - Do: Residue (~25 words) must retain both citations and the pointable-signature claim, e.g. "The flip state of P2X2 (Moffatt and Hume, 2007), later confirmed (Jiang et al., 2012), carried a signature anyone could point to in the record." The delay detail survives in the Introduction.
  - Cross-references: The next sentence's 'Which of the two orderings' keeps its antecedent (unaffected); the closing-paragraph callback 'A mechanism once left a mark anyone could point to' (line 1622) still lands on this residue. Low priority.

- [ ] **COMPRESS**, about 22 words, lines 1467-1469, p.25
  - Passage: "The reading is necessary and not sufficient: memory implies distortion, whiteness certifies nothing, and NR misreports its information by up to a factor of 2217 ..."
  - Why: Duplication (e): the 2217 white-residual example is stated with its mechanism in Results lines 1061-1064 (p.18) and again in the Figure 5-figure-supplement legend (line 3840); the Results copy survives, this Discussion copy is cut back to the principle.
  - Do: Residue: "The reading is necessary and not sufficient: memory implies distortion, whiteness certifies nothing (Results)." Must keep the necessary-not-sufficient asymmetry and the Results pointer.
  - Cross-references: Nothing after this point in the Discussion references the number 2217; the limitations subsection does not use it. Surviving copies: Results 1061-1064 and legend line 3840.

- [ ] **COMPRESS**, about 38 words, lines 1470-1474, p.25
  - Passage: "Once the agonist is removed, a member that has conditioned on the past gains nothing further ... should spend its samples on the jump and its rise." (whole paragraph)
  - Why: Duplicates the Results washout finding at lines 866-876 (which carries the six-interval floor, the 15-16% vs 0.5-0.8% double-counting, and the Figure 3-supplement pointers); the Results copy survives, and only the experimenter-facing design advice is new here.
  - Do: Residue (~35 words) must retain the design recommendation and the pointer, e.g. "Once the agonist is removed, only the closing rate and the unitary current stay informative (Results; Figure 3-figure supplements 1 and 2): an experiment should spend its samples on the jump and its rise."
  - Cross-references: 'spend its samples on the jump and its rise' appears nowhere else in the paper, so the residue must keep it; the 'this is about conditioning' attribution survives in Results 870-872.

- [ ] **COMPRESS**, about 200 words, lines 1475-1492, p.25
  - Passage: "A second result speaks to how a record should be reduced before it is fitted. The acquisition is fast and uniform ... linear in the groups it is given." (whole grouping paragraph)
  - Why: 306-word paragraph re-deriving the record-reduction result that Results lines 1213-1221 already state in full (closure of grouping under the observable, kinetic-vs-amplitude cost asymmetry, window-as-argument legitimacy, instant members describing the wrong random variable); the Results copy survives, the Discussion keeps only the design implication.
  - Do: Residue (~95-110 words) must retain four points: (1) reduction by grouping is a post-acquisition design choice, the same measurement at a larger Delta (Results pointer); (2) orders-of-magnitude reduction for a kinetic question, hardly any for an amplitude one; (3) uneven, transient-following schedules are legitimate because the likelihood takes the window as an argument; (4) fit cost is linear in the groups. E.g. "A second result speaks to how a record should be reduced before it is fitted. Grouping consecutive samples is the same measurement at a larger Delta (Results), so the temporal structure of the record becomes something to design: a record can be reduced by orders of magnitude for a kinetic question and hardly at all for an amplitude one; the groups can grow with the transient they follow, because the likelihood takes the window as an argument; and the cost of the fit follows the design, since the recursion is linear in the groups it is given."
  - Cross-references: Two details exist only here and would leave the paper: 'geometrically' (growth schedule) and 'with a proportionately smaller instrumental variance and no sample discarded'; if wanted, fold 'geometrically' into Results line 1220 at a cost of one word. Nothing else cross-references the Eq. 1 closure-under-grouping claim outside these two passages.
  - **Adjustment (binding)**: Proceeds as the ruled cut copy (Results 1213-1221 survives untouched). Residue must keep the design implication plus the three details existing only here: the geometric growth schedule, 'proportionately smaller instrumental variance and no sample discarded', cost linear in the groups, and the '(Results)' pointer at 1486. Savings 200 to ~180w.

- [ ] **DELETE**, about 40 words, lines 1496-1498, p.25
  - Passage: "It is a concept map, not a phase diagram, its boundaries level sets that a looser criterion moves by up to a decade without changing the layout, and its classical arm is one least-squares fit rather than a fluctuation analysis."
  - Why: Third statement of the Figure 7 caveats: the concept-map/level-set/decade caveat is in the Figure 7 caption (lines 1387-1388) and in Results (1396-1400, 1422-1428), and the ILSE-not-NSFA scope caveat is in Results 1414-1418; those copies survive.
  - Cross-references: The preceding headline sentence (1493-1496, the one-to-three-decades ordering) is retained as the Discussion's one statement of the map result; the not-a-fluctuation-analysis scope also survives in the Discussion NSFA paragraph (1551-1562), so no caveat is lost from the section.

- [ ] **MOVE**, about 30 words, lines 1498-1506, pp.25-26
  - Passage: "The software carries the calibrated member with its score, its Fisher information and the exact simulator ... reproducible only from the engine."
  - Why: Near-duplicate of the Data and code availability section, which already states that the diagnostic returns the distortion matrix, eigenvalue spectrum and first-order bias from one call in either language (line 2174); software inventory belongs there, not in the Discussion.
  - Do: Move the 'cheaper members are not in it, so the comparison between rungs is reproducible only from the engine' disclosure verbatim into Data and code availability (after line 2177); leave in the Discussion one clause: "a reader can run the two identities on their own scheme, at their own open probability, from one call in R or Python (Data and code availability)." Net saving ~30 words after the residue and the moved sentence are counted.
  - Cross-references: The cheaper-members sentence is a reproducibility disclosure and must land in Data availability, not vanish; the 'own scheme, own open probability' clause is the reader-facing point and stays in the Discussion residue. Coherence pass should confirm Data availability's macroir paragraph (lines 2171-2177) absorbs the disclosure without contradiction (macroir carries IR only; the engine carries all rungs).

- [ ] **COMPRESS**, about 18 words, lines 1515-1518, p.26
  - Passage: "the occupancy is held on the probability simplex by damping the update rather than by correcting the state after it, and that damping is smooth rather than a hard minimum, both because an operation that is not differentiable where it binds pumps variance into the score (Methods)"
  - Why: Second statement of the simplex-damping safeguard whose protected disclosure lives in Methods; the Discussion needs only the differentiability rationale that justifies listing it among the paper's specific contributions.
  - Do: Residue (~29 words): "the occupancy is held on the probability simplex by a smooth damping of the update, because an operation not differentiable where it binds pumps variance into the score (Methods)". Must keep the (Methods) pointer.
  - Cross-references: The damping-vs-post-hoc-correction and smooth-vs-hard-minimum contrasts survive in Methods (protected fallback/safeguard disclosure); this compress does not touch that copy. Low priority.
  - **Adjustment (binding)**: Begin the cut AFTER '...and reaches the filter itself:' so 'the derivative of the spectral reconstruction inside it (Appendix 1)' survives - the Appendix Van Loan compress cites it as its retained anchor. Methods keeps the contrasts (enforced on item 69).

- [ ] **COMPRESS**, about 160 words, lines 1521-1541, p.26
  - Passage: "Münch et al. (2022) asked this paper's question of a filter of this family ... Two groups reached one device from different motivations." (whole paragraph)
  - Why: 321-word comparison that can lose half: the three-point contrast survives in compressed form, the placement sentence duplicates the roster (lines 315-317 already place Münch's filter with the recursive instantaneous member), and the 78-word quotation sentence shrinks to its substance with one quote.
  - Do: Residue (~160 words) must retain: (1) Münch asked the same question and answered by counting coverage over repeated simulated data sets; (2) the three separations, each in one tight sentence: their count tests a posterior (prior and sampler inside the verdict) vs identities on the likelihood alone at known parameters; their count returns a coverage probability vs a matrix saying by how much, in which parameter, per interval vs accumulated; one posterior fit per replicate per design point (a line of conditions) vs one pass each (560 cells); (3) one clause placing their filter with the instantaneous recursive member (leaning on the roster); (4) their emission generalization (state-dependent noise, fluorescence) costs the analytical update and obliges a two-moment closure of their own, the same obstruction Eq. 2 meets, keeping the single quote 'the conjugacy property is lost'; (5) the closer 'Two groups reached one device from different motivations.'
  - Cross-references: The two direct quotes appear nowhere else in the paper; the residue keeps one, and if the author wants the 'sufficient statistics' quote preserved it must move to Appendix 2's closure discussion (an Appendix line, 2433, already cites Münch's open-channel-noise device). The roster placement at lines 315-317 must not be cut by the family-section editor, since the residue leans on it.

- [ ] **DELETE**, about 15 words, lines 1621-1622, p.28 (Discussion's closing paragraph)
  - Passage: "What the diagnostic cannot see is the specification the likelihood and the simulator share (Methods)."
  - Why: Duplication (g), third statement: the full shared-specification limitation is the protected third-boundary paragraph of the limitations subsection (lines 1591-1597), which survives as the canonical copy; this closing-paragraph echo is ornamental restatement.
  - Cross-references: Methods lines 1686-1690 carry the second copy and explicitly defer to 'the Discussion': that cross-reference points at the limitations copy (1591-1597), which is protected and survives, so it stays valid; the Methods-span editor should trim their copy against the limitations copy, never the reverse. The closing paragraph still reads cleanly with the sentence removed ('...is the interval that is delivered. A mechanism once left a mark...').

Section editor's note: Span: Discussion, lines 1431-1626 (pp.24-27, spilling to p.28), ~2,680 words; footers at lines 1448, 1502, 1556, 1610 fix the page mapping. The section was already trimmed once, and the ~588 words here (~22% of the span) are what remains extractable without touching protected or first-statement material. Duplication survivorship: (d) Milescu diagnosis survives at the p.10 definition (lines 541-543) and the Results p.14 measurement (lines 861-863), the Discussion copy compresses to a one-line attribution; (e) 2217 survives at Results lines 1061-1064 and the Figure 5-supplement legend (line 3840), the Discussion copy loses the number; (g) shared-specification survives as the limitations third-boundary paragraph (lines 1591-1597, protected), with the closing-paragraph echo deleted here and the Methods copy (lines 1686-1690, outside this span) the remaining candidate for the Methods editor, whose deferral clause 'stated with the others in the Discussion' remains satisfied. Deliberately untouched: the opening calibration paragraph (1432-1436) and the headline map-ordering sentence (1493-1496) as the Discussion's single statements of the two main results; the Del Core cost defense (1542-1550) and the NSFA fairness paragraph (1551-1562), which answer named critics and read as review-demanded; the entire limitations subsection including its hierarchy and hand-over paragraphs (1564-1620, 1622-1624) apart from the one duplication-(g) echo; the Mehra whiteness sentence and the Zadrozny concession, each stated only here. Considered and rejected: compressing the two-algebraic-halves supplement pointer (1440-1443), because the honest residue must keep the not-a-clean-separation caveat and saves under 10 words. Cross-span dependency for the coherence pass: the Münch residue leans on the roster placement at lines 315-317, and the grouping residue leans on Results 1213-1221 surviving in full.

### Results, including main-figure captions

Span: pp. 11-23, lines 571-1430. Section size about 8494 words. Items: 22.

- [ ] **COMPRESS**, about 67 words, lines 576-582, p.11
  - Passage: Throughout, kon and koff are the opening and closing rates, i the unitary current, and Δ̃ and S̃ ...
  - Why: Duplication (a): the natural-units/tilde definition is the second of three copies; the body p.5 copy (line 236) survives, this one becomes a pointer, and the Methods p.36 copy (line 2089) is the Methods span's cut.
  - Do: Residue must retain: a pointer to Table 1 for the eight members; the reminder that Δ̃ and S̃ are the acquisition interval and instrumental noise in the natural units defined above; and that both least-squares arms are run over the whole design plane on the same cells as every other member.
  - Cross-references: kon, koff, i are used throughout Results and are defined in the family section (pp.3-5); S (the PSD) reappears only in Methods; nothing else cites this paragraph.

- [ ] **COMPRESS**, about 47 words, lines 586-594, p.11
  - Passage: Figure 2 reads the ladder at one design cell (Nch = 100, S̃ = 0.01, Δ̃ = 0.1), and each panel shows ...
  - Why: The design-cell numbers duplicate the caption (line 663; keep the caption copy), and the verbal glosses of magnitude and anisotropy are the third statement after the protected Eq. (4) definitions (p.10) and the Figure 2 caption's decode.
  - Do: Residue must retain: cloud away from truth = bias; reported ellipse against empirical = distortion, its two numbers being m and a of Eq. (4), the two orange numbers of each panel; both one when calibrated; read here from a thousand ML fits where every later figure reads them from the score, agreeing only to first order (Appendix 3).
  - Cross-references: Eq. (4) (lines 558-568, protected) and the Figure 2 caption (668-670) must keep their glosses; Figure 4's caption glosses its own colour scale independently.

- [ ] **COMPRESS**, about 34 words, lines 703-707, p.13
  - Passage: they are also the two arms whose analytic Fisher the differenced check flags (Supplementary File 1), and the sizes are consistent ...
  - Why: The arithmetic walk-through (0.71-0.78, squared reciprocal 1.6-2.0) can shrink to one clause; the fuller anchor-check disclosure lives at lines 1006-1013 and is kept.
  - Do: Residue must retain: the sandwich fails for both least-squares arms at a factor 1.6; the failure sits in what the two arms share; one clause that the differenced check (Supplementary File 1) flags their analytic Fisher with consistent size.
  - Cross-references: Lines 1006-1013 independently quote the ~1.3 least-squares re-anchoring factor and cite Supplementary File 1; keep the two consistent.

- [ ] **DELETE**, about 21 words, lines 710-711, p.13
  - Passage: What this figure measures is the separation of the two factors, the harder half and the one the gating variance carries.
  - Why: Ornamental restatement of the identifiability sentence at lines 601-603 (Hines et al., 2014), which is kept.

- [ ] **COMPRESS**, about 42 words, lines 719-726, p.13
  - Passage: The log-likelihood adds up everything a member gets wrong, and it orders the ladder accordingly ...
  - Why: The intermediate nat deltas (eleven nats, sixty-nat gap, instant member below that) are printed on Figure 3 panel A and in Figure 3-source data 1; text and figure should carry the numbers once.
  - Do: Residue must retain: the ladder ordering with the 119-nat top-to-bottom span (pointing to panel A / source data 1); 'the ranking says nothing about how far from the truth its top sits; every other row is read against a value a correct likelihood must give'; the qualification that both least-squares arms fit a noise level to each recording they score.
  - Cross-references: Figure 3 panel A annotations and Figure 3-source data 1 must keep the per-member logL values (they do).

- [ ] **COMPRESS**, about 88 words, lines 727-731 and 835-837, pp.13-15
  - Passage: What remains is measured from three objects and the figure keeps them apart. The residual is the observation against ... (continues after the figure at 835-837: and with variance equal to the information ...)
  - Why: Re-defines the residual, score, and information already defined in the protected diagnostics passage (pp.9-10, lines 489-527), and re-lands 'only the residual can be computed on a real recording' (line 492-493; said again at line 864 where it matters); the Fig 3 caption (818-826) already decodes the rows.
  - Do: Residue must retain: one sentence pointing to the three diagnostics defined in the family section and mapping them to rows (residual: variance (B), memory the black series of (G); score: mean (D), variance against the information (E, F), memory the coloured series of (G); information: (C)).
  - Cross-references: Figure 3 caption's per-row decode (lines 817-826) must survive intact; keep exactly one 'transfers to a real recording' statement in this subsection (line 864's is kept).

- [ ] **DELETE**, about 12 words, lines 841-842, p.15
  - Passage: so averaging the model over the acquisition window is what removes it
  - Why: Third statement of averaging-removes-displacement: the one-cell statement (lines 684-687) and the plane-wide one with numbers (lines 979-982) are both kept.
  - Cross-references: NONE; the sentence remains grammatical ending at '...for the displaced four.'

- [ ] **COMPRESS**, about 50 words, lines 866-876, p.15
  - Passage: The information on its own answers a different question, about the recording rather than about the algorithms ...
  - Why: The mechanism sentence ('A member that has conditioned on the past already carries the number of channels still open ...') is restated in the Discussion (lines 1470-1474), which cites (Results); the measured content stays here.
  - Do: Residue must retain: where each parameter is measured (information floor within six intervals of washout for Nch and kon; koff through the mean and i through the variance of channels still open; Fig 3-supps 1-2); the open-loop double-counting numbers (15-16% of pulse total vs 0.5-0.8%) with 'the same information counted twice'; the noise-level white-score note (supp 3).
  - Cross-references: Discussion 1470-1474 duplicates this paragraph; the coherence pass must keep the mechanism statement in exactly one place - if the Discussion span cuts theirs, this sentence must be restored here.
  - **Adjustment (binding)**: Both-copies conflict resolved: the mechanism sentence STAYS HERE because the Discussion copy (1470-1474) is cut to design advice. Compress shrinks to the framing sentence only; savings 50 to ~18w. The six-interval floor and 15-16% numbers must stay (two supplement legends drop their copies of the six-interval fact).

- [ ] **COMPRESS**, about 30 words, lines 884-890 and 969-971, pp.15-17
  - Passage: Where that crossing sits is worth having in the plane's own coordinates, because it moves ... (continues at 969-971: the longest, which puts it at S̃ = 0.025 ...)
  - Why: The two worked examples (S̃ = 0.025 at ten channels fastest, S̃ ≈ 1400 at ten thousand slowest) restate the law just given in the same sentence.
  - Do: Residue must retain: the crossing law S̃/Nch = 0.0025 at the shortest window to 0.14 at the longest, and that the reference cell of Figures 2-3 sits a factor of 234 in variance below it.
  - Cross-references: Figure 7's boundary paragraph (1406-1410) refers back to this mechanism; the crossing law itself must survive.

- [ ] **COMPRESS**, about 27 words, lines 986-990, p.17
  - Passage: Modelling the correlation across intervals is what removes the distortion, in the amounts above ...
  - Why: The opening sentence restates the subsection title and the numbers already given; the rest of the paragraph is unique.
  - Do: Residue must retain, verbatim or near: 'A member that takes one axis and not the other repairs one failure and keeps the other, so the family cannot be ordered on a single scale, and nothing in that is specific to ion channels: it holds for any state model whose observations are averages over a window.'
  - Cross-references: The generality claim ('holds for any state model whose observations are averages over a window') appears NOWHERE else in the paper (verified by grep) and must survive in the residue.

- [ ] **MOVE**, about 62 words, lines 1014-1019, p.17
  - Passage: The flat floor comes from the split of the distortion into a per-sample part and a correlation part ...
  - Why: The per-member sweep numbers (1.09 to 0.998, 1.10 to 1.007, 1.35 to 1.47, 1.24 to 1.00) belong with the supplement that draws them, not in the body.
  - Do: Figure 4-figure supplement 1 legend. Body keeps a one-sentence residue: the floor comes from the per-sample/correlation split - the members equally faithful to a single interval, parting on the temporal dependence left in the score, which conditioning on both endpoints removes (Figure 4-figure supplement 1).
  - Cross-references: Figure 6 text (1231-1235) and caption (1303-1305) use the same per-sample/correlation decomposition but not these numbers; the supplement legend block (lines 3139-4026) must gain the numbers.
  - **Adjustment (binding)**: Execute by landing the numbers in the Fig 4-supp 1 legend's 'Third' reading and deleting that legend's 'which is the paragraph above'; the captions item on that legend (item 103) is reduced so it does not cut the landing. 62w from the body.

- [ ] **COMPRESS**, about 39 words, lines 1034-1041, p.18
  - Passage: A reader who accepts the distortion will reach for the cheapest repair before changing likelihood ...
  - Why: Duplication (c): the tau_int gloss 'counts how many intervals a likelihood effectively has for every interval it believes it has' is the second of three copies; the p.10 definition (lines 497-499, protected) survives as the only full gloss.
  - Do: Residue must retain: the deflate-the-count setup; tau_int (pointing to its definition above) as the axis Figure 5 carries; the four readings 1.008 (IR), 1.30 (R), 3.93 (NR), 9.26 (ILSE) at S̃ = 0.01; 'about nine times more independent observations than the recording supplies'; computable on a real recording with no ensemble and no known truth.
  - Cross-references: The Figure 5 caption's copy of the gloss is cut in the caption item below; the caption must keep 'reads 1 for a white residual' as its axis decode.

- [ ] **DELETE**, about 24 words, lines 1065-1066, p.18
  - Passage: Where on the design plane a recording leaves that memory is itself a map, drawn for all eight members in Figure 5-figure supplement 3.
  - Why: Bare supplement pointer; the Figure 5 supplement listing (lines 1140-1141) already names supplement 3.
  - Cross-references: No other text cites Figure 5-figure supplement 3; the listing under the caption is its remaining anchor.

- [ ] **COMPRESS**, about 124 words, lines 1116-1136 (caption), p.19
  - Passage: Figure 5. What each likelihood does to the reported error bar, against a number measurable on a single recording. ...
  - Why: The caption restates the subsection's readings at similar length; caption should carry decoding only, and the readings (colour rationale with 8.9-vs-1.6, the ratio law, the reader-holds point, the no-single-correction conclusion) all survive in the Results text at 1042-1074 and 1147-1149.
  - Do: Residue (decoding only) must retain: mark = one grid cell over the same 10,000 recordings as Figure 4, no binning; horizontal axis, shared, a property of the least-squares fit alone = tau_int of its standardized residual, reading 1 for white; vertical axis = distortion of the reported information, point = magnitude (geometric mean of the eigenvalues), bar ends = best- and worst-determined directions; error-bar scale = square root; no threshold drawn, the criterion the reader's; colour = acquisition interval Δ̃; twenty-eight lines per panel, seven visible because the four channel counts land on one another (supp 2); the four members = the family's two axes crossed (Table 1); the least-squares panel reports its own four-parameter fit (Methods), nothing pooled.
  - Cross-references: Text lines 1047-1058 (ratio law and colour rationale) and 1034-1041 must survive as the readings' single home; the Table 1 pointer must remain so the four panel names resolve; dup (c)'s only surviving gloss is p.10.
  - **Adjustment (binding)**: Superseded: execute the captions editor's three finer items on this caption (items 111-113) instead; savings counted there, 0 here.

- [ ] **COMPRESS**, about 38 words, lines 1156-1161, p.19
  - Passage: That reading is on the two invariants; the distortion itself, over the same plane ...
  - Why: Repeats MR-worse-than-R just stated at lines 1154-1155 and frames the supplement as 'saying the same thing'; only the VR reading and the 1.01 endpoint result are new.
  - Do: Residue must retain: Figure 4-figure supplement 3 adds the second partial correction - VR sits on R along the rates and above it on the channel number; neither partial correction buys anything on this plane; the second endpoint is what moves the closing rate to 1.01.

- [ ] **DELETE**, about 19 words, lines 1171-1172, p.19
  - Passage: The contrast is the one the recovery cell reports at a single window, and it holds over two decades of them.
  - Why: Ornamental tie-back to Figure 2; the paragraph's numbers (121.0 vs 121.9, factor 9 to 50, factor 240) already state the contrast.

- **REJECTED, do not execute**: COMPRESS, lines 1213-1221, p.20
  - Passage: One consequence is practical and belongs to the reduction rather than to the acquisition. The sampling rate of the rig is not what this axis is ...
  - Why: Near-verbatim duplicate of the Discussion's grouping paragraph (lines 1475-1492), which is fuller (closure-under-averaging identity, geometric schedule, window-as-argument) and survives; the Results copy shrinks to a link.
  - Do: Residue, one sentence: a fast, uniform record is reduced by grouping before it is fitted, and the numbers above price that grouping - almost nothing for a kinetic question, the whole gain for an amplitude one (Discussion).
  - Cross-references: Discussion 1475-1492 must survive intact - it cites '(Results)' at line 1486 for the finer-groups-buy-fluctuation claim, which lives at lines 1197-1212 and stays; the Discussion span must not also cut its copy.
  - Reason: Both-copies conflict: the discussion editor cut 1475-1492 assuming this copy survives, and this compact copy sits with the sweep numbers it interprets ('What the numbers above say...'). Ruled survivor; the Discussion's larger cut (~180w) proceeds instead.

- [ ] **COMPRESS**, about 25 words, lines 1243-1248, p.21
  - Passage: Figure 7 collects the measurements into one plane, the channel number against the instrumental noise at a single acquisition interval ...
  - Why: The axes and Δ̃ = 0.1 are restated in the caption (lines 1374-1377); keep the numbers in the caption only.
  - Do: Residue must retain: Figure 7 collects the measurements into one plane in which the measured cells are the input and the four boundaries the output.

- [ ] **COMPRESS**, about 41 words, lines 1396-1400, p.23
  - Passage: The two boundaries that say a parameter has stopped being measurable are level sets of one continuous field, the standard error each member actually delivers ...
  - Why: Second verbatim statement of 'delivered standard error = square root of the diagonal of the distortion-corrected covariance'; the first, at lines 1183-1185, survives.
  - Do: Residue must retain: the boundaries are level sets of the delivered standard error (introduced above), drawn for every member and all five parameters in Figure 4-figure supplement 4, so a reader wanting a criterion other than the factor of two can take it off that map.
  - Cross-references: Lines 1183-1185 must survive (they do under these cuts); Figure 4 caption line 958 names the quantity for its line code and stays.
  - **Adjustment (binding)**: Cut only the repeated delivered-SE definition (~30w, not 41); keep the level-set claim, the Fig 4-supp 4 pointer, and the reader-criterion sentence. First statement at 1183-1185 verified untouched.

- [ ] **COMPRESS**, about 23 words, lines 1406-1410, p.23
  - Passage: The boundary at which the classical error bar becomes honest is one the mechanism predicts: the gating variance grows with the channel count ...
  - Why: Third statement of the gating-vs-instrumental-variance mechanism (after lines 884-887 and 993-996; the 993-996 copy is kept as the one clear statement).
  - Do: Residue must retain: the boundary is the one the mechanism above predicts; fitted in logarithms over four measured crossings in each of two parameter directions it runs as instrumental noise proportional to channel number, slopes 1.01 and 1.03, a measured slope rather than an assumed scaling.
  - Cross-references: Item at 884-890 is also trimmed; the surviving mechanism statement is lines 993-996 - coherence pass should confirm it stands.

- [ ] **DELETE**, about 19 words, lines 1427-1428, p.23
  - Passage: ribbons showing where each moves as the interval is swept from 0.01 to 1 in units of τ.
  - Why: A dangling sentence fragment (editing artifact) whose content the Figure 7 caption carries in full ('Ribbons are the same boundary swept over Δ̃ from 0.01 to 1 ...', lines 1381-1382).
  - Cross-references: NONE - but flag to the author that the preceding sentence ('...the ten cells that would pin it have not been run.') is where the paragraph should end.

- [ ] **COMPRESS**, about 30 words, lines 943-967 (caption), p.16
  - Passage: Figure 4. The design plane: how wrong each member's report is, in both moments, on one scale. ...
  - Why: Two small trims in an otherwise decoding-only caption: the theta-pool evaluation rationale ('where the score vanishes by construction so a displaced gradient cannot contaminate it') compresses to '(Methods)', and 'and the difference is not cosmetic' goes.
  - Do: Residue must retain everything else: both halves' definitions, the column nesting, the least-squares four-parameter configuration note, the colour scale with the (A)-factor-on-parameter vs (B)-ratio-of-variances distinction, and the two line codes with both thresholds.
  - Cross-references: Verify Methods states why the distortion is evaluated at the pooled optimum; if it does not, keep the clause and drop this item.

Section editor's note: Span = Results, lines 571-1430 (pp.11-23). The wc count 8,494 includes roughly 700 words of figure furniture (axis labels, tick numbers, in-figure keys); prose is about 7,800, so the 972 proposed is roughly 12% of the section. Results is where the protections concentrate (all seven figures, the m/a summaries, the safeguard disclosures at 1006-1013, 1023-1030, 1236-1241, 1414-1418, 1422-1428, and the readings themselves), so the compressible mass here is caption-text duplication and second statements; the 35-40% overall target must fall mainly on the family section, Discussion, Methods, and appendix-bound material. Duplication survivorship in this span: (a) natural units - body p.5 copy (line 236) survives; Results copy (576-582) compressed to a pointer; Methods copy (line 2089) is the Methods span's cut. (c) tau_int gloss - p.10 definition (lines 497-499) survives as the only full gloss; both my copies (text 1037-1038, caption 1120-1121) are cut, the caption keeping only 'reads 1 for a white residual'. (d) Milescu diagnosis - the Results p.14 copy (lines 861-863) SURVIVES because it anchors the measured lag-one autocorrelations (-0.004 vs 0.191 vs 0.864); the p.10 copy (lines 542-543, inside the correlation-distortion definition) and the Discussion copy (1443-1444) are the other spans' cuts. (e) factor-2217 - the Results p.18 copy (lines 1061-1064) SURVIVES (tied to Fig 5-supp 1); the Discussion copy (1467-1469) is the Discussion span's cut. Cross-span handoffs for the coherence pass: (1) the one-to-three-decades ordering appears at Results 1411-1414 (survives, it is the subsection's result) and Discussion 1493-1496 (their cut); (2) 'concept map, not a phase diagram' appears in the Fig 7 caption 1387-1388 (survives - it travels with the figure) and Discussion 1496-1497 (their cut); (3) the washout/conditioning mechanism is paired between my compressed 866-876 and Discussion 1470-1474 - exactly one must keep the mechanism sentence (I cut mine on the assumption Discussion keeps theirs, which also carries the design advice); (4) the grouping/reduction paragraph is paired between my compressed 1213-1221 and Discussion 1475-1492 - the Discussion copy survives and must not also be cut, since it cites (Results) for lines 1197-1212, which stay. Artifacts noticed for the author: the sentence fragment at 1427-1428 (cut proposed) and the lowercase sentence start 'recursive members are flat instead' at line 1176. Not proposed anywhere: the two-closures material, Table 1, the three diagnostics definitions, the m/a summaries, any figure, the anchor-check and grey-cell disclosures, the NSFA scope statement, or the thin-evidence paragraph.

### The likelihood family and its diagnostics

Span: pp. 3-10, lines 142-570. Section size about 4569 words. Items: 16.

- [ ] **COMPRESS**, about 20 words, lines 157-163, pp.3-4 (first paragraph of 'The observable is an interval average', before Eq. 1)
  - Passage: "A recorded sample is the current averaged over the acquisition window rather than sampled at an instant: the current is low-pass filtered in analogue form before it reaches the converter, so each stored value carries a weighted average of the recent current and not a reading of it at a point."
  - Why: The sentence states the same fact twice (averaged-not-instantaneous, then weighted-average-not-point-reading); one statement plus the analogue-filter reason suffices to set up Eq. 1.
  - Do: Residue must retain: recorded sample = window average, because the current is analogue low-pass filtered before the converter (~32 words).
  - Cross-references: NONE. Eq. 1 (protected) and the post-equation motivation (lines 169-173) carry the argument.

- [ ] **COMPRESS**, about 65 words, lines 174-182, p.4
  - Passage: "Uniform weighting is a first approximation to what the electronics do ... block-average a record down to that scale before analysis." (whole paragraph)
  - Why: The middle sentence (~63 words on the dimensionless group, first-order vs second-order residue, where it lands in the predictive variance) is derivation-level detail that Appendix 1 already carries (its block-average bound is at line ~2255); the body needs only the disclosure and the pointer.
  - Do: Residue must retain three things: (i) uniform window is an approximation to the amplifier's impulse response; (ii) simulator and every member share the uniform window (Methods), so the window does not enter the measurement and the closures are what the diagnostics see; (iii) Appendix 1 bounds the residue and the remedy is to block-average before analysis (~75 words). The cut sentence's content is NOT moved -- it already exists in Appendix 1.
  - Cross-references: This is a safeguard disclosure -- it must survive in compressed form, not be deleted. Clause (ii) is a body-side variant of duplication (g) (simulator-shares-specification): the Methods copy is the surviving full statement, the Discussion copy is the other cut (Discussion span's job); this residue keeps only the one clause. Appendix 1's bound must not be cut by the appendix pass.

- [ ] **COMPRESS**, about 14 words, lines 189-192, p.4 (end of paragraph introducing the two closures)
  - Passage: "An instantaneous member drops that evolution instead, one snapshot standing in for the average the recording holds, which misspecifies the random variable rather than closing it, and in the unhelpful direction, since averaging is what brings the law closer to Gaussian."
  - Why: "One snapshot standing in for the average the recording holds" restates lines 170-173 ("a likelihood built on instantaneous sampling describes the wrong random variable") twenty lines earlier.
  - Do: Residue must retain: an instantaneous member misspecifies the random variable rather than closing it, in the unhelpful direction, since averaging brings the law closer to Gaussian (~28 words).
  - Cross-references: Results use the misspecification-vs-closure distinction for NR's bias (line ~599-600); the distinction survives in the residue.

- [ ] **COMPRESS**, about 28 words, lines 224-227, p.5 (first paragraph after the two-closures box)
  - Passage: "The closure is indifferent to the shape of the law it replaces ... what justifies the closure beyond the two-state case run here."
  - Why: The second sentence's "exact law that degrades with K and with Nch" restates the box's own "closed only at K=2 and proliferating with K" (lines 209-211); only the moments-side of the asymmetry and the beyond-K=2 conclusion are new here.
  - Do: Residue must retain: the moments the closure keeps are closed-form single-channel objects, independent of Nch and of the topology (Appendix 1), which is what justifies the closure beyond the two-state case run here (~30 words).
  - Cross-references: The limitations subsection and Discussion lean on the beyond-K=2 justification; it survives in the residue. The box (protected) keeps the K-degradation half.

- [ ] **DELETE**, about 34 words, lines 233-235, p.5 (final sentence of the three-regimes paragraph)
  - Passage: "Any Gaussian likelihood of this kind must fail somewhere, and this paper's aim is to show that the failures are the predicted degradation of these two closures rather than defects of a particular algorithm."
  - Why: Restates the paper's aim, which the Introduction already carries (lines 111-123) and the Results enact; "must fail somewhere" appears nowhere else, so nothing dangles.
  - Cross-references: The three regime NAMES in the same paragraph (microscopic, telegraphic, Gaussian; lines 228-233) are used by Results at lines 1252 and 1403-1404 and MUST stay; the deletion is the final sentence only.

- [ ] **COMPRESS**, about 30 words, lines 236-243, p.5
  - Passage: "Two design quantities are reported throughout in the natural units of the model, and a tilde marks them ... measured against the single-channel signal power (Methods)." -- duplication (a), body copy
  - Why: Duplication (a): this body copy survives only as the bare first-use definition (Figure 1's caption uses S-tilde on p.8, before Results, so the tilde must be defined by p.5); the interpretive tail "so S-tilde is the noise power accumulated over one channel time constant measured against the single-channel signal power" is verbatim in Methods (lines 2095-2096), which is the surviving full specification.
  - Do: Residue must retain: tilde marks natural units; time in the channel time constant tau = 1/koff (the closing time constant), current in the unitary current i; Delta-tilde = Delta*koff; S-tilde = S*koff/i^2 with S the instrumental noise power spectral density; pointer to Methods (~55 words).
  - Cross-references: Duplication (a) survivor accounting: Methods copy (lines 2089-2099) survives in full; this body copy survives compressed; the third copy at the Results opening (lines 576-579) is the outright cut -- that is the Results span's item, flagged for the coherence pass.
  - **Adjustment (binding)**: Per duplication-(a) ruling this becomes the SOLE full definition: keep the interpretive tail ('noise power accumulated over one channel time constant measured against the single-channel signal power'), since the Methods copy sheds its re-derivation. Savings 30 to ~10w.

- [ ] **MOVE**, about 134 words, lines 244-263, pp.5-6
  - Passage: "Both are per-sample variances, which is what makes them comparable ... which is Delta-tilde/4 for windows short against that correlation time and a quarter for long ones." (the per-sample interpretation, the S-tilde=1 gloss, the gating-variance saturation argument, and the equal-variance display...
  - Why: This derivation of where instrumental and gating variance cross is design-plane machinery, not family definition; the display equation is cited nowhere else in the paper (grep: only hit is line 261), and half its content (the S-tilde=1 interpretation) already exists verbatim in Methods at lines 2098-2099 and 2112-2114, where it will simply dissolve on merge.
  - Do: Methods, the natural-units/design-plane passage at lines 2088-2114, merged so the S-tilde=1 sentences are not duplicated. The BODY must keep the paragraph's final sentence (lines 263-265, ~24 words, not counted as saved): "the channel count and the instrumental noise enter the comparison through their ratio, which is the variable the design plane collapses in (Results)" -- fused onto the compressed tilde-definition paragraph.
  - Cross-references: Results lines 1047-1055 and Figure 5-figure supplement 2 rest on the ratio-collapse claim -- covered by the kept sentence. The equation itself must land intact in Methods because Figure 5's ratio axis is its consequence; coherence pass should confirm the Methods insertion point survives the Methods span's own cuts.
  - **Adjustment (binding)**: The landing (Methods ~2109-2117, verified surviving the Methods cuts) is also main text, so count net-of-landing: ~80w, not 134. Keep the closing ratio-collapse sentence (263-265) in the body; the S-tilde=1 half dissolves against 2113-2114.

- [ ] **COMPRESS**, about 33 words, lines 267-275, p.6 (opening of 'Least squares models no gating...')
  - Passage: "Before any of the choices above comes a prior question: does the method model the gating fluctuations at all? ... because it is the method most readers are using."
  - Why: The clause "so it carries neither an occupancy distribution nor a conductance conditioned on anything, and the temporal structure the two closures approximate is absent from it" elaborates what "fits the deterministic mean and one constant variance" already says and what Table 1's LSE row states.
  - Do: Residue must retain: the prior question (does the method model gating at all), least squares answers no (deterministic mean, one constant variance per residual), it is therefore not a member of the family, and it is the bottom rung and the anchor because it is the method most readers are using (~95 words).
  - Cross-references: Results and Discussion name it "the classical fit" throughout -- naming unaffected; Table 1 rows LSE/ILSE (protected) carry the dropped detail.

- [ ] **COMPRESS**, about 42 words, lines 284-288, p.6
  - Passage: "The window axis is not the family's property alone ... which no pair inside the family does."
  - Why: The second sentence ("Least squares therefore has two arms here, identical but for that decision, so the pair isolates the acquisition average with nothing else moving, which no pair inside the family does") is functionally verbatim in Table 1's caption at lines 475-476 ("LSE and ILSE differ in the mean alone, which isolates the acquisition average with nothing else moving"), which survives because Table 1 is protected.
  - Do: Residue must retain: the window decision exists at every level of the gating description, least squares included, which is why least squares has two arms (~35 words).
  - Cross-references: Results comparisons of LSE against ILSE invoke the isolation point; the surviving copy is Table 1's caption. Coherence pass: Table 1 caption must not lose that clause.

- [ ] **COMPRESS**, about 168 words, lines 289-304, p.6
  - Passage: "The members are named from the two axes and the names are read off them. The suffix is the recursion axis ... and what the least-squares arm discards is the variance (Methods)." (the member-naming walkthrough)
  - Why: The suffix/prefix tutorial spells out in prose exactly what Table 1's layout and caption show (second column = window axis, two blocks = recursion axis), and the "to within 10^-15 pA" LSE/NR equivalence demonstration is stated in Methods at line 1715 with a sharper number (9x10^-16), which is the surviving copy of that fact.
  - Do: Residue (~110 words) must retain, in order: names are read off the two axes (prefix = window: none for an instant, M for the interval mean given the start state, I for the interval mean given both ends; suffix = R recursive, NR non-recursive); five gating members rather than six because without an update the end state is never observed, so MNR coincides with INR; M is the mean conductance, not an abbreviation of macroscopic; VR is MR with the residual interval variance in place of the total one, the one member whose variance form is a choice rather than a consequence (Appendix 2); LSE and ILSE sit off the lattice, each the open-loop member at its own window setting with the gating variance replaced by one fitted constant (Methods).
  - Cross-references: The VR characterization is the Results' negative control and Table 1's caption restates it (protected, survives). Appendix 2 pointer must stay attached to VR. Methods line 1715 must survive the Methods pass, since it becomes the only statement of the LSE=NR mean identity.
  - **Adjustment (binding)**: Residue must keep: the suffix/prefix naming key (one sentence), the five-not-six sentence, the M-not-macroscopic note, VR-as-choice, and the LSE/ILSE off-lattice + variance-discard statement (the Appendix 2 least-squares cut leans on it). Savings 168 to ~150w. The 10^-15 demonstration may drop only because Methods keeps the 9e-16 sentence (enforced on item 63).

- [ ] **DELETE**, about 60 words, lines 305-309, p.6
  - Passage: "Below, a member is named by its code where a number is attached to it and by what it conditions on where the point is structural: NR is the open-loop member ... and the two least-squares arms together the classical fit."
  - Why: This prose-name glossary is Table 1's 'Called, in the Results' column (lines 440-470) restated as a sentence; the column is the surviving definition.
  - Do: Replace with a ~7-word pointer folded into the preceding compressed paragraph: "Table 1 lists the name each goes by in the Results."
  - Cross-references: HIGH-TRAFFIC names: 'boundary-conditioned' appears 19 times, 'open-loop member' 8, 'classical fit' 6, 'one-endpoint' 5, 'variance-corrected' 2, across Results and Discussion. All resolve through Table 1's last column, which is protected; coherence pass must verify no later span cuts that column.

- [ ] **COMPRESS**, about 32 words, lines 310-321, pp.6-7
  - Passage: "Two things follow. IR occupies the top rung ... the ground truth every rung below it is judged against."
  - Why: "which conditions on the full trajectory inside each window and is intractable, which is why the closures are needed at all" repeats the section opening (lines 149-152) and the same sentence's own "the intractable exact likelihood".
  - Do: Residue must retain in full: IR's standing is structural, set before any measurement; the lower rungs are established likelihoods, NR being the independent-interval likelihood of Milescu et al. (2005) and R the macroscopic filter of Moffatt (2007) and the generalized filter of Muench et al. (2022); the exact likelihood sits above all and the exact forward simulation realizes it, so the simulation is the ground truth (~100 words).
  - Cross-references: The literature identifications are review-critical and referenced by the Discussion's Milescu passage (line 1443ff) and Introduction (line 135); they must survive verbatim.

- [ ] **DELETE**, about 11 words, line 341, p.7 (inside the paragraph following Eq. 2)
  - Passage: "It is also the object the interval version has to rebuild."
  - Why: Forward tease whose content the K-squared substitution paragraph (lines 348-352) delivers concretely seven lines later.
  - Cross-references: NONE.

- [ ] **COMPRESS**, about 75 words, lines 344-358, p.7
  - Passage: The three paragraphs after Eq. 2's gain discussion: "What it gives up is not exactness but sufficiency ...", "What stops Eq. 2 from being used on the interval average ...", "In a recursive member the posterior pair becomes the prior ... what the diagnostics below measure."
  - Why: Three paragraphs merge to two; connective scaffolding ("For this reason", "whose only approximation is the Gaussian replacement at the end of it", "which the Gaussian projection does not guarantee") compresses while every distinct claim survives once.
  - Do: Residue (~155 words) must retain all four claims: (i) Eq. 2 consumes exact moments but Var[y-bar | N] depends on N, so (N, y-bar) is not jointly Gaussian and Eq. 2 is not the exact Bayesian update -- a heteroscedastic emission treated as homoscedastic; (ii) the average is not a function of the occupancy at any one time, so no conductance vector on the K states delivers its mean, but on the K^2 pairs of endpoint states one does, and every interval member is that substitution closed by a Gaussian (Appendix 1) -- this is duplication (b)'s SURVIVING copy; (iii) carrying the posterior forward is itself a closure, exact only while the end marginal is a sufficient summary of the boundary posterior; (iv) the recursive members' Gaussian information is therefore a model-based quantity, and how far it departs from the truth is what the diagnostics measure.
  - Cross-references: Duplication (b): this body copy survives; the two Appendix 1 copies are the cuts (appendix span's job). Claim (iv) is the bridge the diagnostics subsection opens from, and the Discussion's honesty claims echo (i); neither may vanish. Duplication (f) sits just above at lines 326-328 ("Chaining intervals leaves a single state variable per junction") -- keep untouched here as the SURVIVOR; the Appendix 1 copy at line 2471 is the cut.

- [ ] **COMPRESS**, about 45 words, lines 416-431, p.8
  - Passage: Figure 1 caption: "One filter step, and what each of the two axes changes. ... No number is quoted from this simulated recording."
  - Why: The sentence "Each stands for more than itself: least squares is the open-loop column of its own window setting with a band of constant height (above), and MR and VR are the recursive column carrying INR's prediction and innovation (Table 1)" (~40 words) restates Table 1's mapping and shrinks to "LSE/ILSE and MR/VR map onto these columns (Table 1)"; the S-tilde sentence tightens.
  - Do: Residue must retain: the setup (Nch=20, green trace, one number per window, 2 ms windows, S-tilde=0.01 so the recorded value sits on the average), the four columns named, the rows (A)-(D) walk, the row-(B) reading of the window axis, the IR dashed-arc-is-covariance-not-path sentence, the colour key, and the closing disclosure "No number is quoted from this simulated recording."
  - Cross-references: Table 1's caption cross-references Figure 1 ("the four that cross the two axes are drawn as one turn of the filter cycle in Figure 1") -- unaffected. Figure 1 itself is kept: it is the only rendering of the filter cycle and Table 1 points to it.

- [ ] **COMPRESS**, about 34 words, lines 552-556, p.10
  - Passage: "The correlation part has a reader-facing translation. Where R is one number times the identity ... no such correction suffices at any factor and the full matrix is needed."
  - Why: The clause "dividing the record's interval count by it gives the count the likelihood effectively has" plus the thousand-intervals-holding-a-hundred illustration echo the tau-int gloss defined 55 lines earlier (lines 497-499) and are re-illustrated concretely in Results at line 1039 ("about nine times more independent observations").
  - Do: Residue must retain: where R is one number times the identity, that number is the factor by which the likelihood over-counts its intervals; where R is anisotropic, or Cs differs from I, no scalar correction suffices and the full matrix is needed (~45 words).
  - Cross-references: The anisotropy caveat is the Discussion's argument for the sandwich over a scalar inflation -- it survives in the residue. Duplication (c): the definitional gloss at lines 497-499 is the SURVIVING copy; this cut removes the in-span echo; the Results p.18 copy (lines 1036-1038) and Figure 5 caption copy (line 1120) are the other spans' cuts (an Introduction copy at lines 124-126 and a figure-supplement-legend copy at lines 4024-4025 also exist -- flag for those spans).

Section editor's note: Span = lines 142-570 (pp.3-10), measured at 4,569 words including Figure 1's caption and Table 1. Roughly 1,700 of those words are protected (two-closures box 198-221, Table 1 with caption 437-483, the three diagnostics definitions with Eqs. 3-4 and the m/a summaries 489-568, Eqs. 1-2); the 825 proposed savings are ~29% of the remaining prose, front-loaded on the member-naming walkthrough (cuts at 289-309 total 228 words) as the span guidance directed. Deliberately untouched and why: the section-opening thesis (149-154, the two-objects framing every member varies on); the two-axes definition with the prequential/open-loop factorizations (276-283, the formulas appear nowhere else); the boundary-state definition (323-328, needed for Eq. 2 and Table 1's IR row, and duplication (f)'s survivor); the theta-sim/theta-pool paragraph (360-366, every figure states which vector it shows); the F_t/J_t accumulation sentence (543-551, verified load-bearing: Figure 3 panels E-F and its source data plot Var(s_t)/F_t and J_T/F_T); the White/Huber/Godambe attributions and the non-circularity and bottom-rung disclosures (513-527). Duplication survivor decisions inside this span, for the coherence pass: (a) natural units -- Methods copy (2089-2099) survives full, body copy survives compressed to the bare definition plus the ratio sentence, the Results-opening copy (576-579) is the outright cut belonging to the Results span; (b) average-not-a-function-of-occupancy -- body copy (348-352) survives, both Appendix 1 copies are cuts; (c) tau-int gloss -- the definition (497-499) survives, the echo at 554 is cut here, the copies at Results 1036-1038, Figure 5 caption 1120, Introduction 124-126, and supplement legend 4024-4025 belong to other spans; (d) Milescu local-time-correlation -- the attribution at the diagnostic's definition (541-543) survives, the Results repeat (862) and Discussion repeat (1443) are other spans' cuts (Discussion can keep 'The diagnosis is Milescu's, now measured' as a clause but not the re-explanation); (f) chaining-single-state-variable -- body copy (326-328) survives, Appendix 1 copy (2471) is cut; (g) simulator-shares-specification -- Methods copy survives, this span retains one clause inside the compressed uniform-weighting residue (cut 2), the Discussion copy is cut. Dependencies the coherence pass must hold: Table 1's 'Called, in the Results' column and its caption clause 'LSE and ILSE differ in the mean alone...' become sole definitions after cuts 9 and 11; Methods line 1715 becomes the sole statement of the LSE=NR mean identity after cut 10; Methods 2088-2114 must absorb the moved equal-variance derivation (cut 7) without duplicating its existing S-tilde=1 sentences.

### Abstract and Introduction

Span: pp. 1-3, lines 1-141. Section size about 1749 words. Items: 10.

- [ ] **COMPRESS**, about 8 words, lines 24-26, p.1 (Abstract)
  - Passage: "the smallest that isolates the approximation's own error," (appositive inside the abstract's sweep sentence "The sweep runs over channel number, noise and interval in a two-state scheme...")
  - Why: The isolation-by-design claim is restated at intro line 102 ("which isolates the error the approximation itself commits") and again in the Discussion limitations (line 1599); the abstract sentence stands without the appositive and every revised claim (two-state, ten thousand exact simulations of the process they approximate, the sweep coordinates) survives.
  - Do: Residue must retain: "The sweep runs over channel number, noise and interval in a two-state scheme, ten thousand exact simulations of the process they approximate at every point." The surviving copy of the isolation claim is intro line 102.
  - Cross-references: Nothing elsewhere quotes the abstract wording; coherence pass should confirm intro line 102's copy is kept by item on lines 98-110 below (it is, in that item's residue).

- [ ] **COMPRESS**, about 35 words, lines 72-75, p.2 (Introduction, prior-likelihoods paragraph)
  - Passage: "differing in how much gating stochasticity they carry and in whether the observable is an instantaneous sample or the interval average a recording actually holds, the most recent, MacroIR, conditioning that average on the channel states at both its ends (Moffatt and Pierdominici-Sottile, 2025)"
  - Why: This 44-word clause pre-explains the family's two axes (what is conditioned on; instantaneous vs interval observable) that the dedicated section states at full length (lines 276-288 plus Table 1, both protected survivors), and MacroIR was already defined with the same citation at lines 44-46.
  - Do: Residue (~8 words) must retain: that the whole-record likelihoods differ in what they condition on, and that MacroIR is the most recent. E.g. "...(Münch et al., 2022), differing in what they condition on, MacroIR the most recent." Survivor of the full two-axes exposition: the family section, lines 276-288 and Table 1.
  - Cross-references: The 'interval average is the observable' idea must still be introduced before line 156; it is, at lines 45-46 (MacroIR 'averages the current over each acquisition interval') and at line 58-59. The Moffatt 2007 expectations at lines 76-79 (unique, tested later) must stay untouched.

- [ ] **COMPRESS**, about 30 words, lines 80-85, pp.2-3 (Introduction, ensemble paragraph)
  - Passage: "No statistic computed on a single recording settles it either. Whether the error bar a method reports matches the spread it delivers is a property of the sampling distribution, and one recording is one draw from it; the generating model is unknown as well, so a departure may belong to the approx...
  - Why: The per-check account of what can and cannot be run on one experiment is given definitively in the protected diagnostics definitions (lines 489-507: 'needing no ensemble and no known truth', 'cannot be run on an experiment at all', 'one recording gives one score'), so the intro needs only the claim, not the sampling-distribution explanation.
  - Do: Compress 92 to ~62 words. Residue must retain: (i) no single-recording statistic settles it; (ii) the confound clause, that with the generating model unknown a departure may belong to the approximation or to the model (stated nowhere else); (iii) the ensemble from known parameters; (iv) "the forward process simulates exactly where its likelihood cannot be evaluated exactly" verbatim.
  - Cross-references: The exact-forward-simulation point recurs at lines 319-321 (family section, kept); Discussion lines 1618-1621 echo the reports-vs-delivers theme but quote nothing from here.

- [ ] **COMPRESS**, about 45 words, lines 86-97, p.3 (Introduction, instruments paragraph)
  - Passage: "The instruments are then classical. For a correctly specified likelihood the score, the gradient of the log-likelihood in the parameters, has mean zero... equals the Fisher information the likelihood reports (Huber, 1967; White, 1982)... Near the maximum that same matrix enters the Bayesian evid...
  - Why: This is the flagged pre-explanation of the diagnostics: the score and information identities are stated three times (abstract lines 20-22, here, and the protected definitions at lines 489-527 with the same Huber/White citations), so the intro copy shrinks to a one-clause pointer while its unique content stays.
  - Do: Compress 176 to ~130 words. Residue must retain: one clause stating both identities with the Huber 1967 / White 1982 citations; "measured... rather than a test to be passed"; the two unique negatives IN FULL (goodness of fit can look good while information is off by an order of magnitude; comparing two approximations answers nothing because both may be wrong); a tightened evidence link (the matrix enters the Bayesian evidence, so misreporting it propagates into every comparison built on it); and the scope-disclaimer sentence VERBATIM ("The correction is derived elsewhere; this paper measures its input and computes no evidence, and it never benchmarks against non-stationary fluctuation analysis (Stepanyuk et al., 2014)"), which reads as review-demanded and is the only statement of the NSFA non-benchmark.
  - Cross-references: Survivor of the full identity definitions: protected section lines 500-527. The evidence-volume point recurs at Discussion line 1616 and Appendix 4 line 3134, so tightening is safe; the Stepanyuk 2014 citation exists only here, keep it.
  - **Adjustment (binding)**: Residue must keep one clause linking misreported information to evidence comparisons (its volume-term mechanism appears nowhere else: Discussion 1616 and Appendix 4 only say evidence is untouched by covariance repair) and the Stepanyuk 2014 citation. Savings unchanged (~45w).

- [ ] **COMPRESS**, about 45 words, lines 98-110, p.3 (Introduction, design-plane paragraph)
  - Passage: "We measured both departures over a plane... The rate constants are held fixed because those coordinates are dimensionless, the closing rate folding into the interval. The scheme is two states at an open probability of one half... most favourable case for a Gaussian treatment of the occupancy; th...
  - Why: The design rationale is restated in the natural-units passage (lines 236-250), Methods (lines 1640-1648: p=0.5 maximizes gating variance, zero skewness) and the Discussion limitations (lines 1599-1616, protected), so the intro keeps the design facts and drops the justifications.
  - Do: Compress 166 to ~120 words. Residue must retain: the plane's three coordinates; rate constants held fixed; two states at open probability one half, single concentration jump; "which isolates the error the approximation itself commits" (this becomes the surviving copy once the abstract appositive is cut); a one-clause conjecture hedge ("whether it bounds anything richer is a conjecture, not a result of this paper"); ten thousand recordings per cell scored by the eight members from least squares to the both-ends filter; and the sentence "An implementation error shows up as a departure like any other, so the same test checks the code" verbatim (its only occurrence). Cut: "because those coordinates are dimensionless, the closing rate folding into the interval" (survives at line 238 and Methods) and "is the most favourable case for a Gaussian treatment of the occupancy" (survives at Methods 1646-1648 and Discussion 1599-1603).
  - Cross-references: Discussion line 1615-1616 restates the conjecture ("that the two-state case bounds anything richer is conjectured") and is the survivor of its full form; no other text quotes "most favourable case".
  - **Adjustment (binding)**: Keep 'which isolates the error the approximation itself commits' (line 102) verbatim, since the abstract's only other early copy is cut by the first item; savings ~40w.

- [ ] **COMPRESS**, about 125 words, lines 111-123, p.3 (Introduction, results-preview paragraph)
  - Passage: "The plane separates two failures usually named as one... On the first the ladder is nearly flat, since averaging the model over the acquisition interval removes the channel-number bias... the discrepancy is carried by the correlation in time of the score... and it carries an ordering independent...
  - Why: This is the third statement of the findings within three pages: every result here is in the abstract (lines 26-34) and stated in full in the Results (bias flatness and its mechanism at 595-603; the score-time-correlation mechanism at 541-543 and Figure 3; IR's failure corner at 1152-1154; the no-region ordering at 1411-1413), so the intro keeps a headline, not the mechanisms.
  - Do: Compress 210 to ~85 words. Residue must retain: the framing sentence "The plane separates two failures usually named as one: what a method gets wrong, and what it says about how wrong it is" (the paper's organizing split, first stated here); one clause that on the first the ladder is nearly flat; one clause that the correlation-blind members misreport their information by an order of magnitude while the member conditioning on both ends is calibrated over almost the whole plane; and the self-audit rationale "a method whose limits are unmapped has been endorsed rather than characterised" (unique to this paragraph). Cut: the averaging-removes-bias mechanism (survivor Results 595-603), "least squares carries no bias in what it fits", the score-time-correlation mechanism and "does not show up one interval at a time" (survivors lines 541-550, protected), the failure location "fewest channels and lowest noise" (survivors abstract 31-32 and Results), and the ordering sentence "no region of the plane gives both..." (survivors: abstract lines 32-34 and Results line 1411-1413; the intro copy is the one cut).
  - Cross-references: Duplication of the no-region ordering is three-way (abstract, here, Results 1411): abstract and Results copies survive. Coherence pass should confirm the Results-section editor keeps line 1411 ("the result the map exists to state") since it becomes the sole full in-text statement.

- [ ] **COMPRESS**, about 31 words, lines 124-128, p.3 (Introduction, takeaway paragraph, first half)
  - Passage: "Making the measurement takes an ensemble. Using it does not: the integrated autocorrelation of the standardized least-squares residual counts how many intervals a fit effectively has for every interval it believes it has, and it costs nothing... and that is a number the reader chose as they chos...
  - Why: This is a fourth copy of the tau_int gloss beyond the three the brief lists (the protected p.10 definition at lines 498-499, Results p.18, Figure 5 caption line 1120), and "the fit was going to be made anyway" also duplicates the Figure 5 caption (line 1121).
  - Do: Compress 76 to ~45 words. Residue must retain: the pivot "Making the measurement takes an ensemble. Using it does not"; that the integrated autocorrelation of the standardized least-squares residual reads the effective interval count off a fit already being made, at no cost (WITHOUT the "counts how many intervals... believes it has" gloss, whose surviving copies are the p.10 definition and the Figure 5 caption); and the analogue-filter pickup claim ("it also picks up the analogue filter on a record sampled faster than its own bandwidth"), which appears NOWHERE else in the paper and must survive here. Cut: the gloss itself and "that is a number the reader chose as they chose the sampling rate" (the reader-holds-it point survives in the Figure 5 caption).
  - Cross-references: Intro copy of known duplication (c): survivors are the p.10 definition (line 498-499, protected) and Figure 5 caption (line 1120); the Results p.18 copy belongs to another editor's span. The filter-pickup claim is load-bearing and unique - flag to coherence pass that it must not be dropped in execution.
  - **Adjustment (binding)**: Cut only the counting gloss and the costs-nothing clause; keep 'Making the measurement takes an ensemble. Using it does not' and the unique analogue-filter-pickup sentence (127-128).

- [ ] **COMPRESS**, about 60 words, lines 128-134, p.3 (Introduction, takeaway paragraph, implementation sentences)
  - Passage: "What is new here is narrow... Its derivatives are. The gradient of the log-likelihood is analytic and comes out of one implementation: the recursion is written once, over a templated C++ type that carries a value and its derivative together... The Fisher information is then wired by hand as the ...
  - Why: The implementation mechanism is restated near-verbatim in Methods p.36 (lines 2137-2145: "one templated C++ type", "wired by hand, as the Gaussian form of Eq. 3", "free of a second pass over the record"), and the no-extra-pass fact also sits at Eq. 3 (lines 517-523); the Methods copy survives, the intro keeps only the novelty claim.
  - Do: Compress 99 to ~40 words. Residue must retain: "What is new here is narrow, and the likelihood's construction is not part of it: that was published with the P2X2 analysis. Its derivatives are" plus one pointer clause, e.g. "the gradient is analytic and the Fisher information follows from the same pass (Methods)". Cut: the templated-C++-type mechanism, the score-returns-with-the-log-likelihood detail, and the wired-by-hand/no-second-pass sentences. Survivor: Methods lines 2137-2145.
  - Cross-references: Discussion lines 1510-1515 ("its differentiation through the recursion", "The derivatives are what the diagnostic could not be run without") lean on derivatives having been introduced - the residue's novelty claim covers that. This is a body-Methods near-duplication of the same kind as the listed pair (g).

- [ ] **COMPRESS**, about 12 words, lines 134-136, p.3 (Introduction, prior-credit sentence)
  - Passage: "Analytic gradients of a likelihood over independent intervals are Milescu et al. (2005); a joint estimate of rates and amplitudes from fluctuation data, carrying standard errors, is Anderson and Stevens (1973)."
  - Why: The Milescu analytic-gradients credit is given more fully in the Discussion (lines 1512-1514, survivor), so the intro's credit sentence can shrink to one clause naming both while the adjacent novelty sentence ("What this work adds... verified rather than assumed", lines 136-138) stays untouched.
  - Do: Compress 31 to ~19 words, e.g. "Standard errors from fluctuation data go back to Anderson and Stevens (1973) and analytic gradients to Milescu et al. (2005)." Must retain both names with both citations (the Anderson-Stevens standard-errors credit exists only here) immediately before the novelty sentence, which must survive verbatim.
  - Cross-references: The novelty claim at lines 136-138 depends on these credits staying adjacent; Discussion 1512-1514 is the surviving full Milescu treatment.

- [ ] **DELETE**, about 11 words, lines 138-139, p.3 (Introduction, final sentence)
  - Passage: "Richer schemes, stationary protocols and filtered recordings remain to be tested."
  - Why: Exact duplicate in function of the protected Discussion limitations (lines 1599-1607: "Richer schemes, the stationary regime, the few-channel boundary, the filtered kernel and experimental data are named companions rather than omissions"), which is the survivor; the intro then closes on the stronger "verified rather than assumed".
  - Cross-references: None found: nothing cross-references this sentence, and the limitations subsection is protected and carries the fuller list.

Section editor's note: Answer to the span question (which explanation survives): the dedicated section pp.3-10 survives as the sole full explanation of both the likelihood family and the diagnostics - its diagnostics definitions and Table 1 are protected - and the Introduction is reduced to claim-level pointers; the items on lines 72-75, 80-85 and 86-97 implement that direction. Span arithmetic: 1,749 words by wc includes the title block and page footers (~60 words); prose is ~1,690, of which the Abstract (297) and the three motivation paragraphs at lines 40-68 (~460) are essentially untouchable - the Abstract by the revised-for-accuracy constraint (only an 8-word appositive comes out; every other sentence is a claim), the motivation because it is the paper's unique hook (flip mark, MacroIR/MacroINR reversal with the factor of six, N*i identifiability, Clerx/IonBench survey - none of it duplicated elsewhere at this level). The 402 proposed words are therefore ~28% of the Introduction proper (lines 39-139, ~1,390 words), concentrated where the intro triple-states things: the results-preview paragraph (lines 111-123, the paper's third statement of its findings after the abstract and before the Results, -125) and the takeaway paragraph (lines 124-139, -114). Two discoveries for the duplication ledger: (1) the intro holds a FOURTH copy of the tau_int gloss (c) at lines 124-126, beyond the three the brief lists - survivors are the p.10 definition and the Figure 5 caption; (2) lines 130-134 vs Methods lines 2137-2145 are an unlisted near-verbatim pair (templated C++ type / Fisher wired by hand / no second pass) - Methods survives. Items for the coherence pass to verify against other editors' spans: intro line 102 becomes the surviving copy of the isolates-the-approximation's-own-error claim once the abstract appositive is cut; Results line 1411-1413 becomes the sole full in-text statement of the no-region ordering (abstract still carries the punchline); the analogue-filter-pickup claim at line 127 is unique in the whole paper and must survive the lines 124-128 compression; and the scope-disclaimer sentence at lines 95-97 (no evidence computed, no NSFA benchmark, Stepanyuk 2014) is unique and kept verbatim.

### Appendices 1-4

Span: pp. 39-55, lines 2205-3138. Section size about 9516 words. Items: 14.

- [ ] **DELETE**, about 170 words, lines 2302-2312 (p. 40, end of 'The interval average, and how the members are built')
  - Passage: "The interval average is not a function of the occupancy at any one time, so no conductance vector on the K states delivers its mean. On the K2 pairs ... carried out in three steps. Run the update on the pairs ... the Gaussian replacement is where the approximation enters."
  - Why: This is the second Appendix-1 copy of duplication (b) plus a full preview of the three-step construction that the subsection 'The three steps that build the interval update' (2456-2467) restates verbatim in function, including the 'legitimate rather than convenient' justification (repeated at 2575-2578) and 'the Gaussian replacement is where the approximation enters' (repeated at 2584-2586); surviving copies: body p. 7 (lines 348-352) and Appendix 1 lines 2456-2467.
  - Cross-references: The next subsection (2318-2321) redefines gamma-bar_i->j and v-bar_i->j in its own words, so no symbol is orphaned; coherence pass should confirm nothing between 2312 and 2456 refers back to 'the three steps' (nothing does) and that body lines 348-352 survive the body-span edits.

- [ ] **COMPRESS**, about 40 words, lines 2468-2473 (p. 43)
  - Passage: "Conditioning on the pair rather than on one end is the same device as static condensation ... the word transition is not used for it anywhere in this paper."
  - Why: Duplication (f): the surviving copy is body p. 7 (line 326-327), where the reader first meets the chaining point; the appendix keeps only what is new here.
  - Do: Residue must retain the static-condensation / spatial-Markov analogy (unique to this passage) and the 'not a transition state in the mechanistic sense' disclaimer; drop only the clause 'chaining intervals leaves a single state variable per junction, since the end state of one interval is the start state of the next' (duplication (f)).
  - Cross-references: Body line 326 must survive in the body span; 'static condensation' appears nowhere else, so the analogy must not be lost in compression.

- [ ] **COMPRESS**, about 60 words, lines 2240-2264 (pp. 39-40, acquisition-kernel subsection)
  - Passage: "Four properties of that table matter more than the numbers in it. From fcD = 1 upward ... The two concessions are one concession."
  - Why: Only the rhetorical scaffolding ('Four properties ... matter more than the numbers', 'The two concessions are one concession') is removable; the numbers themselves are what body line 181 ('Appendix 1 bounds it') leans on.
  - Do: Residue must retain all four numeric facts: the 0.147/(fcD) and 0.075/(fcD) laws with lags >= 2 zero; the variance-sum conservation forced by unit DC gain (misallocation, not loss); pole-count insensitivity and the Gaussian-equivalent width convention with the Colquhoun and Sigworth (1995) citation; the first-order (white) vs second-order (gating) deficit with the factor-27 example and where the residue lands.
  - Cross-references: Body lines 181-182 (bound plus block-average remedy) depend on the table and the four facts; the block-averaging paragraph 2269-2274 backing the 'average the record first' scoping is deliberately untouched.

- [ ] **COMPRESS**, about 85 words, lines 2364-2386 (pp. 41-42, Van Loan alternative)
  - Passage: "An equivalent evaluation avoids diagonalizing Q altogether. F1 and F2 are the strictly upper blocks ... Every number reported in this paper comes from the spectral route."
  - Why: The alternative route needs stating once with its citation and the reason it was not taken, not a 67-word tour of the 3K-block derivative payload.
  - Do: Residue must retain Eq. A7 with the Van Loan (1978) citation, one sentence that it returns the transition matrix and both moments with no spectral machinery and is the construction to use where no reliable eigensolver exists, and one sentence that the spectral route was chosen because it carries the parameter derivatives in closed form by first-order perturbation of the same decomposition, and that every reported number comes from it.
  - Cross-references: Discussion line 1515 leans on 'the derivative of the spectral reconstruction inside it (Appendix 1)': the derivative rationale must survive in the residue.

- [ ] **COMPRESS**, about 80 words, lines 2596-2620 (p. 45, 'The start-conditioned members, stated forward')
  - Passage: "The construction above reaches IR by choosing the boundary pair as the state. The members that stop at the interval's start are usually presented as ... reached without ever forming a pair."
  - Why: The 'usually presented as abbreviations / they are not / worth stating plainly' framing (roughly 4 sentences) restates the subsection's purpose rather than adding content.
  - Do: Residue must retain Eq. A22 with the statement that everything MR does follows from it and needs no object larger than the K states, Eq. A23 identified with A9/A11, and one sentence that both moments are row sums of the tilted-semigroup derivatives (integrals of K-vectors), so MR is reached without ever forming a pair.
  - Cross-references: The falsity paragraph 2621-2628 ('a correctly derived filter for a measurement model that is false', unique at line 2627, with the pointer to Eq. A28) must stay untouched.

- [ ] **COMPRESS**, about 115 words, lines 2679-2709 (p. 46)
  - Passage: "The variance: both members are exact, and their slots differ. By the law of total variance ... which it stops being at the first update."
  - Why: The passage states the MR/IR variance identity a third time ('This is Eq. A24 read from the outside ... which is why the totals coincide') on top of A24 and A37; the surviving anchor is Eq. A24 (line 2646).
  - Do: Residue must retain Eq. A26 with its two underbraces identified once with the two members' slots (total form in w for MR, boundary tilde of A18 for IR), the bold exactness claim (both report the exact variance given mu0, Sigma0, so their agreement measures nothing), and a one-sentence version of the caveat that the conditioning variable is the count vector and a single channel's X0 gives the same expressions only until the first update.
  - Cross-references: Methods Table 2 caption (line 1841) cites Eq. A24 for the state-dependent noise, so A24 and its underbraced equation must remain intact; Results line 692 points to Appendix 2 for the gain and is unaffected.

- [ ] **COMPRESS**, about 85 words, lines 2745-2753 (p. 47)
  - Passage: "Why the same quantity is invisible in one block and decisive in the other. D enters the variance squared and the covariance once ... rather than a defect of the arithmetic that implements it."
  - Why: Pure pedagogy re-deriving in words what the two preceding subsections (A26 and A27/A28) each already state; the phrase and idea appear nowhere else, so one compressed statement suffices.
  - Do: Drop the subsection header and fold two sentences onto the end of the A27/A28 subsection: D enters the variance squared and the covariance once; a second moment can sit in either variance slot (Eq. A24) but a collapsed first moment has nowhere to go, so the D-squared-versus-D asymmetry is the whole MR-to-IR difference and is a consequence of Eq. A22.
  - Cross-references: NONE (grep confirms 'invisible' and 'column profile' occur only here).

- [ ] **COMPRESS**, about 45 words, lines 2756-2762 (p. 47, opening of the linear-filtering-frame subsection)
  - Passage: "Everything above was Gaussian conditioning: a joint law over (state, recorded value) ... the equations that make it checkable are here."
  - Why: The 'forty years of work / take it on faith' ornament restates the Discussion's concession (lines 1507-1509), which survives there.
  - Do: Residue must retain: the derivation above is Gaussian conditioning and borrows nothing from linear filtering; the correspondence with that literature is nevertheless exact and is conceded in the Discussion; the equations that make it checkable follow.
  - Cross-references: Discussion lines 1507-1512 (the Zadrozny concession and the specificity case) must survive in the Discussion span; the sentence 'The relation is conceded in the Discussion' keeps that link.

- [ ] **COMPRESS**, about 90 words, lines 2811-2823 (p. 48, the 'Three things follow' paragraph)
  - Passage: "The augmented-state filter and the boundary-state filter compute the same three blocks ... independent of the old given the occupancy."
  - Why: The content is essential (the reset is operationally load-bearing) but each point is stated with a preamble and an echo that compression removes.
  - Do: Residue must retain all three facts, one sentence each: same three blocks so same update and posterior; the K-squared pair is a device of the derivation, not information the filter carries (the K+1 augmented state delivers the same blocks, the conclusion A14/A15 reached); the per-state variance sits inside Var(z_Delta), which is why the cross block needs F1 only while the variance block needs F2; and the equivalence holds per window with a z = 0 reset because each window's integral is fresh given the occupancy.
  - Cross-references: The reset caveat must survive verbatim in substance; the 'device not information' point is cross-referenced to Eqs. A14/A15, which stay.

- [ ] **COMPRESS**, about 180 words, lines 2837-2855 (p. 48, closing paragraph of the linear-filtering-frame subsection)
  - Passage: "Where the frames genuinely differ is the model rather than the arithmetic. A Kalman derivation presumes a linear-Gaussian state process ... no comparably natural form to the other two members."
  - Why: This is the frame-philosophy passage the span guidance flags: beyond the equivalence statement and A29-A33 the 'two consequences' sentences duplicate points made elsewhere (the Sigma meaning-versus-convention point survives at Methods line 1816; the imported-measurement-equation point is already implicit in A33).
  - Do: Residue (about 95 words) must retain: a Kalman derivation presumes a linear-Gaussian state process, which multinomial occupancy counts with state-dependent transition noise are not; this derivation computes the exact first two moments of the true process and Gaussianizes only at the two named closures, so the routes agree because joint-Gaussian conditioning consumes only two moments; and one sentence on why the pair remains the presentation (linearity in the state is structural, the collapse hands over the next window's prior, the per-state variance slot makes R, MR and IR one construction).
  - Cross-references: Methods line 1816 ('the one moment at which Sigma genuinely is a single channel's covariance') must survive in the Methods span, since the compressed residue no longer carries that point.

- [ ] **COMPRESS**, about 105 words, lines 2892-2903 (p. 50, Appendix 2)
  - Passage: "The tilde operator. Over a bold symbol, the tilde marks a contraction taken on the K2 boundary pairs ... one power of gamma-bar against two."
  - Why: The scalar/natural-units parenthetical is a fourth touch of duplication (a) (surviving copy for my cuts: Methods lines 2089 and 2100-2101, which already carry the two-uses-never-meet disambiguation), and 'the K2 x K2 object is never built' is stated at A1 line 2531 and body line 328.
  - Do: Residue (about 60 words) must retain: tilde over a bold symbol marks a contraction on the K-squared boundary pairs rather than the K states; its two cases are written out in Appendix 1 as Eq. A18 (bilinear) and the gain of Eq. A19, both collapsing to K-by-K arithmetic; the G and H of Eq. A34 tell the two apart, one power of gamma-bar against two.
  - Cross-references: Coherence pass must ensure whichever natural-units tilde copy the Methods/body agents keep (body 236, Results 578, or Methods 2089) survives, since this passage no longer disambiguates the scalar use; Methods 2100-2101 cross-references these boundary-pair tildes and stays valid.
  - **Adjustment (binding)**: KEEP the scalar/natural-units disambiguation parenthetical (2894-2896): the Methods copy of the disambiguation is cut by item 76, making this the sole survivor. Cut only the K-squared-never-built and one-power-vs-two repetitions. Savings 105 to ~85w.

- [ ] **COMPRESS**, about 165 words, lines 2993-3018 (pp. 51-52, Appendix 2, the A37 paragraph and the Results-steps paragraph)
  - Passage: "One consequence of the two switches acting together is worth stating ... The two together are the R-to-IR gap."
  - Why: The on-a-recording caveat (3006-3011) and the 'Neither step is confined to what it flips' passage (3015-3018) repeat, nearly clause for clause, Appendix 1 lines 2666-2671, which survive as the single statement of that caveat (they also carry the VR content and the Results pointer); this resolves the flagged A24/A37 triplication with A24 as anchor, A37 as the switch-level algebra, and the Table 2 caption copy left to the Methods agent.
  - Do: Residue (about 175 words including the equation) must retain: Eq. A37 with the sentence identifying its two extra pieces as what switch 1 moves from the third term of A17 into the second; 'The decompositions differ, the sum does not', and that what changes is the gain, the whole MR-to-IR difference at a common prior; plus a two-sentence statement that MR-to-VR flips switch 3 and VR-to-IR restores the boundary terms, the two steps the Results use, together the R-to-IR gap.
  - Cross-references: Appendix 1 line 2665 '(Appendix 2)' and Methods line 1755 'rearrange into each other (Appendix 2)' both point here, so Eq. A37 and its identification sentence must stay; Appendix 1 lines 2666-2671 must not be cut by any other pass.

- [ ] **COMPRESS**, about 45 words, lines 3020-3027 (p. 52)
  - Passage: "Least squares, which is off the lattice. Classical nonlinear least squares has no setting of the three switches ... the anchor of the comparison rather than a member of the family."
  - Why: Body lines 299-309 already establish that the two least-squares arms sit off the lattice and discard the variance, so the appendix needs only the equation-level restatement.
  - Do: Residue (about 55 words) must retain: no setting of the three switches produces it; it propagates the mean occupancy alone, predicts from Eq. A16, and assigns one constant fitted variance, so there is no gain and no occupancy object is conditioned on the data; it is the anchor, not a member.
  - Cross-references: Body lines 299-309 and the Methods description of the LS arms must survive their spans; the Eq. A16 pointer stays.

- [ ] **COMPRESS**, about 30 words, lines 3064-3067 (p. 53, Appendix 3)
  - Passage: "Both are first order, and the manuscript does not lean on that. The sandwich is verified directly against the empirical covariance ... an expansion cannot certify on its own."
  - Why: Near-verbatim duplicate of body lines 533-534 ('are first-order approximations (Appendix 3). We verify it directly, without leaning ...'), which is the surviving copy.
  - Do: Residue: one clause noting both corrections are first order and the sandwich is verified directly against the empirical covariance and the coverage of Figure 2.
  - Cross-references: Body lines 533-534 must survive the body-span edits; the Figure 2 pointer is kept in the residue.

Section editor's note: Span = Appendices 1-4 (lines 2205-3138, pp. 39-55), ~9,500 words including equations and page footers. Savings ~1,295 words (~14%), deliberately modest since eLife tolerates appendices. Answers to the flagged questions. (1) Linear-filtering frame (2755-2855, 1,175 words): beyond the equivalence statement and Eqs. A29-A33, what earns its keep is the augmentation definition with the Zadrozny citation (2778-2789), the A30-A32 identifications including the F1 check of A32's second term (2790-2810), the per-window z=0 reset, and the A33 MR-misspecification result (2824-2831); the intro ornament, the 'three things' expansion, and the closing frame-philosophy paragraph compress, netting ~315 words. (2) Duplication (b): surviving copies are body p. 7 (lines 348-352) and Appendix 1 lines 2456-2467; the copy at 2302-2312 is deleted. Duplication (f): surviving copy is body p. 7 (line 326); the appendix clause is cut inside the compressed static-condensation sentence. Duplication (a): this span touches it only via Appendix 2's parenthetical (2894-2896), cut on the assumption the Methods copy (2089, cross-noted at 2100-2101) survives. (3) MR/IR variance identity triplication: Eq. A24 (line 2646) survives as anchor (cited by Methods Table 2 caption, line 1841); Eq. A37 survives as the switch-level algebra (pointed to by A1 line 2665 and Methods line 1755); what is cut is the commentary repeating the 'maps from a prior, not recordings' caveat (kept exactly once, at A1 2666-2671) and the A26 subsection's restatement clause. Deliberately untouched as load-bearing or protected: the three-marks notation box (2215-2226); the 0/0 shrinkage safeguard (2387-2404); the exactness/two-closures mapping (2567-2594); the verification tolerances (2857-2868); Appendix 2's table and caption (2971-2987); the Switch-1 trailing-propagator passage with its observed-symptom disclosure (2923-2941); all of Appendix 3's m/a derivation and sampling-floor disclosure (3092-3112); and Appendix 4 in full (target of the Discussion pointer at line 1562). No MOVE items: everything here is already in appendix position, and nothing merits demotion to the repository README.

## After this pass

Re-measure before deciding on more. If 55 pages is still too long, the next tier is structural and riskier: fold what remains of the family section (about 3,700 words after these cuts) to about two pages by pushing the discussion around Eq. 2 into Appendix 1 and the equal-variance derivation into Methods. That would bring the body under 30 pages, at the cost of making Results lean harder on the appendices. Do it only if the journal asks.

Verified total saving after adjustments and rejections: about 7749 words.
---

## Execution log, 2026-08-26

Passes land one section at a time, each with its own commit and a `check.sh` run, in the order the
plan proposes. Word counts below are `check.sh`'s whole-document figure, which is the one comparable
across passes; `wordcount.py`'s main-text figure and the plan's `pdftotext` figure are different
definitions and must not be mixed with it.

| pass | section | items | document words after | commit |
|---|---|---|---|---|
| 1 | Materials and methods, Data availability | 26 | 31,149 (from 33,285) | `8883fbd` |
| 2 | captions and figure-supplement legends | 19 | 30,242 | `f63671a` |
| 3 | Discussion | 10 | 29,753 | `95d2dd1` |
| 4 | Results | 22 | 29,356 | `ecd9617` |
| 6 | Appendices 1 to 3 | 14 | 27,789 | `9d48a22` |
| 5 | the likelihood family, Abstract and Introduction | 26 | 27,910 | `769a8a7` |

**Result.** Whole document 33,285 -> 27,910 words by `check.sh`; main text 13,088 -> 11,313 by
`wordcount.py`, which takes it from 1.14x Munch (the longest comparable article in the 413-article
survey) to 0.98x, and from 1.44x the p95 to 1.25x. Rendered: **69 pages -> 60**, with zero overfull
vertical boxes, where the plan projected ~55. The gap to 55 is the two page-saving items that are
figure work rather than text: merging Figure 3's supplements 1 and 2, and the reflow that follows.
`check.sh` ends at 9 pass, 2 fail (both pre-existing: LINT-SRC and the index), 1 warn (length,
pending D-1).

**What passes 5 and 6 owed and paid.** The family section's MOVE of the equal-variance derivation
had no destination inside that agent's remit and was left preserved-but-commented at the removal
site; it is now landed in Methods' natural-units passage. Had it stayed commented, Figure 5's ratio
axis would have lost the equation it is a consequence of.

**A defect the pass created and a check that now catches it.** Compressing prose by commenting it
out swallowed the opening words of three sentences, which compile silently and print as sentences
with no subject. The same defect already existed twice in HEAD, once from 2026-08-13 and once from
the previous day's correction pass, for five in total. `papers/_program/swallowed_sentences.py` finds
all five signatures and runs as `check.sh` item 12; all five are repaired and the check is green.

**Destinations created by pass 1.** `supplementary_file_1.tex` gains S1.6 (the members as
dispatched: flags, data keys, cells, and the deposit verification of INR = MacroINR), S1.7 (symbols
against implementation names, and the fixed build flags) and S1.8 (the design rationale behind the
three numerical safeguards). `projects/eLife_2025/data_note.md` is new and holds what belongs to the
deposit rather than to the paper: the stale parameter tables, the inert substep argument, the three
data directories and their search order, the parameter-index dictionaries, the bootstrap internals,
the run-ledger provenance and the cross-language check tool.

**Deferred, because it is figure work and not text.** Merging Figure 3's supplements 1 and 2 into a
single four-band page: it needs the producing `.Rmd` rewritten and the page re-rendered, then the
renumbering the item lists. Nothing else in the plan depends on it.

**Executed against the plan's own body text, with reasons recorded at each site.** The INR =
MacroINR identification stayed in Methods rather than travelling with the moved table; Methods item 8
kept the clause the appendix pass depends on; the washout mechanism sentence stayed in the Results;
the m/a prose gloss was not cut, per reconciliation ruling 2.
