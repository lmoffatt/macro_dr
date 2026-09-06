# Audit of the evidence computation (thermo_evidence_dts), 2026-09-01

Request: check whether the evidence machinery still works (unused since
~Dec 2025). Method: read from the CLI inward, no compiling. Real usage
reference: `projects/p2x2/ops/MacroIR/evidence.macroir`.

## The chain

CLI: `make_dts_compiler()` (merged at command_manager.cpp:450) registers
`set_ThermoAlgorithm_dts`, `thermo_evidence_dts`, `thermo_evidence_dts_2`,
`thermo_evidence_dts_continuation`, `_continuation_2`
(legacy/CLI_thermo_evidence_dts.h). The likelihood setter is the LEGACY
`set_Likelihood_algorithm` (CLI_macro_dr_base.h:191, registered in
make_model_compiler), NOT the modern `build_likelihood_function` of the
eLife roster. Core: `thermo_evidence<Adapt_beta=true>` →
`thermo_evidence_loop` (parallel_tempering_linear_regression.h:685):
report → adapt_beta → adjust_beta → reset_statistics → step_stretch →
thermo_jump. Evidence emitted by `save_Evidence::report`
(parallel_tempering_linear_regression.h:106-163) into `*__i_iter.csv`.

## Short verdict

The structure is intact and compiles (the variance_form refactors dragged
it along: d4a45340, 7e1549c5). The exact combo of the P2X2 script is
still the only one accepted and passes the domains. BUT: the dts
continuation is broken by construction, there are two asserts with an
assignment, one UB from a string typo, and the likelihood underneath
(qmodel.h) has 37 commits since December, so the numbers will NOT
reproduce December's.

## Findings (by severity)

1. **dts continuation BROKEN by domains.** `calc_thermo_evidence_dts`
   instantiates domains (adaptive {false}, variance {false}, taylor
   {true}); `_continuation` and `_continuation_2` instantiate (adaptive
   {true}, variance {true}, taylor {true})
   (CLI_thermo_evidence_dts.h:123-131 vs 381-389, 517-525). The
   continuation reloads the ORIGINAL flags from
   `thermo_evidence_dts_<id>.txt` → `to_variant(false)` is not in {true}
   → cryptic stderr error ("0 is not in ...") and no run. It can only
   continue dts_2 runs. Fix: copy the dts domains into `_continuation`
   (or unify into a single {false,true} set).

2. **One single likelihood combo per command** (singleton domains +
   to_variant fails at runtime, variables.h:360-377). dts accepts ONLY
   (adaptive=0, recursive=1, averaging=2, variance=0,
   variance_correction_approximation(taylor)=1), with family_macro,
   variance_form=variance_total and default qdt_method implicit
   (qmodel.h:8300-8342). **taylor=1 FORCED**: the evidence runs with the
   Taylor correction no matter what; the paper's IR is built through a
   different path (build_likelihood_function_with_family). Decide whether
   that is what is wanted for CCO/COC before running.

3. **Assert with assignment**: `assert(beta[beta.size()-1] = 1)` in
   adapt_beta (parallel_tempering.h:1551) and adjust_beta (:1518).
   Release: no-op; Debug: silently forces beta.back()=1. Should be `==`.

4. **Equalizer typo = UB.** If `adapt_beta_equalizer` is not one of the
   recognized strings (parallel_tempering.h:1446-1504),
   calculate_controler_step returns `{}` and the "s" branch reads `d[i]`
   out of range (:1567). No validation at construction. The good value in
   the reference script: "deltaBeta_deltaL_vfm", controler "s".

5. **6-argument calcEvidence with a wrong sign**
   (parallel_tempering.h:422-428): with dL=0 it does not reduce to the
   trapezoid (an extra 2·b1·L2 survives; the first term should be
   `L1*b2 − L2*b1`). TODAY it is dead code for dts (only report_old uses
   it), but it is a landmine if anyone revives the "variance-corrected"
   quadrature.

6. **The live quadrature omits the [0, β_min] tail** (report:141 `if
   (beta0 > 0)`), with β_min ≈ 1/|E₀[logL]| (initial_beta_dts,
   parallel_tempering.h:604-617). Omitted tail = O(1) nat and DIFFERENT
   per model → O(1) bias in Bayes factors. Consistent with the null-test
   prediction "Δ log Z = O(1), does not grow", but do NOT expect exact 0
   even with pushforward priors.

7. **Two evidence columns** in `*__i_iter.csv`: `log_Evidence` uses the
   instantaneous mean over that iteration's walkers
   (mean_logL(thermo_mcmc), parallel_tempering.h:314); the
   `mean_log_Evidence` columns use statistics accumulated since the last
   reset (reset every adapt_beta_every). A single log_Evidence row is
   noise: average over iterations or use mean_log_Evidence.

8. **scheme_10_inact is no longer registered.** Current library =
   scheme_CO, scheme_CCO, scheme_COC, scheme_1, scheme_10_d
   (models_used.h:16). The old P2X2 evidence.macroir fails at get_model.
   For CCO/COC we are fine. The old script's paths
   (data/experiments/..., data/models_Ag/...) do not exist in the current
   layout either.

9. **Continuation, two more fragilities**: (a) extract_parameters_last
   restores betas and ladders that GREW, but if the run ended with fewer
   rungs than beta_size (adjust_beta removes), extract_iter fails and the
   "continuation" silently starts from a prior sample with beta=all-zeros
   → NaNs in adapt_beta (parallel_tempering.h:1875-1925). (b) the
   Maybe_error of the logL recompute for the loaded walkers is discarded
   (parallel_tempering_linear_regression.h:832).

10. Minor: `number_trials_until_give_up` and
    `save_every_param_size_factor` are requested in the DSL and go
    nowhere; the legacy setter PERMUTES variance_approximation ↔
    variance_correction_approximation between signature and tuple
    (CLI_macro_dr_base.h:74-80) and the consumer compensates, so the DSL
    NAMES are what counts: variance_correction_approximation → taylor,
    variance_approximation → variance; init_seed=0 = random
    (random_device), fix the seed for the null test.

## What was NOT verified

- Numbers against December: qmodel.h has 37 commits since 2025-12
  (canaries, shrinkage at Qdt sources, "restore interval conductance
  variance in the non-recursive path", SymmetricMatrix full-storage...).
  The likelihood changed ⇒ the evidence changes even though dts is
  intact. Before trusting: a short smoke run with scheme_CO and, if any
  December outputs were kept, compare log_Evidence.
- The inside of step_stretch/thermo_jump at the β=0 rung.
- Verified positive: the sampler is stretch-move, it does NOT use
  Lapack_SymmPosDef_inv (which is broken); signatures of
  new_thermo_Model_by_max_iter_dts (qmodel.h:9605) and of the two
  function table makers are OK.

## dts vs dts_2, settled 2026-09-02

The one in use (and the one that survives) is **dts**. dts_2 differs in
three things and all three are now obsolete:

1. Domains: dts_2 pins adaptive {true}, variance {true}. Adaptive is
   DEPRECATED (Luciano: its purpose was to prevent negative Pmean from
   the recursive update, now solved by another mechanism, the
   canaries/simplex machinery). The modern builder
   (include/macrodr/cmd/likelihood.h:48,72) already pins adaptive {false}
   in every branch.
2. Memoization: dts memoizes Calc_Qdt_step, dts_2 does not. The obvious
   suspicion (key = (Agonist_step, fs) WITHOUT the parameters, plus
   ff[i_th] shared across walkers) was chased to the bottom and is
   UNFOUNDED: each log_Likelihood does `f_local = f.create("_lik")`
   (qmodel.h:6781) and create() yields empty memoizers → cache lives
   within one evaluation, correct and faster. Intent documented in commit
   b0724eca ("caches within a log_Likelihood call").
3. The continuations carry the dts_2 domains → today they would only
   continue dts_2 runs.

**Orphan combo**: dts forces (recursive=1, av=2, variance=0, taylor=1).
The modern builder does not offer variance=0+taylor=1 (non-taylor branch:
taylor {false}; taylor branch: variance {true}, av {1,2}). The December
evidence ran a likelihood the canonical roster no longer contains.
(Resolved below: the combo was not an orphan member, it was IRTV
mislabeled, since variance is dead weight under taylor=true.)

Three-level plan: A = launch with dts as-is (equalizer spelled right,
fixed seed, no continuations). B = surgical cleanup: delete dts_2 and
continuation_2, copy the dts domains into the remaining continuation,
`==` in the asserts, validate strings with Maybe_error. C =
consolidation: evidence on the modern likelihood_algorithm_type, killing
the legacy 6-tuple and the permuted setter.

## taylor (variance correction, MRT/IRT) also dead, 2026-09-02

Luciano's verdict, 2026-09-02, verbatim (translated): IRT "is an
algorithm that never worked well; it is implemented, but every time I
tested it, it gave wrong results" ("nunca funcionó bien, está
implementado pero las veces que lo testeé me dio mal"). So NOT just an
abandoned branch: a known-defective member. Convergent with four
independent sources:

1. The paper excludes it explicitly: supplementary_file_1.tex:310 ("the
   Taylor-corrected variants MRT and IRT are not part of the study",
   taylor off in every reported run; and "variance approximation ON").
2. papers/_program/axes.md:72: canonical roster = NR, INR, R, MR, IR with
   taylor=false.
3. Production: taylor=true ran ONLY in the macro_IRTV lane of the April
   2026 sweeps; absent from the last ~60 runs.
4. The Bessel plan does not inherit it (macroir_bessel_plan; memory
   archive: Taylor/IRT closed).

Live code that FORCES it to true: ALL FOUR legacy evidence commands
(CLI_thermo_evidence.h:217,347; CLI_thermo_evidence_dts.h;
CLI_thermo_evidence_fraction_dts.h:143,283;
CLI_thermo_levenberg_evidence.h:124,255, note the Levenberg that is
slated for reactivation too) + CLI_function_table.h:788,909 + the taylor
branch of the modern builder (likelihood.cpp:1942). In other words: **the
entire evidence surface is built exclusively on the dead branch**, and on
top of that with variance=0, while the paper runs variance=1. The
December evidence used a doubly non-canonical member.

HOMONYMY (Luciano 2026-09-02): there are TWO distinct taylors and only
one dies.

- **Taylor of the LOGLIKELIHOOD** = `uses_taylor_variance_correction_
  aproximation`. Gates the σ correction terms in y_mean/y_var
  (qmodel.h:3764-3799: sSg, sSs, gvar_i, delta_emu...) and the
  MacroR2/MRT/IRT rank-2 Woodbury recursion (~4700-5100). THIS one is
  dead (tested, gave wrong results). It is the one the evidence commands
  force to true, and the only one the retarget touches.
- **Taylor of Qdt** = `calc_Qdt_taylor` (qmodel.h:2271): numerical method
  for the propagator/moments. ALIVE and load-bearing: universal fallback
  when the eigendecomposition fails (calc_Qdt:3165) and the derivative
  path delegate of calc_Qdt_schur (:2502-2512). Orthogonal to the
  likelihood member.

Naming trap feeding the homonymy: the modern `qdt_method` axis (0=eig,
2=schur) is labeled "taylor_qdt" in the model report (qmodel.h:8289) and
the script argument `taylor_qdt_approximation=true` means SCHUR, not
Taylor (likelihood.h:35-38). Three names carrying "taylor" for two
mechanisms, and one of the three is not even Taylor.

DECISION (Luciano 2026-09-02): taylor_variance_correction = zombie code,
keep it so it can eventually be revived to study whether the theory is
valid and why it fails, but off in all upcoming experiments. Taylor of
Qdt = valid option (fallback when eig fails, probably less efficient),
alive.

RETARGET DONE 2026-09-02 (not compiled; Luciano compiles):
- CLI_thermo_evidence_dts.h: dts and BOTH continuations →
  (adaptive {false}, recursive {true}, averaging {2}, variance {true},
  taylor {false}) = the paper's canonical IR. This also FIXES the broken
  continuation (it now pairs with dts). dts_2 keeps its condemned
  adaptive {true} (deletion candidate, untouched).
- parallel_tempering.h: both asserts `=` → `==` (adapt_beta,
  adjust_beta); `d.empty()` guard in adapt_beta, so an equalizer typo now
  warns on stderr and leaves the ladder alone instead of UB.
- SCRIPTS: flags INVERTED relative to the old P2X2 script:
  variance_approximation=1, variance_correction_approximation=0.
  Everything else unchanged (adaptive 0, recursive 1, averaging 2).
- PENDING in the same vein: the Levenberg
  (CLI_thermo_levenberg_evidence.h:124,255) and thermo_evidence /
  fraction_dts still force taylor=true; same retarget if reactivated.
  Validating equalizer strings in the setter with Maybe_error was not
  done (the runtime guard suffices for now).
- The December outputs stop being a regression reference (different
  member + qmodel drift). Validation of the new pipeline: the twin null
  test + the exact conjugate case ALREADY in the code
  (bayesian_linear_regression.h, evidence(conjugate{},...) used by
  report_model).
- With taylor dead and adaptive dead, the modern builder collapses to the
  non-taylor macro branch + micro + nonlinearsqr → plan C (evidence on
  likelihood_algorithm_type) got smaller than it was.
- The MacroR2/MRT/IRT branch of qmodel.h becomes an archiving candidate;
  do not start yet (deprecate in stages).

## "Says vs does" investigation of taylor_variance_correction, HEAD 2026-09-02

Trigger: Luciano recalled a discrepancy between what the correction says
it does and what it does. The archived memory (May-June 2026) already
documented three bugs (missing ½ factor in the σ² direction; β using
V_pred where the derivation requires V_obs; the supplement §6
substitution (gS+vS)·Σ_post mathematically wrong, verified numerically 4×
too small) and the PSD issue. State at today's HEAD:

**Fixed and consistent:**
- The three May fixes ARE in the code (Newton loop with self-consistent V
  at the iterate, ½ factors present, the corrected form
  δ/(2V)·[(gS+vS) − sm·(b′+vSv)·vS]).
- The supplement (theory/macroir/docs/Macro_IRT/macroirt_supplement.tex)
  WAS corrected: eq. mu_post_IRT (line ~370) is the correct form and the
  pseudo-code (~551) warns "the SM expansion is explicit, NOT
  (gS+vS)·Sigma_post". The §6 misdirection is gone.

**LIVE discrepancies (says vs does), verified at HEAD:**
1. **The variance flag (V) is DEAD LETTER under taylor=true.** In the
   observation block, variance_correction wins before variance::value is
   consulted (qmodel.h:3764 vs 3849/4029/4060); in the whole rank-2
   update branch (4675-5120) variance::value does not appear ONCE.
   ⇒ IRT ≡ IRTV exactly. This corrects the "orphan combo" above: the
   December evidence (variance=0, taylor=1) was not an orphan member, it
   was IRTV mislabeled.
2. **taylor=true silently inert** at av=0 (the branch requires
   averaging>0, qmodel.h:4675 → falls through to the standard Kalman with
   no notice) and at recursive=false (safely_calculate_Algo_State:5969,
   "Non-recursive path doesn't apply variance_correction", commented but
   no warning).
3. **Stale header comment**: qmodel.h:4658 says "rank-1 quasi-Laplace...
   Newton step in direction (γ̃ᵀΣ+ṽᵀΣ)" (the old scheme); the code is the
   exact rank-2 with the 2×2 Woodbury.
4. **The det guard admits indefinite K**: |det|≤1e-30 rejects only
   near-singular; det<0 (m22 = c − 2V²/N < 0 at small N) passes silently
   (4725-27, 4789-91). And α_σ is still DIAGONAL (4810/4817), which the
   June audit measured as false-passing ~100% of unsafe steps. The guard
   claims to protect PSD and does not.
5. **Zero tests**: no unit test exercises MRT/IRT; the σ̄²=0 → IR collapse
   test recommended in May was never written.

**History of the empirical verdict**: Luciano's "IRT gave wrong results"
matches what was measured 2026-05-13 (logL systematically worse than IR,
Information Distortion 10-30×, posterior covariance massively
underestimated) on the pre-fix code; IRT's earlier apparent advantage was
an artifact of the gvar_i bug ("it was never the T"). Post-fix IRT was
NEVER empirically re-validated (the IRTV sweeps are from April, pre-fix).
So the zombie's state is: fixes applied + doc corrected + no validation +
incomplete PSD guard.

**Revival checklist**: (a) det-sign/eigen check on the 2×2 K instead of
the diagonal α_σ (or Joseph form); (b) unit test of the σ̄²=0 → IR
collapse (one step, catches all three bugs at once); (c) decide the
variance-flag semantics under taylor (ignored today) and make the inert
combos loud (av=0, recursive=false); (d) re-validate on figure_2
post-fix; (e) update the comment at 4658.

## Telescopic evidence (stepping stone), 2026-09-02

Source: ~/Projects/"Intercambio de replicas"/ (August 2026 talk;
derivation in gemini/Z cociente.md; inventory in material_establecido.md,
which already anchors the prior art: stepping stone, harmonic mean of
Newton & Raftery 1994, Neal). The identity: ln Z = Σₖ ln E_{p_βk}[L^{Δβk}],
exact for any finite grid; the naive estimator and the harmonic mean are
the K=1 case in the two directions; a factor's variance is finite iff
Z_{β+2Δβ} is finite (which is why the harmonic mean blows up and the
stepping stone does not). Literature names: stepping-stone sampling (Xie,
Lewis, Fan, Kuo, Chen 2011 Syst Biol); the sequential version is annealed
importance sampling (Neal 2001).

**Fit with thermo_evidence_dts: it comes FREE from the runs:**

1. `*__i_beta__i_walker.csv` (save_likelihood) already emits (iter, beta,
   i_walker, logLik) per rung. The telescopic estimator is computable
   OFFLINE: per adaptation window (ladder constant between
   adapt_beta_every), factor_k = logmeanexp(Δβₖ·logLik at rung k),
   ln Ẑ = Σ factors. Zero C++ changes; one R script.
2. **Fixes the [0, β_min] tail**: the first factor E_prior[L^{β_min}]
   uses the β=0 rung the dts ladder ALREADY carries and evaluates. The
   trapezoid's O(1) tail-omission bias disappears.
3. **No discretization bias**: the identity is exact on the grid; what
   remains is Jensen ~1/(2·ESS) per segment + MCMC correlation (batch
   means). The trapezoid stays as the SECOND estimator of the same run →
   the evidence-method comparison wanted this week, with no extra runs.
4. The adaptive dts ladder equalizes dβ·dlogL, which is exactly what
   controls each factor's ESS → the ladder is already near optimal for
   stepping stone. Per-segment diagnostic: ESS of w = exp(Δβ·(logL−max)).
5. Relatives: cuevi = the same telescoping over data fractions (the
   machinery exists: CLI_thermo_evidence_fraction_dts, taylor retarget
   pending). Optional future: swap statistics → Bennett/BAR, would
   require logging the proposed ΔlogL, today only counts exist.

R-script care: log-sum-exp always; group by window between adaptations
(the ladder CHANGES with iter, use the csv's beta column, never assume it
fixed); ESS per factor as the validity criterion; Jensen bias estimable
as −V̂ar(factor)/2 per segment.

**IMPLEMENTED ONLINE 2026-09-02** (Luciano's idea: have save_Evidence
compute trapezoid AND telescopic): save_Evidence now carries streaming
log-sum-exp accumulators PER SEGMENT, in TWO lifetimes (mirroring the
log_Evidence vs mean_log_Evidence convention) and TWO directions:

- instantaneous per event (plog_Evidence_ss / log_Evidence_ss and _dn):
  only THIS report's walkers. Luciano's question "is it pointless? the
  rare ones carry the weight": correct, with weights exp(Δβ·logL) the
  rare high-logL draws dominate, so the instantaneous estimate from ~32
  samples is biased LOW precisely by missing them → NOT an estimator, a
  DIAGNOSTIC: its series shows burn-in/drift, and the gap
  mean(instantaneous) vs window is the Jensen bias (the weight of the
  rare draws) made visible. With a well-adapted ladder
  (Δβ·sd(logL)~O(1)) the gap must be small; large = unresolved segment.
- window (mean_plog_Evidence_ss / mean_log_Evidence_ss and _dn):
  accumulated across report events, automatic reset when the ladder
  changes. THIS is the estimator. Alongside: ss_ess_up/dn (weight ESS per
  segment) and ss_count.

The cumulative INCLUDES the [0, β_min] segment (fixes the trapezoid's
tail). GUARD against "looks converged": with heavy tails the ESS can lie
too; the up-down bracket is the honest check because the two directions
have OPPOSITE tail sensitivity: if it does not close, the segment is not
resolved even if the ESS looks healthy. Negligible cost. Not compiled
(Luciano compiles). Assumes the ascending (dts) ladder. The offline R
script remains for old outputs and batch means.

Note from the talk material: "[descartado] ... en el paper solo usamos
stepping stone" (only stepping stone is used in the paper), Luciano's
prior decision to make stepping stone primary over the trapezoid;
consistent with making the telescopic column the main estimator and TI
the cross-check.

## DEO (non-reversible exchange, Syed et al. 2022) on macro_dr, analysis 2026-09-02

Luciano's question: can the deterministic even/odd scheme (Okabe 2001;
Syed, Bouchard-Côté, Deligiannidis, Doucet 2022 JRSS-B) be applied easily
to MacroIR? Answer: yes, ~15 lines in ONE function, with one design
decision due to the ensemble.

**How it swaps today** (`thermo_jump_mcmc`, parallel_tempering.h): every
thermo_jumps_every iterations it attempts ALL adjacent pairs (ib, ib+1)
at once; within each pair, n_walkers/2 random matchings (disjoint shuffle
halves: "lower" role = shuffled[0..n/2), "upper" role = shuffled[n/2..n),
hence no data race between pairs sharing a rung). The kernel is a product
of Metropolis moves on disjoint components = reversible; each walker's
temperature index (tracked by i_walkers/id_walker) diffuses → round trips
O(N²).

**The DEO change**: (1) parity = (jump counter) % 2; attempt only pairs
with ib%2==parity; (2) match ALL walkers one-to-one between the two rungs
(not n/2): if only half participate, directional persistence breaks half
the time and the lift is diluted. Each parity kernel remains reversible,
the alternating composition is not, π invariant, marginals identical →
NOTHING changes for the evidence (TI and telescopic only gain ESS).
calc_logA, swap and stats stay the same.

**Interactions**: thermo_jumps_every=1 (the reference script default) is
exactly what DEO wants. adapt_beta: the existing equalizers
(Acceptance_vfm / deltaBeta_deltaL_vfm) equalize exactly what DEO theory
asks to equalize (per-pair rejection rates) → the adaptive ladder is
already DEO-optimal in spirit; only note that with alternating parity
each pair accumulates stats half as often → perhaps double
adapt_beta_every. adjust_beta (insert/remove rungs) coexists fine, parity
is by index at each sweep.

**Measurable BEFORE implementing, from existing outputs**: Λ (the
communication barrier) = Σ per-pair rejection rates, computable from the
thermo_jump_stat columns already emitted → predicts the DEO round-trip
rate 1/(2+2Λ). And current round trips can be counted offline by
following id_walker in `*__i_beta__i_walker.csv`.

**Recommended sequence**: run the null test with the current sampler
(validated); from the smoke run measure round trips and Λ; if index
diffusion is the bottleneck, enable DEO as a separate change and
re-measure. Do not mix the sampler change with the evidence validation in
the same batch.

## save_Score: IMPLEMENTED 2026-09-03 (not compiled; Luciano compiles)

- Class: legacy/parallel_tempering.h:1991 (beside save_Predictions), with
  report_title and the report_model no-op.
- Report overload: legacy/qmodel.h:9297-9390 (beside the save_Predictions
  one). Direct `dlogLikelihood(ff[i_th], lik, data.get_Parameter(...), y,
  x)` per walker, fork+omp, beta==1 every event / full ladder every
  n_beta·interval events, extraction idiom copied verbatim from
  evaluate_likelihood_as_dlogPs_impl. Failed Maybe → stderr line + row
  skipped.
- Wiring: qmodel.h:9701, appended to the dts save_mcmc tuple with the
  Save_Evidence_every pair (the bandwidth budget self-throttles it,
  ~every 88 iterations at reference settings, ~+10% compute).
- Files: `*__i_beta__i_walker__i_par_score.csv` (iter, iter_time, i_beta,
  num_beta, beta, i_walker, id_walker, logL, i_par, dlogL) and
  `*__i_beta__i_walker__i_par_j_par_fim.csv` (..., i_par, j_par, gfi;
  lower triangle).
- Validation on first run: (1) logL column must equal save_likelihood's
  logLik row-by-row (~1e-12); (2) FD spot-check of dlogL via
  calc_likelihood at theta ± h·e_j. NOTE: with max_iter ~1000 and the
  self-throttled interval (~88) the FULL-ladder events are rare
  (~every 700 iterations → maybe one); raise max_iter or lower
  sampling_interval for Bartlett-per-temperature coverage; beta=1 rows
  arrive every ~88 regardless.

### First real measurement, figure_0 local 2026-09-05 (commit 5ff1280a-dirty)

Reference job: scheme_CCO stationary, 310 samples, p=8, 16 rungs, 32
walkers, max_values=128. As shipped the event does the FULL ladder (512
dlogLikelihood evaluations), so the event fires every 180 iterations
(num_values 55 × 512 / 128) and costs 30.4 s at 2 threads, 18.7 at 4,
~14 at 8. That is 24-30% of total wall time, not the ~10% estimated
above, and the event scales poorly with threads (2→8 gives only 2.1×).

Why dlogL/logL ≈ 60× and not ~10×: measured logL = 1.9 ms, dlogL =
119 ms (2 threads). It is NOT a hidden numeric path — selfDerivative
is analytic (qmodel.h:8678), the Gaussian FIM is the per-step analytic
formula XXT(d_y_mean)/var + XXT(d_y_var)/(2 var²) (qmodel.h:6152), and
d(Qdt) is memoized alongside Qdt (CLI_function_table.h:117, the
Memoiza_overload holds both value types). The correct analytic floor
for forward-mode MATRIX products is 1+2p ≈ 17× (d(AB)=dA·B+A·dB), not
p+1 ≈ 9×; the remaining ~3.5× is allocation overhead: with k=3 the
matrices are 3×3, flops are negligible, and every Derivative op plus
the per-step GFI accumulation allocates fresh objects. Evidence that it
is memory/allocator-bound: thread-seconds per dlogL WORSEN from 119 ms
(2 threads) to ~227 ms (8 threads) inside the event.

TODO (future, needs recompile, expect ≤2-3×, decided 2026-09-05):
- qmodel.h:6158-6162: accumulate the per-step GFI in place (+=) instead
  of `current_GFI = current_GFI + t_GFI` (a fresh SymPosDef allocation
  per timestep, 310 per evaluation).
- qmodel.h:9329-9338 (save_Score report): preallocate `all_scores`
  (resize + index assignment) instead of push_back of large objects
  inside the omp region; optionally collapse the walker×beta loops with
  a dynamic schedule.
- qmodel.h / qmodel_types.h range canaries (check_g*_in_range,
  to_Probability*): rate-limit the [warn] band's std::cerr line (print
  the first N per process, then count + max excursion, summary at exit).
  Measured 2026-09-05: ~2300 warn lines/iteration under MCMC (walkers at
  low beta live in the warn band), 267 MB of stderr in a 4-minute run,
  unbuffered cerr from inside omp regions (thread contention + ~15
  flushes per line). Until then the figure_0/figure_1 runners filter
  these lines out of the logs (filter_warns in ops/slurm/run_figure_0.sh,
  ops/slurm/run_figure_1.sh, ops/local/dispatch_figure_1_local.sh); the
  binary still pays the write cost, the disk does not.
The big lever stays the CADENCE (how many score samples per temperature
are actually needed), which is a design decision, not code.

SECOND, DEEPER FORM OF THE SAME DEFECT, ALSO FIXED 2026-09-06 (same rebuild):
the savers derived their cadences independently as point_size/max_values,
yielding mutually incommensurate integers (measured in the figure_1
campaign: parameters every 128, score every 180; other savers land on
220/310/4960 depending on settings), so ANY cross-saver join only met at
rare common multiples. Fix:
aligned_sampling_interval() (parallel_tempering.h, used by every saver in
the live dts path) rounds each budget UP to the next power of two, making
every pair of cadences nested by construction and aligned with the
adapt_beta windows (adapt_beta_every is a power of two in the lanes). At
campaign settings score and parameters both land on 256. cuevi/levenberg/
fraction savers deliberately left untouched (outside the live path).

DESIGN DEFECT FOUND AND FIXED IN CODE 2026-09-06 (needs the next rebuild):
the score csv carried dlogL but not theta, so the tempered-target identity
test Var(beta*dlogL + dlogprior) = E[beta*GFI + I_prior] could only use the
iterations where save_Parameter's cadence (128) happened to coincide with
save_Score's (180): 5 of 166 events, wasting ~97% of the expensive dlogL
data. Fixed by appending a par_value column to each score row
(parallel_tempering.h report_title + qmodel.h report): every score event is
now self-sufficient. Campaigns run before the rebuild (figure_1 v1,
figure_2 v1) remain limited to the coincident events for this test.

Also measured 2026-09-05: the trapezoid-vs-telescopic bracket closes on
first real data (windowed ss_up −352.90 vs ss_dn −352.82 at iter 888 of
the discovery run).

THIRD DEFECT OF THE SAME FAMILY, FOUND AND FIXED 2026-09-06 — COLUMN
LABELS PERMUTED. `Moment_statistics` is `Vector_Space<count, mean,
variance>` (moment_statistics.h:1090), so every statistics triple written
through `.sep()` emits count, mean, variance; both save_Evidence and
save_likelihood titled those triples mean, var, count. Consequence: in
every csv written before this fix, each statistics triple must be read
as count, mean, variance. This RETRACTS the "open check" noted on
2026-09-05 (that the windowed trapezoid `mean_log_Evidence` printed 0):
it was not a broken accumulation, it was the count column under the
mean's name. Verified on the campaign (rep2 s910121, cold rung, last
iteration): the column labeled mean_logL holds 32 (the walker count),
var_logL holds −319.3 (the mean logL), count_logL holds 2.0 (its
variance); the windowed trapezoid is −343.00, sitting in the column
labeled var_log_Evidence, 3.4 nats below the telescopic bracket
(−339.9/−339.3) exactly as its O(dbeta^2) discretization bias predicts.
The telescopic columns added on 2026-09-02 are plain scalars emitted in
title order and ARE correctly labeled, so the figure_1 estimates (built
on mean_log_Evidence_ss / _ss_dn) are unaffected. Fix: titles now say
count, mean, variance; the windowed logL triple is renamed *_logL_w to
end its pre-existing name clash with the varLik triple.

## save_Score design (Bartlett/FIM at save time), original plan 2026-09-02

Luciano's idea: the saving hook is the one moment where paying dlogL makes
sense. Design, with every anchor verified:

- **Template to copy**: `save_Predictions`. Class beside it in
  parallel_tempering.h (two ofstreams, own interval pair, title +
  variadic no-op fallbacks); the real `report` overload in qmodel.h
  beside the save_Predictions one (qmodel.h:9188), because it needs the
  likelihood machinery. Signature names what it needs and absorbs the
  rest: `report(FunctionTable& f, std::size_t iter, const Duration&,
  save_Score<var::Parameters_transformed>& s, thermo_mcmc<...> const&
  data, Prior const& prior, t_logLikelihood const& lik, const Data& y,
  const Variables& x, ...)`.
- **The compute, RE-CORRECTED (Luciano: skip the stages, use
  selfDerivative)**: call `dlogLikelihood(ff[i_th], lik,
  data.get_Parameter(i_walker, i_b), y, x)` DIRECTLY on the incoming
  lik. Why this is right in the tempering world and wrong in the modern
  one, so nobody re-imports the detour: there are TWO worlds. The modern
  likelihood_algorithm wraps the type-erased single-signature
  `interface::IModel<Parameters_values>` (likelihood.h:40), whose virtual
  dispatch cannot take Derivative types, hence the paper commands'
  load_dmodel + rewrap dance (calculate_mdlikelihood_impl,
  likelihood.cpp:626-642). The dts/tempering lik instead carries the
  CONCRETE legacy model (decltype(model0) template parameter, templated
  call operator), so `dlogLikelihood`'s internal `selfDerivative(p)`
  (qmodel.h:8678) lifts the parameters and the concrete model's
  templated operator() propagates Derivative types with no reload.
  Precedent for the same direct call: the Levenberg tempering
  (parallel_levenberg_tempering.h:291, 1048). CAUTION (Luciano
  2026-09-03): Levenberg was functional at some point, no guarantee now,
  and in fact CLI_thermo_levenberg_evidence.h is NOT included in
  command_manager.cpp, so nothing in the current build instantiates
  those calls; the precedent is parse-level only. Return is the same
  `dMacro_State_Hessian_minimal` = {Derivative<logL,
  Parameters_transformed>, Gaussian_Fisher_Information}; extraction:
  primitive() = logL, derivative()() = score vector, get<GFI>() = the
  Gauss-Newton FIM.

  **Risk split and why direct is still safe**: the NUMERIC core is shared
  with the paper path. make_model_interface virtualizes into the SAME
  concrete model's templated operator() with Derivative types, so the
  concrete-lambda Derivative instantiations and the
  log_Likelihood<dMacro_State_Hessian_minimal> engine are exactly what
  the August dlik/fisher lanes validated. The direct path's delta is one
  template instantiation (Model = concrete type instead of the IModel
  wrapper), currently instantiated NOWHERE in the build; its dominant
  failure mode is a compile error at the first build, not silent wrong
  numbers.

  **Self-validation built into the output, independent of Levenberg's
  health**: (1) primitive(Derivative<logL>) goes into the score CSV and
  must equal the sampler's logLik in save_likelihood's CSV for the same
  (iter, i_beta, i_walker); those come from the non-derivative path, so
  row-by-row equality (~1e-12) validates the value channel on every
  event, for free. (2) Offline finite-difference spot check of the score:
  evaluate logL at θ ± h·e_j for a few saved walkers with the existing,
  validated `calc_likelihood` CLI command and compare against the emitted
  ∇logL. Zero new C++. If maximum current-validation is ever preferred
  over simplicity, the fallback is the load_dmodel+rewrap route (the
  paper's exact instantiation), at the cost of the cmd-layer dependency.
- **Placement, simplified**: with no cmd-layer dependency, follow the
  save_Predictions precedent exactly: class beside save_Predictions in
  parallel_tempering.h, report overload in qmodel.h beside the
  save_Predictions one (:9188). Only compile-check pending: the
  dlogLikelihood instantiation for the dts member (recursive=1, av=2,
  variance=1, taylor=0) over the concrete legacy model; the Levenberg
  header instantiates its own combo the same way.
- **Files** (long format, house style): `__i_beta__i_walker__i_par_score.csv`
  (iter, i_beta, beta, i_walker, id_walker, logL, i_par, dlogL_dpar) and
  `__i_beta__i_walker__i_par_j_par_fim.csv` (lower triangle of GFI).
- **Cost control, no new knob**: pass the SAME (sampling_interval,
  max_number_of_values_per_iteration) pair as the other savers; the
  existing per-saver throttle interval = max(interval,
  point_size/max_values) self-scales because score's point_size
  (rungs·walkers·(k + k(k+1)/2)) is 20-40× evidence's, so its events are
  automatically 20-40× rarer, and compute only happens on events. With
  the reference settings (interval 1, max_values 128, 8 rungs, 32
  walkers, k=8): one score event every ~88 iterations ≈ +10% total cost.
- **β=1 priority** (save_Predictions' trick at qmodel.h:9225): evaluate
  the β=1 rung every event (F/J and posterior score, cheap), the full
  ladder every n_beta·interval events (Bartlett per temperature).
- **Wiring**: add `save_Score<var::Parameters_transformed>` to the
  save_mcmc tuple in new_thermo_Model_by_max_iter_dts (qmodel.h:9605-9620)
  with the same pair. Only the dts maker; Saving_intervals type untouched.
- **Offline (R)**: prior score ∇logP = −Σ⁻¹(θ−μ) from the prior csv and
  saved θ; Bartlett-1 per β (mean of β·∇logL + ∇logP = 0, ESS-aware);
  Cov(score) vs E[β·GFI] + Σ_prior⁻¹ (the F/J diagnostic per rung, GFI is
  the Gauss-Newton proxy so deviations = misspecification or proxy error,
  not necessarily non-convergence); scores as control variates for
  E_β[logL] (variance reduction for BOTH evidence estimators).
- **Paper anchor (checked in the eLife .macroir files)**: the paper's
  score/FIM producer is the DSL command `calc_dlikelihood_predictions`
  (figure_3_mle.macroir Stage 4, "score / dlikelihood at theta_sim AND
  theta_pool") on the modern likelihood_algorithm, routing to
  calculate_mdlikelihood_impl. Two flavors exist: the _predictions one
  emits per-interval evolution (the 286 MB dlik dumps); save_Score wants
  the MINIMAL flavor (recording-total dMacro_State_Hessian_minimal),
  same stack, small output.
- **Feasibility note**: the retarget made this safe; dts now instantiates
  the canonical IR member whose dlogLikelihood path is the well-trodden
  MLE one. Failed Maybe per walker: skip row, count, one stderr note.
- Decisions open: β=1 priority on/off; FIM full triangle vs diagonal;
  all walkers vs subsample (recommended: all, the throttle already pays).
- Estimated size: ~120-line class + ~90-line report overload + 1-line
  tuple insertion.

## Measurement-level diagnostic (plogL variance), analysis 2026-09-03

Luciano: the variance of the per-measurement log-likelihoods (plogL)
looks like an interesting test; saving everything fills the disk, but
saving MOMENTS could be the option; not now, a next step. Analysis:

**What plogL is, and the name.** Per-step conditional logL increment:
plogL_t = −½log(2πv_t) − ½χ²_t (qmodel.h:3533), with eplogL = its
model-conditional mean −½log(2πv_t) − ½ (there is even a Poisson-noise
variant, calculate_elogL, qmodel.h:4102-4119) and vplogL = ½. The totals
logL/elogL/vlogL accumulate these per step and are LIVE
(update_macro_state, qmodel.h:5999-6003) and EMITTED per walker
(save_likelihood: logLik, elogLik, vlogLik). The name IS equivocal,
three different 'p' meanings coexist: plogL (per-step), p_P_mean/p_y
(predicted vs observed, same function), plog_Evidence (per-β-segment).
Rename to logL_step/elogL_step/vlogL_step when this area is next touched
(one meaning, one type). The class defs at qmodel.h:215-222 are
commented out; live ones elsewhere.

**What the test is.** plogL_t − eplogL_t = −½(χ²_t − 1) EXACTLY (the
log-det cancels), so the variance-of-plogL test is a pure dispersion
test on standardized residuals: Var_emp(χ²_t) vs 2, i.e. the FOURTH
moment of z_t. It complements the paper's r̄²_std (second moment,
figure_7 family) with tail sensitivity.

**Free TODAY, zero changes**: the total-level version. At β=1,
z = (logLik − elogLik)/√vlogLik per walker (columns already in
save_likelihood's csv; increments are martingale differences so
Var(total) = Σ½ = T/2). z at posterior-typical θ should be O(1);
|z| ≫ 3 = misfit. On the twin-null (well-specified simulated data) this
is another guaranteed-O(1) channel. Use it in the null-test R analysis
immediately. Caveat: a goodness-of-fit check, not a convergence one; and
it is a β=1 story (data fixed; tempered rungs are deliberately
mis-weighted, the z-profile across β reads as "where the ladder unlearns
the data", not as misfit).

**What the total washes out, and the moments design (next step, not
now).** A 6% overdispersion concentrated in 100 steps vanishes inside
T/2. Two aggregation directions, both disk-safe:
(a) across-t per evaluation (scalars per walker): m2 = Σ(plogL−eplogL)²,
optionally m3, m4 → the dispersion test proper; ~4 extra columns.
(b) across-draws per t (the WAIC direction + misfit map): running
mean/var of plogL_t per t at β=1, dumped rarely; T×4 doubles ≈ 160 KB
per window at T=5000, vs the 286 MB per-(t,walker,iter) raw dumps which
stay off. Gives lppd_t/p_eff (WAIC ingredients) and the
which-millisecond-breaks-the-model map.

**Where to compute WITHOUT touching the hot path**: not in the step loop
(hot path, and extending Patch_State/logLs types has wide blast radius).
Mirror save_Score: a low-cadence saver calling `logLikelihoodPredictions`
(alive, paper-exercised, save_Predictions uses it at qmodel.h:9231; its
Evolution already materializes per-t plogL) and REDUCING Evolution to
moments before writing (precedent for reduce-from-Evolution:
reduce_micro_gradient_all_to_hessian_minimal, likelihood.cpp:766). β=1
priority, own throttle; cost ≈ one predictions pass on the cold chain
per event. The always-on in-step m2 accumulator stays as a later option
(hot path, Luciano writes it) only if low cadence proves insufficient.

**Sequence**: now, nothing new in C++; use the free total-z test in the
null analysis. After the 2×1 smoke + save_Score: save_MeasurementMoments
per the design above. Rename plogL when touching.

## Inline Experiment for the dts command, 2026-09-04

Luciano's requirement: the evidence lanes will SWEEP the experiment the
way the eLife_2025 figure lanes did (create_experiment in the .macroir,
segments injected by the dispatch script), and the flow must be the same
as the likelihood commands (calc_dlikelihood_predictions): the Experiment
object goes straight into the run, no intermediate file.

Done (legacy/CLI_thermo_evidence_dts.h, not compiled):
- Core factored out: `run_thermo_evidence_dts(filename, model, prior,
  likelihood, recording, const Experiment&, thermo_algo, ...)` holds the
  whole former body.
- Two entry points, SAME DSL name `thermo_evidence_dts`, dispatched by
  argument type: the historical one with experiment_file_type (loads the
  file, writes the thermo_evidence_dts_<id>.txt restart file, then calls
  the core) and the new one with `const Experiment&` (stateless, no
  restart txt, calls the core directly).
- Registrations disambiguated with static_cast on the function pointer
  (the overload made the plain `&calc_thermo_evidence_dts` ambiguous).
- Consequence: `thermo_evidence_dts_continuation` applies to file-based
  runs only. dts_2 and the continuations untouched.
- projects/macroir_next: data/experiments/ no longer needed for the
  evidence lane; the stationary protocol lives in the .macroir.

## Launching CCO/COC this week

Use `thermo_evidence_dts` with the post-retarget flags (INVERTED relative
to the old P2X2 script): adaptive_aproximation=0,
recursive_approximation=1, averaging_approximation=2,
variance_correction_approximation=0, variance_approximation=1. Equalizer
"deltaBeta_deltaL_vfm", controler "s", fixed seed ≠ 0, models scheme_CCO
and scheme_COC (both registered), COC prior centered on the twin
(tmp/scheme_COC_twin_par.csv), and a long stationary protocol for the
null test (≥10τ of equilibration at the matched agonist concentration).
Continuations work again after the retarget, but keep single runs for the
validation batch. Read the evidence off mean_log_Evidence (trapezoid) and
mean_log_Evidence_ss with its up-down bracket (telescopic); a single
log_Evidence row is noise.

## beta_min from the prior logL variance, 2026-09-05

Luciano's design question: start with beta_0 = 0 and set beta_min from the
VARIANCE of logL at beta = 0. Implemented in initial_beta_dts
(parallel_tempering.h): beta_min = 1/sd_0(logL) over the prior samples the
initializer already has (all rungs start at beta=0: n_beta x n_walkers
prior draws feed the same statistics the old mean rule used), fallback to
the old 1/|E_0[logL]| when the variance is degenerate, capped at 0.5.

Why variance and not mean: the two consumers of the bottom rung are
governed by beta_min*sd_0(logL), not by beta_min*|E_0[logL]|. (1) The
first stepping-stone factor E_0[L^beta_min] is an importance average with
weights exp(beta_min*logL) over prior draws; its ESS degrades like
exp(beta^2*Var_0) - 1, so beta_min*sd_0 ~ 1 keeps it estimable. (2) The
acceptance of the bottom thermo jump goes as exp(-dbeta*dlogL), same
scale. And dbeta*dlogL ~ const is exactly what the dts equalizer
(deltaBeta_deltaL_vfm) enforces on the rungs above, so the initial bottom
rung now obeys the same criterion as the adapted ladder. The rules differ
most when the prior is informative (truth-centered priors of the Figure 1
lanes: sd_0 << |E_0|), where the mean rule wastes decades of ladder.

Same day: Current_Baseline settled (Luciano): log-scale like everything,
mean 1 pA (log10 = 0). The six seeded csvs in
projects/macroir_next/data/models were fixed from the root copies'
0/-inf to 1/0. The -inf prior question from the audit is closed.
