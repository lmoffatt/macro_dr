# VR (residual-variance recursive) — implementation plan

> Opened 2026-07-20. Telegraphic. Seams verified against LIVE code by workflow scout `a38f613bd828dbb98`
> (the `gvar_i_overcount_audit.md` line numbers are STALE — live code moved to ~4457-4620/8124-8295;
> the `family` int flag from `nonlinearsqr_lse_plan.md` has ALREADY landed). Re-confirm before editing.
> Paper context: paper 1's mechanism (`papers/1_method/decisions.md`, `[[project_vr_variance_form]]`).

## 0. What VR is
`VR` = `MR`'s mean and gain (av=1, start-conditioned, no boundary cross-cov) + `IR`'s **residual**
predictive variance (instead of `MR`'s total). It differs from MR only in the observation variance,
from IR only in the mean/gain. Roster step: `R → MR → VR → IR`, where **MR→VR is the variance step**
and **VR→IR is the gain step**. Paper 1's thesis needs VR to come out **over-confident** (real variance
removed without the boundary gain); if VR is calibrated, the "MR's problem is the gain" thesis falls.

## 1. THE LINCHPIN — the residual is already computed on the MR path
The predictive-variance block runs for **both** av=1 and av=2 (`qmodel.h:4532`,
`if constexpr (variance::value && averaging::value>0)`). Inside it, the `ms` lambda already has both
forms; only an inner `if constexpr (averaging::value == 2)` picks residual vs total:

```cpp
// legacy/qmodel.h:4544-4557  (ms lambda)
auto& t_gsqr_i = get<gsqr_i>(t_Qdt);
if constexpr (averaging::value == 2) {                       // IR: residual
    auto& t_gtotal_ij = get<gtotal_ij>(t_Qdt);
    auto& t_gmean_ij  = get<gmean_ij>(t_Qdt);
    Matrix<double> u(p_P_mean().size(), 1, 1.0);
    auto gvar_i_residual = t_gsqr_i() - elemMult(t_gtotal_ij(), t_gmean_ij()) * u;
    return getvalue(p_P_mean() * gvar_i_residual);
} else {                                                     // MR: total
    auto gvar_i_total = t_gsqr_i() - elemMult(t_gmean_i(), t_gmean_i());
    return getvalue(p_P_mean() * gvar_i_total);
}
```

All three residual inputs (`gsqr_i`, `gtotal_ij`, `gmean_ij`) are in the `Qdtm` MR already receives
(`qmodel_types.h:1222-1223`). **So VR's compute delta is ONE selector flip in this lambda. No new
field, no new Qdt term, no touch to `calc_Qdt*`.** VR's mean and gain come free from the av=1 `else`
branches: `gSg` @`qmodel.h:4524`, `gS` @`qmodel.h:4587`. (The non-recursive/gap path
`safely_calculate_Algo_State_non_recursive` @4359 adds no `ms` term, so the variance form is moot
there.)

This changes the effort calculus entirely: **the whole implementation cost is the flag surface**, not
the math.

## 2. DECIDED (2026-07-20): option (a), the orthogonal flag, end to end

New orthogonal `uses_variance_form_aproximation<int>` {total=0, residual=1}, named `constexpr int`
constants (not magic ints), mirroring the just-landed `family`/`qdt_method` idiom. Compiler-checked, no
overloading of `averaging`. VR = (recursive, av=1, variance_form=residual); IR stays (recursive, av=2);
MR stays (recursive, av=1, variance_form=total).

### The ONE semantic decision, at the compute site — do NOT let IR regress
Today the `ms` residual-vs-total split keys on `averaging::value == 2` (`qmodel.h:4546`). IR gets
residual *for free* from that test. If you replace that test with `variance_form::value == residual`,
then **IR must be built with variance_form=residual or it silently reverts to the total variance** — a
regression in the published algorithm, and the worst kind because it compiles and runs.

Two ways to write the `ms` selector; pick with eyes open:

- **(A) OR-guard (recommended for the first landing, zero IR-regression risk):**
  `if constexpr (averaging::value == 2 || variance_form::value == residual)` → residual; else total.
  IR (av=2) is unchanged whatever its variance_form; VR (av=1, residual) adds the residual at av=1; MR
  (av=1, total) unchanged. Default variance_form = total everywhere is then safe for every existing
  build. Slightly impure (two flags can express residual) but cannot break IR.
- **(B) fully orthogonal:** `if constexpr (variance_form::value == residual)` → residual; else total,
  AND change IR's construction to carry variance_form=residual explicitly. Cleaner axis, but you must
  find and fix every site that builds IR, or IR regresses. More surface, higher risk.

**Recommend (A) now.** It gives the clean orthogonal flag everywhere in the type system while keeping
the one compute predicate regression-proof. A later cleanup can tighten to (B). The gSg (@4524) and gS
(@4587) boundary branches key on `averaging::value == 2` and **must not change** — VR at av=1 correctly
takes the non-boundary mean and gain. Only the `ms` term moves.

### ToString
`qmodel.h:640` renders only av — add a variance_form suffix (e.g. `_res`) so VR's output label ≠ MR's,
or the CSV/label collision silently merges the two.

**Naming:** display `VR`, data key `macro_VR`, verbatim end-to-end (axis label → CSV `algorithm` →
R `ALGOS`+`ALGO_LAB`); a mismatch silently drops rows to NA. NOT a Taylor variant — the engine flag is
`taylor_variance_correction`; Methods must say the `V` of `VR` is unrelated. Do not fold the `V`/Taylor
`MRV`/`IRV` cleanup into this change.

## 3. Ordered tasks (file:symbol — done-check). macro_VR = one more algorithm, end to end.

1. **`legacy/qmodel_types.h`** — ADD `uses_variance_form_aproximation<int>` {variance_total=0,
   variance_residual=1} + named `constexpr int` constants + `_value` class + `_c` concept. Mirror
   `uses_averaging_aproximation` (@94-98) + concept (@145-147, NOT the micro concept @169-170); cf.
   `uses_family_aproximation` @110-115. — *compiles; the two constants resolve.*
2. **`legacy/qmodel.h` compute** — thread the new param through, mirroring `variance_correction`:
   - `safely_calculate_Algo_State_recursive` template params + `requires` (@4457-4464).
   - The `ms` selector (@4544-4557): OR-guard per §2 (A): `averaging::value==2 || variance_form::value
     ==variance_residual` → residual; else total. **gSg @4524 and gS @4587 untouched.**
   - `safely_calculate_Algo_State` (@5900-5919); `Macror` (@6316,6335); `MacroR2::operator()` (@8117);
     `log_Likelihood` template (@6682-6684) + the ~9 `.log_Likelihood<…>` call sites
     (@8513/8538/8662/8696/8732/8750/8776/8798/8816).
   - `ToString` (@640): append a `variance_form` suffix (residual→`_res`). — *compiles; a residual
     build carries a distinct label; verification (0) below passes.*
3. **`legacy/qmodel.h` model classes** — `Likelihood_Model_constexpr` param + `requires` + `using` +
   ctor (@8124-8155); `Likelihood_Model_regular` member + `range_*` + ctor + **`cartesian`** (@8223-8228)
   + **`Likelihood_Model_variant_impl`** (@8230-8255) + **`get_variant`** (@8259-8294). The bulk; the
   variant cartesian grows by ×2 on the macro-non-vc slot only. — *variant instantiates; arity right.*
4. **`include/macrodr/cmd/likelihood.h`** — add `variance_form_approximation` param to
   `build_likelihood_function_with_family` (@29-33) + a `constexpr_Var_domain<…uses_variance_form…,
   variance_total, variance_residual>` slot in the **macro-non-vc branch ONLY** (@46-65; do NOT widen
   the vc/micro/LSE branches @70/95/120 — cartesian blow-up) + declval (@154-155). The bool wrapper
   `build_likelihood_function` (@141-151) passes `variance_total` as default. — *the VR combo is in the
   variant; existing 8-arg bool scripts still compile.*
5. **`src/cli/command_manager.cpp:1024-1032`** — bump `build_likelihood_function_with_family`
   arity/signature + named arg `variance_form_approximation`. The bool `build_likelihood_function`
   registration unchanged. — *family builder callable with the new named arg; old scripts unaffected.*
6. **`legacy/qmodel.h` visit / raw builders** — the raw-ModelPtr branch builders that enumerate the
   macro variant must include the new slot or they fail to compile; default them to `variance_total`.
   (Same class of edit the family flag needed; grep the sites that construct the macro `Likelihood_Model`
   variant.) — *whole tree compiles.*
7. **`projects/eLife_2025/ops/local/figure_3_mle_G.macroir`** (@81-89) — switch the builder from the
   bool `build_likelihood_function` to `build_likelihood_function_with_family`: replace
   `micro_approximation = algo_micro_approximation` with `family_approximation = algo_family_approximation`
   and ADD `variance_form_approximation = algo_variance_form_approximation`. (Template: the LSE macroir
   already calls `_with_family`.) — *macroir parses; macro_* run unchanged with family=0, vf=total.*
8. **`projects/eLife_2025/ops/slurm/dispatch_figure_3_G.sh`** (@130-141 case block, @165-169 injection) —
   migrate from the `micro` bool column to `family`+`variance_form` int columns (mirror
   `dispatch_figure_3_LSE.sh` @134-145/175):
   - every existing case gains `family=0` (macro) / `family=1` (micro) and `variance_form=0`;
   - **ADD** `macro_VR) recursive=true; averaging=1; taylor=false; family=0; variance_form=1 ;;`
   - replace the `--algo_micro_approximation` injection with `--algo_family_approximation = indexed_int_by(...)`
     and ADD `--algo_variance_form_approximation = indexed_int_by(axis=algorithm_axis, values=[$variance_form])`;
   - assemble both into the run args (@195-196). — *`N_ALGO="macro_VR" ...G.sh` dispatches; CSV shows a
     `macro_VR` `algorithm` column.*
9. **(beyond "hasta dispatch", but data will not plot without it) R figures** — append `macro_VR` to
   `ALGOS` + `ALGO_LAB` (display `VR`) in each target `.Rmd`; bump `wrap_plots` ncol for the
   filename-token figures. Verbatim `macro_VR` or rows drop to NA. — *the VR series renders.*

**Compile gate:** do NOT `[[feedback_no_compile]]` build in this environment; hand the tree to Luciano
to compile. The done-checks above are the acceptance oracle.

## 4. Config (mirror nonlinearsqr §7, but VR is a full likelihood so fewer constraints)
VR uses the full emission model (mean + variance + gain), so it does NOT need the LSE amplitude/pure-
variance Fixes. Same `create_parameters` as MR/IR. Only the algorithm flag differs.

## 5. Verification ladder
- **(0) same μ as MR:** VR's predictive MEAN must equal MR's exactly (same av=1 mean path). Diff = 0 at
  fixed θ, one recording. Catches a wrong mean/gain branch.
- **(1) variance below MR, at/above IR:** VR's predicted `y_var` = MR's minus `N·μ·Var_j[gmean_ij|i]`
  (the residual removes the end-state spread). Check `y_var(VR) < y_var(MR)` per interval, and that
  VR's `y_var` equals IR's `ms` term but with MR's `gSg` (no boundary cross-cov in the mean part).
- **(2) score/Fisher:** consistent Newton pair, as for the other macro members (no special handling).
- **(3) the paper's question:** the parameter-space distortion `C = H^{-1/2} J H^{-1/2}`. Prediction:
  VR over-confident (C_ii > 1), MORE than MR if the gain is what MR was implicitly relying on. If
  VR ≈ calibrated, STOP and re-open the mechanism thesis (`papers/1_method/decisions.md`).

## 6. Run
Once building: add `macro_VR` to `dispatch_figure_3_G.sh` and run the same band-A grid as R/MR/IR
(N_ch × noise, interval inside). VR does NOT need the band-C fill (paper 2). Same n_sims as the
band-A cells (10⁴) so it pools with them (`../../../papers/_program/machinery.md` §6.1).
