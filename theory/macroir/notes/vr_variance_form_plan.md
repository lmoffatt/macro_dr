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

## 2. Two options for the flag surface

### (a) New orthogonal flag `uses_variance_form_aproximation` (total=0 / residual=1)
The clean, compiler-checked, durable path. Mirrors the just-landed `family`/`qdt_method` idiom.
Seam sites (from the scout):
- **Type machinery** `legacy/qmodel_types.h`: class + `_value` (mirror `94-98`), `_c` concept (mirror
  `145-147`, NOT the micro concept @169-170). Cf. `uses_family_aproximation` @110-115.
- **Compute thread** (mirror `variance_correction`, a proven ~144-site thread): template params +
  `requires` on `safely_calculate_Algo_State_recursive` (`qmodel.h:4457-4464`) and the ONE compute
  site (the `ms` selector @4544-4557); `safely_calculate_Algo_State` (@5900-5919); `Macror`
  (@6316,6335); `MacroR2::operator()` (@8117); `log_Likelihood` template + ~9 call sites
  (@6682-6684, calls @8513/8538/8662/8696/8732/8750/8776/8798/8816).
- **Model classes** `legacy/qmodel.h`: `Likelihood_Model_constexpr` (@8124-8155); `Likelihood_Model_regular`
  member + `range_*` + ctor + `cartesian` (@8223-8228) + `Likelihood_Model_variant_impl` (@8230-8255)
  + `get_variant` (@8259-8294). **This is the bulk.**
- **Builder** `include/macrodr/cmd/likelihood.h`: param on `build_likelihood_function_with_family`
  (@29-33) + a `constexpr_Var_domain<…uses_variance_form…>` slot in the macro-non-vc branch ONLY
  (@46-65; do NOT widen the vc/micro/LSE branches @70/95/120 — cartesian blow-up) + declval @154-155.
  Keep the bool `build_likelihood_function` wrapper passing a default (total).
- **Command** `src/cli/command_manager.cpp:1024-1032`: arity/signature + named arg.
- **Dispatch** new `--algo_variance_form_approximation = indexed_int_by(...)`; `.macroir` →
  `build_likelihood_function_with_family(..., variance_form_approximation = ...)`. Template already
  proven by `dispatch_figure_3_LSE.sh` (family injection).

### (b) 4th averaging value `av=3` ("av=1 mean + residual variance")
The fast throwaway path. One compute line + domain + ToString; **reuses `algo_averaging_approximation`,
no new command arity, `.macroir` untouched**.
- Compute: `ms` selector @`qmodel.h:4546` → `== 2 || == 3`. Every OTHER averaging branch already routes
  `av=3` (>0, ≠0, ≠2) onto the av=1 side automatically (gSg @4515 else, gS @4583 else, dynamic @5518) —
  which is what VR wants.
- Domain: `likelihood.h:49` add `3` to the macro-non-vc branch ONLY.
- ToString `qmodel.h:640` mis-labels av≠2 as `"__"` → av=3 needs a fix or it collides with MR in output.
- Dispatch: `macro_VR) recursive=true; averaging=3; taylor=false; micro=false ;;` and nothing else.
- **HAZARD (the gvar_i audit's own root-cause lesson):** this overloads `averaging` = "endpoint count"
  with a variance-form meaning; `av=3` is not "3 endpoints". Every `== 2`/`else` averaging branch (~20
  sites) must be re-audited to confirm av=3 lands on the av=1 side except `ms` — fragile, NOT
  compiler-enforced. A missed branch = silently wrong covariance.

## 3. Recommendation: (b) to get the number, (a) if VR ships
The science question comes FIRST and is cheap: **does VR come out over-confident?** That answer decides
whether VR is a durable member at all. If VR is calibrated, the mechanism thesis falls and VR may not
survive — in which case the full option-(a) plumbing was wasted.

So:
1. **Now, to answer the science:** option (b), restricted to the macro-non-vc domain, `ms` selector +
   ToString fix + one dispatcher line. Fastest path to a running `macro_VR` column and a first
   over-confidence measurement. Treat as a throwaway experiment.
2. **If VR survives its own test** (comes out over-confident, earns a figure and a published name):
   migrate to option (a) before it ships, for the compiler-checked hygiene and to stop overloading
   `averaging`. This is the scout's recommendation for anything durable, and it matches the nomenclature
   decision that VR opens a THIRD axis (variance form), orthogonal to `av` (`../../../papers/_program/nomenclature.md`).

**Naming:** display `VR`, data key `macro_VR`, verbatim end-to-end (axis label → CSV `algorithm` →
R `ALGOS`+`ALGO_LAB`); a mismatch silently drops rows to NA. NOT a Taylor variant — the engine flag is
`taylor_variance_correction`; Methods must say the `V` of `VR` is unrelated. Do not fold the `V`/Taylor
`MRV`/`IRV` cleanup into this change.

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
