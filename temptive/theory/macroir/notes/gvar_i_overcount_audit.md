---
date: 2026-05-12
status: finding
scope: MacroIR predictive variance, MacroIRT dispatch
affects: published eLife 2025 figures, current legacy/qmodel.h
---

# Finding: MacroIR over-counts predictive variance via gvar_i confusion

## TL;DR

1. **`gvar_i` is mathematically ambiguous** in this codebase. It denotes two
   different quantities:
   - **Total variance** per starting state: `Var(Ā | X₀=i) = gsqr_i − gmean_i²`
   - **Expected residual conditional variance**: `Σⱼ Pᵢⱼ · gvar_ij = (gtotal_var_ij)·𝟏`
   These differ by `Var_j[gmean_ij | i]` — the variance across end states of
   the boundary-conditioned mean current.
2. **MacroIR (av=2) double-counts** that variance-of-conditional-mean term:
   it appears once in the IR tilde scalar `gSg` and again inside `gvar_i`
   when `gvar_i` is the *total* form (Qdtm-style).
3. **The eLife 2025 paper's MacroIR** uses the *total* form
   (`gsqr_i − gmean_i²`) almost always, so the published figure_2 results are
   subject to this over-count. In ill-conditioned regimes where the eigen
   path fails and the code falls back to `calc_Qdtm_taylor`, the
   *residual* form is used silently — so the same field name carries two
   different numerical values depending on which dispatch path succeeded.
4. **MacroIRT in the current code is not actually applying the Taylor
   correction** at all — a refactor since the submission rewired `Macror`
   to call a different kernel that drops the `variance_correction` template
   parameter, leaving the inline Taylor block at
   [legacy/qmodel.h:3362](../../../legacy/qmodel.h#L3362) as dead code.
5. **Cheap fix** for IR: compute `gvar_i_residual = gsqr_i − (gtotal_ij ∘ gmean_ij)·𝟏`
   from existing Qdtm fields, O(K²). Removes the double-count without
   touching `calc_Qdt`.

---

## How we found it

### Step 1 — The investigation started from MRT scaffolding

The original ask was to scaffold MacroMRT (Macro_MR + Taylor σ² correction)
and add it to figure_2. Adding the dispatch entry was trivial — MRT is
already in the type-domain of `build_likelihood_function`'s
"taylor-variance-correction branch" ([include/macrodr/cmd/likelihood.h:62-81](../../../include/macrodr/cmd/likelihood.h#L62-L81)),
which admits `averaging ∈ {1, 2}` with `variance_correction=true`. So at the
type level MRT is reachable.

### Step 2 — But the kernel didn't seem to handle it

Tracing the runtime call graph for any `(averaging, variance_correction)`
combination through `Macror`:

```
log_Likelihood (qmodel.h:5135)
  → MacroR2{}(...) (line 5437 or 5594)
  → MacroR2::operator() (line 6204-6210)
  → Macror (line 4769-4798)
  → safely_calculate_Algo_State<dynamic, recursive, averaging, variance>(...)
                                                      ^^^^^^^^^
                                              variance_correction dropped here
```

`safely_calculate_Algo_State_recursive` ([qmodel.h:3995-4380](../../../legacy/qmodel.h#L3995-L4380))
never sees `variance_correction`, so the inline Taylor block at
[qmodel.h:3362-3431](../../../legacy/qmodel.h#L3362-L3431) (inside
`safely_calculate_y_mean_yvar_Pmean_PCov`) is never reached. A repo-wide grep
for `safely_calculate_y_mean_yvar_Pmean_PCov` returns only self-references.

**Initial conclusion (later revised):** the Taylor block is dead code; IRT
silently equals IR.

### Step 3 — But IR ≠ IRT empirically (per Luciano)

The user's experimental observation contradicted the dead-code reading: macro_IR
and macro_IRT do produce different figure_2 results. So either the dead-code
analysis was wrong, or the difference comes from somewhere upstream of the
kernel.

### Step 4 — Looking for a flag swap in the dispatch

Audited every layer for a swap between `variance` and `variance_correction`:

| Layer | Result |
|---|---|
| CLI binding ([command_manager.cpp:719-727](../../../src/cli/command_manager.cpp#L719-L727)) | No swap — positions 4 and 5 map correctly |
| `build_likelihood_function` | No swap |
| `Likelihood_Model_regular` constructor and `get_variant()` | No swap |
| `Likelihood_Model_constexpr` template params | No swap |
| `log_Likelihood` dispatch on `variance_correction::value` | No swap |
| `MacroR2 → Macror → safely_calculate_Algo_State` | `variance_correction` dropped, not swapped |

No swap. But this confirmed: *the dispatch on `variance_correction::value`
in `log_Likelihood` (line 5239 vs 5445) selects between `calc_Qdtm` (vc=false)
and `calc_Qdt` (vc=true)* — and that's the only runtime difference between IR
and IRT in current code.

### Step 5 — Following gvar_i into the formulas

Both `Qdtm` and `Qdt` carry a field called `gvar_i`, but the formulas for
populating it are different:

- **`calc_Qdtm_eig`** ([qmodel.h:1554-1572](../../../legacy/qmodel.h#L1554-L1572)) computes
  ```
  gsqr_i = 2·Σ_{k₁,k₂} V(i,k₁) · v_WgV(k₁,k₂) · E₂_sqr(k₁,k₂) · v_Wg(k₂)
  gvar_i = gsqr_i − gmean_i²
  ```
  This is `Var(Ā | X₀=i)` — the *total* variance per starting state. Cost: O(K²) on top of E₂.
- **`calc_Qdt_eig`** ([qmodel.h:1644-1674](../../../legacy/qmodel.h#L1644-L1674)) computes
  ```
  WgV_E3(n₁,n₃) = Σ_n₂ v_WgV(n₁,n₂) · v_WgV(n₂,n₃) · E₃(λ_{n₁}, λ_{n₂}, λ_{n₃})
  gtotal_sqr_ij = V · WgV_E3 · W · 2
  gtotal_var_ij = gtotal_sqr_ij − gtotal_ij ∘ gmean_ij
  gvar_ij       = gtotal_var_ij / Pᵢⱼ            (boundary-conditioned variance)
  gvar_i        = gtotal_var_ij · 𝟏              (sum over end states)
  ```
  This is `Σⱼ Pᵢⱼ · gvar_ij = E_j[Var(Ā | i, j) | i]` — the *expected
  residual conditional variance*. Cost: O(K³) for the triple sum.

By the law of total variance:
```
Var(Ā | X₀=i)  =  E_j[Var(Ā | i, j) | i]  +  Var_j[E(Ā | i, j) | i]
gvar_i(Qdtm)   =  gvar_i(Qdt)             +  (Σⱼ Pᵢⱼ · gmean_ij² − gmean_i²)
```

The "extra" term in Qdtm is `Var_j[gmean_ij | i]` — the spread of the
conditional mean across end states.

### Step 6 — Where the over-count appears

For `averaging=2` (the IR/IRT family), the kernel computes the predictive
variance using the IR tilde scalar at [qmodel.h:4045-4057](../../../legacy/qmodel.h#L4045-L4057):

```cpp
gSg = TranspMult(t_gmean_i, SmD) * t_gmean_i +
      p_P_mean * (elemMult(t_gtotal_ij, t_gmean_ij) * u);
```

The second term `μ · Σⱼ gtotal_ij · gmean_ij = μ · Σⱼ Pᵢⱼ · gmean_ij²`
**already encodes the variance-of-conditional-mean contribution**. Then on
[line 4063](../../../legacy/qmodel.h#L4063):

```cpp
if constexpr (variance::value && averaging::value > 0) {
    auto ms = getvalue(p_P_mean() * get<gvar_i>(t_Qdt)());
    r_y_var() = r_y_var() + N * ms;
}
```

If `gvar_i` is the Qdtm-style *total* form, then `μ · gvar_i` adds back the
`Σⱼ Pᵢⱼ · gmean_ij² − gmean_i²` piece a second time. That's the over-count.

For `averaging=1` (the MR family), `gSg` has only the first term — no
boundary cross-cov — so no double-count, and the Qdtm-style total variance
is exactly what's needed. **MR is correct as-is.**

---

## What the paper actually used (and the silent inconsistency)

In `macro_dr_submission/qmodel.h` the dispatch went through
`safely_calculate_y_mean_yvar_Pmean_PCov` (the same function that now holds
the dead Taylor block in current code). The submission's `Macror`
([submission qmodel.h:3576-3586](file:///home/lmoffatt/Code/macro_dr_submission/macro_dr_submission/qmodel.h)) did
forward `variance_correction` correctly:

```cpp
auto Maybe_all =
    safely_calculate_y_mean_yvar_Pmean_PCov<recursive, averaging, variance,
                                            variance_correction>(
        t_prior, t_Qdt, m, Nch, p_y, fs);
```

And inside that function, the av=2 + vc=false branch used
`+N·μ·gvar_i` from `t_Qdt`, where `t_Qdt` was a `Qdtm` (because the upstream
dispatch `if constexpr (!variance_correction.value)` selected `calc_Qdtm`).

**So the published paper's MacroIR was already double-counting.** This is
not a regression introduced by recent refactors; it has been present since
the submission code.

### The eigen-vs-Taylor inconsistency

`calc_Qdtm` dispatches to either `calc_Qdtm_eig` (preferred) or
`calc_Qdtm_taylor` (fallback when eigen fails) at
[submission qmodel.h:2510-2522](file:///home/lmoffatt/Code/macro_dr_submission/macro_dr_submission/qmodel.h):

```cpp
auto t_Qeig = ...;
if (t_Qeig) {
    auto Maybe_Qdt = calc_Qdtm_eig(f, m, r_Qeig, ...);
    if (Maybe_Qdt) return Maybe_Qdt.value();
}
return calc_Qdtm_taylor(m, t_Qx, ...);
```

These two paths populate the same `gvar_i` field with **different formulas**:

| Path | Formula | Quantity |
|---|---|---|
| `calc_Qdtm_eig` (line 2046-2086) | `gsqr_i − gmean_i²` | Total `Var(Ā\|i)` |
| `calc_Qdtm_taylor → Qn_to_Qdtm` (line 2453-2483) | `gtotal_var_ij · 𝟏` | Residual `Σⱼ Pᵢⱼ · gvar_ij` |

So **the same field name carries different numerical values depending on
which numerical path succeeded**. In the submission:

- **Well-conditioned regimes** (eigen succeeded, ~all of figure_2):
  Qdtm-style total → MacroIR over-counts.
- **Ill-conditioned regimes** (eigen failed → Taylor fallback):
  Qdt-style residual → MacroIR is correct.

This is a silent regime-dependent change in semantics. The current
`legacy/qmodel.h` preserves the same two formulas in the same two functions,
so this inconsistency is still present.

---

## Implications for the eLife 2025 paper

1. **MacroIR's predictive variance is biased upward** by
   `N·μ·(Σⱼ Pᵢⱼ · gmean_ij² − gmean_i²)` per interval, in the parameter
   regimes where the eigen path succeeds (which is essentially all of the
   published figures). The bias is a genuine probabilistic miscount, not a
   numerical artifact.
2. **The size of the bias depends on the model**: the larger
   `Var_j[gmean_ij | i]` is — i.e., the more the conditional-mean current
   spreads across end states — the larger the over-count. For
   `scheme_CO`-like models where boundary states give similar
   `gmean_ij`, it can be small; for richer schemes it can be substantial.
3. **MacroIRT (Taylor correction) was in fact functional in the submission
   code**, computing the rank-2 Laplace correction. So any IRT-vs-IR
   contrast in the published results reflects:
   (a) the over-count present in IR but corrected by the Taylor block (which
       uses Qdt's gvar_ij, gtotal_var_ij directly), and
   (b) the genuine rank-2 Laplace contribution.
   These two effects are confounded in the published comparison.
4. **MacroIRT in current code does not apply the Taylor correction at all**
   (regression since submission). The IRT-vs-IR gap in current figure_2 is
   only the silent gvar_i numerical-path inconsistency, which is small.
5. **MacroMR is correct** in both submission and current code — the M
   family wants Qdtm-style total variance, which is what `gvar_i` provides
   on the eigen path.

### Honest accounting for any followup paper

The over-count was a long-standing latent bug in the submission code. Three
options for any followup work:

- **Acknowledge and fix**: re-run figure_2 with the corrected
  `gvar_i_residual` formula (cheap, see below) and report the cleaner
  results. Note the original over-count as an erratum or supplementary
  discussion.
- **Maintain submission semantics**: keep the over-count for
  reproducibility, document the bias explicitly, and quantify its
  magnitude in ranges of practical interest.
- **Restore IRT and re-evaluate**: revert `Macror` to call the function
  with the Taylor block, fix the IR over-count, and re-run figure_2. Use the
  diamond comparison (MR / MRT / IR / IRT) to decompose the gain into the
  boundary-state lift contribution and the Taylor σ² contribution
  separately — this is what makes the IRT method scientifically defensible
  as distinct from "IR with the over-count fixed".

---

## Cheap fix (no `calc_Qdt` needed)

From the law of total variance:

```
gvar_i(residual)[i]  =  gsqr_i[i]  −  Σⱼ gtotal_ij[i,j] · gmean_ij[i,j]
                     =  gsqr_i[i]  −  ((gtotal_ij ∘ gmean_ij) · 𝟏)[i]
```

All three inputs are already in `Qdtm`. Cost: O(K²) on top of existing Qdtm —
no need for the O(K³) `WgV_E3` triple-sum that `calc_Qdt_eig` requires.

The numerical hazard is dividing by small `Pᵢⱼ` when `gmean_ij` is computed
internally; `Qn_to_Qdtm` already uses `elemDivSafe(..., min_P)` for the
existing `gmean_ij` field, so re-using that field rather than re-deriving is
safe.

Implementation choice (recommended): compute on the fly inside the av=2 path
of the kernel, with a comment pointing back to this document.

---

## Action items

1. **Document the finding in the eLife 2025 decision log**
   ([papers/macroir-elife-2025/02_decision_log.md](../../../papers/macroir-elife-2025/02_decision_log.md))
   with a pointer to this note.
2. **Decide on submission-fidelity vs corrected-results** for any followup
   paper.
3. **If correcting**: implement `gvar_i_residual` in the av=2 + vc=false
   branch of `safely_calculate_Algo_State_recursive`
   ([legacy/qmodel.h:4063](../../../legacy/qmodel.h#L4063)).
4. **Optional, if pursuing IRT seriously**: revert `Macror`
   ([legacy/qmodel.h:4769-4798](../../../legacy/qmodel.h#L4769-L4798)) to
   call `safely_calculate_y_mean_yvar_Pmean_PCov` so the Taylor block
   becomes live again. Without this, IRT in current code is just IR with
   slightly different gvar_i.
5. **Disambiguate `gvar_i` permanently**: either rename one of the two
   formulas (e.g. `gvar_i_total` vs `gvar_i_residual`) or constrain each
   `Qdt`-flavor type to one definition only. The current state where the
   same field name carries different mathematical meanings depending on
   numerical dispatch is the underlying root cause and will keep biting.
6. **Cross-check `Qn_to_Qdtm` consistency**: confirm whether the eigen-path
   total-variance formula and the Taylor-fallback residual-variance formula
   are intentional (and which is "correct" by the submission's own
   reckoning). If unintentional, fix the inconsistency before any new
   results are produced.

---

## Files / lines referenced

- `legacy/qmodel.h:1554-1572` — `calc_Qdtm_eig` gvar_i formula (total)
- `legacy/qmodel.h:1644-1674` — `calc_Qdt_eig` gvar_i formula (residual)
- `legacy/qmodel.h:3319-3494` — `safely_calculate_y_mean_yvar_Pmean_PCov` (the once-live, now-dead kernel)
- `legacy/qmodel.h:3362-3431` — inline rank-2 Taylor block (dead in current code)
- `legacy/qmodel.h:4045-4057` — av=2 gSg formula (the boundary cross-cov term that double-counts with Qdtm-style gvar_i)
- `legacy/qmodel.h:4063` — the `+N·μ·gvar_i` line where the double-count happens
- `legacy/qmodel.h:4769-4798` — current `Macror`, calls a function that drops `variance_correction`
- `legacy/qmodel.h:5135-5601` — `log_Likelihood`, dispatches `vc=false→Qdtm`, `vc=true→Qdt`
- `legacy/qmodel_types.h:991-995` — `Qdtm` and `Qdt` field lists
- `macro_dr_submission/qmodel.h:3576-3586` — submission's `Macror`, calls the function that propagates `variance_correction`
- `macro_dr_submission/qmodel.h:2453-2483` — `Qn_to_Qdtm` uses residual-variance formula (silent inconsistency with `calc_Qdtm_eig`)
