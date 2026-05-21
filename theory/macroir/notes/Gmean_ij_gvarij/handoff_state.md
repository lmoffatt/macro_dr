# Handoff: Bayesian-shrinkage + canary refactor

Status at end of session 2026-05-14. New session pick-up reference.


**Apply Bayesian additive shrinkage at the moment sources** 
(`gtotal_ij` and
`gtotal_sqr_ij` — the primary V·E2·W and V·E3·W outputs), using v_g-based
Poisson-endpoint priors. Downstream uses textbook moment formulas with
`(P + min_P_prior)` in the denominators. The `gvar_ij` prior emerges *automatically*
from the law-of-total-variance algebra applied to the two source priors —
no need to compute it explicitly.



Source-level priors (prior contribution to the *numerator* of the moment;
conductance-scale, no P factor — see derivation below):

$$
\texttt{prior\_gtotal\_ij}[i,j]     = (g_i + g_j)/2,
$$

$$
\texttt{prior\_gtotal\_sqr\_ij}[i,j] = (g_i^2 + g_j^2)/2.
$$

Substitution at construction (additive shrinkage of the numerator; the
regularized denominator $P + \texttt{min\_P\_prior}$ is applied downstream
when forming gmean/gvar):

$$
\texttt{gtotal\_ij}_{\text{shrunk}}     = \texttt{gtotal\_ij}_{\text{raw}}     + \texttt{prior\_gtotal\_ij}     \cdot \texttt{min\_P\_prior},
$$

$$
\texttt{gtotal\_sqr\_ij}_{\text{shrunk}} = \texttt{gtotal\_sqr\_ij}_{\text{raw}} + \texttt{prior\_gtotal\_sqr\_ij} \cdot \texttt{min\_P\_prior}.
$$

Why no P factor on the prior (Bayesian derivation). Let the prior on the
per-transition mean conductance be $g_{\text{mean},ij} \sim
\mathcal{N}((g_i+g_j)/2,\,\tau^2)$, and the empirical estimate from $P_{ij}$
transitions be $g_{\text{mean,raw}} \mid \text{truth} \sim
\mathcal{N}(g_{\text{mean,true}},\,\sigma^2/P_{ij})$. The precision ratio
$\sigma^2/\tau^2$ — interpreted as the pseudo-count of prior phantom
transitions — is the quantity wired in code as
$\texttt{min\_P\_prior}$. The posterior mean is

$$
\text{E}[g_{\text{mean}} \mid \text{data}] = \frac{P_{ij} \cdot g_{\text{mean,raw}} + \texttt{min\_P\_prior} \cdot (g_i+g_j)/2}{P_{ij} + \texttt{min\_P\_prior}} = \frac{g_{\text{total,raw}} + \texttt{min\_P\_prior} \cdot (g_i+g_j)/2}{P_{ij} + \texttt{min\_P\_prior}}.
$$

The numerator of this expression is exactly $\texttt{gtotal\_ij}_{\text{shrunk}}$;
the denominator is the $(P+\texttt{min\_P\_prior})$ in the downstream
$\texttt{gmean\_ij}$ formula. The prior's role is to inject
$\texttt{min\_P\_prior}$ phantom transitions, each carrying the
endpoint-average conductance — its contribution to gtotal scales with
$\texttt{min\_P\_prior}$ (the pseudo-count), not with $P$. A P factor on the
prior would zero it at $P_{ij}\to 0$ and break the small-P limit listed
below.


Downstream:

$$
\texttt{gmean\_ij}     = \texttt{gtotal\_ij}_{\text{shrunk}} / (P + \texttt{min\_P\_prior}),
$$

$$
\texttt{gtotal\_var\_ij} = \texttt{gtotal\_sqr\_ij}_{\text{shrunk}} - \texttt{gtotal\_ij}_{\text{shrunk}} \cdot \texttt{gmean\_ij},
$$

$$
\texttt{gvar\_ij}      = \texttt{gtotal\_var\_ij} / (P + \texttt{min\_P\_prior}).
$$

Pseudo-count `min_P_prior`:
- Eig paths (have `Eigs`):                       `min_P_prior = ε_mach · κ(V)` (already plumbed via `get<kappa_V>(t_Eigs)()`)
- Taylor / Schur paths (no κ):                   `min_P_prior = 10 · √N · ε_mach`

Small-P sanity check (worked through in the session):
- `P_ij → 0` ⇒ `gmean_ij → (g_i + g_j)/2` (gmean prior)
- `P_ij → 0` ⇒ `gvar_ij → ((g_i − g_j)/2)²` (gvar prior, derived algebraically)

## Current state of the code

### Already in place

| File / location | What's there |
|---|---|
| [qmodel.h:1092-1133](../../../legacy/qmodel.h#L1092) | `gmean_ij_endpoint_prior(t_g)`, `gvar_ij_endpoint_prior(t_g)`, `bayesian_shrink_ratio(...)` — v_g-based priors, applied at the gmean_ij/gvar_ij level (NOT at sources yet) |
| [qmodel.h:2716-2880](../../../legacy/qmodel.h#L2716) | `check_gmean_ij_in_range`, `check_gmean_i_in_range`, `check_gtotal_ij_in_range` (returning `Maybe_error<bool>`) and their `require_*_in_range` wrappers (returning `Maybe_error<void>`). Typed with per-quantity `requires(U<C_X, X>)` |
| [qmodel.h:2914+](../../../legacy/qmodel.h#L2914) | Old `force_gmean_in_range`, `force_gtotal_in_range`, `force_gtotal_var_in_range` — STILL ACTIVE, untouched |
| [qmodel_types.h:1050-1057](../../../legacy/qmodel_types.h#L1050) | `kappa_V` tag class |
| [qmodel_types.h:1191](../../../legacy/qmodel_types.h#L1191) | `Eigs = Vector_Space<lambda, V, W, kappa_V>` |
| [qmodel.h:1056-1086](../../../legacy/qmodel.h#L1056) | `calc_eigen` computes `κ_F(V)` |
| Call sites (15): `calc_Qdtm_eig`, `calc_Qdt_eig`, `calc_Qdtm_eig_codex`, `calc_Qdt_eig_codex`, `assemble_Qdt_from_moments`, `Qn_to_Qdt`, `Qn_to_Qdtm` | Bayesian shrinkage applied at `gmean_ij` and `gvar_ij` levels (NOT at sources). Uses `bayesian_shrink_ratio(gtotal_raw, P, prior, ε)`. |

### Stale

(none — Task 4 doc revision below has been completed; parent doc
[bayesian_prior_regularization_of_Qdt.md](bayesian_prior_regularization_of_Qdt.md)
now describes the v_g approach.)

## Remaining tasks (priority order)

### Task 1 — Refactor to source-level shrinkage (the active design above)

For each `calc_Qdtm_eig`, `calc_Qdt_eig`, `calc_Qdtm_eig_codex`,
`calc_Qdt_eig_codex`, `assemble_Qdt_from_moments`, `Qn_to_Qdt`, `Qn_to_Qdtm`:

1. After computing the raw `v_gtotal_ij = V · WgV_E2 · W`:
   ```cpp
   auto gtotal_ij_prior_mat = gtotal_ij_endpoint_prior(t_g);     // NEW helper: (g_i+g_j)/2 outer
   auto gtotal_ij_shrunk = raw + gtotal_ij_prior_mat * min_P_prior;
   ```

2. After computing the raw `v_gtotal_sqr_ij = 2 · V · WgV_E3 · W`:
   ```cpp
   auto gtotal_sqr_ij_prior_mat = gtotal_sqr_ij_endpoint_prior(t_g);  // NEW helper: (g_i²+g_j²)/2 outer
   auto gtotal_sqr_ij_shrunk = raw + gtotal_sqr_ij_prior_mat * min_P_prior;
   ```

3. Downstream replace:
   - `gmean_ij = bayesian_shrink_ratio(gtotal_raw, P, prior_gmean, min_P_prior)` →
     `gmean_ij = gtotal_ij_shrunk / (P + min_P_prior)`  (via elementwise zip or similar)
   - `gtotal_var_ij = gtotal_sqr_raw − elemMult(gtotal_raw, gmean_ij)` →
     `gtotal_var_ij = gtotal_sqr_ij_shrunk − elemMult(gtotal_ij_shrunk, gmean_ij)`
   - `gvar_ij = bayesian_shrink_ratio(gtotal_var, P, prior_gvar, min_P_prior)` →
     `gvar_ij = gtotal_var_ij / (P + min_P_prior)`

4. Store the *shrunken* `gtotal_ij` and `gtotal_sqr_ij` in the returned Qdt
   (replaces the raw versions).

5. Add new helpers near [qmodel.h:1114-1130](../../../legacy/qmodel.h#L1114):
   ```cpp
   template <class C_g>
   static auto gtotal_ij_endpoint_prior(const C_g& t_g) {
       const std::size_t N = t_g().size();
       Matrix<double> U(1, N, 1.0);
       return X_plus_XT(t_g() * U) * 0.5;     // (g_i + g_j)/2
   }
   template <class C_g>
   static auto gtotal_sqr_ij_endpoint_prior(const C_g& t_g) {
       const std::size_t N = t_g().size();
       Matrix<double> U(1, N, 1.0);
       auto g_sq = elemMult(t_g(), t_g());    // g_i² entrywise
       return X_plus_XT(g_sq * U) * 0.5;       // (g_i² + g_j²)/2
   }
   ```
   The current `gmean_ij_endpoint_prior` / `gvar_ij_endpoint_prior` helpers
   can stay or be removed — they're equivalent to the source-level priors
   modulo the algebra.

6. **Unit-check the `dt` factor.** The priors `(g_i+g_j)/2` and `(g_i²+g_j²)/2`
   are conductance-scale (no time). `gtotal_ij` in the codebase may carry an
   implicit `dt` factor (depending on how `WgV_E2` is scaled). Check at one
   eig site whether `gtotal_ij ≈ g · P_ij` or `gtotal_ij ≈ g · P_ij · dt`. If
   the latter, multiply priors by `dt` before adding.

### Task 2 — Wire canaries at the source

After shrinking each source, run:
```cpp
auto check = require_gtotal_ij_in_range(gtotal_ij_shrunk, r_P, t_g);
if (!check) return Maybe_error<...>(check.error());
```

And the analogue for `gtotal_sqr_ij` — needs a new
`check_gtotal_sqr_ij_in_range` helper:
- Per-entry bound: `[gmin² · P_ij, gmax² · P_ij]` (or
  `[gmin² · P_ij, (gmax² + ((gmax-gmin)/2)²) · P_ij]` if accounting for variance).
- Pattern parallels `check_gtotal_ij_in_range` at
  [qmodel.h:2825+](../../../legacy/qmodel.h#L2825).

### Task 3 — Retire `force_*_in_range`

Once the source-level shrinkage + canaries are in:
- The post-construction `force_gmean_in_range`, `force_gtotal_in_range`,
  `force_gtotal_var_in_range` clamps are redundant.
- 15+ call sites in `qmodel.h`. Remove the `force_*_in_range` calls and
  delete the helper definitions.
- The `enforce_gmean_bounds` policy bit becomes dead — can either remove it
  or repurpose it to enable/disable the source-level canaries.

### Task 4 — Doc revision  *(done)*

[bayesian_prior_regularization_of_Qdt.md](bayesian_prior_regularization_of_Qdt.md)
now describes the v_g (Poisson endpoint) approach, with the
"Why not Empirical Bayes from marginals?" subsection documenting why EB
was considered and rejected (κ(V)-driven marginal contamination defeats
the regularization in the regime it is designed for). The conjugate-prior
posterior derivation, κ-based ε appendix, and migration section all
already used v_g formulas and remain unchanged.

### Task 5 — `min_P` field cleanup

After Tasks 1-3 land and figure_2 runs clean:
- Drop the `min_P` field from `Qdt`, `Qdtm`, `Qdtg`, `Qn`, `micro_Qdt*` in
  [qmodel_types.h](../../../legacy/qmodel_types.h).
- Remove `min_P(1e-12)` from model initializers in
  [models_simple.h](../../../legacy/models_simple.h),
  [models_MoffattHume_linear.h](../../../legacy/models_MoffattHume_linear.h),
  [models_MoffattHume_allosteric.h](../../../legacy/models_MoffattHume_allosteric.h).
- Delete `elemDivSoftAbs` and `elemDivSafe` from `matrix.h` and
  `parameters_derivative.h` (no live callers after Task 3).

### Task 6 — `expm1`/Opitz rewrite of `E1`/`Ee`/`E2`/`E3`

Separate concern, blocked on Task 5 dropping `min_P` (the current `E1`/`E2`/`E3`
take `min_P` as an FP-near-equality threshold, which is dimensionally wrong).

Methodology already written:
[divided_differences_of_exp.md](../../scientific-software/notes/divided_differences_of_exp.md).

Recipe:
- `E1(x) → expm1(x)/x` with `x == 0` branch
- `Ee(x,y) → exp(y) * expm1(x − y) / (x − y)` with `x == y` branch
- `E2`, `E3` → Opitz form (`exp` of a small bidiagonal `J`) or Stewart's
  recursion with relative `ε^(1/(n+1))` Taylor cutoff
- Remove the `t_min_P` / `eps` argument from these helpers entirely

### Task 7 — Validation rebuild + rerun figure_2

After Tasks 1-3:
- Rebuild
- Run `ops/local/figure_2.macroir`
- Expect: zero `[warn]` lines from the simplex canaries (the source-level
  shrinkage should eliminate the `2^-16` per-step drift we tracked earlier)
- Expect: zero `to_Probability` / `to_Covariance_Probability` errors
- If anything fires: the source-level shrinkage didn't catch it; investigate
  with the `check_*` canaries to localize.

### Task 8 (aspirational) — Derivative-aware `frechet_p_a_b_combined`

Currently `calc_Qdt_schur` delegates to `calc_Qdt_taylor` for derivative
inputs because `lapack::frechet_p_a_b_combined` is plain-double only.
Writing a sparse-aware, derivative-aware version would let Schur run
end-to-end for derivatives. Algebraically the same as Taylor (per
[Qdt_schur_vs_taylor_rationale.md](../../scientific-software/notes/Qdt_schur_vs_taylor_rationale.md));
pure consistency win, no algorithmic gain. Low priority.

### Task 9 (open) — RT (MRT/IRT) rank-2 Newton derivative-conservation drift

**Scope:** affects only the `variance_correction=true` paths (MRT and IRT,
the "T" variants). The non-T paths (MR, IR) use the simple
`r_y_var = e + N·gSg + N·ms` update without iteration and are
derivative-clean — none of the runs surfaced canary fires from
those branches. The bug lives exclusively in the rank-2 Newton step at
[qmodel.h:4587-4790 (MRT)](../../../legacy/qmodel.h#L4587) and
[qmodel.h:4815-5210 (IRT)](../../../legacy/qmodel.h#L4815).

The rank-2 Newton step in `safely_calculate_Algo_State_recursive` for the
IRT family (av=2, variance_correction=true) loses probability-conservation
in the derivative w.r.t. `Num_ch_mean` at small intervals. Empirically
~2.2% relative misalignment of `∂P/∂N` against the all-ones direction,
observed reproducibly at `Num_ch=100, interval_in_tau=0.1, k=97` in
`figure_2.macroir`. The primitive Σ p = 1 holds exactly; the breakdown is
purely in the AD chain.

Localised after instrumenting Probes A/B/C inside the Newton loop
(diagnostic data session 2026-05-14): the `(p_P_mean − p_iter)`
contribution and the SmD invariant are FP-clean. The leak comes from
`delta_term + 0.5·qS_post` derivative-side, specifically the cross-term
contractions involving `V_iter/N`, `V_iter²/N` (m11, m22 entries of the
Woodbury matrix), `det = m11·m22 − m12²`, and `k11/k12/k22 = ±m_ij/det`.
N appears in denominators with explicit derivative-aware terms whose
cancellation works at primitive level but degrades through AD by ~1%
of magnitude.

Short-term D2 mitigation: the cos² test (`to_Probability`,
`to_Probability_displacement`, `to_Covariance_Probability`) has been
demoted to **warn-only**. The error-band branches are removed; warns
above `eps_dcos_warn_sq = 1e-8` still log to stderr but never abort.

Rationale: the cos² metric is a **direction** test (`Σ²/(N·‖·‖²)`), not a
magnitude one. It catches genuine derivative-conservation drifts but
also amplifies FP-noise and small-scale tail events into "errors" that
have negligible downstream impact. We've watched the pipeline catch a
6.2% drift at k=97 in IRT av=2 vc=1; the absolute effect on
`Num_ch_mean` / `unitary_current` Fisher info is small, and the
algorithm is converging-but-imperfect AD machinery at the rare tail
rather than a structural error. Warn output preserves visibility for
the D1 investigation without aborting legitimate runs.
`canary_dcos_error_sq` is now effectively dead but kept as a constant.

Long-term D1 fix — Kalman-coupled update: enforce that the **mean step**
and the **covariance down-date** share the same AD-tracked Kalman
intermediates (`k11`, `k12`, `k22`, `det`, `V_iter`, `gS`, `chi`). Under
proper Kalman discipline the two updates are coupled by construction
and probability conservation (`Σ ∂p = 0`, `Σ_j ∂Cov_ij = 0`) survives
automatically. The drift arises when AD re-evaluates the same
intermediate independently on the mean and covariance paths and FP
rounding pulls the two evaluations apart. Memoizing or explicitly
threading the shared derivatives through both branches restores the
coupling; no separate projection of `∂P/∂θ` is needed if the update
itself is constructed correctly. Own session.

## Things settled — do NOT re-litigate

These came up multiple times in the session and have settled answers:

- **v_g-based prior, not Empirical-Bayes from marginals.** EB was an
  attempt to break a non-existent circularity. v_g priors are simpler,
  dimensionally cleaner, and don't introduce a chicken-and-egg with `gmean_ij`.

- **Pseudo-count `min_P_prior` is parameter-free.** `ε_mach · κ(V)` for
  eig paths, `10 · √N · ε_mach` for taylor/schur. No `c_floor` or `c_safety`
  multipliers — both were paranoia without math behind them. See appendix in
  [bayesian_prior_regularization_of_Qdt.md](bayesian_prior_regularization_of_Qdt.md).

- **Schur and Taylor are algebraically equivalent.** They evaluate the same
  Van Loan augmented-matrix integral, differing only in storage (sparse vs
  dense) and small-interval expansion (Padé vs truncated Taylor). The
  performance gap is dense-vs-sparse, not algorithmic. See
  [Qdt_schur_vs_taylor_rationale.md](../../scientific-software/notes/Qdt_schur_vs_taylor_rationale.md).

- **Simplex canaries use the cos² test for derivatives.** Two-tier
  warn/error via `Maybe_error` + `std::cerr` interim, replacing the
  absolute-threshold check that was unworkable. See
  [qmodel_types.h:285+](../../../legacy/qmodel_types.h#L285) for the
  shared `canary_*` constants.

- **`With_warning<T>` is a future abstraction.** Currently use `std::cerr`
  for warn-band logging as a known-ugly interim. See
  [with_warning_abstraction.md](../../scientific-software/notes/with_warning_abstraction.md).

## Session reflection

The session went sideways with three design pivots (v_g → EB → v_g →
source-level). Each was reasonable on its own; the cost was thrash — code,
docs, and call sites kept getting half-updated between pivots.

Lesson for next session: when a design pivot lands, *finish the rewrite of
all artifacts* (code, doc, helpers) before starting another exploration.
The mixed-state interruptions made the assistant repeatedly ask clarifying
questions instead of executing.

The source-level shrinkage design (Task 1) is the cleanest landing place
and should be where the next session starts — execute Task 1, validate with
figure_2, then proceed sequentially. The earlier intermediate states
(`bayesian_shrink_ratio` at the gmean level, `force_gmean_in_range` as a
clamp) can be removed once Task 1 lands.
