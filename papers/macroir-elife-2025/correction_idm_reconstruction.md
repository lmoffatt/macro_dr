# Correction: the IDM reconstruction uses the wrong square root

> Opened 2026-07-14. Actionable note: what is wrong, why, what to change, and what does **not** change.
> Closes the item recorded in `00_master_plan_v2.md` §5 (Supplement): *"two-path reconstruction of the distortion matrix (~1 but not exactly, flagged open question)"*. It is not a numerical mystery; it is a missing orthogonal factor, and the fix is exact.

## The finding, in one line

`IDM = SDM^(1/2) · CDM · SDM^(1/2)` is **false**. The exact identity is `IDM = K · CDM · Kᵀ` with `K = H^(-1/2) · J_sample^(1/2)`, which satisfies `K Kᵀ = SDM` but is **not** the symmetric principal square root of SDM.

## Why

The definitions (both in the supplement and in the code) are, with H the Fisher anchor (F_b or G_b), J_s = J_sample, J = J_total:

  IDM = H^(−1/2) J H^(−1/2)  SDM = H^(−1/2) J_s H^(−1/2)  CDM = J_s^(−1/2) J J_s^(−1/2)

The chain H → J_s → J composes exactly:

  IDM = H^(−1/2) J H^(−1/2)
      = H^(−1/2) J_s^(1/2) · [ J_s^(−1/2) J J_s^(−1/2) ] · J_s^(1/2) H^(−1/2)
      = **K CDM Kᵀ**,  K := H^(−1/2) J_s^(1/2).

`K Kᵀ = H^(−1/2) J_s H^(−1/2) = SDM`, so K is a legitimate factor of SDM in the congruence sense. But `K·K ≠ SDM`, and **K is not symmetric**: it is a product of two SPD matrices, which is symmetric only if they commute. The symmetric principal root SDM^(1/2) and K differ by an orthogonal matrix Q (polar decomposition, K = SDM^(1/2) Q), and Q = I **iff H and J_sample commute**.

Being SPD is not the issue; every matrix in the chain is SPD. Non-commutativity is. The printed identity silently sets Q = I.

## What is affected, and what is not

**NOT affected. No published or plotted number changes.**

- **The three matrices themselves.** IDM, SDM and CDM are each computed **directly** from (H, J, J_s), never through one another (`src/core/likelihood.cpp:3398`, `:3474`, `:3488`; and the Gaussian twins at `:3575`, `:3610`). Their values, spectra and interpretations stand.
- **The volume decomposition, which is what the maps plot.** det is multiplicative, so
    **log det IDM = log det SDM + log det CDM**
  holds identically, with no assumption. Verified on the data (IR, N_ch = 10⁴, noise 0.1, pool): 0.0661 vs 0.0662.
- **The physical story.** Two Gaussian approximations, two failure modes: per-sample non-Gaussianity (SDM) and temporal correlation (CDM). Unchanged.
- Consequently: **no re-run, no re-plot.**

**Affected.**

1. The matrix identity as printed in `theory/macroir/docs/Likelihood_Information_Distortion/supplement_information_distortion_main.tex` §10.
2. The observable `Likelihood_Information_Distortion_Reconstituted` (`idm_reconstituted` in the CSVs), which implements the wrong version.
3. Any statement that composes **trace, eigenvalues, Frobenius norm or the per-parameter diagonal** across the decomposition. Only the determinant composes.

## Empirical size of the error (433ed13, battery_pool, noise 0.1, bootstrap median)

Measured by comparing the emitted `idm` against the emitted `idm_reconstituted`, with ‖Q − I‖_F recomputed from the same CSVs:

| algo | N_ch | ‖IDM − I‖_F | ‖REC − IDM‖ / ‖IDM‖ | ‖Q − I‖_F |
|---|---|---|---|---|
| IR | 10 | 1.00 | 3.1% | 0.30 |
| IR | 10000 | 0.13 | 0.06% | 0.03 |
| R | 10000 | 0.73 | 0.34% | 0.03 |
| MR | 10000 | 1.76 | 0.64% | 0.04 |
| NMR | 10000 | 262.6 | 3.7% | 0.13 |
| NR | 10000 | 224.4 | 9.1% | 0.20 |

The reconstruction error tracks ‖Q − I‖ row by row, which is the prediction: it vanishes when the H and J_sample frames align, and it is largest for the algorithms whose distortion is largest. The residual is **not** numerical noise and it will not go away with more bootstrap replicates.

## The changes

### C1 — Code: fix the reconstruction operator

**Site:** `src/core/likelihood.cpp:3501`

```cpp
auto idm2 = Likelihood_Information_Distortion_Reconstituted(
    lapack::apply_sqrt_congruence(W_SDM, cdm().value(),
                                  "reconstituted distortion matrix", k_psd_rtol, k_psd_atol)
        .value_or(f_nan_spd()),
    sdm().parameters_ptr());
```

`apply_sqrt_congruence(W_SDM, ·)` computes SDM^(1/2) · CDM · SDM^(1/2). Replace with a congruence by **K = H^(−1/2) J_s^(1/2)**, both factors already available as `PSDDecomposition`s at that point in the function: `W_F_b` (whose `inv_sqrt` gives H^(−1/2)) and `W_Js` (whose `sqrt_vals` gives J_s^(1/2)).

**New helper** in `legacy/lapack_headers.h`, next to the existing congruence family (`apply_normalized_congruence` at :2844, `apply_inverse_congruence` at :2865, `apply_sqrt_congruence` at :2886):

```
// result = K · B · Kᵀ  with  K = A^{-1/2} · C^{1/2},  built from the two
// decompositions W_A (inv_sqrt) and W_C (sqrt_vals).
// This is the exact reconstruction IDM = K · CDM · Kᵀ, K Kᵀ = SDM.
apply_factor_congruence(const PSDDecomposition& W_A,   // H  → A^{-1/2}
                        const PSDDecomposition& W_C,   // J_s → C^{1/2}
                        const SymPosDefMatrix<double>& B,   // CDM
                        ...)
```

Both are already computed in scope. Watch the retained-subspace convention: K mixes two different retained bases, so the result must be embedded back through the H basis (the frame IDM lives in), and if either subspace is empty the result should stay undefined rather than zero-filled, per the §8 convention.

After the fix, `idm_reconstituted == idm` to machine precision, which turns the observable from an unexplained near-miss into a genuine numerical self-check.

**Also flag as wrong**, though it is not on the active path: `Lapack_C_h_R_C_h` (`legacy/lapack_headers.h:2634`) uses the **Cholesky** factor of SDM (L · CDM · Lᵀ). L also satisfies L Lᵀ = SDM and is also not K, so it is wrong in the same way, by a different orthogonal. Fix or delete it; do not leave a third variant of the same identity in the tree.

### C2 — Theory: fix the printed identity

**Site:** `theory/macroir/docs/Likelihood_Information_Distortion/supplement_information_distortion_main.tex` §10.

Replace `C = C_sample^(1/2) R C_sample^(1/2)` with the exact statement, and add the one-sentence reason (K Kᵀ = C_sample, K not symmetric, the two coincide iff H and J_sample commute). Add the determinant corollary explicitly, because it is the version the figures use and it is unconditional.

### C3 — Manuscript: what may and may not be claimed

- **May:** the information *volume* decomposes exactly, log det IDM = log det SDM + log det CDM, so the distortion of the information volume splits, with no residual, into a per-sample part and a temporal part. This is the sentence the maps support.
- **May:** IDM = K CDM Kᵀ exactly, as the structural statement of the chain.
- **May not:** compose the diagonals, traces, eigenvalues or Frobenius norms across the decomposition.
- Kill the "additivity residual = 0.000" line in `figure_5_master_STATUS.md` as evidence of anything. log-det additivity is an identity; the check only detects the indefinite NR/NMR cells where the log-dets are ill-defined. It is a numerical guard, not a finding.

### C4 — Planning docs

- `00_master_plan_v2.md` §5: the supplement item "two-path reconstruction (~1 but not exactly, flagged open question)" is **closed**. Either drop the supplement panel or repurpose it (see below).
- `diagnostics_plan.md` §"Two things in this section are not what they look like", item 1: already carries the analysis; its recommended fix (redefining the correlation factor in the H frame, R_H) is **superseded** by the K-congruence, which is strictly better because it keeps CDM frame-independent (the property `src/core/likelihood.cpp:3570` deliberately exploits so that one CDM serves both the F and G anchors).
- `03_metrics_diagnostics.md` §2: the word "decomposition" is fine, but say *volume* decomposition where the claim is quantitative.

## A by-product worth keeping

**‖Q − I‖ is a meaningful observable in its own right.** Q = SDM^(−1/2) · K is the rotation between the Fisher frame and the J_sample frame; it is the identity exactly when H and J_sample commute, that is, when the per-sample distortion is aligned with the information geometry of the model. In the table above it orders the algorithms the same way the distortion does. If the supplement panel for the two-path reconstruction is repurposed rather than dropped, this is what it should plot.

## Checklist

- [ ] `apply_factor_congruence` helper in `legacy/lapack_headers.h`, with the retained-subspace convention handled.
- [ ] Call site `src/core/likelihood.cpp:3501` switched to it (and the Gaussian twin, if one exists on the G path).
- [ ] `Lapack_C_h_R_C_h` fixed or deleted.
- [ ] Regression test: `idm_reconstituted == idm` to ~1e-12 on a cell with large distortion (NR, N_ch = 10⁴), which is exactly where the current code is 9% off. Put it next to `tests/math/test_idm_matrix.cpp`.
- [ ] Supplement §10 corrected.
- [ ] `figure_5_master_STATUS.md` additivity line reworded.
- [ ] `00_master_plan_v2.md` §5 supplement item closed.
</content>
