# Components — MAP (hierarchy / index)

> **STATUS: PARKED (2026-07-15).** This scheme is a proof-of-concept, not live. The paper's live index is
> `../00_master_list.md`. The three demo nodes that proved it (MR/sandwich/fact) are in
> `../archive/components_poc/`. This scheme earns propagation only when `_LOG.md` records that master_list
> failed (the D-0 rerun-refresh test). Below = the plan IF revived; do not build from it now.
>
> The hierarchy lives HERE, not in the nodes (DITA map). `✓` = built (archived), `·` = pending.
> Governed by `_SPEC.md`. Whether the scheme works: `_LOG.md`.

## algorithms/ — type `reference`
- ✓ `[[MR]]` — Mean-averaged Recursive (av=1). Cautionary intermediate; over-confident 1.5–2.1.
- · `[[NR]]` — Non-recursive instantaneous (av=0). ≈ Milescu 2005. Overconfident ×10–16, cond ∝N_ch².
- · `[[NMR]]` — Non-recursive mean (av=1). Unbiased, over-confident; possible speed niche.
- · `[[R]]` — Recursive instantaneous (av=0). ≈ Moffatt 2007 / Münch 2022. ×1.3.
- · `[[IR]]` — Interval/boundary-conditioned Recursive (av=2) = MacroIR. Sole calibrated survivor.

## concepts/ — type `concept`
- · `[[recursion-axis]]` — the occupancy-covariance-propagated axis (N vs R).
- · `[[window-axis]]` — the conductance/endpoint axis (0/1/2 endpoints); NON-monotone (MR worse than R).
- · `[[boundary-state]]` — the (i₀,iₜ) pair; static condensation on the time axis. (scoped, not "transition state")
- · `[[measurement-not-test]]` — exact-sim ground truth ⇒ measure a known misspecification, no null.

## claims/ — type `claim`
- · `[[novelty]]` — owner: `decisions/D-3`. Recursion-thesis + measurement + physical-identity; concede algorithm+apparatus.
- · `[[MR-sign]]` — owner: `decisions/D-4`. Observable-variance vs parameter-covariance; MR over-confident.
- findings (type `claim`): · `[[fisher-to-zero]]` · `[[cond-Nch2]]` (built evidence 2026-07-15) · `[[MR-worse-than-R]]`

## params/ — type `reference`
- units owner: `decisions/D-2`. · `[[Current_Baseline]]` (log-scale inheritance; NR degeneracy) · `[[Current_Noise]]` (dimensionless noise_in_conductance_tau) · others.

## axes/ — design; type `concept`
- owner: memory `project_design_axes_nondimensional`. · `[[N_ch]]` · `[[interval]]` (internal to cells) · `[[noise]]` · `[[P_open-deferred]]`.
