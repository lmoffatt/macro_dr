# fact-MR-overconfident — claim

> type: `claim` (finding) · group: `facts` · governed by `../_SPEC.md`
> A fact touching ≥2 dimensions ⇒ its own node (`_SPEC` Rule 1). Dimensions it references list it; it does not live inside them.

## Claim
MR under-reports parameter uncertainty (over-confident): empirical/Fisher parameter-covariance ratio **1.5–2.1** (geom-mean 1.80), same side as NMR.

## Grounds (evidence → sources)
`[[data-1c2ae6f]]` MR battery cells, recomputed 2026-07-15. `computed-via: [[sandwich]]` (emp cov vs Gaussian-Fisher cov).

## Warrant
ratio > 1 ⇔ the reported Fisher covariance is smaller than the empirical ⇔ over-confident. Standard sandwich reading.

## Qualifier (scope/strength)
noise 0.1, canonical N_ch. `[OPEN]` full-grid/other-noise pending Gaussian fill. Sign vs the observable-variance statement: see `[[MR-sign]]`.

## Rebuttal (attack surface)
Do not conflate with the *observable* variance direction (`nomenclature.md:52` "overestimates variance") — different object; that is `[[MR-sign]]`.

## Status / Provenance
`live` · decided `decisions/D-4` (pending Luciano approval).

## Dimensions-referenced
`[[MR]]` · `[[NMR]]` · `[[window-axis]]` · calibration.
