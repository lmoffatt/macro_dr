# Figure 3—figure supplement 2 caption

**Figure 3—figure supplement 2. The checks the body figure has no room for.**

All seven members of the lattice, in cost-ladder order, at the same cell and over the same 1,000 recordings as Figure 3, carrying the checks that figure does not draw: (A) the mean standardized residual, the one of the three properties of r_t that Figure 3 leaves out; (B) the score bias for the three parameters Figure 3 does not colour, the opening rate k_on, the unitary current i and the instrumental noise level σ_noise; and (C) their score autocorrelation. Rows, scales, bands, agonist shading and the calibrated reference line are as in Figure 3. Between the two figures each check is drawn once and none twice.

**(A) The first residual moment separates nobody, which is the row's claim.** Over the whole record the mean standardized residual runs from 0.0149 for least squares down to 0.0037 for the boundary-conditioned filter, monotonically along the ladder, and every member is within 0.015 of zero. Read with Figure 3's residual-variance row and its residual autocorrelation, this is the ordering of the three properties of r_t: the mean carries almost no verdict, the variance carries one, and the autocorrelation carries the most.

**(B) The bias in the three parameters the body does not draw, and the non-monotonicity along the ladder.** Counting informative intervals whose interval excludes zero against a chance level of 0.05: IR reads 0.067, 0.050 and 0.050, **INR reads 0.100, 0.100 and 0.000**, LSE 0.188 and 0.163, NR 0.475 to 0.510, and the three recursive members between them 0.409 and 0.975. The cheapest member that models the interval average is second only to the boundary-conditioned filter and ahead of all three recursive members, which cost more. On the noise direction it never fails at all.

That is the same ordering the per-interval information ratio gives (Figure 3—figure supplement 1), and together the two say that everything checkable one interval at a time is in order for INR. Its accumulated information ratio is nevertheless 20.8 and its lag-one score autocorrelation 0.87, so what it lacks is the recursion and the failure is entirely in the accumulation.

**The noise level is the control direction, and it behaves unlike the rest.** It is the only parameter measured before the agonist arrives, where there is no gating and the recorded variance is instrumental alone: it carries information at all 100 intervals, 20 of them pre-agonist, while k_on and i carry none before the pulse. In (C) its score is white even for the members whose kinetic directions are not, reading −0.004 at lag one for NR and −0.008 for R, both intervals covering zero, against 0.35 and 0.21 for those same members on k_on. The temporal memory that carries the accumulated information error therefore lives in the gating directions and not in the instrumental one, which is the conclusion Figure 4 reaches from the design plane, where the noise direction departs from calibration in neither recursive member.

**Least squares carries no noise curve, and the reason is structural.** Its noise level is not a free direction of the likelihood but a plug-in, estimated as the residual sum of squares over the same record, so the score in that direction vanishes identically: zero informative intervals, and an autocorrelation that is undefined rather than small. This is the same fact as its mean squared standardized residual being pinned at exactly 1.000, seen from the score side instead of the residual side. The curve is absent by construction and not by omission.

**(C) The autocorrelations.** At lag one LSE reads 0.833 and 0.836 on k_on and i, NR 0.347 and 0.478, INR about 0.87, R 0.207 and 0.195, and IR −0.013, −0.011 and −0.010, whose intervals narrowly exclude zero at this ensemble size. Figure 3's row G carries the same statistic for k_off and N_ch beside the residual autocorrelation, which reads 0.846, 0.510, 0.118 and −0.010 across the four members.

<!-- Source: projects/eLife_2025/figures/paper_both/figure_3_S3.Rmd (three rows, seven
     columns, three parameters, 7.0 x 5.8 in; widened from four columns on 2026-07-31 once
     macro_INR had run). Supersedes the residual triptych built earlier on
     2026-07-31, which was three quarters contained in the parent; what was new in it went two ways,
     the residual autocorrelation into Figure 3's row G and the residual mean into row A here. The
     superseded notebook is kept at tmp/triptych_superseded.Rmd until the thread closes.
     Data: figures/data/digest/figure_3_digest_{LSE,NR,INR,R,MR,VR,IR}.rds, re-extracted 2026-07-31 with all
     six param_index values so that the noise direction (3) is available; engine 0ffbda7, seed
     20260722. B = 400. Numbers printed by the caption-numbers chunk on each knit. -->
