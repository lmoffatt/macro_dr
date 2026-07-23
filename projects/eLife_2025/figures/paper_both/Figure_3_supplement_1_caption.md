# Figure 3—figure supplement 1 caption

**Figure 3—figure supplement 1. Where in time each parameter is measured, and the per-step calibration of that information, along the recursive ladder.**

For the same 1,000 recordings as Figure 3 (a patch of N_ch = 100 two-state channels, instrumental noise 1e-4, measurement interval 0.1 in units of the model's characteristic time τ = 1/k_off), the per-step information about each parameter is resolved in time. Rows come in pairs, one pair per parameter (the opening and closing rates k_on and k_off, the unitary current i, and the channel number N_ch); the columns are the four rungs of the recursive ladder, R (Recursive), MR (Mean Recursive), VR (Variance Recursive) and IR (Interval Recursive). All four are recursive, so every difference below is attributable to the interval conditioning alone. The grey band marks the agonist pulse and time is in milliseconds.

The taller row of each pair is the information profile, on a log10 axis: the expected per-step Fisher information F_t = E[I_t] (blue), computed as I_t = (∂ y_mean/∂θ)² / y_var + ½ (∂ y_var/∂θ)² / y_var² from the predictive mean y_mean and variance y_var of the interval, and the variance of the per-step score J_t = Var(s_t) over the recordings (orange). The two coincide wherever the filter is calibrated. The shorter row is their per-step ratio log10(J_t / F_t), on a linear axis with a green line at zero (the calibrated value); the shaded band is a 90% bootstrap interval over the recordings, and a point is drawn only where the parameter carries information (F_t above a small floor).

**When each parameter is measured.** The information about k_on, i and N_ch drops to zero the moment the agonist is removed: the filter has already extracted it during the pulse and, conditioning on the past, gains nothing further from watching the channels relax. The closing rate k_off is the exception, staying informative through the decay, which the closing kinetics govern, and the unitary current keeps a conductance-noise floor throughout. This is the analytic signature of a conditional filter: it registers when it has finished measuring a quantity. It holds for all four rungs, so it is a property of the recursion rather than of the interval conditioning.

**The per-step ratio is mildly but systematically below one for every rung except IR.** Median log10(J_t / F_t) over the informative steps is −0.18 (k_on), −0.14 (k_off), −0.19 (i) and −0.18 (N_ch) for R; −0.25, −0.23, −0.28 and −0.25 for MR; −0.15, −0.13, −0.17 and −0.15 for VR; and +0.006, +0.033, +0.014 and +0.005 for IR. Only IR's interval covers zero at most steps (64 to 91% of them, against 0 to 28% for the other three). So per step the score varies **less** than the information predicts, by up to a factor of about 1.5 for MR.

**The sign flips on accumulation, and that is the point.** Per step J_t < F_t, yet over the whole recording J_T / F_T rises above one (Figure 3E: MR 1.19, VR 1.20, against IR 1.08). A discrepancy that changes sign between the per-step and the accumulated statistic cannot be a per-sample modelling error; it is what positive cross-time correlation of the score does, and Figure 3F shows that correlation directly. This figure is the premise of that argument: the per-step magnitudes are mildly off and bounded, far too small to account for the accumulated behaviour.

The two non-recursive members of the family, NR and MNR, are computed in the same run and stored with these data, but are measured in no figure of this paper; their much larger accumulated mismatch belongs to the companion paper.

<!-- Source: projects/eLife_2025/figures/paper/figure_3_supplement_1.Rmd (renamed 2026-07-22 from
     figure_4.Rmd on renumbering; this is now a supplement of Figure 3, not a body figure).
     Roster ALGOS <- c("R","MR","VR","IR").
     Every number above is printed by the notebook's caption-numbers chunk on each knit; re-read it
     after any re-render rather than copying these values forward.
     Data: figures/data/figure_3_time_dlik_{R,MR,VR,IR}.csv, engine 0ffbda7, regenerated 2026-07-22.
     NOTE: the previous caption said the per-step identity "holds for all five approximations (the
     ratio sits on zero within the bootstrap interval)". The numbers do not support that and never
     did; the defensible claim is the one now written above, mild and bounded, with IR the only rung
     statistically indistinguishable from one. -->
