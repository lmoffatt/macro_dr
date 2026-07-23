# Figure 4—figure supplement 1 caption

**Figure 4—figure supplement 1. The distortion-induced bias across all five parameters: diffuse for R, confined for IR.**

Figure 4 shows the bias in the two directions that carry it. This supplement shows all five parameters, because which directions are affected is itself a result. The plane and the conventions are those of Figure 4: columns are the channel number N_ch across the four counts 10 to 10000, the horizontal axis inside each panel is the acquisition interval scaled by the closing rate, Δ·k_off, over the seven swept values, the vertical axis is the instrumental noise label, and rows are the parameter (outer) then the algorithm (inner), R above IR, so each parameter's two rows are adjacent. The quantity is the distortion-induced bias in log10 units, evaluated at the true parameter (a bias evaluated at the optimum is zero by construction), shrunk to the bound of its 95% bootstrap interval nearest zero, so a cell whose interval covers zero is white.

R carries a bias in every direction. The fraction of cells whose interval excludes zero is 0.79 for the opening rate k_on, 0.57 for the closing rate k_off, 0.88 for the unitary current i, 0.87 for the channel number N_ch, and 0.75 for the instrumental noise σ_noise. IR sits at a floor of about 0.13 in four of the five and rises to 0.24 in the instrumental-noise direction. The interval conditioning removes the bias from four of the five directions rather than only reducing it. The instrumental-noise row is the shared reference: it is the one direction where IR too carries an appreciable bias, so both algorithms depart from the null there.

Six cells are drawn grey rather than in the scale colour, all R at N_ch 10000 and Δ·k_off = 1: three in the i direction and three in the N_ch direction, one per noise row. These are the near-unidentified direction described in the Figure 4 legend, where the Gaussian Fisher covariance has a condition number of about 10^5. A cell is greyed only where its value is off the scale and the condition number exceeds 3×10^4, so a genuinely large but well-conditioned departure would still be shown; the sign is meaningful and the magnitude in those cells is not.

<!-- Source: projects/eLife_2025/figures/paper/figure_4_supplement_1.Rmd (shares figure_4_supp_common.R).
     Data: battery_sim_G, anchored at theta_sim, Gaussian anchor, figures/data/1c2ae6f and 87889e6 by
     search path; 4 N_ch x 3 noise, 7 intervals inside each. Numbers printed by the caption-numbers
     chunk on each knit; re-read after a re-render. Open across Figure 4 and its supplements: whether
     the kappa ~ 1e5 cells are greyed as unidentified or shown as bias/SE. -->
