# Figure 4—figure supplement 4 caption

**Figure 4—figure supplement 4. Sample and correlation distortion against channel count, resolved by noise: R and IR agree on the per-sample part and split on the correlation part, and raising the noise cures R only at few channels.**

The mechanism of Figure 4—figure supplement 3, generalised to the two carrier directions and to the noise axis, and shown as curves because it is a dependence on the channel count. The horizontal axis of every panel is the channel number N_ch on a log scale across the four counts 10 to 10000, and the vertical axis is the distortion, one for a calibrated filter. The recursive filter R (Recursive) and the interval-conditioned filter IR (Interval Recursive) are both drawn in each panel; the band is the interquartile range over the seven acquisition intervals. Rows are the parameter (the closing rate k_off, the channel number N_ch) each split into its two components (the per-sample part, the correlation part); columns are the instrumental noise label, 0.1, 1 and 10.

Two readings sit in one figure. First, the split. In the sample rows the R and IR curves lie on each other and both approach one as the channel count grows, so the two algorithms do not differ in per-sample fidelity. In the correlation rows R stays elevated and rises with the channel count while IR descends to one, so the difference between the two algorithms is the correlation between successive intervals. In k_off at noise 0.1 R's correlation distortion runs 1.35, 1.36, 1.43, 1.47 across the four channel counts while IR runs 1.24, 1.02, 1.02, 1.00; in N_ch it runs 0.99, 1.26, 1.34, 1.34 for R against IR near one. More channels push the per-sample part toward one and the correlation part up by a matching amount, which is why the total stays flat and more channels do not help R.

Second, the effect of noise. Reading a correlation row left to right, R's distortion starts lower at higher instrumental noise, because at few channels the instrumental noise dominates the gating variance and there is little to distort, but it re-emerges toward 1.3 to 1.5 by N_ch 10000 at every noise level. For k_off it runs 1.35 to 1.47 at noise 0.1, 1.08 to 1.43 at noise 1, and 1.03 to 1.37 at noise 10. Raising the instrumental noise cures R only where the channels are few; because the gating variance grows with the channel count, the correlation distortion returns at high channel count regardless of the noise. This is the quantitative form of the boundary between the regime where a likelihood is needed and the regime where instrumental noise removes the need.

<!-- Source: projects/eLife_2025/figures/paper/figure_4_supplement_4.Rmd (shares figure_4_supp_common.R).
     Data: battery_pool_G, theta_pool, Gaussian anchor, figures/data/1c2ae6f and 87889e6 by search
     path; components Gaussian_Sample_Distortion and Likelihood_Correlation_Distortion; k_off and
     N_ch; 4 N_ch x 3 noise, median and IQR over the 7 intervals. Numbers printed by the
     caption-numbers chunk on each knit; re-read after a re-render. The delta cut was left out
     deliberately: it is non-monotone and carries the kappa ~ 1e5 corner, so the interval dependence
     stays in supplement 2's maps. -->
