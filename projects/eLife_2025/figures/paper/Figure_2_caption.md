# Figure 2 caption

**Figure 2. Uncertainty calibration and bias along the interval-likelihood ladder.**

Each simulated recording is a patch of N_ch = 100 two-state (closed ⇌ open) channels. From it the four model parameters are estimated by maximum likelihood: the opening and closing rates k_on and k_off, the unitary current i, and the channel number N_ch, all on a base-10 log scale. The estimation is repeated over 10,000 independent recordings under each rung of the recursive ladder: R (Recursive), MR (Mean Recursive), VR (Variance Recursive) and IR (Interval Recursive). The four differ only in what the likelihood conditions the interval-averaged conductance on, so the columns read left to right as one added piece of interval structure at a time. Every panel is centred on the true parameter value, marked by a cross.

(A) The cloud of estimates (grey points) for two parameter pairs, one per row: the kinetic pair (log10 k_on, log10 k_off) and the amplitude pair (log10 N_ch, log10 i). Columns are the four approximations. In each panel three 95% ellipses are drawn at the cloud mean: the empirical covariance of the cloud (grey), which is the estimator's true spread; the covariance the approximation reports from its Fisher information (orange, dashed); and the distortion-corrected sandwich covariance (blue, dotted). The two coloured numbers are ellipse-area ratios. The orange number is empirical over Fisher, the factor by which the Fisher information underestimates the true spread (the overconfidence factor). The blue number is empirical over corrected, near one wherever the sandwich correction recovers the spread. Two further markers report the bias: the filled circle is the empirical mean of the cloud, and the open circle is the truth plus the bias predicted by the distortion theory. The distance from the cross to the filled circle is the realised bias, and the open circle sitting on the filled circle means the theory predicts that bias.

The overconfidence factor is 1.32 for R, 1.97 for MR, 2.18 for VR and 1.02 for IR on the kinetic pair, and 1.09, 1.53, 1.77 and 1.00 on the amplitude pair. The sandwich correction returns every one of them to within a tenth of one. **VR is the most over-confident member of the ladder, more so than MR**, which is the panel's central result: VR removes the part of the predicted observable variance that conditioning on the interval endpoint would explain, without taking the boundary term in the gain that would justify removing it, and the reported uncertainty degrades rather than improves. The step that recovers calibration is therefore the gain, not the variance.

(B) The full correlation structure of IR, the approximation that calibrates, over every pair of the four parameters, with the same ellipses, markers and numbers as in (A). The three ellipses coincide and the numbers read 1× in every pair (the second moment is calibrated), and the three centre markers coincide on the cross (the first moment is unbiased).

The true values are k_on = 10, k_off = 100, i = 1 and N_ch = 100. Fisher and corrected covariances are the Gaussian-Fisher family, evaluated at the pooled estimate and scaled by the reciprocal of the group size (10 recordings), so they are comparable to the empirical spread of the cloud. The instrumental noise is 1e-4 (Current_Noise; the "0.1" that labels the data is the noise on the conductance-τ scale, a separate parametrisation) and the measurement interval is 0.1 in units of the model's characteristic time τ.

The two non-recursive members of the family, NR and MNR, are named in the Theory section and measured in no figure of this paper.

<!-- Source: projects/eLife_2025/figures/paper/figure_2.Rmd on the GAUSSIAN anchor
     (ANCHOR <- "gaussian"): R, MR, IR from figures/data/1c2ae6f, VR from figures/data/0ffbda7,
     all at nsim 10000, so no n_sim mismatch inside any panel.
     Ratios recomputed from the battery_pool_G files 2026-07-22; if the figure is re-rendered on
     other data these five numbers move and must be re-read, not copied. -->
