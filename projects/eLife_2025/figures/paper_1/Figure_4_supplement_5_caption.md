# Figure 4—figure supplement 5 caption

**Figure 4—figure supplement 5. The distortion-corrected covariance against the empirical multivariate distribution of the estimator: a Mahalanobis quantile-quantile plot.**

Figure 4's block A compares, one parameter axis at a time, the empirical spread of the estimates against the covariance the algorithm reports and against the distortion-corrected covariance. A per-axis ratio is a projection of the joint distribution; this supplement asks whether the reported covariance matches the full multivariate distribution of the estimator. The standard way to reduce a p-dimensional distribution to one scalar with a known reference is the squared Mahalanobis distance of each estimate to the cloud mean, d² = (θ̂ − θ̄)ᵀ Σ⁻¹ (θ̂ − θ̄), which is distributed as χ² with p degrees of freedom when θ̂ is Gaussian with covariance Σ. The distance is taken to the empirical mean rather than the truth, so it tests the second moment alone, free of the bias.

Each panel is a quantile-quantile plot of the sorted empirical d² against the χ²_p quantiles, for the two candidate covariances: the covariance the algorithm reports from its Gaussian Fisher information (orange), and the distortion-corrected sandwich covariance (blue). A curve on the diagonal means the covariance matches the empirical distribution. The reported Fisher curve rises above the diagonal in the upper tail, so the estimates are more dispersed than the reported covariance claims, and its nominal 95% region covers only 91%. The sandwich curve lies on the diagonal and its 95% region covers 95%. This is the multivariate form of the over-confidence Figure 4 shows per axis, and it identifies the sandwich, not the Fisher information, as the covariance the estimator actually has.

The effect is mild in six dimensions, and this is the same marginalisation as block A's ellipse. The joint Mahalanobis distance averages over all directions, so the strong per-axis over-confidence in the closing rate (a factor of about 1.5 in Figure 4) is diluted against the conservative channel-number direction, and the joint under-coverage is 91% against the nominal 95% rather than a larger gap. The per-direction magnitudes are the subject of Figure 4; this supplement shows that the joint distribution is the one the corrected covariance predicts. The framework is the misspecified-likelihood sandwich covariance and its information matrix test (Huber 1967; White 1982); the Mahalanobis reduction to χ²_p is the standard validation of a robust covariance against a simulated estimator distribution.

Both parameters and both intervals give the same reading. The panels are the two acquisition intervals, Δ·k_off = 1 and 0.05.

<!-- Source: projects/eLife_2025/figures/paper/figure_4_supplement_5.Rmd.
     Cloud: figures/data/0ffbda7/figure_3_G_nch_10_nsim_100000_macro_IR_noise_0.05_mle_cloud_runs.csv,
     group size 100 (1000 fits), the only high-fits cloud on disk; NOT block A's noise-0.1 design
     point, because 1000 fits are needed for a clean Q-Q at p = 6 and that run exists only here.
     Covariances (Gaussian_Fisher_Covariance and Gaussian_Distortion_Corrected_Covariance) from the
     matching nsim-10000 battery in 1c2ae6f, scaled by the group size. p = 6 (on, off,
     unitary_current, Current_Noise, Current_Baseline, Num_ch_mean). Coverage: Fisher 0.91, sandwich
     0.95 at both intervals; median d^2 sandwich 5.24-5.28 against chi2_6 median 5.35. Numbers printed
     by the caption-numbers chunk on each knit; re-read after a re-render.
     Citations (Huber 1967, White 1982) are given from memory and must be verified against the PDFs
     before the manuscript quotes them. -->
