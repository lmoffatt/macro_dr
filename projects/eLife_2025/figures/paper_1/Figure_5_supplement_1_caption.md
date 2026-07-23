# Figure 5—figure supplement 1 caption

**Figure 5—figure supplement 1. The per-recording precision behind the design trade-off: the achievable variance of each parameter in a single recording.**

The same measurement as Figure 5, with the channel-count scaling removed. The vertical axis is the distortion-corrected variance of the parameter directly, on a log scale, which is the variance of the estimate in one recording of N_ch channels, that is the precision the experiment achieves before any pooling. The plane is otherwise identical: the horizontal axis is the acquisition interval Δ·k_off, colour is the channel count N_ch over the four counts 10 to 10000, rows are the five parameters and columns are the instrumental noise label from 0.05 to 10. This and Figure 5 are the same data; on the log axis they differ only by a vertical shift of log(N_ch) per line, so each curve keeps its shape and only the vertical ordering of the channel-count lines changes meaning, from pooled variance at fixed budget in Figure 5 to per-recording variance here.

Two things read directly here that Figure 5's budget scaling hides. First, the ordering. More channels give a smaller per-recording variance in every direction, so the N_ch = 10000 line sits at the bottom of each panel, which is the intuitive statement that a larger patch measures better. For the rates the lines are evenly spaced by about a decade per decade of channels, the signature of information linear in channels; for the amplitude and count directions the lines bunch together at high channel count, the signature of the sub-linear information that makes concentrating costly in Figure 5. Second, the instrumental-noise row nearly collapses: the per-recording variance of σ_noise spans only a factor of 3.2 across the four channel counts (0.0086 to 0.028 at noise 0.1), so the noise is measured about equally well at any channel count. Figure 5 renders that same fact as a misleading factor of 3200, because multiplying a near-constant per-recording variance by N_ch inflates it with the channel count.

The vertical bars are the 95% bootstrap interval. As in Figure 5, only IR is shown, and the covariance is the distortion-corrected one at the pooled optimum.

<!-- Source: projects/eLife_2025/figures/paper/figure_5.Rmd, the Figure_5_supplement_1.pdf output
     (direct corrected covariance, no N_ch scaling). Same data and filters as Figure 5.
     Quantity: Probit_statistics_Gaussian_Distortion_Corrected_Covariance diagonal at theta_pool.
     Numbers printed by the caption-numbers chunk on each knit; re-read after a re-render. -->
