# Figure 6 caption

<!-- 2026-08-26: the slope cited below was −0.56, which is the midpoint of two candidates and
     matches neither. figure_7.html reports −0.50 for BH, the boundary Figure 7 actually draws, and
     −0.62 for `floor`/k_off in the geometry chunk's table, which is a different quantity. Corrected
     to −0.50 here and in the Figure 7 subsection of 04_results.tex in the same pass; the arithmetic
     of the sentence works better with it, one decade of channels being worth half a decade of noise
     against the 0.65 and 0.3 decades of the sampling it defends. -->

<!-- RENUMBERED 2026-08-10: was Figure 5. The tau_int decision panel entered the body ahead of
     this figure as Figure 5 and the usage map moved to Figure 7. -->

**Figure 6. Where the boundary-conditioned filter stops reporting an honest error bar, and which of the two discarded moments is responsible.**

Figure 4 leaves IR's panels almost empty, which is the result but also a question: the filter is calibrated across most of the design plane, so where does it stop being so, by how much, and why. This figure answers on the corner where the departure lives. It is the same measurement as Figure 4, the Gaussian information distortion evaluated at the pooled optimum, on the same runs; only the resolution differs, because on Figure 4's scale a departure of a factor 1.5 is a faint tint and here it needs an axis.

The three blocks are the distortion and the two parts it decomposes into: (A) the total, (B) the part carried by the per-sample non-Gaussianity of the score, and (C) the part carried by its correlation across intervals. Each block is a grid of the design plane in the same orientation as Figures 4 and 7: the channel count N_ch increases to the right, the dimensionless instrumental noise increases upward, and inside every cell the horizontal axis is the dimensionless sampling interval Δ·k_off, from 0.01 to 1, shared by every panel and labelled once per block. The vertical axis is the distortion factor on a log scale, shared by all twelve columns so that they can be compared, with the green line at 1 marking a report that is neither over- nor under-confident. Colour is the parameter, the closing rate k_off and the channel number N_ch, the two directions that carry the departure. Bands are the 95% percentile bootstrap interval over recordings; their median half-width is 0.03, so the departures below are many band widths from 1. The channel count is sampled in steps of about 0.65 decades against 0.3 decades of noise, which matches the −0.50 slope of the region-0 boundary of Figure 7: noise is roughly twice as potent per decade, so it needs twice the resolution.

**The departure is two-sided within a single cell.** At the fewest channels and the lowest noise, the reported error bar of k_off is too small by a factor of 1.62 at the shortest interval and 1.52 at the longest, while that of N_ch is too large by a factor of 0.65 at the shortest interval and is back at 1.07 by the longest. The same likelihood is therefore over-confident about the rate and conservative about the count at the same time, which is why a single scalar summary of miscalibration would report almost nothing here. Both departures grow monotonically as channels are removed and as noise falls, and at a dimensionless noise of 1 the whole grid is flat, which is the instrumental noise making the emission more Gaussian and relaxing the approximation, as the Theory section predicts.

**The two parts respond to the sampling interval in opposite directions.** At the worst cell the correlation part of k_off falls from 1.52 at the shortest interval to 1.27 at the longest, while the per-sample part rises from 1.06 to 1.19. Measured as the ratio of the departure at the shortest interval to that at the longest, the correlation part is 1.9 for k_off and 11.9 for N_ch, while the per-sample part is 0.31 for k_off. Sampling more coarsely therefore buys back the correlation part and not the per-sample one. The two parts compose: their product reproduces the total to within 1% over the whole grid (median 1.000, 5th to 95th percentile 0.998 to 1.008), so panel A can be read as the product of B and C.

**What this figure cannot decide.** The macroscopic family applies both of its Gaussian approximations at every step and the score sees only their product, so a departure measured here cannot be assigned to the occupancy approximation or to the interval-signal one; separating them needs the microscopic recursion, which is a companion paper's subject. And the grid is truncated on the side where the effect grows: at the lowest simulated noise the departure of k_off is still accelerating (1.22, 1.37, 1.52 over the last three noise steps at ten channels), so the figure locates the onset of the departure but does not close its far side, and no channel count can be quoted as the threshold below which the approximation fails.

<!-- Source: projects/eLife_2025/figures/paper_both/figure_6.Rmd (was figure_5.Rmd). Data: figures/data/1c2ae6f, one
     run for all 24 drawn cells (checked: no provenance seam). Component paths are figure_4_common.R's
     COMPS, i.e. Probit_statistics_Likelihood_Gaussian_Information_Distortion for A,
     Probit_statistics_Gaussian_Sample_Distortion for B and
     Probit_statistics_Likelihood_Correlation_Distortion for C, all at battery_pool_G (theta_pool).
     Lines are the bootstrap mean `m`, NOT `Dconf` (figure_4_common.R:200), which is the conservative
     CI edge collapsed to 1 and belongs to the maps only.
     Numbers: worst cell = N_ch 10, dimensionless noise 0.005. Noise displayed = swept label / 10
     (NOISE_AXIS_UNITS.md). Interval ratios and the composition check are in
     tmp/fig5_layouts/analysis.R and critique.R. -->
