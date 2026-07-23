# Figure 4—figure supplement 3 caption

**Figure 4—figure supplement 3. The information distortion split into a per-sample part and a correlation part, on the closing rate: R's floor is the correlation.**

This supplement explains why more channels do not rescue R (Figure 4C). The information distortion in the closing rate k_off is decomposed into two factors: a per-sample part, the distortion of the single-interval Gaussian approximation, and a correlation part, the distortion contributed by the statistical dependence between successive intervals. The layout follows Figure 4: columns are the channel number N_ch across the four counts 10 to 10000, the horizontal axis inside each panel is the acquisition interval scaled by the closing rate, Δ·k_off, the vertical axis is the instrumental noise label. Rows are the component (total, sample, correlation) then the algorithm (R, IR). Colour is the distortion, shrunk to the 95% bootstrap bound nearest one, red above one and blue below.

Read down the R block. Its per-sample part goes to one as the channel count grows (median over the intervals 1.09, 1.07, 1.02, 1.00 across N_ch 10, 100, 1000, 10000 at noise 0.1), so the single-interval Gaussian approximation becomes exact, exactly as it does for IR (1.10, 1.07, 1.04, 1.01). What stays elevated for R is the correlation part (1.35, 1.36, 1.43, 1.47), which is IR's only where the channels are few and then falls to one (1.24, 1.02, 1.02, 1.00). R and IR therefore agree on per-sample fidelity and differ on the correlation between intervals, which is what conditioning on the interval endpoints removes. R's total distortion is its correlation part, and that is the flat floor of Figure 4C. The values quoted are raw medians of the estimate, matching Figure 4—figure supplement 4; the map itself colours the conservative bound of the 95% bootstrap interval, so a cell near one reads slightly closer to one in colour than in the number.

The three components are read multiplicatively, total close to sample times correlation. This is an empirical near-identity on this data, not an algebraic one (the exact relation is a matrix congruence), and it holds to a median relative error of 0.03% and a maximum of 1.7% over the volume shown. The exactly additive statement is that the log-determinant of the total distortion equals the sum of the log-determinants of the two parts.

<!-- Source: projects/eLife_2025/figures/paper/figure_4_supplement_3.Rmd.
     Data: battery_pool_G, theta_pool, Gaussian anchor, figures/data/1c2ae6f and 87889e6 by search
     path; k_off only; 4 N_ch x 3 noise, 7 intervals. Components:
     Likelihood_Gaussian_Information_Distortion (total), Gaussian_Sample_Distortion (sample),
     Likelihood_Correlation_Distortion (correlation). Do NOT use
     Likelihood_Information_Distortion_Reconstituted (it implements the false product identity and is
     NaN here). Numbers printed by the caption-numbers chunk on each knit; re-read after a re-render. -->
