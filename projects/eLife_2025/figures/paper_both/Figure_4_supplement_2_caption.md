# Figure 4—figure supplement 2 caption

**Figure 4—figure supplement 2. The information distortion across all five parameters: diffuse for R, confined for IR.**

The second-moment companion to Figure 4—figure supplement 1, in the same layout: columns are the channel number N_ch across the four counts 10 to 10000, the horizontal axis inside each panel is the acquisition interval scaled by the closing rate, Δ·k_off, the vertical axis is the instrumental noise label, and rows are the parameter (outer) then the algorithm (inner), R above IR. The quantity is the information distortion, the ratio of the variance the score actually has to the variance the algorithm's Gaussian Fisher information claims, evaluated at the pooled optimum, shrunk to the bound of its 95% bootstrap interval nearest one, so a cell whose interval covers one is white. Above one (red) the reported interval is too narrow, below one (blue) too wide.

R's distortion is diffuse. The fraction of cells departing by at least 15% is 0.52 for the opening rate k_on, 0.68 for the closing rate k_off, 0.27 for the unitary current i, and 0.38 for the channel number N_ch. IR's is confined: 0.02, 0.11, 0.00 and 0.04 in the same four directions, so its largest departures, in k_off and N_ch, fade to zero with channel count while k_on stays negligible. The instrumental-noise direction σ_noise departs in neither algorithm, the shared reference that shows the shrinkage rule is not painting every cell. The interval conditioning localises the distortion to the kinetic and count directions and removes it from the rest, and this is what justifies Figure 4 showing k_off and N_ch. Figure 4—figure supplement 4 shows that what survives for R in each of these directions is the correlation part.

Six cells (R's k_on and N_ch at the ill-conditioned N_ch 10000, Δ·k_off = 1 corner, three noise levels each) are drawn grey rather than in the scale colour: their value falls below the 0.5 floor and the condition number there exceeds 3×10^4, so the direction is not identified and the magnitude is not meaningful. They are the same near-unidentified corner noted in the Figure 4 legend.

<!-- Source: projects/eLife_2025/figures/paper/figure_4_supplement_2.Rmd (shares figure_4_supp_common.R).
     Data: battery_pool_G, anchored at theta_pool, Gaussian anchor, figures/data/1c2ae6f and 87889e6
     by search path; 4 N_ch x 3 noise, 7 intervals inside each. Numbers printed by the caption-numbers
     chunk on each knit; re-read after a re-render. -->
