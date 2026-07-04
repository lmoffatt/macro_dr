# Figure 4 (alt) caption

**Figure 4 (alternative). The per-step information profiles alone, without the per-step calibration rows.**

The information rows of Figure 4 shown on their own, one row per parameter (k_on, k_off, i, N_ch) across the five approximations (NR, MNR, R, MR, IR). Each panel carries the expected per-step Fisher information F_t = E[I_t] (blue) and the score variance J_t = Var(s_t) (orange) on a log10 axis, over the same 1,000 recordings and the same cell as Figure 4; the grey band marks the agonist pulse and time is in milliseconds. Dropping the per-step ratio rows (which sit on the calibrated value for all five approximations) leaves the information geometry over time on its own: the recursive filters stop gaining information about k_on, i and N_ch when the agonist is removed, having conditioned on the past, while the closing rate k_off stays informative through the decay. See the Figure 4 legend for the full description.
