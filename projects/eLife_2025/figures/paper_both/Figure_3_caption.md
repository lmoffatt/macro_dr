# Figure 3 caption

**Figure 3. The calibration cascade in time, over the whole implemented family.**

One common set of 1,000 simulated recordings (a patch of N_ch = 100 two-state, closed ⇌ open, channels, instrumental noise 1e-4, acquisition interval 0.1 τ with τ = 1/k_off) is scored by each of the seven members, with no fitting: least squares (LSE), the open-loop gating likelihood (NR), the interval non-recursive member (INR), the recursive filter on instantaneous samples (R), the two one-endpoint recursive members (MR, VR) and the boundary-conditioned filter (IR), in cost-ladder order. Each row is a time-resolved test of one calibration property. Unless stated otherwise the quantity is a mean over the 1,000 recordings and the band or vertical bar is its 95% nonparametric bootstrap interval, B = 400; the green line is the value a calibrated member would give and the grey band marks the agonist pulse. Rows A and B are properties of the data and are drawn neutral; rows C to G are per parameter, k_off in blue and N_ch in vermillion. Time is in milliseconds; row G is resolved in lag.

**Read across the columns, the figure is a factorial on the two switches of the lattice**, interval averaging against recursion, rather than a list. Interval averaging alone (INR) repairs the data level and does nothing for the accumulation: mean squared standardized residual 0.9995 on row B, against an accumulated information ratio of 20.8 on row F. Recursion alone (R) does the opposite, 0.782 against 1.07. Only the two together repair both, 0.999 and 1.08 for IR, while the member with neither (NR) fails both, 1.153 and 13.6. Least squares sits outside the lattice as the classical baseline.

The two one-endpoint members are the other half of that reading. They do not lie between R and IR: on the ellipse-area measure of Figure 2 the recursive instantaneous filter reads 1.32 on the kinetic pair while MR reads 1.97 and VR 2.19, and only then does IR return to 1.02. Conditioning the interval-averaged conductance on one endpoint therefore costs calibration rather than buying it, and what recovers it is the second endpoint entering the gain. VR, which removes the between-endpoint part of the variance without restoring it in the gain, is worse than MR, which localises the whole R-to-IR distance in the update rather than in the variance.

(A) The output. One example recording, the same one in every column: the predicted current (mean, with the member's own predictive ±1 s.d. band) over the observed current, on a scale shared across columns. The annotation is the ensemble-mean total log-likelihood over the 1,000 recordings with its standard error, and read across the row it is the whole cost ladder: −262.3 ± 0.7 for LSE, −230.2 ± 0.8 for NR, −215.4 ± 0.6 for INR, −147.2 ± 0.2 for R, −153.9 ± 0.2 for MR, −148.9 ± 0.2 for VR and −143.4 ± 0.2 for IR.

(B) Residual variance: the mean squared standardized residual, r² with r = (observed − predicted mean)/predicted s.d., one for a calibrated member, on a log axis.

(C) Per-interval Fisher information F_t, as log10 on a linear axis and against no reference: this row says where in the recording each parameter is measured at all, and for how long.

(D) Score bias, standardized: the mean per-interval score over √F_t, zero when the score is unbiased, dimensionless so that two parameters with different units share an axis.

(E) Per-interval information ratio Var(s_t)/F_t, one when the score variance matches the information at each interval.

(F) Accumulated information ratio Var(Σ s_t)/Σ F_t = J_T/F_T, one when the same equality holds over the whole recording.

(G) Autocorrelation against lag, three series: the standardized residual in black and the score for each parameter in its colour. The score's is why (E) can hold where (F) does not, since a correlated score makes the variance of the sum exceed the sum of the per-interval informations. The residual's is the same memory at the level of the data, and it is the only quantity in this figure an experimentalist can compute on a real recording, needing no ensemble and no known truth. The two agree in order and not in magnitude: at lag one the score reads about 1.6 to 1.7 times the residual for the members in which either is large (NR 0.864 against 0.510, R 0.191 against 0.118), and the two coincide only for least squares, 0.856 against 0.846, whose score is proportional to its residual by construction.

**Scale conventions.** Any non-negative quantity is on a log10 axis with the calibrated value on the ticks; signed quantities are on a linear axis with zero on the ticks. Row C is the exception, log10(F_t) on a linear axis, because it never passes near one. Rows B, D, F and G use a separate scale per block, the three open-loop columns against the four recursive ones, because the first fail by a factor of twenty and the second by twenty percent and one axis cannot show both.

**The channel number is not drawn in the least-squares column**, by exclusion and not by absence of data. Least squares does not model the gating fluctuations, the only thing that separates the channel count from the unitary current, so with them unmodelled the two enter only through their product and its information about N_ch would be an artefact of the parameterisation rather than something it delivers. F_t over the pulse is 1076.55 for both, to every digit. figure_4_common.R applies the same exclusion through EXCLUDE_ROWS.

<!-- Source: projects/eLife_2025/figures/paper_both/figure_3.Rmd. Promoted 2026-07-31 from what had
     been Figure 3—figure supplement 1: the body walked four members and the supplement seven, and
     holding two figures apart for two columns stopped making sense once every column carried a
     claim the paper already makes in prose. The four-column predecessor, the only notebook that
     builds this figure straight from the ~1 GB dumps, is archived at
     figures/archive/figure_3_4col_superseded_20260731/ and stays runnable as a check on the digests.
     Data: figures/data/digest/figure_3_digest_{LSE,NR,INR,R,MR,VR,IR}.rds, written by
     figure_3_digest.R from engine 0ffbda7, regenerated 2026-07-22 with seed = 20260722; the INR
     dump is the 2026-07-31 run of macro_INR, and NMR, the build that lost the N*ms interval-variance
     term, is on disk and drawn nowhere. Every number above is printed by the notebook's
     caption-numbers chunk on each knit; re-read it after a re-render rather than copying forward. -->
