# Figure 5 caption

**Figure 5. The design trade-off: how the precision of the interval-conditioned filter trades against channel count and acquisition interval, at a fixed channel budget.**

Figures 2 to 4 establish which likelihood to use, IR (Interval Recursive), and that its reported uncertainty is honest. This figure uses that honesty to advise the experiment. Only IR is shown, because the quantity below is a covariance the algorithm reports, and only a calibrated algorithm's reported covariance can guide design. The covariance used is the distortion-corrected one, and for IR it is close to what the algorithm reports directly, so it is the precision the experiment actually achieves rather than an optimistic claim.

The quantity is a design trade-off, not an atlas of precision. The vertical axis is N_ch times the distortion-corrected variance of the parameter, on a log scale. At a fixed total channel budget split into recordings of N_ch channels each, the variance of the pooled estimate is N_ch times the single-recording variance divided by the budget, so this quantity is that pooled variance up to the fixed budget constant. A line that is flat across N_ch means the information is linear in channels and concentrating the budget into few large recordings costs nothing; a line that rises means concentrating loses information, so the same budget spent as many small recordings gives a better estimate.

The horizontal axis of each panel is the acquisition interval scaled by the closing rate, Δ·k_off, on a log scale. Colour is the channel count N_ch, the design choice, over the four counts 10 to 10000. Rows are the parameter (the opening and closing rates k_on and k_off, the unitary current i, the channel number N_ch, and the instrumental noise σ_noise). Columns are the instrumental noise, given as the dimensionless label used in the run, from 0.05 to 10 (the physical Current_Noise is that label divided by 1000). Noise is on columns because it is set by the recording rig rather than chosen, so the columns show how the advice shifts with the rig while the design choice itself is the within-panel spread over N_ch. Vertical bars are the 95% bootstrap interval.

**The rates and the amplitude directions give opposite design advice.** The closing rate k_off is nearly flat across channel count (its pooled variance at noise 0.1 runs 0.32, 0.38, 0.40 and 0.40 across N_ch 10, 100, 1000 and 10000, a factor of 1.24 end to end), so its information is close to linear in channels and concentrating is close to neutral. The opening rate k_on is similar (a factor of 3.1). The amplitude, count and noise directions rise steeply: over the same channel range the pooled variance grows by a factor of 288 for i, 60 for N_ch and 3200 for σ_noise. For these directions, concentrating the channel budget into few recordings loses information heavily, and spreading it across many recordings is the better design. The per-recording precision that lies behind these numbers, which orders the opposite way, is in Figure 5—figure supplement 1.

**A noisier rig penalises few channels.** Reading a row across the noise columns, k_off stays flat until the highest noise, where at few channels its pooled variance rises (to 0.81 at N_ch 10, noise 10, against 0.32 at noise 0.1), because there the instrumental noise is large next to the small-N_ch gating variance. The amplitude directions keep their steep rise across the whole noise range.

The interval trade-off runs the same way in every panel: the pooled variance is smallest at short intervals and grows toward the coarse end, mildly for the rates and more for the amplitude directions, so a short acquisition interval is the better design for precision, consistent with the calibration trade-off of Figure 4.

<!-- Source: projects/eLife_2025/figures/paper/figure_5.Rmd (migrated from figure_6_precision.Rmd,
     which was on the numeric anchor 433ed13 with five algorithms). This is the Gaussian anchor, IR
     only, all six noise levels as columns. Body = Figure_5.pdf, the N_ch x cov (pooled-variance)
     version; the direct-covariance version is Figure_5_supplement_1.pdf.
     Quantity: Num_ch x Probit_statistics_Gaussian_Distortion_Corrected_Covariance diagonal, at the
     pooled optimum theta_pool. Data figures/data/1c2ae6f (and 87889e6 for noise 1 and 10) by a search
     path; each CSV stamps its own engine hash. Cells with a non-finite corrected covariance or a
     bootstrap CI spanning more than a factor of 8 drop out of the lines. Every number is printed by
     the notebook's caption-numbers chunk on each knit; re-read after a re-render.
     For IR the corrected covariance equals the Fisher covariance closely for the rates and sigma_noise
     and is up to about 24% smaller for i and N_ch at few channels (IR is conservative there), so
     using the corrected one does not overstate the precision. -->
