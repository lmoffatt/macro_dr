# Figure 3 caption

**Figure 3. The calibration of the filter over time, from its output through the residual to the accumulated score, along the recursive ladder.**

One common set of 1,000 simulated recordings (a patch of N_ch = 100 two-state, closed ⇌ open, channels, instrumental noise 1e-4, measurement interval 0.1 τ with τ the model's characteristic time) is scored by each rung of the ladder: R (Recursive), MR (Mean Recursive), VR (Variance Recursive) and IR (Interval Recursive). All four are recursive and differ only in what the likelihood conditions the interval-averaged conductance on, so every difference below is attributable to that conditioning alone. Each row is a time-resolved test of one calibration property, and the columns are the four approximations. Unless stated otherwise the quantity is a mean over the 1,000 recordings, the shaded band or vertical bar is its 95% interval, and the green line is the value a calibrated filter would give. The grey band marks the agonist pulse. Time is in milliseconds.

(A) The output. One example recording: the predicted current (blue, mean ± one standard deviation) over the observed current (grey points), on a scale shared by the four columns. The number in each column is the ensemble-mean total log-likelihood ± its standard error over the 1,000 recordings, the quantity the filter maximizes: R −147.2 ± 0.2, MR −153.9 ± 0.2, VR −148.9 ± 0.2, IR −143.4 ± 0.2. The ladder spans 10.4 nats and IR ends highest, so the data are about e^10 times more probable under the fully conditioned filter than under the mean-recursive one.

(B) Residual variance: the mean squared standardized residual, r² with r = (observed current − predicted mean) / predicted standard deviation, which is one for a calibrated filter (the data-level second moment, on a log10 axis). IR holds one throughout; R and MR sit below it for the whole pulse, over-predicting the spread; VR holds one except for a sharp excursion at the concentration jump, where removing the between-endpoint part of the variance costs most because the population is redistributing fastest.

(C) Score bias for k_off: the mean per-interval score s_t = ∂ logL_t / ∂ k_off, which is zero when the score is unbiased (the inference-level first moment).

(D) Per-interval information ratio: Var(s_t) / F_t, with F_t the per-step Fisher information, which is one when the score variance matches the information at each interval.

(E) Accumulated information ratio: Var(Σ s_t) / Σ F_t = J_T / F_T, which is one when the same equality holds over the whole recording (log10 axis). At the end of the recording it reads R 1.07 (95% CI 0.99 to 1.17), MR 1.19 (1.10 to 1.29), VR 1.20 (1.09 to 1.30) and IR 1.08 (1.00 to 1.17). **MR and VR are the pair whose interval falls clear of one**, and VR, which removes part of the predicted variance without taking the boundary term in the gain, is no better than MR.

(F) Score autocorrelation for k_off against lag, zero for an uncorrelated score. At lag one it reads IR −0.004 ± 0.004, indistinguishable from zero, against R 0.191 ± 0.004, VR 0.256 ± 0.004 and MR 0.357 ± 0.004. **Only the fully conditioned filter produces a white score.** The temporal correlation of the score is why the per-interval identity (D, near one) can fail once the score is summed (E): a correlated score makes the variance of the sum exceed the sum of the per-step informations.

Row F is resolved in lag; rows A to E are resolved in time and carry the agonist shading. The score is shown for k_off as a representative parameter; all four parameters appear in the figure supplement. The two non-recursive members of the family, NR and MNR, are computed in the same run and stored with these data, but are measured in no figure of this paper.

**On reading (E) as Figure 2's overconfidence factor:** it carries the same direction and identifies the same over-confident pair, but it is a single-parameter scalar and does not reproduce Figure 2's ellipse-area factors numerically. R in particular sits near one here while its kinetic-pair ratio in Figure 2 is 1.32. State the correspondence qualitatively, not as an equality.

<!-- Source: projects/eLife_2025/figures/paper/figure_3.Rmd, roster ALGOS <- c("R","MR","VR","IR").
     Every number above is printed by the notebook's caption-numbers chunk on each knit; re-read it
     after any re-render rather than copying these values forward.
     Data: figures/data/figure_3_time_dlik_{R,MR,VR,IR}.csv, engine 0ffbda7, regenerated 2026-07-22
     with seed = 20260722 (previously seed = 0, i.e. irreproducible). -->
