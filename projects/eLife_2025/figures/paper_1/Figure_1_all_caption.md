# Figure 1 (all algorithms) caption

**Figure 1, extended. One filter step across the whole approximation family.**

This is Figure 1 with the two non-recursive members of the family restored, and it is the version to
read when the question is how the recursion axis compares with the conditioning axis. It is not the
paper's Figure 1: paper 1 holds recursion fixed and measures only the recursive ladder, so the
non-recursive columns appear here for orientation and carry no measured claim.

The current from a patch of N_ch = 20 identical two-state (closed ⇌ open) channels is recorded as a time-average over successive measurement intervals, with added instrumental noise. Each column is one Gaussian (moment-closure) approximation to the likelihood of that recording: NR (Non-Recursive), MNR (Mean Non-Recursive), R (Recursive), MR (Mean Recursive), VR (Variance Recursive), and IR (Interval Recursive). They differ along two axes and one further distinction. Recursion: whether the estimate from one interval is carried into the next as its prior (R, MR, VR, IR) or each interval is predicted from the same fixed start (NR, MNR). Averaging: whether the predicted current uses the channel conductance at a single instant (NR, R) or its time-average over the interval (MNR, MR, VR, IR). Endpoints: whether the quantity carried across the update is a single state (NR, MNR, R, MR, VR) or the joint distribution of the state at the two interval endpoints (IR). MR and VR sit at the same position on all three: they differ only in the predictive variance, MR taking the total per-start-state conductance variance and VR the residual part that remains after conditioning on where the interval ends. Two consecutive measurement windows are shown (6 to 10 ms) so that the per-interval operations are visible; the first ten milliseconds of the twelve-millisecond recording are given as a companion panel.

Rows follow one turn of the filter cycle. (A) Prior: the open probability P_open predicted by the Markov step, before the current of the interval is measured. (B) Observation: the predicted current and its spread (orange), the measured current, and the innovation (blue), which is the measured minus the predicted current. The spread is drawn as an error bar for the instantaneous approximations (NR, R) and as a band for the averaging ones (MNR, MR, VR, IR), so the averaging axis reads directly off the observation row. (C) Posterior: P_open after the Bayes update. NR and MNR run open-loop and have no update step, marked "no update". (D) The per-step log-likelihood accumulated over the window, Σ logL, reset to zero at the window start so the six columns begin together.

Colour encodes the filter phase and is shared across the rows. Green is the channel current, the resolved simulated current, drawn in A and C as the realised open fraction (channel current divided by N_ch), the ground truth the estimate tracks. Orange is the predict phase: the Markov step, the prior, and the predicted current. Blue is the observation and the innovation. Purple is the Bayes update and the posterior. Black is the log-likelihood. Arrows mark operations. A state carried in from the neighbouring step is drawn hollow (points) or dashed (segments) to set it apart from the state computed within the panel.

<!-- Source: projects/eLife_2025/figures/paper/figure_1_all.Rmd (roster) + figure_1_panels.R (panels).
     Writes Figure_1_all.pdf and Figure_S1_all.pdf, never the paper's Figure_1.pdf.
     The data key for the mean-non-recursive algorithm is NMR and it prints as MNR. -->
