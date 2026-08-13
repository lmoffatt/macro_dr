# Figure 3 caption

**Figure 3. Where the bias and the distortion are produced, over the whole implemented family.**

One common set of 1,000 simulated recordings (a patch of N_ch = 100 two-state, closed ⇌ open,
channels, instrumental noise 1e-4, acquisition interval 0.1 τ with τ = 1/k_off) is scored by each of
the eight members, with no fitting: the two least-squares arms (LSE on the instantaneous mean, ILSE
on the interval-averaged one), the open-loop gating likelihood (NR), the interval non-recursive
member (INR), the recursive filter on instantaneous samples (R), the two one-endpoint recursive
members (MR, VR) and the boundary-conditioned filter (IR), in cost-ladder order. Each row is a
time-resolved test of one calibration property. Unless stated otherwise the quantity is a mean over
the 1,000 recordings and the band or vertical bar is its 95% nonparametric bootstrap interval,
B = 400; the green line is the value a calibrated member would give and the grey band marks the
agonist pulse. Rows A and B are properties of the data and are drawn neutral; rows C to G are per
parameter, k_off in blue and N_ch in vermillion. Time is in milliseconds; row G is resolved in lag.

**Read across the columns, the figure is a factorial on the two switches of the lattice**, interval
averaging against recursion, and the two switches repair two different failures. Interval averaging
removes the bias of row D: the fraction of informative steps whose interval excludes zero falls from
0.45 to 0.10 on k_off from NR to INR, and from 0.50 to 0.04 from R to IR. Recursion removes the
distortion of row F: 13.63 to 1.07 from NR to R, and 20.84 to 1.08 from INR to IR. Neither switch
does the other's work, so the family cannot be ordered on one scale. Least squares sits outside the
lattice as the classical baseline.

**The accumulated ratio of row F is not a third measurement.** At every point it factorises exactly
into the two the figure already draws, the per-interval ratio of row E weighted by the information
and the factor by which the correlation of row G inflates the variance of a sum,

    Var(Σ s_t) / Σ F_t  =  [ Σ Var(s_t) / Σ F_t ] × [ Var(Σ s_t) / Σ Var(s_t) ] .

Read that way the two recursive members are doing different things. IR carries both factors at one,
1.07 and 1.01 on k_off. R reaches a comparable 1.07 as the product of a per-interval variance it
under-reports by 28% and a correlation that inflates by 48%, and MR of 0.58 against 2.05: the
accumulated calibration of recursion alone is a cancellation of two errors, which is also why
per-interval calibration is not monotone along the cost ladder (Figure 3—figure supplement 1). For
the members that do not condition on the past the failure is almost all correlation, INR reading
1.01 and 11.76 on N_ch.

**The two one-endpoint members are the other half of that reading.** They do not lie between R and
IR: on the accumulated ratio R reads 1.07, MR 1.19 and VR 1.20 before IR returns to 1.08, and on the
ellipse-area measure of Figure 2 the same order is 1.32, 1.97 and 2.18 before IR's 1.02.
Conditioning the interval-averaged conductance on one endpoint costs calibration rather than buying
it, and what recovers it is the second endpoint entering the gain. Which of the two partial
corrections is worse depends on the measurement: VR is further out on the ellipse areas and on the
accumulated ratio for N_ch (1.43 against MR's 1.19), while MR leaves the more correlated score
(lag-one 0.36 against 0.26) and the larger accumulated displacement of row D.

(A) The output. One example recording, the same one in every column: the predicted current (mean,
with the member's own predictive ±1 s.d. band) over the observed current, on a scale shared across
columns. The annotation is the ensemble-mean total log-likelihood over the 1,000 recordings with its
standard error, and read across the row it is the whole cost ladder: −262.30 ± 0.73 for both
least-squares arms, −230.15 ± 0.83 for NR, −215.37 ± 0.64 for INR, −147.17 ± 0.20 for R,
−153.87 ± 0.20 for MR, −148.94 ± 0.22 for VR and −143.45 ± 0.24 for IR, a span of 118.9 nats.

(B) Residual variance: the mean squared standardized residual, r² with r = (observed − predicted
mean)/predicted s.d., one for a calibrated member, on a log axis. The row states its result twice
and in opposite directions for least squares, whose mean over the recording is 1.0000 while 98% of
the individual steps exclude one: the variance it reports is right on average and wrong at almost
every step. The means are 1.1526 for NR, 0.9995 for INR, 0.7816 for R, 0.6904 for MR, 0.8118 for VR
and 0.9989 for IR.

(C) Per-interval Fisher information F_t, as log10 on a linear axis and against no reference: this
row says where in the recording each parameter is measured at all, and for how long. The count of
informative steps carries the row's result on its own: N_ch is informative at all 80 intervals for
NR and INR and at 43 to 46 for the four recursive members, since a member that has conditioned on
the past already carries the number of channels still open and learns nothing further from it once
the agonist is gone.

(D) Score bias, standardized: the mean per-interval score over √F_t, zero when the score is
unbiased, dimensionless so that two parameters with different units share an axis. Reading the row
by summing it is the one thing it does not support: what governs the estimate is the accumulated
score divided by the whole information matrix and not by the information in each direction, since
the parameters are correlated, and the first-order displacement is b = F⁻¹·E[Σ s_t]. In units of
the marginal standard error, and over all four identified parameters rather than the two drawn, the
two interval-averaged members cover zero everywhere and no other member does: R is displaced −1.38
on the unitary current and +0.82 on the channel number, MR −1.96 and +1.14, VR −1.09 and +0.57, NR
−0.95 and +0.92, and IR at most 0.036 on any of the four. Dividing instead by each parameter's own
Fisher diagonal, which is what summing the drawn row amounts to, gives −0.34 for R on the channel
number where the marginal displacement is +0.82, the sign included.

(E) Per-interval information ratio Var(s_t)/F_t, one when the score variance matches the information
at each interval.

(F) Accumulated information ratio Var(Σ s_t)/Σ F_t = J_T/F_T, one when the same equality holds over
the whole recording.

(G) Autocorrelation against lag, three series: the standardized residual in black and the score for
each parameter in its colour. The score's is why (E) can hold where (F) does not, since a correlated
score makes the variance of the sum exceed the sum of the per-interval informations. The residual's
is the same memory at the level of the data, and it is the only quantity in this figure an
experimentalist can compute on a real recording, needing no ensemble and no known truth. The two
agree in order and not in magnitude: at lag one the score reads about 1.6 to 1.7 times the residual
for the members in which either is large (NR 0.864 against 0.510, R 0.191 against 0.118), and the
two coincide only for least squares, 0.856 against 0.846, whose score is proportional to its
residual by construction.

**Magnitude and significance are both reported, because they answer different questions.** Interval
by interval the departures are small (the median of log10(J_t/F_t) over the informative steps lies
within ±0.19 of zero, a factor of 1.6 at the widest) and almost always resolvable (row E excludes
one at 95% of the steps for least squares, 72% for R and every step for VR on N_ch). The fractions
quoted for rows B, D and E are over 1,000 recordings and the steps are correlated, so they measure
what that many recordings resolve rather than how large a departure is, and they carry no p-value.

**Scale conventions.** Any non-negative quantity is on a log10 axis with the calibrated value on the
ticks; signed quantities are on a linear axis with zero on the ticks. Row C is the exception,
log10(F_t) on a linear axis, because it never passes near one. Rows B, D, F and G use a separate
scale per block, the four open-loop columns against the four recursive ones, because the first fail
by a factor of twenty and the second by twenty percent and one axis cannot show both.

**The channel number is not drawn in either least-squares column**, by exclusion and not by absence
of data. Least squares does not model the gating fluctuations, the only thing that separates the
channel count from the unitary current, so with them unmodelled the two enter only through their
product and its information about N_ch would be an artefact of the parameterisation rather than
something it delivers. F_t over the pulse is 1076.55 for both, to every digit. figure_4_common.R
applies the same exclusion through EXCLUDE_ROWS.

<!-- Source: projects/eLife_2025/figures/paper_both/figure_3.Rmd. REWRITTEN 2026-08-12 on the guion
     at papers/1_method/docs/manuscript-drafts/sections/figure_3_guion.md. What changed and why:
     the file had been left at the seven-member roster (no ILSE) since the 2026-08-05 run, and its
     factorial paragraph read the two switches off r_std^2 and J/F, which are a data-level witness
     and the right quantity; it now reads them off rows D and F, which ARE the bias and the
     distortion. Three things are new and none was drawn differently: the accumulated score bias
     (row D summed, in standard errors), the exact factorisation of row F into rows E and G, and
     the fractions-excluding-the-null for rows B, D and E, which the notebook has always printed
     and no caption reported. The MR-against-VR ordering is now stated as measurement-dependent:
     the old text asserted VR worse than MR with a mechanistic story, and row G of this same figure
     reverses that order.
     Promoted 2026-07-31 from what had been Figure 3-figure supplement 1: the body walked four
     members and the supplement seven, and holding two figures apart for two columns stopped making
     sense once every column carried a claim the paper already makes in prose. The four-column
     predecessor, the only notebook that builds this figure straight from the ~1 GB dumps, is
     archived at figures/archive/figure_3_4col_superseded_20260731/ and stays runnable as a check on
     the digests.
     Data: figures/data/digest/figure_3_digest_{LSE,LSE_av0,NR,INR,R,MR,VR,IR}.rds, written by
     figure_3_digest.R from engine 0ffbda7, regenerated 2026-07-22 with seed = 20260722; the INR
     dump is the 2026-07-31 run of macro_INR, and NMR, the build that lost the N*ms interval-variance
     term, is on disk and drawn nowhere.
     RE-KNIT PENDING. Every number above is printed by the notebook's caption-numbers chunk on each
     knit, and the accumulated-bias and factorisation block was added to that chunk on 2026-08-12,
     so the .html on disk is one knit behind. Those numbers were measured from the same digests with
     the notebook's own mask and estimators (tmp/fig3_check.R); the point estimates are
     deterministic and the bootstrap interval edges will move in the last digit under the chunk's
     own seed. Re-read the chunk output after the next knit rather than copying forward. -->
