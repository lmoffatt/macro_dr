# Figure 7 caption

<!-- 2026-08-26: the short form (Figure_7_caption.tex) was PASTED into the manuscript on this
     date, replacing a caption that described the pre-2026-07-31 drawing (four boundaries, colour =
     parameter, line type = question, ribbons that the `map` chunk does not draw, vertical strips).
     This long form and that short form now agree with the panel; if the drawing changes, change
     this file first. -->

<!-- RENUMBERED 2026-08-10: was Figure 6. Archive directory and archived variant filenames keep
     their old number, which is dated and correct for what they are. -->

**Figure 7. What a macroscopic recording can be asked for, and which method can deliver it: a usage map over channel number and instrumental noise.**

The plane is the design space: the number of channels in the patch, N_ch, against the dimensionless instrumental noise, both on log scales. Every cell of it is filled, and the fill says which parameters survive there, counting a parameter as recovered when its distortion-corrected standard error is inside a factor of two. Five regions partition the plane. The key states each one as a recommendation, because a reader of this figure is choosing a method and a preparation rather than reading a diagnostic, and it lists them in their vertical order on the panel, so that it reads as a section through the map. In the near-white region nothing is measurable. In the amber one the rates alone survive and classical least squares (ILSE, the interval-averaged arm, which is the classical arm this map is built on) is calibrated there, so the cheaper method is the one to use. In the pale green one the rates alone survive but least squares reports an interval that is too narrow, so the likelihood is needed. In the dark green one the unitary conductance i is recovered as well, together with the closing rate k_off and usually with the opening rate and the channel count; this is the region only the recursive likelihood can reach, for a reason given below. The hatched corner is the one in which the likelihood's own interval is not calibrated: the hatching is drawn over the yield colour rather than replacing it, so the corner keeps saying that the parameters are recovered there and that it is the report which has to be corrected. Everything is at one acquisition interval, Δ·k_off = 0.1.

Two classes of line cross the plane and they are drawn differently because they answer different questions. A white line cased in dark marks where a method's reported error bar becomes calibrated, drawn where its information distortion crosses 1.15, a 7% error on the reported standard deviation, and there are two of them, one per method. A thin black line marks where a parameter stops being measurable, drawn where its distortion-corrected standard error crosses a factor of two, and there are four, one per parameter. Each line is a power law fitted through the measured crossings, which are the open dots on it. Two of the four black lines carry the partition, the one for i and the one for k_off; the other two subdivide the dark green region without partitioning it, marking where the opening rate and the channel count individually drop out. The fitted slopes are +1.11 for i, +1.47 for N_ch, +2.09 for k_off and +2.20 for k_on, +1.03 for the calibration of least squares and −0.50 for the calibration of the likelihood. That last one is the only boundary with a negative slope: the region below it closes as instrumental noise rises, because Gaussian noise makes the emission more Gaussian and so repairs the closure the likelihood assumes.

**The partition is measured rather than imposed.** Only six subsets of the four parameters occur anywhere in the 80 cells: all four in 51 cells, {k_off, i, N_ch} in 1, {k_off, i} in 11, {k_on, k_off} in 5, {k_off} alone in 6, and none in 6. Those six nest at exactly the resolution the figure draws, with no counterexample: no cell recovers the amplitude without also recovering a rate, the channel count never survives without the conductance, and no subset carries the opening rate without the closing one. Counting the survivors rather than naming them would not separate the cases that matter here, since the 11 cells holding {k_off, i} and the 5 holding {k_on, k_off} both have two parameters and sit on opposite sides of the question of which method is usable. Classifying every cell twice, once from its own measurements and once from the four fitted lines, 78 of the 80 agree. Both disagreements fall outside the support of the fit that misplaces them, at N_ch = 5 for the boundary of i, which is fitted from 10 upward, and at N_ch = 200 for the calibration boundary of the likelihood, which is fitted only to 100.

**Why the dark green region belongs to one method.** Least squares holds the unitary current fixed and cannot separate the channel count from it, so wherever the fluctuations still carry the amplitude, least squares has nothing to report about it. This is a property of the method and not a finding about its error bar, and it holds over the 63 cells that recover the conductance. Where the two methods do compete, in the 11 cells that recover the rates alone, the interval least squares reports is calibrated in 3 of them and too narrow in the rest. Reading upward at a fixed channel number, then, a recording loses the amplitude first and the rates second, and least squares becomes usable only after the amplitude is already gone.

**Where it is calibrated, least squares is also no less precise on the rates**, which is why the key recommends it rather than merely permitting it. Over the 32 cells in which both methods were run, the ratio of the distortion-corrected standard error, least squares over the likelihood, has a median of 1.04 on the closing rate and 1.01 on the opening rate. The recursion buys resolution on the rates only where channels are few: the median ratio is 1.89 at five channels and 1.11 at ten, and falls to 1.03 at a hundred, 1.01 at a thousand and 0.99 at ten thousand. This is the information budget seen from the design plane, since the rates are carried by the shape of the deterministic transient, which both methods fit. The one region where the recursion does buy a great deal on the rates is the hatched corner, at a median ratio of 1.64 and a maximum of 3.87, and that is exactly the region where its own reported interval has to be corrected before the extra precision can be claimed. In the amber cells the residual cost of the cheaper method is a closing-rate interval about a quarter wider in two of the three, which is far inside the resolution of the map itself.

**What is measured and what is not.** The dots are the 80 simulated cells, which run from 5 to 10,000 channels and from 0.005 to 10⁶ in dimensionless noise. Beyond the last simulated column the fills are veiled with white stripes and the veil is labelled, so the 56% of the width that carries no cell at all is visible as such. The veil begins half a grid step past the last column, at N_ch = 3.16 × 10⁴, because the columns that carry the sweep up to high noise are one per decade, at 10, 100, 1000 and 10⁴, so the column that would have come next is 10⁵. Each boundary is solid only over the range of channel numbers where its own crossings were found and dashed outside it, and those ranges differ: 10 to 10⁴ for the conductance, 10 to 10³ for k_off, 5 to 10⁴ for the opening rate and the channel count, 10 to 10⁴ for the calibration of least squares, and 5 to 100 for the calibration of the likelihood. That last boundary is fitted over one decade and drawn over seven, which is why it thins where it is extended, and it is also the one responsible for one of the two misclassified cells.

**The noise axis is dimensionless by construction.** It is S̃ = σ²/(B τ i²), with σ the RMS current noise measured in a bandwidth B, τ = 1/k_off the relaxation time of the channel and i the unitary current. A whole-cell rig at 1 pA RMS in a 1 kHz band, recording 1 pA channels with a 10 ms relaxation, sits at S̃ = 0.1. No second axis in picoamps is drawn, because the conversion goes as the square of the unitary current: an axis placed at an assumed 1 pA would be two decades wrong for a 10 pA channel, which is more than the resolution of the map itself, and the error would be silent. The three preparation footprints carry the conversion instead, each of them placed from its own reported ranges of current, bandwidth and relaxation time rather than from one value assumed for all three.

**Three real recording configurations are placed on the same plane**, drawn as heavy outlines with reversed tabs so they cannot be taken for either class of boundary: an excised outside-out patch at 10 to 300 channels, whole-cell recording at 20 to 10⁴, and an oocyte under two-electrode voltage clamp at 10⁵ to 5 × 10⁷. Their channel ranges are derived from the smallest and largest currents each configuration can hold rather than assumed, and their noise ranges come from the rig-to-rig spread reported for each. Whole cell ends at 10⁴ channels, which is where the simulations end, so that configuration is covered by measured cells over its whole width. The oocyte box lies entirely inside the veiled half.

**What the map cannot decide.** It is a concept map and not a phase diagram. Its boundaries are level sets of continuous diagnostics, so the criterion moves them: loosening the measurability threshold from a factor of two to a factor of ten lifts the four parameter boundaries by a median factor of 9.1 to 21.2 in noise (9.13 on N_ch, 9.37 on i, 21.19 on k_off, 18.75 on k_on), and loosening the calibration threshold from a distortion of 1.15 to 1.30 lowers the least-squares boundary by a factor of 0.47 to 0.50. The `floor` rows of that table are the geometry chunk's k_off floor and not a boundary this map draws, so the sensitivity of the likelihood's own calibration boundary is not among the numbers the knit prints. The layout survives all of that and the positions do not, so no cell should be read to better than about a decade in noise. The amber region is the thinnest claim in the figure, resting on 3 cells, and it opens only above about 100 channels. All of it is the two-state scheme at a single open probability and a single concentration-jump protocol; under stationarity the mean is flat and the balance between the methods inverts.

<!-- Source: projects/eLife_2025/figures/paper_both/figure_7.Rmd (was figure_6.Rmd), chunk `map` (revision AM, adopted
     2026-07-31; the AF-AL series, which filled cells by HOW MANY parameters survive instead of by
     WHICH, is in tmp/fig6/). The three superseded drawings of the same measurement are in
     figures/archive/figure_6_superseded_20260731/ and the notebook writes them straight there:
     Figure_6_lines.pdf (chunk `render`), which is the version for checking rather than for reading,
     the five-region predecessor Figure_6_regions.pdf, and Figure_6_frontiers.pdf. The sibling
     Figure_6_hatch.pdf (written as Figure_7_hatch since the renumbering) hatches the uncalibrated
     corner over the yield colour instead of replacing it with purple, so that corner keeps saying
     that the amplitude is recoverable there and only the interval is wrong.
     CHOICE CLOSED 2026-08-26 IN FAVOUR OF THE HATCH, and the manuscript's Figure 7 float now
     includes Figure_7_hatch. Three things decided it: the key's own row already says "use MacroIR,
     but not its reported interval", which is the hatch's claim and not the purple one's; Figure 6
     resolves that same corner at about a factor of 1.5, correctable, while the purple fill erases
     what is recoverable there; and the corner is where the recursion buys the MOST precision on the
     rates (median ratio 1.64), so a fill that reads as "nothing to see here" is backwards. The two
     textures do not collide: the veil is white, wide and on the right-hand columns, the hatch black,
     tight and in the left-hand corner, and no cell carries both. The purple decision of 2026-07-31
     is not erased, it is superseded; reverting is one token in the includegraphics plus "violet" in
     the caption.
     Data: battery_pool_G at theta_pool, nsim 10000, figures/data/{1c2ae6f,87889e6,0ffbda7} by search
     path; IR = macro_IR, LSE = nonlinearsqr. N_ch is NOT capped to four columns: the calibration
     boundary uses the 5/20/50/200/500 columns and is invisible without them.
     ALL BOUNDARIES ARE IR-SOURCED. The predecessor took the rate boundaries from least squares,
     which holds the unitary current fixed and cannot measure k_on at all, and so understated what
     the likelihood delivers.
     CALIBRATED, NOT "HONEST", from 2026-07-31: an honest confidence interval is one whose coverage
     holds uniformly over a class of parameters (Li 1989; Baraud; Genovese and Wasserman), which is
     a stronger and different property than the one measured here, and the manuscript's own term is
     already calibration. Still to propagate (2026-07-31): 10 uses of "honest" in 00_abstract.tex, 12
     in 04_results.tex, 3 in 01_introduction.tex, 1 in 05_discussion.tex, and a handful in the
     captions of Figures 3-supplement-1, 4, 5 and 5-budget.
     PROPAGATION CLOSED 2026-08-26, and the count above was already stale when it was read back. In
     live text, comments stripped and counted after the Figure 7 subsection was rewritten: the
     abstract had 1, 04_results.tex 5, the introduction 0, the discussion 2, the appendix 1, against
     40 uses of "calibrated" across the paper. Most of the
     sweep had been done and nobody recorded it. Everything that names the measured property now
     says calibrated; two uses survive on purpose, in the discussion (a sweep-level bootstrap on
     non-stationary fluctuation analysis "returns honest intervals with no model", a resampled
     interval from a route this paper does not measure) and in the derivation appendix ("the honest
     recipe discards or flags the first interval", the colloquial sense).
     Numbers printed by the chunk on every render, including the cell-versus-line confusion table;
     re-read after a re-render. Subset counts, region counts and the criterion sensitivities are
     measured, not fitted. Noise displayed = swept label / 10 (NOISE_AXIS_UNITS.md).
     CROSS-REFERENCE TO FIX: Figure_6_caption.md (was Figure_5_caption.md) cites "the -0.56 slope of the region-0 boundary of
     Figure 6". That was the LSE-mixed k_off floor of the predecessor; this figure's region-0
     boundary is IR-sourced and has slope -0.50.
     Open: whether the lower vertex, where the map closes below N_ch 10, is pinned by running N_ch 2
     and 5 at noise 0.1-10, ten cells. -->
