# Figure (number pending) caption — the log-likelihood plane

Written 2026-08-12 against `figure_4_logL.Rmd` / `Figure_4_logL.pdf`, the merged two-side version.
It supersedes nothing: this figure did not exist before. Every number below is printed by the
notebook's own chunks (`logl-normalisation`, `sub-nat-signal`, `what-the-likelihood-ranks`) and was
read off the 2026-08-12 render, not carried from a draft.

The slot is undecided. The case for the body, ahead of Figure 4, is that a reader who has just been
shown which approximation is the better density then meets a figure saying that ranking does not
predict the inference. The case for a supplement is page budget. Whoever decides should note the
figure is 7.4 x 6.6 in and leaves about 3 in of a page for the caption, which is what merging the
two halves bought: nine columns alone came to 11.1 in and seven alone to 14.0, against about 9.6 in
of usable text height.

---

**Figure N. Δ log-likelihood over the design plane: which approximation is closer to the truth, and
what each step of the cost ladder buys.**

Ten thousand recordings are simulated at every cell of the design grid and scored by every member of
the ladder. Each panel is a difference of log-likelihoods between two members, summed over the
recording and averaged over the ten thousand, with both members evaluated at the simulation truth
θ_sim. Evaluating them at one parameter is what makes the difference interpretable: it is a
difference of Kullback-Leibler divergence rates between the two approximations of the same density,
so a positive value says the member named on the upper line of the column is the closer density to
the process that generated the data, by that many nats. No estimator enters, which is what separates
this page from Figure 4: there the question is what a member's report gets wrong, here it is what
its density gets wrong, and the second is the property a likelihood ratio, an information criterion
or a Bayes factor would rank on.

(**A**) the edges of the family's lattice, one step at a time, grouped by which axis of the family
the step moves: modelling the gating variance, recursion, averaging the model over the acquisition
window, and conditioning on one endpoint of the interval instead of two. Every difference is taken
in the direction more machinery minus less, so a positive value always means the extra step bought
density and no panel is read backwards. (**B**) every member against `IR`, grouped by family and
cheapest first inside each group, which answers the other question the ladder raises: not what a
step buys, but how far behind each member ends up.

**Layout.** Columns carry the subtraction on two lines, minuend above and subtrahend below. Rows are
the channel count from 10 to 10⁴, each given a height proportional to the range of instrumental noise
the sweep reaches there, so a decade of noise is the same distance in every row. Inside a panel the
horizontal axis is the acquisition interval in units of the closing time constant and the vertical
axis the dimensionless instrumental noise, and one decade is the same physical distance on both axes
everywhere in the figure, so the shape of a field can be read as a shape.

**Units.** Nats per recording, and the recording is 10 τ in every cell of this grid. The summed
log-likelihood is therefore already a quantity per unit time and is not divided by the number of
intervals: doing so would multiply by Δ̃/10 and impose a hundredfold ramp along the horizontal axis
that is nothing but the sample count. Measured both ways, the raw form is the flatter one along that
axis in eight of the nine contrasts of (A).

**Colour** is shared by the two halves, which is what allows them to be compared rather than merely
placed side by side, and is signed about zero. The two thresholds are band edges: the solid contour
is 3 nats, where the recording tells the two densities apart at all, and the dashed one 100 nats,
where it does so decisively. Grey is the terminal band, reached by 44 of 3360 cells, 35 of them the
ill-conditioned corner Figure 4 also greys.

**What white means, and it is measured rather than assumed.** Each member is scored on its own
ensemble of ten thousand recordings, so a difference between two of them carries a standard error,
and a cell whose difference falls inside twice that error is drawn at exactly zero. That floor is not
uniform inside a panel: the spread of a summed log-likelihood grows as the square root of the
interval count, so twice the standard error runs 0.64 nats at Δ̃ = 0.01 and 0.064 at Δ̃ = 1. The
scale reaches down to 0.1 nats for that reason and no further: no cell anywhere on the grid can
resolve 0.03. A sub-nat band therefore appears only where the local floor sits below it, which is why
the faint structure fades toward the fine end of every panel. The collapse moves 1349 of 3360 cells
and changes the colour band of ten of them, since the innermost band is already ±1 nat.

**The readings are in the text.** In brief: the recursion is worth 1.93 nats at the median of the
plane and modelling the gating variance −0.46, while averaging the model over the window is worth
0.04 to 0.27 and the two partial rungs are worth nothing, `MR` giving up more than 3 nats against `R`
on 18 per cent of the plane and as much as 83. And in (B) there is not one cell of 1470 where a
cheaper member is more than 3 nats ahead of `IR`; the most any of them manages is 1.13.

---

## Sources

- `figure_4_logL.Rmd` (layout, the two sides on one scale, the CI collapse, the calibration of
  `per_row`, and the numbers quoted above)
- `figure_4_logL_data.R` (the extraction, the θ_sim anchor, and the probit/statistic trap)
- `../figure_4_source_data/figure_4_source_data_logL.csv` (2086 rows, with its provenance stamp)
- `figure_4_common.R`, `figure_4_layout.R` (geometry, unmodified: this notebook shadows `ggsave`
  locally rather than editing them)
- `../../NOISE_AXIS_UNITS.md` (the dimensionless noise is one tenth of the swept label; the axis
  applies the conversion on display only)

## Open

The number and the slot. Also, whether the two redundant columns of (A) are kept: the two
four-member squares are drawn with all their edges, so (ILSE−LSE) + (INR−ILSE) = (NR−LSE) + (INR−NR)
holds by construction. They are a visible closure check and they cost two columns, and columns are
what keep this figure short, so the trade is not obviously worth reversing.
