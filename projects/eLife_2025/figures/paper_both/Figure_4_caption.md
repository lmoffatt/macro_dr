# Figure 4 caption

Written 2026-07-31 against the merged figure (`figure_4.Rmd`, `Figure_4.pdf`). It replaces the
caption for the superseded R-vs-IR map, whose notebook is archived at
`../archive/figure_4_SUPERSEDED_20260730.Rmd`. Nothing from that version survives here, including
its (A) (B) (C), which now collide with the two letters this figure uses. The numbers it quoted
were printed by chunks that no longer exist; do not carry any of them forward without re-measuring.

**Figure 4. The design plane: how wrong each member's report is, in both moments, on one scale.**

Ten thousand recordings are simulated at every cell of the design grid and scored by each member of
the cost ladder. (**A**) The distortion-induced bias, evaluated at the simulation truth
$\theta_{\mathrm{sim}}$: how far the estimate is displaced. (**B**) The information distortion,
evaluated at the pooled optimum $\theta_{\mathrm{pool}}$, where the score vanishes by construction
so that the comparison is not contaminated by a displaced gradient: how far the uncertainty each
member reports is from the uncertainty it delivers.

**Layout.** Columns nest twice, member and then parameter, the closing rate $k_{\mathrm{off}}$ and
the channel number $N_{\mathrm{ch}}$. Least squares holds one column in each half: its
four-parameter configuration fixes the unitary current and the noise level at the simulated values
(Methods), so a channel-number comparison against it would not be like for like. Rows are the
channel count from $10$ to $10^4$, each given a height proportional to the range of instrumental
noise the sweep reaches at that count, so a decade of noise is the same distance in every row.
Inside a panel the horizontal axis is the sampling interval in units of the closing time constant,
$\Delta/\tau$, and the vertical axis the dimensionless instrumental noise. One decade is the same
physical distance on both axes everywhere in the figure, so the shape of a field can be read as a
shape.

**Colour** is shared by the two halves: the factor by which the member's report is wrong, on a
scale centred on one, so pale is right and saturated is wrong in either direction. What a factor
means differs between the halves, and the difference is not cosmetic. In **A** it is a factor on
the parameter, so $2$ says the estimate is out by a factor of two. In **B** it is a ratio of
variances, so $2$ says the reported error bar is $\sqrt 2$, about $40\%$, too narrow. Blue is the
conservative direction in both. The scale is clipped at $10^{\pm 3}$, which touches $22$ of $3276$
cells, twenty of them the open-loop member, where being wrong by $10^3$ and by $10^{12}$ say the
same thing.

**Lines carry two codes, and neither borrows from the other.** Colour says which quantity the line
cuts: dark hairlines cut the colour field itself, white lines on a dark casing cut the
distortion-corrected standard error, which answers a different question, whether the parameter can
be measured at all rather than whether the error bar is honest. White rather than a second hue
because the field is red and blue, so hue would fail for red-green colour blindness exactly where
the lines matter most, while luminance survives every kind of colour vision and greyscale. Dash
says which criterion: solid is a factor $1.15$ and dashed a factor $2$, the two thresholds used
throughout, so the four lines are two criteria on two quantities rather than four things to learn.
Both criteria are edges of the colour bands by construction, which is why the dark pair needs no
key: the bar names them where they sit.

**Reading it.** The left half is pale over most of the plane and the right half is not, and because
both are on one scale that comparison is legitimate rather than an impression: the estimates land
close to the truth in places where the uncertainty reported with them is badly wrong. Within each
half the members are ordered by cost and the saturation falls with it. The boundary-conditioned
member is the only one pale across most of the plane in both moments, and where it is not is the
few-channel, low-noise corner in which the occupancy Gaussian is misspecified.

## Sources

- `figure_4.Rmd` (layout, the shared factor scale, the two line codes, the clipping count, and the
  isotropy constraint that fixes the height from the width)
- `figure_4_common.R` (`EXCLUDE_ROWS` for the least-squares column, `nch_noise_span` for the row
  heights, `CRIT_SOL` / `CRIT_DSH` for the two criteria)
- `../../NOISE_AXIS_UNITS.md` (the dimensionless noise is one tenth of the swept label; the axis
  applies the conversion on display only, the data are untouched)
- `papers/1_method/docs/manuscript-drafts/sections/06_methods.tex` (the two least-squares
  configurations and the two anchors)

## Open, and it belongs to the figure rather than to the caption

The two halves are put on one scale and asked to be compared, but a factor on the parameter and a
ratio of variances are not the same kind of quantity, which is why the colour paragraph has to
spend three sentences saying so. Plotting the square root of the distortion would put all three
quantities in the figure, including the standard error the white lines cut, into one unit and
retire those sentences. The notebook header carries the cell counts that would change. Until it is
decided, the caption states the difference in words rather than letting the shared scale imply it
away.
