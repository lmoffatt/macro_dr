## figure_4_layout.R — the Figure 4 plane, as machinery, so that the figure and its supplements
## cannot drift apart.
##
## Extracted verbatim from `figure_4.Rmd` on 2026-07-31, when the second and third notebook on the
## same layout appeared. What was `moment` there is a SIDE here: each of the two halves is a spec,
##   list(lab = <the label on the bar>, letter = <"A"/"B">, data = <a data frame>, zf = <function>)
## and everything else -- the nesting, the rows pinned to the union, the shared factor scale, the
## two line families, the gaps, the colour bar and the key -- is identical by construction.
##
## A notebook using this declares, in order:
##   FIG4_ALGOS / FIG4_ALGO_LAB / FIG4_DATA_DIRS   (optional, before sourcing figure_4_common.R)
##   source("figure_4_common.R"); source("figure_4_layout.R")
##   SIDES <- list(<side>, <side>)
##   fig4_render(SIDES, per_row = <nominal in per decade>, out = "<file>.pdf")
##
## THE COLOUR SCALE IS SHARED BY CONSTRUCTION, not by convention. SH_FAC / SH_L / SH_PAL / SH_LAB
## and the two criteria are built HERE, once, out of DIST_SOL and DIST_DSH in figure_4_common.R,
## and fig4_render() draws the bar from the same objects the blocks fill from. So every figure on
## this layout carries the identical bands, the identical 1.15 and 2, and the identical clip at a
## factor 10^3, and a colour means the same thing across the whole set. Changing the scale for one
## figure is not possible without changing it for all of them, which is the point.
##
## PER_ROW is NOT shared: it is solved per figure, because the null-unit split depends on how many
## columns the composition carries. See the note above fig4_render().
##
## Sourced by: figure_4.Rmd, figure_4_supplement_sample_corr.Rmd,
##             figure_4_supplement_recursive_ladder.Rmd

source("noise_units.R")  # the one definition of S_tilde and of when to convert

suppressMessages({ library(tidyverse); library(patchwork) })

# ---- typesetting and furniture ------------------------------------------------------------
# Typeset numbers: 10^4, never 1e+04. Returns strings to be parsed.
fmtnum <- function(v) vapply(v, function(x) {
  if (!is.finite(x)) return("")
  if (x == 0) return("0")
  sg <- if (x < 0) "-" else ""; a <- abs(x)
  e  <- round(log10(a)); ispow <- isTRUE(all.equal(a, 10^e, tolerance = 0.02))
  if (a >= 1000 || a <= 0.001) {
    if (ispow) sprintf("%s10^%d", sg, e)
    else sprintf("%s%g %%*%% 10^%d", sg, signif(a / 10^e, 2), e)
  } else if (a > 0.7 && a < 1.45) {
    # NEAR 1 the label must round-trip, not round: a fixed three significant figures printed the
    # standard-error scale's 1.005 as "1" and its 1.015 as "1.02", which collided with the real
    # 1.02 edge two bands along. Use the FEWEST decimals that reproduce the number exactly, so
    # 1.005 prints 1.005, 1.015 prints 1.015, 1.02 prints 1.02 and 0.971 prints 0.971.
    # capped at three decimals with a fallback: the diverging scale's blue arm is 1/1.35, 1/1.15
    # and friends, which never round-trip, and an uncapped rule printed 0.970874 where the figure
    # has always said 0.971.
    dec <- 0
    rt  <- function(k) isTRUE(all.equal(round(a, k), a, tolerance = 1e-12))
    while (dec < 3 && !rt(dec)) dec <- dec + 1
    paste0(sg, if (rt(dec)) sprintf("%.*f", dec, a)
               else format(signif(a, 3), trim = TRUE, scientific = FALSE))
  } else {
    paste0(sg, format(signif(a, 2), trim = TRUE, scientific = FALSE))
  }
}, "")

# Colour bar with the tick labels alternating below / above, which doubles the room each one has
# and is what lets the crowded tail hold 7 pt type. Its title is NOT an axis title: see cbtitle.
cbar <- function(pal, edge_lab, lin, sol = numeric(0), dsh = numeric(0)) {
  n <- length(pal)
  # Where a value sits along the bar. This USED to be `match(v, lin) - 0.5`, which demands that
  # every mark be an exact band edge and drops it silently otherwise: figure_4_common.R warns
  # about it in as many words, and it still bit, because this figure feeds the marks declared
  # against DIST_LIN to the merged scale SH_L, where log10(0.87) is not an edge. Both criteria
  # vanished without a message. Interpolating puts a mark wherever its value falls.
  pos <- function(v) stats::approx(lin, seq(0.5, n + 0.5, 1), v, rule = 2)$y
  mkx <- c(
    if (length(dsh)) list(annotate("segment", x = pos(dsh), xend = pos(dsh), y = .5, yend = 1.5,
                                   colour = "grey25", linewidth = .35, linetype = "dashed")),
    if (length(sol)) list(annotate("segment", x = pos(sol), xend = pos(sol), y = .5, yend = 1.5,
                                   colour = "black", linewidth = .5)))
  ed <- data.frame(pos = seq(0.5, n + 0.5, 1), lab = as.character(edge_lab),
                   stringsAsFactors = FALSE)
  ed$below <- seq_len(nrow(ed)) %% 2 == 1
  PT <- 7 * 0.3527
  ggplot(data.frame(x = seq_len(n)), aes(x, 1, fill = factor(x))) +
    geom_tile(colour = "grey92", linewidth = .15) + mkx +
    scale_fill_manual(values = pal, guide = "none") +
    annotate("segment", x = ed$pos[ed$below], xend = ed$pos[ed$below], y = .5, yend = .32,
             colour = "grey55", linewidth = .3) +
    annotate("segment", x = ed$pos[!ed$below], xend = ed$pos[!ed$below], y = 1.5, yend = 1.68,
             colour = "grey55", linewidth = .3) +
    annotate("text", x = ed$pos[ed$below], y = 0.10, label = ed$lab[ed$below],
             size = PT, vjust = 1, family = "Helvetica", parse = TRUE) +
    annotate("text", x = ed$pos[!ed$below], y = 1.90, label = ed$lab[!ed$below],
             size = PT, vjust = 0, family = "Helvetica", parse = TRUE) +
    scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(-0.35, 2.35), clip = "off") +
    theme_void(base_family = "Helvetica") +
    theme(plot.margin = margin(2, 4, 1, 0))
}

# The colour bar's title, as its own plot rather than an axis title. An axis title is a FIXED
# width, and patchwork hands every row in a composition the widest fixed left offset in it: this
# one label, plus the 46 pt margin that used to hold the bar clear of it, was reserving 1.8 of
# the 7.4 inches on every row of maps, for a legend three rows down. Priced in the same column
# units as everything else it costs what it is worth. The right margin is what keeps the leftmost
# tick of the bar, which overhangs its own left edge by half a label, off the last word.
cbtitle <- function(s) ggplot() + theme_void(base_family = "Helvetica") +
  annotate("text", x = 1, y = .5, label = s, size = 7 * 0.3527, hjust = 1, lineheight = .9) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme(plot.margin = margin(0, 16, 0, 0))

# Key for the lines drawn on the maps. TWO families, and they are not interchangeable: the
# iso-distortion contours are plain dark hairlines on the fill, the standard-error contours are
# white cored on a dark casing so they survive the saturated end of the scale. Each sample takes
# the linewidths and the dash pattern blk() actually draws, and sits on a chip of the tone its
# family is meant to be read against: the white core would vanish into the page otherwise, and
# it is the casing, not the core, that a reader sees where the fill is pale.
linekey <- function(labs, cased, chip, dash, x2 = .32, sz = 7) {
  ent <- function(x0, lab, dashed) {
    lt <- if (dashed) dash else "solid"
    c(list(annotate("rect", xmin = x0, xmax = x0 + .05, ymin = .12, ymax = .88,
                    fill = chip, colour = NA)),
      if (cased)
        list(annotate("segment", x = x0, xend = x0 + .05, y = .5, yend = .5, colour = "grey10",
                      linewidth = if (dashed) .70 else .85, linetype = lt),
             annotate("segment", x = x0, xend = x0 + .05, y = .5, yend = .5, colour = "white",
                      linewidth = if (dashed) .35 else .45, linetype = lt))
      else
        list(annotate("segment", x = x0, xend = x0 + .05, y = .5, yend = .5,
                      colour = if (dashed) "grey25" else "black",
                      linewidth = if (dashed) .22 else .30, linetype = lt)),
      list(annotate("text", x = x0 + .062, y = .5, label = lab, size = sz * 0.3527,
                    hjust = 0, family = "Helvetica")))
  }
  ggplot() + ent(0, labs[1], FALSE) + ent(x2, labs[2], TRUE) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void(base_family = "Helvetica") + theme(plot.margin = margin(4, 0, 0, 0))
}

# `gap` is the space below the box, and it is the ONLY thing setting the distance to the row
# under it: the rows sit flush otherwise (see hsp).
boxed <- function(lab, sz, fill, panel = NULL, gap = 0, parse = FALSE) ggplot() +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = fill, colour = "grey20",
           linewidth = .28) +
  annotate("text", x = .5, y = .5, label = lab, size = sz * 0.3527, fontface = "bold",
           family = "Helvetica", parse = parse) +
  (if (is.null(panel)) NULL else
     annotate("text", x = .012, y = .5, label = panel, size = (sz + 1.5) * 0.3527, hjust = 0,
              fontface = "bold", family = "Helvetica")) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  theme_void() + theme(plot.margin = margin(0, 0, gap, 0))

# THE X AXIS, DRAWN ONCE. A ggplot scale belongs to every facet column of its block, so tick
# labels asked for on a block with five columns print five times and, with panel.spacing at zero,
# the "1" of one column lands exactly on the "0.01" of the next. That never showed while the
# labelled block was LSE, which has a single column; it appeared the moment the parameter became
# the outer level. Since all 22 columns carry the identical range, the honest fix is one axis for
# the figure: the blocks draw ticks and no numbers, and this strip goes in the x-title row, one
# column wide, aligned under the first column of the first block.
xaxis_strip <- function(sz = 7) ggplot() +
  annotate("segment", x = 0, xend = 1, y = 1, yend = 1, colour = "grey20", linewidth = .3) +
  annotate("segment", x = c(0, .5, 1), xend = c(0, .5, 1), y = 1, yend = .78,
           colour = "grey20", linewidth = .3) +
  annotate("text", x = c(0, 1), y = .62, label = c("0.01", "1"), size = sz * 0.3527,
           vjust = 1, family = "Helvetica") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1.05), clip = "off") +
  theme_void(base_family = "Helvetica") + theme(plot.margin = margin(0, 0, 0, 0))

txt <- function(s, sz = 8, top = 0, bot = 0) ggplot() + theme_void(base_family = "Helvetica") +
  annotate("text", x = .5, y = .5, label = s, size = sz * 0.3527) +
  theme(plot.margin = margin(top, 0, bot, 0)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off")

# A spacer whose width is written down, and whose HEIGHT costs nothing. patchwork's own
# plot_spacer() keeps the default 5.5 pt plot margin on every side, so a spacer asked for zero
# width still costs 11 pt across; six of them between the members were spending 0.9 in of the 7.4
# on white, and no knob in this file reached it. The vertical half of the same margin was worse,
# because a row's margin is the widest one among its plots: the single spacer holding the two
# halves apart was therefore putting 11 pt of white under EVERY row it appeared in, which is what
# stood the moment labels off the member labels and those off the parameter strips. The rows are
# flush now and the space between them is whatever `boxed(gap =)` asks for.
hsp <- function(pt) ggplot() + theme_void() +
  theme(plot.margin = margin(0, pt / 2, 0, pt / 2))

# ---- the shared scale, the criteria, and the standard error ---------------------------------
# Distortion-corrected standard error: same file and same anchor as the distortion that colours
# the map, so the white contours and the colour field are commensurable.
sedat <- dplyr::filter(.fig4_src$sedat, algo %in% ALGOS)

# SHARED LIMITS, one range per N_ch row for EVERY block, taken as the union over the members.
# Without this the axis is quantitatively wrong for one of them. `facet_grid(space = "free_y")`
# sizes each row from the range of the block it is in, and patchwork then gives every block the
# same total height, so a member whose sweep reaches a noise the others never ran gets its own
# noise scale. IR runs down to 0.05 where LSE, NR and R start at 0.1: 27.204 decades against
# 26.000, a 4.6 per cent difference, and only ONE y axis is drawn for the eight blocks. Measured
# on the PDF before the fix, IR's row boundaries sat 2 px below everyone else's, and the tick
# printed 0.05 landed on 0.0275 of IR's data, a factor 1.8 out, growing to a factor 2 at the foot
# of each row. Pinned to the union, a member with no cell at 0.05 shows empty there, which is
# what "not run" should look like, instead of stretching its rows to hide it.
#
# It also makes the isotropy calibration true for the whole figure rather than for IR alone:
# every block now carries the 27.204 decades that `nch_noise_span` already assumed, so PER_ROW
# needs no re-solving.
XLIM <- range(c(bias$lx, stat$lx))
YLIM <- dplyr::bind_rows(dplyr::select(bias, Num_ch, ly),
                         dplyr::select(dplyr::filter(stat, comp == "total"), Num_ch, ly)) %>%
  dplyr::filter(Num_ch %in% NCHS) %>%
  dplyr::group_by(Num_ch) %>%
  dplyr::summarise(lo = min(ly), hi = max(ly), .groups = "drop")

# Y_BRK's lowest break IS the shared floor, so once the limits are pinned it prints on the very
# edge of every row, hard against the top label of the row below. Drop it: a tick on a flush panel
# edge names a value the reader can already see and costs a collision. It also costs width, since
# "0.005" is the widest label on the axis and the axis is a fixed left offset for every row.
YB <- Y_BRK[Y_BRK > min(YLIM$lo)]

# The scale is declared as FACTORS and the two validity criteria are MEMBERS of it, so an
# iso-line falls exactly on a colour boundary instead of cutting a band in half. That is also
# why the criteria need no key of their own: the bar labels them, 1.15 and 2, and the line sits
# on the label. The ladder is geometric in log10, each step about 2.2 times the last, so the fine
# structure near 1 that resolves IR survives and the coarse tail to 10^3 is what it was.
#
# The magnitudes come from figure_4_common.R so the bar, the maps and the supplements cannot
# drift apart, but the pair is MIRRORED here rather than taken as declared: DIST_SOL is
# c(0.87, 1.15) and 0.87 is a rounded 1/1.15, off by 0.0002 in log10. That rounding is free on a
# scale whose edges it never has to hit, and not free here, where the low criterion IS an edge.
CRIT_SOL <- max(DIST_SOL)   # 1.15, strict, ALWAYS solid
CRIT_DSH <- max(DIST_DSH)   # 2,    loose,  ALWAYS dashed
# 1.01 ADDED 2026-08-06, and it is the innermost cut the scale can honestly carry. The white band
# used to run 0.971 to 1.03, i.e. everything within three per cent of the null was one colour.
# Counted on the six-member roster of Figure 4, that band held 975 of 2149 cells, of which 843 sit
# at EXACTLY 1 because their bootstrap interval brackets it and Dconf collapses them there. The
# split at one per cent divides the remaining 132 almost in half, 67 inside and 65 between one and
# three per cent, and 96 of those 132 are R and IR: it is the pair the fine end of the scale exists
# to separate, and until now they shared a colour. So the new band is neither empty nor a noise
# floor dressed as structure — the CI collapse is what guarantees the second, since a cell whose
# interval reaches the null cannot land in it at all.
SH_FAC <- c(1.01, 1.03, 1.07, CRIT_SOL, 1.35, CRIT_DSH, 5, 10, 100, 1000)
SH_L   <- c(rev(-log10(SH_FAC)), log10(SH_FAC))
SH_PAL <- bandcols(SH_L)
SH_LAB <- fmtnum(10^SH_L)
# cutting the fill, so mirrored: wrong by 1.15 either way
ISO_SOL <- c(-1, 1) * log10(CRIT_SOL)
ISO_DSH <- c(-1, 1) * log10(CRIT_DSH)

# The standard error takes THE SAME two criteria and the same style code, so the two readings are
# comparable number for number: "the report is wrong by more than 1.15" against "the error bar is
# wider than 1.15". Until 2026-07-30 this pair was 2 solid and 1.1 dashed, which inverted the code
# against the fill's lines and put a number, 1.1, that matched nothing else in the figure.
# CONSEQUENCE, deliberate: the factor-2 line is the frontier Figure 6 draws, and here it is now
# the dashed one. The two figures still agree on the value; only the weight differs.
SE_SOL <- log10(CRIT_SOL)
SE_DSH <- log10(CRIT_DSH)

# ---- the scales a side can carry ---------------------------------------------------------------
# A side draws with ONE of these. The default is the diverging factor scale above, and every plane
# figure built before 2026-08-01 uses it, which is what keeps a colour meaning the same thing
# across the set. A quantity that is one-SIDED (a standard error, an autocorrelation time, a
# condition number: all of them >= their null and none of them meaningful below it) wastes half of
# a diverging scale and gets its own sequential one instead. The two criteria stay wherever they
# still mean something.
#
# Sequential and not diverging also means the luminance is monotone, so these read in greyscale and
# under any colour vision, which a two-armed scale only does on each arm separately.
SEQ_RAMP <- c("#FFFFF0", "#FDE725", "#7AD151", "#22A884", "#2A788E", "#414487", "#440154")
seqcols  <- function(breaks) colorRampPalette(SEQ_RAMP, space = "Lab")(length(breaks) - 1)

SCALE_FACTOR <- list(brk = SH_L, pal = SH_PAL, lab = SH_LAB, sol = ISO_SOL, dsh = ISO_DSH,
                     title = "factor by which\nthe report is wrong")

# Standard error and integrated autocorrelation share a scale, because they are the same KIND of
# number: a factor on the parameter and a factor on the effective sample size, both with a null of
# 1 and both with 1.15 and 2 as the criteria the whole set uses. Fine from 1 to 2, coarse after,
# clipped at 1000 — the grid's upper quartile runs away (measured: q75 of the corrected standard
# error is 1.9 for LSE and 200 for NR) and a tail that chased it would flatten everything below.
# Two one-sided scales, and they differ ONLY in how finely they cut the low end.
#
# SCALE_ONESIDED is the general one: fine from 1.02 to 2, coarse after, clipped at 1000. It is what
# a quantity whose interesting range starts around a few per cent takes, e.g. the integrated
# autocorrelation of the residual, which runs 1.00 to 9.3.
SEQ_FAC <- c(1, 1.02, 1.05, 1.10, CRIT_SOL, 1.35, CRIT_DSH, 3, 5, 10, 100, 1000)
SEQ_L   <- log10(SEQ_FAC)
SCALE_ONESIDED <- list(brk = SEQ_L, pal = seqcols(SEQ_L), lab = fmtnum(SEQ_FAC),
                       sol = log10(CRIT_SOL), dsh = log10(CRIT_DSH), title = "")

# SCALE_SE is for the standard error and for nothing else, because the standard error is the one
# quantity here whose useful range reaches down to a few PER MIL: the closing rate's error at ten
# thousand channels is 0.0016 to 0.0071 in log10 units, a factor 1.004 to 1.016, so on the scale
# above that whole column came out one flat band. Three extra cuts below 1.02 open it without
# touching anything from 1.02 up, so the two scales agree wherever they overlap and the criteria
# are the same members, 1.15 and 2.
# The cuts below 1.02 were first placed at 1.002 and 1.005, and measuring the source data showed
# that was the wrong place: only 11 of 5725 cells fall below 1.01 at all. What IS populated is the
# band those cuts were meant to open, (1.01, 1.02], which holds 217 of the 247 k_off cells at ten
# thousand channels. The cut that splits it is 1.015 — below it 102 cells, below 1.02 231 — and it
# is the one that separates the members there, since k_off at 1e4 runs 1.003 to 1.019 across the
# seven. So: 1.005 kept as a floor marker, 1.01 and 1.015 doing the work.
SE_FAC <- c(1, 1.005, 1.01, 1.015, 1.02, 1.05, 1.10, CRIT_SOL, 1.35, CRIT_DSH,
            3, 5, 10, 100, 1000)
SE_L   <- log10(SE_FAC)
SCALE_SE <- list(brk = SE_L, pal = seqcols(SE_L), lab = fmtnum(SE_FAC),
                 sol = log10(CRIT_SOL), dsh = log10(CRIT_DSH),
                 title = "standard error, as a\nfactor on the parameter")

# The condition number needs its own, and it is the reason it could not share: measured over this
# grid it runs from about 400 at ten channels to 1e5 and past it at ten thousand, so on a scale
# whose fine structure sits between 1.02 and 1.35 every cell would saturate. Decade bands, and the
# one mark it has is 3e4, which is the ill-conditioning threshold this project already used when
# the maps still greyed cells out. No dashed criterion: there is only one threshold, and inventing
# a second to fill the style code would be inventing a number.
# Floor at 100 and ceiling at 1e8 because that is where the data is: measured over the grid at
# noise 0.1 the smallest median is 359 (VR at a hundred channels) and the largest 4.9e7 (LSE at ten
# thousand). Bands below the floor would be empty and a ceiling at 1e7 would clip LSE alone.
KAP_FAC <- c(100, 300, 1000, 3000, 1e4, 3e4, 1e5, 1e6, 1e7, 1e8)
KAP_L   <- log10(KAP_FAC)
SCALE_KAPPA <- list(brk = KAP_L, pal = seqcols(KAP_L), lab = fmtnum(KAP_FAC),
                    sol = log10(3e4), dsh = numeric(0),
                    title = "condition number of\nthe Gaussian Fisher")

# The INNER column level. Two parameters by default, which is what Figure 4 and the plane
# supplements carry; a notebook declares FIG4_PARS / FIG4_PSHORT before sourcing to put something
# else there. Two uses so far: the remaining parameters (k_on, i, sigma_noise), and the residual
# diagnostics, which have no parameter index at all and put the QUANTITY on this level instead.
# The names are the values of the data frame's `param` column; the labels are parsed as plotmath.
PARS   <- if (exists("FIG4_PARS")) FIG4_PARS else c("off", "Num_ch_mean")
PSHORT <- if (exists("FIG4_PSHORT")) FIG4_PSHORT else c(off = "k[off]", Num_ch_mean = "N[ch]")

# the row variable is named once, on the top row, and carries only the number below it
nchstrip <- function(nchs) {
  lab <- vapply(sort(nchs), function(n) {
    e <- round(log10(n)); if (e == 1) "10" else sprintf("10^%d", e) }, "")
  lab[1] <- sprintf('atop("N"[ch], %s)', lab[1]); lab
}

# ---- one block per (side, member) -----------------------------------------------------------
# One block per (moment, member): its own ggplot, faceted parameter x N_ch. The member's own box
# is drawn separately and composed above it, because facet_grid with two column variables does
# not span and patchwork gives exact control.
# `side` is one half of the figure, declared by the notebook: which data frame it draws, how the
# fill is read out of it, and what the bar over it says. Everything a half CAN differ in is in
# there, which is what stops the three notebooks from re-implementing the block.
# `outer` names which level owns the BLOCK, and the other one becomes the facet columns inside it.
# "member" is Figure 4's nesting, member outermost and parameter inside. "param" transposes it, so
# the algorithms sit side by side within one parameter, which is the arrangement to use when the
# question is how the members COMPARE on the same quantity rather than what one member reports.
# The column count is identical either way; only the adjacency changes.
blk <- function(a, side, nchs, ylab, ystrip, xlab, se_lines = TRUE, outer = "member") {
  sc <- if (is.null(side$scale)) SCALE_FACTOR else side$scale
  if (outer == "member") {
    keep <- setdiff(PARS, EXCLUDE_ROWS$param[EXCLUDE_ROWS$algo == a])
    d    <- side$data %>% dplyr::filter(algo == a, param %in% keep, Num_ch %in% nchs)
    have <- PARS[PARS %in% unique(d$param)]
    inner_lab <- PSHORT
    inner_of  <- function(x) x$param
  } else {
    keep <- setdiff(ALGOS, EXCLUDE_ROWS$algo[EXCLUDE_ROWS$param == a])
    d    <- side$data %>% dplyr::filter(param == a, algo %in% keep, Num_ch %in% nchs)
    have <- keep[keep %in% unique(d$algo)]
    # quoted, so label_parsed prints them upright like the member boxes and not as italic symbols
    inner_lab <- setNames(sprintf("'%s'", ALGO_LAB[ALGOS]), ALGOS)
    inner_of  <- function(x) x$algo
  }
  if (!length(have)) return(NULL)
  zf <- side$zf
  nlev <- sort(unique(d$Num_ch)); strip <- nchstrip(nlev)
  mk <- function(x) x %>% dplyr::mutate(
    nchf = factor(vapply(Num_ch, function(n) strip[match(n, nlev)], ""), levels = strip),
    pf   = factor(inner_lab[inner_of(.)], levels = inner_lab[have]))
  # OFF THE CHART, GREY, and grey as a TERMINAL BAND of the scale rather than as a layer of tiles
  # stamped over the field. A cell past the top edge of a one-sided scale is not "a very large
  # value", it is a quantity that ran away: for the standard error it is a parameter the design
  # does not determine, where the sandwich is unbounded and the number reaches 10^19115. Painting
  # it the darkest band says "very big", which is a claim; grey says "not on this scale".
  #
  # Doing it as one more band is what makes the grey CONTINUOUS: the same geom_contour_filled draws
  # it, so its boundary is exactly the contour where the value crosses the top edge, and it joins
  # up the way every other band does. Tiles gave a lattice of little squares with their own edges,
  # which read as data. It also puts the grey in the COLOUR BAR, where it explains itself.
  # Opt-in per side, so no figure built before 2026-08-01 changes.
  if (isTRUE(side$grey_offscale)) {
    sc$brk <- c(sc$brk, 30)          # 10^30: past anything, and finite so contouring works
    # grey40 and not a light grey: the panel background is grey85 and MEANS "not run", so a pale
    # terminal band read as the same absence. This one has to read as the END OF THE RAMP, past the
    # darkest colour, which is a result and not a gap. The two greys sit side by side wherever a
    # member with a short sweep neighbours a parameter that runs off scale.
    sc$pal <- c(sc$pal, "grey40")
    sc$lab <- c(sc$lab, "")          # the sentinel edge carries no number
  }
  d   <- mk(d) %>% dplyr::mutate(zz = pmin(pmax(zf(.data), min(sc$brk)), max(sc$brk)))
  sed <- mk(if (outer == "member")
              sedat %>% dplyr::filter(algo == a, param %in% have, Num_ch %in% nchs)
            else sedat %>% dplyr::filter(param == a, algo %in% have, Num_ch %in% nchs))
  # the union range of every member, pinned into this block; see YLIM
  lim <- mk(YLIM %>% dplyr::filter(Num_ch %in% nchs) %>%
              dplyr::reframe(Num_ch = rep(Num_ch, 2), ly = c(lo, hi)) %>%
              (\(z) if (outer == "member") tidyr::crossing(z, param = have)
                    else tidyr::crossing(z, algo = have))() %>%
              dplyr::mutate(lx = rep(XLIM, length.out = dplyr::n())))
  # A member swept at ONE noise level cannot be contoured: geom_contour_filled needs two distinct
  # y and silently draws nothing, which would read as "ran and came out white" rather than "not
  # run yet". Those cells go down as tiles one decade tall, the device figure_4_common.R already
  # uses for the unidentified marker, cut on the SAME breaks so the colours are the scale's own.
  # The branch stops firing by itself the moment a second noise level lands. INR is the only
  # member in this state today.
  # width 0.40: the seven intervals are log-spaced 0.3 and 0.4 decades apart alternately, so a
  # 0.30 tile leaves a white slot wherever the step is 0.4 and the strip reads as broken
  fill_lyr <- if (length(unique(d$ly)) < MIN_CONTOUR_LEVELS)
    geom_tile(aes(fill = cut(zz, sc$brk, include.lowest = TRUE)), width = 0.40, height = 1.0)
  else geom_contour_filled(breaks = sc$brk)
  g <- ggplot(d, aes(lx, ly, z = zz)) +
    fill_lyr +
    geom_blank(data = lim, aes(lx, ly), inherit.aes = FALSE) +
    scale_fill_manual(values = sc$pal, drop = FALSE, guide = "none") +
    # The two validity criteria as iso-lines of the fill itself, the same pair the colour bar
    # ticks, so a reader can trace where a panel crosses them without counting bands. They were
    # drawn by `figure_4_bias.Rmd` and `figure_4_distortion.Rmd` through mapblock(iso_sol=,
    # iso_dsh=); merging those two files into this one rewrote the block from scratch and left
    # them behind. Put back 2026-07-30. ONE pair serves both halves, where the stacked pair
    # needed two (BIAS_SOL in log10 against DIST_SOL as a ratio), because here both halves are
    # already the same log10 factor: that is what the shared scale buys. Both levels are band
    # EDGES (see SH_FAC), so each line runs along its own colour boundary. The dashed one carries
    # more weight than the solid: a broken line of the same width reads as fainter than it is.
    (if (length(sc$dsh)) geom_contour(breaks = sc$dsh, colour = "grey25", linewidth = .38,
                                     linetype = "dashed") else NULL) +
    (if (length(sc$sol)) geom_contour(breaks = sc$sol, colour = "black", linewidth = .3)
     else NULL) +
    geom_point(size = .04, colour = "grey35", alpha = .3) +
    # The standard error, same two criteria and same style code as the lines above: solid is
    # 1.15, dashed is 2. What separates the two families is COLOUR, not dash. Dark cuts the fill,
    # white cuts the error, and a reader who has learnt "dashed means a factor 2" carries it from
    # one to the other instead of having to unlearn it.
    # se_lines = FALSE for a figure whose FILL already is the standard error: there the white pair
    # would trace the boundary of the field it sits on, which reads as an error rather than as a
    # second reading. `sed` is still computed, so the branch costs nothing but the layer.
    (if (se_lines) list(
      geom_contour(data = sed, aes(lx, ly, z = sev), breaks = SE_SOL,
                   colour = "grey10", linewidth = .85, inherit.aes = FALSE),
      geom_contour(data = sed, aes(lx, ly, z = sev), breaks = SE_SOL,
                   colour = "white", linewidth = .45, inherit.aes = FALSE),
      geom_contour(data = sed, aes(lx, ly, z = sev), breaks = SE_DSH,
                   colour = "grey10", linewidth = .8, linetype = "22", inherit.aes = FALSE),
      geom_contour(data = sed, aes(lx, ly, z = sev), breaks = SE_DSH,
                   colour = "white", linewidth = .42, inherit.aes = FALSE, linetype = "22"))
     else NULL) +
    facet_grid(nchf ~ pf, labeller = label_parsed, scales = "free_y", space = "free_y") +
    # expand = 0 is what actually closes the gap between adjacent maps: the default 5% per side
    # is painted panel background, so two flush panels still show 10% of empty between the data.
    scale_x_continuous(breaks = log10(c(.01, .1, 1)), expand = c(0, 0),
                       labels = if (xlab) c("0.01", "", "1") else c("", "", "")) +
    # the swept label is ten times the dimensionless noise (NOISE_AXIS_UNITS.md); display only
    scale_y_continuous(breaks = YB, expand = c(0, 0),
                       labels = parse(text = fmtnum(dimensionless_noise(10^YB)))) +
    labs(x = NULL, y = if (ylab) "dimensionless instrumental noise" else NULL) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(panel.background = element_rect(fill = "grey85", colour = NA),
          panel.grid.minor = element_blank(), panel.spacing = unit(0, "pt"),
          legend.position = "none",
          axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7),
          axis.title.y = if (ylab) element_text(size = 8, margin = margin(r = 1))
                         else element_blank(),
          strip.text.x = element_text(size = 7, margin = margin(2.5, 0, 2.5, 0)),
          strip.text.y = if (ystrip) element_text(size = 7, angle = 0, lineheight = .9,
                                                  margin = margin(0, 2.5, 0, 2.5))
                         else element_blank(),
          plot.margin = margin(0, 0, 0, 0)) +
    # blanking the ticks is not enough: axis.ticks.length is still reserved, 2.75 pt of dead
    # space at the left edge of every block but the first
    (if (!ylab) theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                      axis.ticks.length.y = unit(0, "pt")) else NULL)
  list(g = g, n = length(have), sc = sc)
}

# ---- assembly -------------------------------------------------------------------------------
# One SIDE of the composition: the member blocks of one quantity, left to right, with their label
# boxes. `first_overall` carries the y axis (only the leftmost block of the leftmost side has it,
# because an axis title is a fixed width that patchwork charges to every row) and `last_overall`
# carries the N_ch strips on the right. The name is historical: there were two sides when this was
# Figure 4's two moments, and there can now be one, two or three.
half <- function(side, nchs, first_overall, last_overall, se_lines = TRUE,
                 outer = "member", drawn_axis = FALSE) {
  keys <- if (outer == "member") ALGOS else PARS
  bs <- lapply(seq_along(keys), function(j)
    blk(keys[j], side, nchs,
        ylab   = first_overall && j == 1,
        ystrip = last_overall  && j == length(keys),
        xlab   = if (drawn_axis) FALSE else j == 1,
        se_lines = se_lines, outer = outer))
  ok <- !vapply(bs, is.null, logical(1)); bs <- bs[ok]
  gs <- list(); ws <- numeric(0); hs <- list()
  for (i in seq_along(bs)) {
    if (i > 1) { gs <- c(gs, list(hsp(GAP_PT))); ws <- c(ws, GAPA)
                 hs <- c(hs, list(hsp(GAP_PT))) }
    gs <- c(gs, list(bs[[i]]$g)); ws <- c(ws, bs[[i]]$n)
    k  <- keys[ok][i]
    hs <- c(hs, list(if (outer == "member")
                       boxed(ALGO_LAB[k], 7, "grey85", gap = BOX_GAP)
                     else boxed(PSHORT[[k]], 7, "grey85", gap = BOX_GAP, parse = TRUE)))
  }
  list(g = gs, w = ws, h = hs, sc = lapply(bs, function(b) b$sc))
}

# ---- the layout constants, all of them fixed except the one that is not ------------------------
GAPA    <- 0.0     # gap between members, in column units (grows with the panels)
GAP_PT  <- 5.5     # gap between members, fixed part, in points: half of what plot_spacer() cost
HALF_PT <- 11      # gap between the two sides, in points
BOX_GAP <- 0.66 * GAP_PT   # under a label box: side to member, member to the parameter strips.
                           # Two thirds of the distance between members, so the vertical breaks
                           # read as subordinate to the horizontal ones rather than competing
CB_TITLE_W  <- 2.6 # the colour bar's title, in the same column units as the maps
# Below this many DISTINCT noise levels a block is drawn as tiles instead of a filled contour.
# Three, not two: `figures_build_plan.md` section 6 bans the two-level contour, which is a linear
# ramp between two rows of data wearing the clothes of a field, and one level draws nothing at all
# and says nothing about it. A member that has not been swept should look like it. The threshold
# stops applying to a member the moment its third noise level lands, with no edit anywhere.
MIN_CONTOUR_LEVELS <- 3
ANC <- c(ttl = 0.22, hrow = 0.20, xrow = 0.10, cb = 0.34, krow = 0.22)   # the fixed rows
OVERHEAD_IN <- 0.990  # the y axis title, its tick labels and the N_ch strips on the right,
                      # MEASURED off two renders (0.993 at 14 columns, 0.987 at 18). Not a
                      # function of the roster, and not the 0.850 a first pass assumed.

# PER_ROW is a nominal knob, NOT the realized scale: the six layout heights are null units, so the
# panels get their share of what is left after every fixed row (strips, tick labels, the colour
# bar) has taken its inches. It therefore has to be re-solved for EVERY column count, which is why
# it is an argument and not a constant. Measured on the PDF, the map row resolves to
#   panel_height = (H - F) * hC/(hC + 0.86) - 0.185        [hC = per_row * summed noise decades]
# with F the composition's total fixed height and 0.185 the parameter strip inside the map row.
# F is NOT a constant of the design: 0.919 while the rows still carried the spacer's margins,
# 0.389 once they came out, 0.476 refit against the four-member render. Any row, strip or margin
# that changes moves it, so the formula is a starting point and never the answer.
#
# DO NOT SOLVE IT ON PAPER. The geometry depends only on the column count and on the noise span of
# each row, NOT on any value inside the maps, so a figure drawn from synthetic cells has the same
# layout as the real one and renders in seconds instead of after a gigabyte of CSV.
#   Rscript tmp/calib_fig4.R                  renders candidates for each roster
#   python3 tmp/measure_fig4.py tmp/calib_*.pdf   measures panel width and band off the raster
# Measured that way, 2026-07-31: 18 columns want per_row 0.179 (x 0.1567, y 0.1562, 99.7 per cent)
# and 16 columns want 0.204 (x 0.1817, y 0.1807, 99.5 per cent). The response is very nearly
# linear, y ~ 0.873 * per_row, so a new column count needs one round of three candidates.
fig4_render <- function(sides, per_row, out, nchs = NCHS,
                        # NULL: take it from the scale the sides carry, which for the default
                        # diverging one is the same string this argument used to hold
                        cb_title = NULL,
                        se_title = "distortion-corrected\nstandard error",
                        xtitle   = "dimensionless sampling interval",
                        se_lines = TRUE, key = TRUE, outer = "member",
                        drawn_axis = FALSE) {
  stopifnot(length(sides) >= 1, length(sides) <= 3, outer %in% c("member", "param"))
  hs2 <- lapply(seq_along(sides), function(i)
    half(sides[[i]], nchs, first_overall = i == 1, last_overall = i == length(sides),
         se_lines = se_lines, outer = outer, drawn_axis = drawn_axis))
  # Interleave the sides with the HALF_PT hole. One side and there is no hole; three and there are
  # two. Everything downstream reads the assembled widths, so nothing else knows how many there are.
  W <- HW <- numeric(0); Gs <- Hs <- Ts <- Xs <- list()
  for (i in seq_along(hs2)) {
    if (i > 1) {
      W  <- c(W, 0.30);            HW <- c(HW, 0.30)
      Gs <- c(Gs, list(hsp(HALF_PT))); Hs <- c(Hs, list(hsp(HALF_PT)))
      Ts <- c(Ts, list(hsp(HALF_PT))); Xs <- c(Xs, list(hsp(HALF_PT)))
    }
    W  <- c(W, hs2[[i]]$w);        HW <- c(HW, sum(hs2[[i]]$w))
    Gs <- c(Gs, hs2[[i]]$g);       Hs <- c(Hs, hs2[[i]]$h)
    Ts <- c(Ts, list(boxed(sides[[i]]$lab, 8, "grey92", sides[[i]]$letter, BOX_GAP)))
    Xs <- c(Xs, list(NULL))   # filled below: the x row is built once, across the whole width
  }

  ttl  <- wrap_plots(Ts, nrow = 1, widths = HW)
  hrow <- wrap_plots(Hs, nrow = 1, widths = W)
  row  <- wrap_plots(Gs, nrow = 1, widths = W)
  # THE X ROW, ONCE ACROSS THE FIGURE. It used to be one title per side, which is fine while the
  # sides are the same width and breaks the moment they are not: a five-column side got the same
  # sentence as a thirteen-column one and wore it. Now the strip sits under the first column and
  # the title is centred on everything to its right, whatever the sides are doing.
  xrow <- if (drawn_axis)
            wrap_plots(list(xaxis_strip(), txt(xtitle, top = 2, bot = 4)),
                       nrow = 1, widths = c(1, sum(W) - 1))
          else wrap_plots(list(txt(xtitle, top = 2, bot = 4)), nrow = 1, widths = sum(W))
  # ONE BAR PER SCALE, not one bar per figure. Sides that share a scale share a bar spanning the
  # whole width, which is the case for every figure whose quantities are commensurable and is what
  # keeps a colour meaning one thing across the set. Sides that do NOT (an autocorrelation time
  # next to a condition number: 1 to 13 against 400 to 1e7) get a bar each, under their own
  # columns, because putting them on one scale would leave one of the two saturated in every cell.
  # the scale as the BLOCKS ended up drawing it, terminal grey band included
  scs    <- lapply(seq_along(sides), function(i) {
    b <- Filter(Negate(is.null), hs2[[i]]$sc)
    if (length(b)) b[[1]] else if (is.null(sides[[i]]$scale)) SCALE_FACTOR else sides[[i]]$scale
  })
  shared <- length(scs) == 1 || all(vapply(scs[-1], function(s) identical(s, scs[[1]]), TRUE))
  mkbar  <- function(sc, tw, ww, ttl_override = NULL) {
    lab <- if (!is.null(ttl_override) && nzchar(ttl_override)) ttl_override else sc$title
    wrap_plots(list(cbtitle(lab), cbar(sc$pal, sc$lab, sc$brk, sol = sc$sol, dsh = sc$dsh)),
               nrow = 1, widths = c(tw, ww - tw))
  }
  cb <- if (shared) {
    mkbar(scs[[1]], CB_TITLE_W, sum(W), cb_title)
  } else {
    # a narrower title per bar, since each one now has a side's width and not the page's
    parts <- list()
    for (i in seq_along(scs)) {
      if (i > 1) parts <- c(parts, list(hsp(HALF_PT)))
      # the title cannot take a fixed 2.2 columns of a side that only has five: cap it at a third
      # of the side, so a narrow bar keeps two thirds of its width for the bands
      wi <- HW[[if (i == 1) 1 else 2 * i - 1]]
      parts <- c(parts, list(mkbar(scs[[i]], min(2.2, wi / 3), wi)))
    }
    wrap_plots(parts, nrow = 1, widths = HW)
  }
  # Only the standard error gets a key. The iso lines on the fill do not need one: their levels
  # are edges of the colour scale, so the bar names them where the line sits.
  fac  <- function(v) sprintf("a factor %g", v)
  krow <- wrap_plots(list(cbtitle(se_title),
                          linekey(fac(c(CRIT_SOL, CRIT_DSH)),
                                  cased = TRUE, chip = "grey62", dash = "22")),
                     nrow = 1, widths = c(CB_TITLE_W, sum(W) - CB_TITLE_W))

  dec <- sum(vapply(nchs, nch_noise_span, numeric(1)))
  hC  <- per_row * dec
  # The key row goes only where there is a key to explain. Dropping it changes the null-unit split
  # and therefore per_row, so a figure without it is calibrated separately; see the note above.
  g   <- if (key) ttl / hrow / row / xrow / cb / krow +
           plot_layout(heights = c(ANC[["ttl"]], ANC[["hrow"]], hC,
                                   ANC[["xrow"]], ANC[["cb"]], ANC[["krow"]]))
         else ttl / hrow / row / xrow / cb +
           plot_layout(heights = c(ANC[["ttl"]], ANC[["hrow"]], hC,
                                   ANC[["xrow"]], ANC[["cb"]]))
  H   <- hC + sum(ANC) - if (key) 0 else ANC[["krow"]]
  ggsave(out, g, width = 7.4, height = H)

  # The x scale is DERIVED: the page less OVERHEAD_IN less the member gaps, which are the only
  # part that grows with the roster. The y scale is the formula above. Both are printed next to
  # what they claim so the number to check against the PDF is on the log.
  ncol   <- sum(vapply(hs2, function(h) sum(h$w), numeric(1)))
  gaps   <- sum(vapply(hs2, function(h) (length(h$w) - 1) / 2, numeric(1))) * GAP_PT / 72 +
            (length(hs2) - 1) * HALF_PT / 72
  colstr <- paste(vapply(hs2, function(h) sprintf("%d", round(sum(h$w))), ""), collapse = " + ")
  wpan  <- 7.4 - OVERHEAD_IN - gaps
  pan   <- (H - 0.476) * hC / (hC + 0.86) - 0.185
  cat(sprintf(paste0("%s: 7.40 x %.2f in, %s columns, per_row %.4f nominal\n",
                     "  x derived: %.4f in per decade, %.2f mm per cell over 7 intervals\n",
                     "  y expected: %.4f in per decade (%.2f in of panel over %.3f decades)\n",
                     "  isotropy claimed %.1f%% -> MEASURE THE PDF, this is the formula\n"),
              out, H, colstr, per_row,
              wpan / (2 * ncol), 25.4 * wpan / (7 * ncol),
              pan / dec, pan, dec, 100 * (pan / dec) / (wpan / (2 * ncol))))
  invisible(g)
}
