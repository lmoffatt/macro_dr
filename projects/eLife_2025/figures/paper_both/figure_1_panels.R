## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
library(patchwork)


## -----------------------------------------------------------------------------
Nch= 20

fs = 50e3
nsamples= 100
noise=1e-4 *fs/nsamples


## -----------------------------------------------------------------------------
# The two vectors are zipped BY POSITION (see the read loop below): row i of `filenames`
# is labelled with row i of `algorithms`. A wrong filename errors at read.csv; a wrong
# ORDER silently mislabels a whole column. Edit them together, always.
filenames=c(
  "../data/figure_1_likelihood_diagnostic_LSE_av0",
  "../data/figure_1_likelihood_diagnostic_LSE",
  "../data/figure_1_likelihood_diagnostic_NR",
  "../data/figure_1_likelihood_diagnostic_INR",
  "../data/figure_1_likelihood_diagnostic_R",
  "../data/figure_1_likelihood_diagnostic_MR",
  "../data/figure_1_likelihood_diagnostic_VR",
  "../data/figure_1_likelihood_diagnostic_IR",
  "../data/figure_1_simulation"

)

algorithms=c("LSE"
             ,"ILSE"
             ,"NR"
             ,"INR"
             ,"R"
             ,"MR"
             ,"VR"
             ,"IR"
             ,"SIM")



## -----------------------------------------------------------------------------

d=data.frame()
for (i in 1:length(filenames)){
  f <- paste0(filenames[i], ".csv")
  if (!file.exists(f)) { message("figure_1_panels: skipping missing dump ", f); next }  # LSE dump may not exist yet
  d0=read.csv(f, check.names = FALSE)
  d0$algo=algorithms[i]
  d=rbind(d,d0)
}



## -----------------------------------------------------------------------------
d_s=d%>%filter(algo=="SIM", scope=="simulation_sub")%>%select(sub_index,sample_index,step_start,step_middle,step_end,component_path,value)%>%pivot_wider(names_from = component_path, values_from = value)

d_s_mean<-d_s%>%group_by(sample_index)%>%summarize(step_start=min(step_start), step_end=max(step_end), patch_current= mean(patch_current))

d_IR=d%>%filter(algo=="IR", scope=="evolution", value_row==1| is.na(value_row),  is.na(value_col), !is.na(step_start))%>%mutate(component= component_path %>%
           str_replace("^Algo_State_Dynamic\\.", "") %>%
           str_replace("^Patch_State\\.", ""))%>%
  select(sample_index,step_start,step_middle,step_end,component,patch_current,value)%>%
  pivot_wider(names_from = component, values_from = value)


## -----------------------------------------------------------------------------
d_NR=d%>%filter(algo=="NR", scope=="evolution", value_row==1| is.na(value_row),  is.na(value_col), !is.na(step_start))%>%mutate(component= component_path %>%
           str_replace("^Algo_State_Dynamic\\.", "") %>%
           str_replace("^Patch_State\\.", ""))%>%
  select(sample_index,step_start,step_middle,step_end,component,patch_current,value)%>%
  pivot_wider(names_from = component, values_from = value)

# LSE (least squares on the mean): non-recursive like NR. Its predictive variance is a single CONSTANT
# sigma_hat^2 = SSE/n (the marginalised noise scale that makes the logL well-defined), so its row-B
# band is FLAT. Dump = figure_1.macroir's LSE block (family = 2). Columns used are VERIFIED present in
# figure_1_likelihood_diagnostic_LSE.csv: P_mean_t20_y1, y_mean, y_var, patch_current, logL.
# Empty until that dump exists (tolerant read).
d_LSE=d%>%filter(algo=="LSE", scope=="evolution", value_row==1| is.na(value_row),  is.na(value_col), !is.na(step_start))%>%mutate(component= component_path %>%
           str_replace("^Algo_State_Dynamic\\.", "") %>%
           str_replace("^Patch_State\\.", ""))%>%
  select(sample_index,step_start,step_middle,step_end,component,patch_current,value)%>%
  pivot_wider(names_from = component, values_from = value)

# ILSE: least squares WITH the acquisition average (av = 1). Same extraction as LSE above;
# the family emits the same keys for both arms, verified against both dumps on 2026-08-05.
d_ILSE=d%>%filter(algo=="ILSE", scope=="evolution", value_row==1| is.na(value_row),  is.na(value_col), !is.na(step_start))%>%mutate(component= component_path %>%
           str_replace("^Algo_State_Dynamic\\.", "") %>%
           str_replace("^Patch_State\\.", ""))%>%
  select(sample_index,step_start,step_middle,step_end,component,patch_current,value)%>%
  pivot_wider(names_from = component, values_from = value)

d_INR=d%>%filter(algo=="INR", scope=="evolution", value_row==1| is.na(value_row),  is.na(value_col), !is.na(step_start))%>%mutate(component= component_path %>%
           str_replace("^Algo_State_Dynamic\\.", "") %>%
           str_replace("^Patch_State\\.", ""))%>%
  select(sample_index,step_start,step_middle,step_end,component,patch_current,value)%>%
  pivot_wider(names_from = component, values_from = value)

d_R=d%>%filter(algo=="R", scope=="evolution", value_row==1| is.na(value_row),  is.na(value_col), !is.na(step_start))%>%mutate(component= component_path %>%
           str_replace("^Algo_State_Dynamic\\.", "") %>%
           str_replace("^Patch_State\\.", ""))%>%
  select(sample_index,step_start,step_middle,step_end,component,patch_current,value)%>%
  pivot_wider(names_from = component, values_from = value)

d_MR=d%>%filter(algo=="MR", scope=="evolution", value_row==1| is.na(value_row),  is.na(value_col), !is.na(step_start))%>%mutate(component= component_path %>%
           str_replace("^Algo_State_Dynamic\\.", "") %>%
           str_replace("^Patch_State\\.", ""))%>%
  select(sample_index,step_start,step_middle,step_end,component,patch_current,value)%>%
  pivot_wider(names_from = component, values_from = value)

# VR is built from MR's formulas with the RESIDUAL predictive variance instead of the
# total one, so its panels below are MR's CODE on d_VR — but not MR's numbers. y_var
# divides the innovation (chi = dy / y_var) and the covariance down-date ((N/y_var)*XTX(gS)),
# so the whole recursion diverges. Measured on this recording: y_var 26% below MR at the
# median; sample 0 identical (no agonist) and sample 1 identical in mean and covariance
# only because that prior is deterministic (gS = 0, so the update is null whatever the
# variance); posterior apart from sample 2, propagated mean from sample 3.
d_VR=d%>%filter(algo=="VR", scope=="evolution", value_row==1| is.na(value_row),  is.na(value_col), !is.na(step_start))%>%mutate(component= component_path %>%
           str_replace("^Algo_State_Dynamic\\.", "") %>%
           str_replace("^Patch_State\\.", ""))%>%
  select(sample_index,step_start,step_middle,step_end,component,patch_current,value)%>%
  pivot_wider(names_from = component, values_from = value)






## -----------------------------------------------------------------------------
common_theme = theme_classic()+ theme( legend.position = "bottom")

# Okabe-Ito colour-blind-safe semantic palette, shared across panels so patchwork guides="collect"
# merges everything into ONE legend. Colour = FILTER PHASE: predict (orange) groups the Markov step
# with the prior it produces; update (purple) groups the Bayes step with the posterior it produces
# (this operation->state link holds for MacroR/MR/VR/IR alike). SHAPE separates operation (arrow) from
# state (point/segment); ALPHA is the carried-over channel: an arrow's ORIGIN (the state carried in
# from the previous step) is drawn faded (ORIG_ALPHA) while the state the panel PRODUCES stays at full
# alpha. One mechanism for points and segments (the prior produced by the Markov panel reappears faded
# as the origin of the Bayes panel's arrow, and the posterior the other way). Hue and alpha never collide.
SEM <- c("channel current"          = "#7FCDBB",   # light green — the true resolved current (background reality)
         "Markov & prior"            = "#E69F00",   # orange      — PREDICT: Markov step, prior, AND predicted current
         "observation & innovation"  = "#0072B2",   # blue        — OBSERVE: observed current AND the innovation (obs - pred)
         "Bayes & posterior"         = "#CC79A7",   # purple      — UPDATE : Bayes step AND posterior
         "logLikelihood"             = "#000000")   # black       — the step OUTPUT/score; achromatic + colorblind-safe. Plot as a BOLD line so it does not blend with the (black) axes
sem_scale <- ggplot2::scale_colour_manual(values = SEM, limits = names(SEM), name = NULL, drop = FALSE,
                                          breaks = setdiff(names(SEM), "logLikelihood"),   # logL row is self-labelled -> drop from legend
                                          labels = function(x) ifelse(x == "Markov & prior", "Markov, prior & prediction", x))  # orange spans prior (state) + prediction (measurement)

# common x (time) range over the FULL recording (anchored at 0); each figure applies its OWN crop
# later via build_figure(), so the arrows are always the full-graph arrows, just zoomed.
XLIM <- range(c(0, d_s$step_middle, d_IR$step_start, d_IR$step_end), na.rm = TRUE)
xcommon <- coord_cartesian(xlim = XLIM)

# --- common y-scales, recomputed FROM THE DATA so they adapt if the data changes (no hardcoding) ---
.dP <- function(x) unlist(x[, grepl("^P_mean", names(x)), drop = FALSE])          # every P(open) column
.dI <- function(x) c(-x$patch_current+1, -x$y_mean+1, -x$y_mean+1 - sqrt(x$y_var), -x$y_mean+1 + sqrt(x$y_var))  # current (with +1 baseline compensation)
.algos <- list(d_NR, d_INR, d_R, d_MR, d_VR, d_IR)
YP  <- range(unlist(lapply(.algos, .dP)), na.rm = TRUE)                           # prior & posterior share this
YI  <- range(c(-d_s$patch_current, unlist(lapply(.algos, .dI))), na.rm = TRUE)    # observation row

compact_theme <- common_theme +theme(
  axis.title.x = element_blank(),
  axis.text.x  = element_blank(),
  axis.ticks.x = element_blank(),
  legend.position = "right",
  plot.margin = margin(1, 5, 1, 5)
)
x_only_ticks=theme(axis.title.x = element_blank(), axis.text.x = element_blank())
y_only_ticks=theme(axis.title.y = element_blank(), axis.text.y = element_blank())
no_y_axis= theme( axis.text.y = element_blank(), # Removes y-axis text
         axis.title.y = element_blank(), # Removes y-axis title
          axis.ticks.y = element_blank(),
        axis.line.y = element_blank())
no_x_axis= theme( axis.text.x = element_blank(), # Removes y-axis text
         axis.title.x = element_blank(), # Removes y-axis title
          axis.ticks.x = element_blank(),
        axis.line.x = element_blank())
y_no_title=theme(axis.title.y = element_blank())


## -----------------------------------------------------------------------------
markov_arrow = arrow(length = unit(0.07, "in"), type="closed")
bayes_arrow= arrow(length = unit(0.07, "in"), type="closed")
pred_arrow = arrow(length = unit(0.07, "in"), type="closed")

# BOW OF THE MARKOV ARC, 2026-08-12. geom_curve places its control point at a perpendicular offset of
# about (curvature/2) x (segment length) in INCHES, so the bow deepens when the panel widens while
# the data range does not. At 0.6 and six columns the arc already dipped to the axis in the open-loop
# columns, whose priors sit at the bottom of row A's shared range; at four columns the panel went
# from ~1.0 to ~1.55 in and the arc left the panel, taking the arrow with it (Luciano, on the first
# four-column render). Two things fix it together and neither alone: this value, and the asymmetric
# bottom pad on yP in build_figure. RE-CHECK BOTH if the roster length or the figure width changes.
MK_CURV <- 0.3

# de-emphasis for an arrow's ORIGIN: the state carried in from the previous step is drawn faded, so it
# recedes to context while the full-colour mark is the state the panel produces. One value for points
# and segments alike; kept high enough that the hue still reads (not grey) after PDF downsample.
ORIG_ALPHA = 0.45




## -----------------------------------------------------------------------------
fA<-ggplot(d_IR)+geom_line(data=d_s, aes(x=step_middle, y= -patch_current, color="channel current"),
                linewidth = 0.5) +
  geom_segment(aes(x=step_start,xend=step_end,y=-patch_current+1,color="observation & innovation"), linewidth = 1)+
  geom_line(data=data.frame(role=names(SEM), x=NA_real_, y=NA_real_), aes(x=x, y=y, colour=role), na.rm=TRUE)+
  sem_scale + ylab("current (pA)") + common_theme+x_only_ticks + xcommon
fA


## ----fB_IR--------------------------------------------------------------------
## IR's state is a PAIR of endpoint occupancies plus the covariance between them, not a trajectory
## across the interval, so it is drawn as two marks tied by a dashed arc and not as a segment: a
## segment claims a path inside the interval that the algorithm never computes, and a straight tie
## keeps tracing that path even when dashed, which is why the tie bows.
## The tie is a quadratic Bezier built in INCHES. The panel box is not square and it changes shape
## between the 7.0 in body figure and the 5.6 in supplement, so a bow defined in normalised units
## comes out a different shape in each. Being parametric is what lets the operation arrows land ON the
## arcs instead of on the chord, a straight line that is not drawn: the Bayes arrow is anchored where
## each arc crosses the interval midpoint (ir_y_at), exactly, at both of its ends.
##
## WHERE THE MARKOV STEP STARTS. At the shared disc, not at the middle of the previous pair: only ONE
## marginal crosses an interval boundary. The propagation takes the occupancy at t_j and builds the
## joint over (t_j, t_j+1); the previous interval's start state is marginalised out and does not
## enter, which is the Markov property. An arrow leaving the middle of the previous arc would claim
## that the whole pair feeds the step. Where it ENDS is unchanged, the middle of the arc it produces,
## because what the step delivers is the pair and not either of its endpoints: the right endpoint on
## its own is one more marginal, and the covariance would be left out of the picture.
##
## THE SHARED INSTANT. The right end of the posterior pair of interval j IS the left end of the
## prior pair of interval j+1: one instant with two roles, which is why a single-coloured mark reads
## wrong there. It is a disc split down the middle, the left half in the colour of the posterior it
## closes and the right half in the colour of the prior it opens, the split following the time axis
## so that left is before and right is after. The alpha convention decides which half is faded (full
## for the state the panel produces, faded for the one it carries in), so the two halves swap
## strength between the prior row and the posterior row. Both rows draw the same marks; they differ
## only in that swap and in the operation arrow.
##
## COST: these two panels are the only ones in this file that are NOT window-agnostic. They need the
## crop and the panel box, so build_figure() computes both and builds them per call instead of
## building once and cropping.
IR_N_ARC   <- 61          # odd, so index (n+1)/2 IS t = 0.5
IR_TIE_LWD <- 0.35        # thinner than the operation arrows (0.5): the tie must not read as one
IR_TIE_LTY <- "dashed"
IR_RHO     <- 0.033       # disc radius, inches; ~20% over a size-1.5 point, the only two-role mark
IR_PT_R    <- 0.027       # radius of a size-1.5 point, inches (0.375 * its fontsize)
IR_GAP     <- 0.006       # air between an arrowhead or tail and the mark it points at

# data <-> inches, for one panel box (pw, ph in inches)
ir_gmap <- function(xr, yr, pw, ph) list(
  x2i = function(x) (x - xr[1]) / diff(xr) * pw, y2i = function(y) (y - yr[1]) / diff(yr) * ph,
  i2x = function(u) xr[1] + u / pw * diff(xr),   i2y = function(v) yr[1] + v / ph * diff(yr), ph = ph)

ir_arc <- function(x1, y1, x2, y2, g, n = IR_N_ARC) {
  u1 <- g$x2i(x1); u2 <- g$x2i(x2); v1 <- g$y2i(y1); v2 <- g$y2i(y2)
  nu <- -(v2 - v1); nv <- u2 - u1; L <- sqrt(nu^2 + nv^2)      # chord rotated +90 deg, |n| = chord
  off <- min(0.15 * L, 0.12 * g$ph)                            # apex offset in inches, capped so a
  cu <- (u1 + u2)/2 + 2*off*nu/L                               # wide panel does not get a balloon
  cv <- (v1 + v2)/2 + 2*off*nv/L
  t <- seq(0, 1, length.out = n)
  data.frame(x = g$i2x((1-t)^2*u1 + 2*t*(1-t)*cu + t^2*u2),
             y = g$i2y((1-t)^2*v1 + 2*t*(1-t)*cv + t^2*v2))
}
ir_ok <- function(dd) which(complete.cases(dd[, c("x1","y1","x2","y2")]))

ir_tie <- function(dd, g, alpha) {
  i <- ir_ok(dd); if (!length(i)) return(NULL)
  p <- do.call(rbind, lapply(i, function(k) {
    a <- ir_arc(dd$x1[k], dd$y1[k], dd$x2[k], dd$y2[k], g); a$id <- k; a$role_ <- dd$role_[k]; a }))
  geom_path(data = p, aes(x = x, y = y, group = id, colour = role_),
            linetype = IR_TIE_LTY, linewidth = IR_TIE_LWD, alpha = alpha, na.rm = TRUE)
}
ir_mid <- function(dd, g) {                                    # exact arc midpoints, one per pair
  m <- (IR_N_ARC + 1)/2; out <- data.frame(mx = rep(NA_real_, nrow(dd)), my = rep(NA_real_, nrow(dd)))
  for (k in ir_ok(dd)) { a <- ir_arc(dd$x1[k], dd$y1[k], dd$x2[k], dd$y2[k], g)
    out$mx[k] <- a$x[m]; out$my[k] <- a$y[m] }
  out
}
# pull the two ends of a chord in by d1 and d2 INCHES, so an arrow stops short of the marks it runs
# between instead of emerging from under them
ir_shrink <- function(x1, y1, x2, y2, g, d1, d2) {
  u1 <- g$x2i(x1); u2 <- g$x2i(x2); v1 <- g$y2i(y1); v2 <- g$y2i(y2)
  L <- sqrt((u2 - u1)^2 + (v2 - v1)^2); L[!is.na(L) & L == 0] <- NA
  fu <- (u2 - u1)/L; fv <- (v2 - v1)/L
  data.frame(x    = g$i2x(u1 + d1*fu), y    = g$i2y(v1 + d1*fv),
             xend = g$i2x(u2 - d2*fu), yend = g$i2y(v2 - d2*fv))
}
ir_y_at <- function(dd, g, x0) {                               # where an arc crosses a vertical line
  out <- rep(NA_real_, nrow(dd))
  for (k in intersect(ir_ok(dd), which(!is.na(x0)))) {
    a <- ir_arc(dd$x1[k], dd$y1[k], dd$x2[k], dd$y2[k], g, n = 801)
    out[k] <- a$y[which.min(abs(a$x - x0[k]))] }
  out
}
# Half of the shared-instant disc, as a grob in ABSOLUTE units placed at a data coordinate. Built as
# a polygon in data space it came out an ellipse, because that needs the true panel box and the
# estimate below is off by ~15% in aspect; in inches it is round whatever the box turns out to be,
# and the radius is exact. Only the position is data-driven.
ir_disc <- function(pts, side, role, g, alpha, n = 28) {
  pts <- pts[complete.cases(pts), , drop = FALSE]; if (!nrow(pts)) return(NULL)
  th <- if (side == "left") seq(pi/2, 3*pi/2, length.out = n) else seq(-pi/2, pi/2, length.out = n)
  gb <- grid::polygonGrob(x = grid::unit(0.5, "npc") + grid::unit(IR_RHO * cos(th), "in"),
                          y = grid::unit(0.5, "npc") + grid::unit(IR_RHO * sin(th), "in"),
                          gp = grid::gpar(fill = adjustcolor(SEM[[role]], alpha.f = alpha), col = NA))
  lapply(seq_len(nrow(pts)), function(i)
    annotation_custom(gb, xmin = pts$x[i], xmax = pts$x[i], ymin = pts$y[i], ymax = pts$y[i]))
}

## the state pairs, one per acquisition interval j:
##   posterior j : (start_j, P_mean_t10_y1_j)     -> (end_j, P_mean_t20_y1_j)
##   prior     j : (start_j, P_mean_t20_y1_{j-1}) -> (end_j, P_mean_t11_y0_j)
## so the shared instant is (start_j, P_mean_t20_y1_{j-1}). A posterior's right end is shared with
## the next prior's left end EXCEPT in the last interval, which has no next prior and keeps a circle.
ir_pairs <- function(d) {
  v <- d %>% transmute(xs = step_start, xm = step_middle, xe = step_end,
                       post_y1 = P_mean_t10_y1, post_y2 = P_mean_t20_y1,
                       pri_y1  = lag(P_mean_t20_y1), pri_y2 = P_mean_t11_y0)
  list(v = v,
       post = with(v, data.frame(x1=xs, y1=post_y1, x2=xe, y2=post_y2, role_="Bayes & posterior")),
       pri  = with(v, data.frame(x1=xs, y1=pri_y1,  x2=xe, y2=pri_y2,  role_="Markov & prior")),
       share  = data.frame(x = v$xs, y = v$pri_y1),
       post_l = data.frame(x = v$xs, y = v$post_y1),
       pri_r  = data.frame(x = v$xe, y = v$pri_y2),
       post_r = data.frame(x = v$xe, y = v$post_y2)[is.na(lead(v$pri_y1)), , drop = FALSE])
}
IRP <- ir_pairs(d_IR)

## the marks, identical in both state rows; a_post / a_pri carry the produced-vs-carried-in swap
ir_marks <- function(g, P, a_post, a_pri) c(
  list(ir_tie(P$post, g, a_post), ir_tie(P$pri, g, a_pri),
       geom_point(data = P$post_l, aes(x, y, colour = "Bayes & posterior"), alpha = a_post, na.rm = TRUE),
       geom_point(data = P$post_r, aes(x, y, colour = "Bayes & posterior"), alpha = a_post, na.rm = TRUE),
       geom_point(data = P$pri_r,  aes(x, y, colour = "Markov & prior"),    alpha = a_pri,  na.rm = TRUE)),
  ir_disc(P$share, "left",  "Bayes & posterior", g, a_post),   # one layer per disc: annotation_custom
  ir_disc(P$share, "right", "Markov & prior",    g, a_pri))    # does not vectorise over positions

ir_prior_panel <- function(g, P = IRP) {
  # ORIGIN at the shared disc, the only thing that crosses the boundary; destination unchanged, the
  # middle of the arc, because what the step produces is the pair and not either of its ends
  m <- ir_mid(P$pri, g)
  ar <- ir_shrink(P$v$xs, P$v$pri_y1, m$mx, m$my, g, IR_RHO + IR_GAP, 0)
  ggplot() + ir_marks(g, P, ORIG_ALPHA, 1) +
    geom_curve(data = ar, aes(x = x, y = y, xend = xend, yend = yend, colour = "Markov & prior"),
               curvature = MK_CURV, angle = 90, ncp = 10, linewidth = 0.5,
               arrow = markov_arrow, na.rm = TRUE) +
    scale_fill_manual(values = SEM, guide = "none") +
    sem_scale + ylab("P(open)") + common_theme + x_only_ticks + guides(colour="none") + xcommon
}
ir_post_panel <- function(g, P = IRP) {
  x0 <- P$v$xm                                                  # anchored on both arcs at the same
  ar <- data.frame(x = x0, y = ir_y_at(P$pri, g, x0),           # x, so the update arrow is vertical
                   xend = x0, yend = ir_y_at(P$post, g, x0))
  ggplot() + ir_marks(g, P, 1, ORIG_ALPHA) +
    geom_segment(data = ar, aes(x = x, y = y, xend = xend, yend = yend, colour = "Bayes & posterior"),
                 linewidth = 0.5, arrow = bayes_arrow, na.rm = TRUE) +
    scale_fill_manual(values = SEM, guide = "none") +
    sem_scale + ylab("P(open)") + common_theme + x_only_ticks + guides(colour="none") + xcommon
}

# registry default, for the whole recording; build_figure() rebuilds both with the real crop and box
fB_IR <- ir_prior_panel(ir_gmap(XLIM, YP, 1.06, 1.18))


## -----------------------------------------------------------------------------
fC_IR<-ggplot(d_IR)+
  geom_rect(aes(xmin=step_start, xmax=step_end, ymin=-y_mean+1-sqrt(y_var), ymax=-y_mean+1+sqrt(y_var), fill="Markov & prior"), alpha=0.3)+
  geom_segment(aes(x=step_start, xend= step_end, y=-patch_current+1, color="observation & innovation"), linewidth=1, alpha=1)+
  geom_segment(aes(x=step_start, xend= step_end, y=-y_mean+1, color="Markov & prior"), linewidth=1)+
  geom_segment(aes(x=(step_start+step_end)/2, xend=(step_start+step_end)/2, y=-y_mean+1, yend=-patch_current+1, color="observation & innovation"),
               arrow=pred_arrow, linewidth=0.5, na.rm=TRUE)+
  scale_fill_manual(values = SEM, guide = "none") +
  sem_scale + ylab("predicted current (pA)") + common_theme+x_only_ticks + guides(colour="none") + xcommon
fC_IR



## -----------------------------------------------------------------------------

# Same marks as the prior row with the produced/carried-in strengths swapped (see the block above);
# registry default only, build_figure() rebuilds it with the real crop and panel box.
fD_IR <- ir_post_panel(ir_gmap(XLIM, YP, 1.06, 1.18))





## -----------------------------------------------------------------------------
fE_IR<-ggplot(d_IR)+
 # geom_segment(aes(x=step_start, xend=step_end, y=logL))+
  geom_line(aes(x=step_end,  y=cumsum(logL), color="logLikelihood"))+
  geom_point(aes(x=step_end,  y=cumsum(logL), color="logLikelihood"))+
  

  sem_scale + ylab("logL") + common_theme+guides(colour="none") + xcommon+xlab("time (ms)")
fE_IR



## -----------------------------------------------------------------------------
# ROW A, THE INSTANT EACH MEMBER PREDICTS AT (2026-08-05). The x is set by the averaging axis and
# NOT by which P_mean key the member happens to emit. A member with av = 0 predicts the current at
# the SAMPLE, which sits at the middle of the acquisition window, so its prior is drawn at
# step_middle; a member with av = 1 predicts the average OVER the window, whose prior is the
# occupancy the window opens with, so it is drawn at step_start.
#   step_middle   LSE, NR, R      (av = 0, instantaneous)
#   step_start    ILSE, INR, IR   (av = 1, interval)
# Aligning by the key's own meaning instead (t2_y0 at step_end, t20_y1 at step_start) was tried and
# is wrong: it puts NR half an interval late, because what row A shows is the state that produced
# the prediction drawn in row B, not wherever the emitted key happens to sit.
# NR reads P_mean_t15_y0, the mid-window state, and not P_mean_t2_y0. The engine emits both since
# 2026-08-05: t2_y0 is the occupancy the window CLOSES with, which is the prior threaded into the
# next window, and t15_y0 is the one the prediction is taken FROM, half an interval in. NR is av=0,
# so those differ, and row A wants the second. Verified against the rebuilt dump: t15_y0 equals the
# nonlinearsqr av=0 arm's emitted series to seven decimals (0.0906346, 0.2255942, 0.3160603,
# 0.3767015, 0.4173506), which is the same occupancy recovered independently from NR's own y_mean.
fB_NR<-ggplot(d_NR)+
  geom_point(aes(x=step_middle, y=P_mean_t15_y0, color="Markov & prior"), linewidth = 1, alpha=1)+
geom_curve(aes(x=lag(step_middle), xend = step_middle, y=lag(P_mean_t15_y0),
   yend = P_mean_t15_y0, color="Markov & prior"),curvature = MK_CURV,  angle = 90, ncp=10,linetype = 1, linewidth = 0.5, arrow=markov_arrow)+
  sem_scale + ylab("P(open)") + common_theme+x_only_ticks + guides(colour="none") + xcommon

fB_NR


## -----------------------------------------------------------------------------
fC_NR<-ggplot(d_NR)+
  geom_errorbar(aes(x=step_middle, ymin=-y_mean+1-sqrt(y_var), ymax=-y_mean+1+sqrt(y_var),width=0.0005, color="Markov & prior"), alpha=1)+
  geom_segment(aes(x=step_start, xend= step_end, y=-patch_current+1, color="observation & innovation"), linewidth=1, alpha=1)+
  geom_point(aes(x=step_middle,  y=-y_mean+1, color="Markov & prior"), linewidth=1)+
  geom_segment(aes(x=(step_start+step_end)/2, xend=(step_start+step_end)/2, y=-y_mean+1, yend=-patch_current+1, color="observation & innovation"),
               arrow=pred_arrow, linewidth=0.5, na.rm=TRUE)+
  scale_fill_manual(values = SEM, guide = "none") +
  sem_scale + ylab("predicted current (pA)") + common_theme+x_only_ticks + guides(colour="none") + xcommon
fC_NR


## ----LSE panels: like NR (non-recursive), with a CONSTANT least-squares variance (Luciano 2026-07-23)
## LSE is least squares on the mean, noise scale marginalised to a single plug-in sigma_hat^2 = SSE/n
## (mean(r_std^2)=1 confirms y_var IS that). Drawn in NR's grammar — POINT mean + ERRORBAR — NOT
## macroir's band + interval-spanning segment. Its errorbars are all the SAME width (the flat
## sigma_hat^2), against NR/R's width that varies with the per-interval gating variance.
## fB_LSE (row A, P(open)) uses P_mean_t20_y1 — VERIFIED present in the LSE dump (its ONLY P_mean key;
## LSE is built averaging=1 so it follows IR's naming, NOT NR's P_mean_t2_y0). y_mean, y_var,
## patch_current, logL also verified present in figure_1_likelihood_diagnostic_LSE.csv.
fB_LSE<-ggplot(d_LSE)+
  geom_point(aes(x=step_middle, y=P_mean_t20_y1, color="Markov & prior"), linewidth = 1, alpha=1)+
geom_curve(aes(x=lag(step_middle), xend = step_middle, y=lag(P_mean_t20_y1),
   yend = P_mean_t20_y1, color="Markov & prior"),curvature = MK_CURV,  angle = 90, ncp=10,linetype = 1, linewidth = 0.5, arrow=markov_arrow)+
  sem_scale + ylab("P(open)") + common_theme+x_only_ticks + guides(colour="none") + xcommon

# THE LEAST-SQUARES PAIR, 2026-08-05. LSE is av = 0 and ILSE is av = 1, and they differ in the data:
# at the first interval y_var reads 8.39 for the un-averaged arm against 0.94 for the averaged one.
# The pair is drawn so that the ONE thing that changes is the thing the pair is about. Both bands are
# a single constant sigma_hat^2, so neither breathes with the gating; what differs is whether the
# prediction is a POINT at the sample (av = 0, NR's grammar) or a SEGMENT across the interval
# (av = 1, the grammar the macro interval members use). This is the reason the pair earns a column:
# on NR/INR and R/IR the band width moves as well, so two things move at once, and here only one does.
fC_LSE<-ggplot(d_LSE)+
  geom_errorbar(aes(x=step_middle, ymin=-y_mean+1-sqrt(y_var), ymax=-y_mean+1+sqrt(y_var),width=0.0005, color="Markov & prior"), alpha=1)+
  geom_segment(aes(x=step_start, xend= step_end, y=-patch_current+1, color="observation & innovation"), linewidth=1, alpha=1)+
  geom_point(aes(x=step_middle,  y=-y_mean+1, color="Markov & prior"), linewidth=1)+
  geom_segment(aes(x=(step_start+step_end)/2, xend=(step_start+step_end)/2, y=-y_mean+1, yend=-patch_current+1, color="observation & innovation"),
               arrow=pred_arrow, linewidth=0.5, na.rm=TRUE)+
  scale_fill_manual(values = SEM, guide = "none") +
  sem_scale + ylab("predicted current (pA)") + common_theme+x_only_ticks + guides(colour="none") + xcommon

fB_ILSE<-ggplot(d_ILSE)+
  geom_point(aes(x=step_start, y=P_mean_t20_y1, color="Markov & prior"), linewidth = 1, alpha=1)+
geom_curve(aes(x=lag(step_start), xend = step_start, y=lag(P_mean_t20_y1),
   yend = P_mean_t20_y1, color="Markov & prior"),curvature = MK_CURV,  angle = 90, ncp=10,linetype = 1, linewidth = 0.5, arrow=markov_arrow)+
  sem_scale + ylab("P(open)") + common_theme+x_only_ticks + guides(colour="none") + xcommon

# fC_ILSE = fC_INR's interval grammar (rect over the window, mean drawn as a segment across it) on the
# CONSTANT least-squares variance, so every rectangle has the same height and only its span says that
# the prediction is an interval average.
fC_ILSE<-ggplot(d_ILSE)+
  geom_rect(aes(xmin=step_start, xmax=step_end, ymin=-y_mean+1-sqrt(y_var), ymax=-y_mean+1+sqrt(y_var), fill="Markov & prior"), alpha=0.3)+
  geom_segment(aes(x=step_start, xend= step_end, y=-patch_current+1, color="observation & innovation"), linewidth=1, alpha=1)+
  geom_segment(aes(x=step_start, xend= step_end, y=-y_mean+1, color="Markov & prior"), linewidth=1)+
  geom_segment(aes(x=(step_start+step_end)/2, xend=(step_start+step_end)/2, y=-y_mean+1, yend=-patch_current+1, color="observation & innovation"),
               arrow=pred_arrow, linewidth=0.5, na.rm=TRUE)+
  scale_fill_manual(values = SEM, guide = "none") +
  sem_scale + ylab("predicted current (pA)") + common_theme+x_only_ticks + guides(colour="none") + xcommon



## -----------------------------------------------------------------------------
fE_NR<-ggplot(d_NR)+
 # geom_segment(aes(x=step_start, xend=step_end, y=logL))+
  geom_line(aes(x=step_end,  y=cumsum(logL), color="logLikelihood"))+
  geom_point(aes(x=step_end,  y=cumsum(logL), color="logLikelihood"))+
  

  sem_scale + ylab("logL") + common_theme+guides(colour="none") + xcommon+xlab("time (ms)")
fE_NR



## -----------------------------------------------------------------------------
fB_INR<-ggplot(d_INR)+
  geom_point(aes(x=step_start, y=lag(P_mean_t2_y0), color="Markov & prior"), linewidth = 1, alpha=1)+
geom_curve(aes(x=lag(step_start), xend = step_start, y=lag(lag(P_mean_t2_y0)),
   yend = lag(P_mean_t2_y0), color="Markov & prior"),curvature = MK_CURV,  angle = 90, ncp=10,linetype = 1, linewidth = 0.5, arrow=markov_arrow)+
  sem_scale + ylab("P(open)") + common_theme+x_only_ticks + guides(colour="none") + xcommon

fB_INR


## -----------------------------------------------------------------------------
fC_INR<-ggplot(d_INR)+
  geom_rect(aes(xmin=step_start, xmax=step_end, ymin=-y_mean+1-sqrt(y_var), ymax=-y_mean+1+sqrt(y_var), fill="Markov & prior"), alpha=0.3)+
  geom_segment(aes(x=step_start, xend= step_end, y=-patch_current+1, color="observation & innovation"), linewidth=1, alpha=1)+
  geom_segment(aes(x=step_start, xend= step_end, y=-y_mean+1, color="Markov & prior"), linewidth=1)+
  geom_segment(aes(x=(step_start+step_end)/2, xend=(step_start+step_end)/2, y=-y_mean+1, yend=-patch_current+1, color="observation & innovation"),
               arrow=pred_arrow, linewidth=0.5, na.rm=TRUE)+
  scale_fill_manual(values = SEM, guide = "none") +
  sem_scale + ylab("predicted current (pA)") + common_theme+x_only_ticks + guides(colour="none") + xcommon
fC_INR



## -----------------------------------------------------------------------------
fE_NMR<-ggplot(d_INR)+
 # geom_segment(aes(x=step_start, xend=step_end, y=logL))+
  geom_line(aes(x=step_end,  y=cumsum(logL), color="logLikelihood"))+
  geom_point(aes(x=step_end,  y=cumsum(logL), color="logLikelihood"))+
  

  sem_scale + ylab("logL") + common_theme+guides(colour="none") + xcommon+xlab("time (ms)")
fE_NR



## ----fB_R---------------------------------------------------------------------
fB_R<-ggplot(d_R)+
  geom_point(aes(x=lag(step_middle), y=lag(P_mean_t15_y1), color="Bayes & posterior"), alpha=ORIG_ALPHA)+
  geom_point(aes(x=lead(step_middle), y=lead(P_mean_t15_y0), color="Markov & prior"),  alpha=1)+
geom_curve(aes(x=lag(step_middle), xend = step_middle, y=lag(P_mean_t15_y1),
   yend = P_mean_t15_y0, color="Markov & prior"),curvature = MK_CURV,  angle = 90, ncp=10,linetype = 1, linewidth = 0.5, arrow=markov_arrow)+
  sem_scale + ylab("P(open)") + common_theme+x_only_ticks + guides(colour="none") + xcommon

fB_R


## -----------------------------------------------------------------------------
fC_R<-ggplot(d_R)+
  geom_errorbar(aes(x=step_middle, ymin=-y_mean+1-sqrt(y_var), ymax=-y_mean+1+sqrt(y_var), width=0.0005,color="Markov & prior"), alpha=1)+
  geom_segment(aes(x=step_start, xend= step_end, y=-patch_current+1, color="observation & innovation"), linewidth=1, alpha=1)+
  geom_point(aes(x=step_middle,  y=-y_mean+1, color="Markov & prior"), linewidth=1)+
  geom_segment(aes(x=(step_start+step_end)/2, xend=(step_start+step_end)/2, y=-y_mean+1, yend=-patch_current+1, color="observation & innovation"),
               arrow=pred_arrow, linewidth=0.5, na.rm=TRUE)+
  scale_fill_manual(values = SEM, guide = "none") +
  sem_scale + ylab("predicted current (pA)") + common_theme+x_only_ticks + guides(colour="none") + xcommon
fC_R



## -----------------------------------------------------------------------------

# Pre-calculate to avoid NA/length errors in the plot call
fD_R <- d_R %>% filter(step_start > 0) %>% ggplot() +
  geom_point(aes(x = step_middle, y = P_mean_t15_y0, color = "Markov & prior"), alpha = ORIG_ALPHA) +
  geom_point(aes(x = step_middle, y = P_mean_t15_y1, color = "Bayes & posterior")) +
  geom_curve(aes(x = step_middle, xend = step_middle, y = P_mean_t15_y0, yend = P_mean_t15_y1,
                 color = "Bayes & posterior"),
             linewidth = 0.5, curvature = -0.9, angle = 180, ncp = 10, arrow = bayes_arrow) +
  sem_scale + common_theme +
  ylab("P(open)") + guides(colour = "none") + xcommon + x_only_ticks
fD_R





## -----------------------------------------------------------------------------
fE_R<-ggplot(d_R)+
 # geom_segment(aes(x=step_start, xend=step_end, y=logL))+
  geom_line(aes(x=step_end,  y=cumsum(logL), color="logLikelihood"))+
  geom_point(aes(x=step_end,  y=cumsum(logL), color="logLikelihood"))+
  

  sem_scale + ylab("logL") + common_theme+guides(colour="none") + xcommon+xlab("time (ms)")
fE_NR



## ----fB_MR--------------------------------------------------------------------
fB_MR<-d_MR%>%filter(sample_index>0)%>%ggplot()+
  geom_point(aes(x=step_start, y=P_mean_t1_y1, color="Bayes & posterior"), alpha=ORIG_ALPHA)+
  geom_point(aes(x=step_end, y=P_mean_t2_y1, color="Markov & prior"),  alpha=1)+
geom_curve(aes(x=step_start, xend = step_end, y=P_mean_t1_y1,
   yend = P_mean_t2_y1, color="Markov & prior"),curvature = MK_CURV,  angle = 90, ncp=10,linetype = 1, linewidth = 0.5, arrow=markov_arrow)+
  sem_scale + ylab("P(open)") + common_theme + x_only_ticks + guides(colour="none") + xcommon

fB_MR


## -----------------------------------------------------------------------------
fC_MR<-ggplot(d_MR)+
  geom_rect(aes(xmin=step_start, xmax=step_end, ymin=-y_mean+1-sqrt(y_var), ymax=-y_mean+1+sqrt(y_var), fill="Markov & prior"), alpha=0.3)+
  geom_segment(aes(x=step_start, xend= step_end, y=-patch_current+1, color="observation & innovation"), linewidth=1, alpha=1)+
  geom_segment(aes(x=step_start, xend= step_end, y=-y_mean+1, color="Markov & prior"), linewidth=1)+
  geom_segment(aes(x=(step_start+step_end)/2, xend=(step_start+step_end)/2, y=-y_mean+1, yend=-patch_current+1, color="observation & innovation"),
               arrow=pred_arrow, linewidth=0.5, na.rm=TRUE)+
  scale_fill_manual(values = SEM, guide = "none") +
  sem_scale + ylab("predicted current (pA)") + common_theme+x_only_ticks + guides(colour="none") + xcommon
fC_MR



## -----------------------------------------------------------------------------

# Pre-calculate to avoid NA/length errors in the plot call
fD_MR <- d_MR %>% filter(step_start > 0) %>% ggplot() +
  geom_point(aes(x=step_start, y=P_mean_t1_y1, color="Bayes & posterior"), alpha=1)+
  geom_point(aes(x=lead(step_end), y=lead(P_mean_t2_y1), color="Markov & prior"),  alpha=ORIG_ALPHA)+
geom_curve(aes(x=step_end, xend = step_end, y=P_mean_t2_y1,
   yend = lead(P_mean_t1_y1), 
                 color = "Bayes & posterior"),
             linewidth = 0.5, curvature = -0.9, angle = 180, ncp = 10, arrow = bayes_arrow) +
  sem_scale + common_theme +
  ylab("P(open)") + guides(colour = "none") + xcommon + x_only_ticks
fD_MR





## ----fB_VR--------------------------------------------------------------------
# VR's three panels are MR's code, verbatim, on d_VR. The visible difference to look for
# is the predictive band in fC_VR (narrower: it is drawn from the residual y_var), and
# then, downstream of it, a posterior and a mean that drift away from MR's because the
# variance divides the gain.
fB_VR<-d_VR%>%filter(sample_index>0)%>%ggplot()+
  geom_point(aes(x=step_start, y=P_mean_t1_y1, color="Bayes & posterior"), alpha=ORIG_ALPHA)+
  geom_point(aes(x=step_end, y=P_mean_t2_y1, color="Markov & prior"),  alpha=1)+
geom_curve(aes(x=step_start, xend = step_end, y=P_mean_t1_y1,
   yend = P_mean_t2_y1, color="Markov & prior"),curvature = MK_CURV,  angle = 90, ncp=10,linetype = 1, linewidth = 0.5, arrow=markov_arrow)+
  sem_scale + ylab("P(open)") + common_theme + x_only_ticks + guides(colour="none") + xcommon

fB_VR


## ----fC_VR--------------------------------------------------------------------
fC_VR<-ggplot(d_VR)+
  geom_rect(aes(xmin=step_start, xmax=step_end, ymin=-y_mean+1-sqrt(y_var), ymax=-y_mean+1+sqrt(y_var), fill="Markov & prior"), alpha=0.3)+
  geom_segment(aes(x=step_start, xend= step_end, y=-patch_current+1, color="observation & innovation"), linewidth=1, alpha=1)+
  geom_segment(aes(x=step_start, xend= step_end, y=-y_mean+1, color="Markov & prior"), linewidth=1)+
  geom_segment(aes(x=(step_start+step_end)/2, xend=(step_start+step_end)/2, y=-y_mean+1, yend=-patch_current+1, color="observation & innovation"),
               arrow=pred_arrow, linewidth=0.5, na.rm=TRUE)+
  scale_fill_manual(values = SEM, guide = "none") +
  sem_scale + ylab("predicted current (pA)") + common_theme+x_only_ticks + guides(colour="none") + xcommon
fC_VR


## ----fD_VR--------------------------------------------------------------------
fD_VR <- d_VR %>% filter(step_start > 0) %>% ggplot() +
  geom_point(aes(x=step_start, y=P_mean_t1_y1, color="Bayes & posterior"), alpha=1)+
  geom_point(aes(x=lead(step_end), y=lead(P_mean_t2_y1), color="Markov & prior"),  alpha=ORIG_ALPHA)+
geom_curve(aes(x=step_end, xend = step_end, y=P_mean_t2_y1,
   yend = lead(P_mean_t1_y1),
                 color = "Bayes & posterior"),
             linewidth = 0.5, curvature = -0.9, angle = 180, ncp = 10, arrow = bayes_arrow) +
  sem_scale + common_theme +
  ylab("P(open)") + guides(colour = "none") + xcommon + x_only_ticks
fD_VR


## -----------------------------------------------------------------------------
fE_MR<-ggplot(d_MR)+
 # geom_segment(aes(x=step_start, xend=step_end, y=logL))+
  geom_line(aes(x=step_end,  y=cumsum(logL), color="logLikelihood"))+
  geom_point(aes(x=step_end,  y=cumsum(logL), color="logLikelihood"))+
  

  sem_scale + ylab("logL") + common_theme+guides(colour="none") + xcommon+xlab("time (ms)")
fE_MR



## ----fig.width = 7.0, fig.height=9--------------------------------------------
# channel (truth) added to the observation row; ONE panel carries the shared legend
chan <- geom_line(data = d_s, aes(x = step_middle, y = -patch_current, color = "channel current"), linewidth = 0.5)
# realised open fraction = channel current / N_ch: the ground truth the P(open) belief is tracking
chan_P <- geom_line(data = d_s, aes(x = step_middle, y = -patch_current / Nch, color = "channel current"), linewidth = 0.5)
# legend keys (NA coords -> invisible, legend only; logLikelihood excluded via the scale breaks):
# the three PHASES are operations -> arrow keys; channel current is the truth -> a plain line key.
.ops <- c("Markov & prior", "observation & innovation", "Bayes & posterior")
legend_arrow <- geom_segment(data = data.frame(role = .ops, x = NA_real_, xend = NA_real_, y = NA_real_, yend = NA_real_),
                             aes(x = x, xend = xend, y = y, yend = yend, colour = role), na.rm = TRUE,
                             arrow = arrow(length = unit(0.05, "in"), type = "closed"))
legend_line  <- geom_line(data = data.frame(role = "channel current", x = NA_real_, y = NA_real_),
                          aes(x = x, y = y, colour = role), na.rm = TRUE)
# algorithm names spelled out (defined at first use), as a subtitle under each abbreviation
# VR's spelled-out name follows the roster's own logic: the prefix letter says WHAT the
# conditioning is about (Mean, Variance, Interval), so V = Variance. It is one string,
# changed here and nowhere else if the paper settles on another wording.
# Keyed by the PRINTED abbreviation, which is what colhead() receives through .DISP.
# The interval prefix reads the same way throughout: I means the quantity is conditioned on the
# acquisition interval rather than on an instant, so ILSE is least squares on the interval-averaged
# mean and INR is the open-loop gating likelihood on the same average. MNR was this member's printed
# name until 2026-08-01; the data key, the file token and the header now all read INR.
FULL <- c(LSE = "Least Squares", ILSE = "Interval Least Squares",
          NR = "Non-Recursive", INR = "Interval Non-Recursive", R = "Recursive",
          MR = "Mean Recursive", VR = "Variance Recursive", IR = "Interval Recursive")
# column headers live in their OWN strip row (colhead) so the top-left panel is free to carry the
# row's filter-phase title, facet_grid style (algorithms on top, phases on the left of each row).
colhead <- function(ab) ggplot() + theme_void() + ggtitle(ab, subtitle = FULL[[ab]]) +
  theme(plot.title    = element_text(hjust = 0.5, size = 9),
        plot.subtitle = element_text(hjust = 0.5, size = 7))
# (the legend is attached inside build_figure, to whichever column comes first)
scP  <- scale_y_continuous(limits = YP)     # prior & posterior (full-data range, window-independent)
scI  <- scale_y_continuous(limits = YI)     # observation current (full-data range, window-independent)

# CAPTION NOTE: in the observation row the predictive spread is a BAND for gmean_i / averaged-
# conductance algorithms (NMR, MR, VR, IR) and an ERRORBAR for plain-g algorithms (NR, R) — the
# band-vs-errorbar distinction encodes the averaging (av) axis, on purpose.

# CAPTION NOTE: in the Markov (prior) and Bayes (posterior) panels of the recursive algorithms
# (R, MR, VR, IR) each arrow starts from a faded mark, the state carried in from the previous step, and
# points to the full-colour state that the panel computes; the non-recursive algorithms (NR, NMR) have
# no such carried-in state.

# --- build ONE figure for a given time window and save it. -----------------------------
# COLUMNS ARE A PARAMETER (`cols`): this file is the single source of the figure, and the two
# notebooks that use it differ only in the roster they pass.
#   figure_1.Rmd      -> c("R","MR","VR","IR")            the paper figure, the recursive ladder
#   figure_1_all.Rmd  -> c("NR","MNR","R","MR","VR","IR")  every algorithm, for exploration
# Nothing below hardcodes a column name: the FIRST column of `cols` carries the y quantity and
# the legend, the middle one carries the x title, and the "no update (open loop)" placeholder is
# used for whichever of NR/MNR are present.
#
# sel = NULL is the whole recording; sel = c(i, j) crops the TIME axis to those windows (full data
# kept, so every lag/lead arrow is exactly the full-graph arrow, just zoomed). logL is the exception:
# it re-accumulates from 0 at the window start, and when zooming all columns share ONE logL axis
# (unified). The full recording keeps two logL blocks (the non-recursive pair crashes to ~-95) when
# any non-recursive column is present. y-NUMBERS live only on the first column (+ stage-label title);
# every other column keeps tick MARKS only. x-NUMBERS live on the bottom row (D). ---

# panel registries, so a roster is just a vector of names
.DAT <- list(LSE = d_LSE, ILSE = d_ILSE, NR = d_NR, INR = d_INR, R = d_R, MR = d_MR, VR = d_VR, IR = d_IR)
.FB  <- list(LSE = fB_LSE, ILSE = fB_ILSE, NR = fB_NR, INR = fB_INR, R = fB_R, MR = fB_MR, VR = fB_VR, IR = fB_IR)
.FC  <- list(LSE = fC_LSE, ILSE = fC_ILSE, NR = fC_NR, INR = fC_INR, R = fC_R, MR = fC_MR, VR = fC_VR, IR = fC_IR)
.FD  <- list(R = fD_R, MR = fD_MR, VR = fD_VR, IR = fD_IR)     # non-recursive members (LSE, NR, NMR) have no posterior update
.NAIVE <- c("LSE", "ILSE", "NR", "INR")
# A roster is written in DATA keys. The mean-non-recursive algorithm is a pre-existing wart: its data
# key is NMR, its file token and its printed header are MNR (see the filenames/algorithms vectors and
# the FULL map). Keep the two apart here rather than propagating the wart into the roster.
.DISP <- c(LSE = "LSE", ILSE = "ILSE", NR = "NR", INR = "INR", R = "R", MR = "MR", VR = "VR", IR = "IR")

# hgt (2026-08-05): the body figure and its supplement need different heights. At 7.5 in the
# Figure 1 float was 102.5 pt taller than a page could hold with its 287-word legend, and LaTeX
# dropped the last 31 words plus the supplement line off the bottom. 6.05 in puts the art at
# 435 pt, the height Figure 2 already fits at. Font sizes are in points and do not scale with the
# device, so the panels get shorter and the text stays where eLife's floor needs it.
build_figure <- function(sel, outfile, cols, hgt = 7.5, wdt = 7.0) {
  stopifnot(length(cols) >= 2, all(cols %in% names(.DAT)))
  unify <- !is.null(sel)
  first <- cols[1]                                  # carries the y title and the legend
  xcol  <- cols[ceiling(length(cols) / 2)]          # carries the single x title
  win <- function(df) if (is.null(sel)) df else dplyr::filter(df, sample_index %in% sel)
  dW  <- lapply(.DAT[cols], win)

  # time crop
  if (is.null(sel)) {
    xlim <- range(c(0, d_s$step_middle, d_IR$step_start, d_IR$step_end), na.rm = TRUE)
  } else {
    xlim <- range(c(win(d_s)$step_start, win(d_s)$step_end), na.rm = TRUE)
  }

  # X AXIS = THE ACQUISITION INTERVAL, NOT THE CLOCK (2026-08-12, Luciano). The figure's subject is
  # what one interval contains and what the recording keeps of it, so the ticks are the interval
  # BOUNDARIES and the label is the index of the interval that starts there (the last tick carries
  # the next index, since it is the boundary that closes the crop). The dotted verticals are the
  # same boundaries drawn inside every panel: without them the reader has to infer where a window
  # begins from the marks themselves, which is what the figure is trying to teach. Seconds are gone
  # from the axis; the interval duration is a caption fact.
  bnd <- win(d_s) %>% dplyr::group_by(sample_index) %>%
    dplyr::summarise(s = min(step_start), e = max(step_end), .groups = "drop") %>%
    dplyr::arrange(sample_index)
  xbnd <- c(bnd$s, max(bnd$e))                       # boundaries: the dotted verticals
  xbrk <- (bnd$s + bnd$e) / 2                        # labels sit at the MIDDLE of the interval they
  xlbl <- bnd$sample_index                           # name, not at its edge, so two intervals read
  vgrid <- geom_vline(xintercept = xbnd, linetype = "dotted",   # as two rather than as three ticks
                      linewidth = 0.25, colour = "grey70")

  # AUTO y-ranges from the (windowed) data, SHARED across the roster per row. The full recording
  # keeps the data range (which reaches 0); a zoom uses the window's own range (non-zero) so the
  # arrows are large, padded a touch so arrowheads at the extremes are not clipped.
  pad <- function(r, lo = 0.06, hi = 0.06) r + c(-lo, hi) * diff(r)
  yP <- range(unlist(lapply(dW, .dP)), na.rm = TRUE)
  yI <- range(c(-win(d_s)$patch_current, unlist(lapply(dW, .dI))), na.rm = TRUE)
  # row A's bottom pad is larger than its top: the Markov arc bows DOWNWARD (MK_CURV) and the
  # open-loop columns' priors are the lowest points in the row, so a symmetric pad clips the arc and
  # the arrowhead with it. The pad and MK_CURV are one fix in two places; change them together.
  if (!is.null(sel)) { yP <- pad(yP, lo = 0.20); yI <- pad(yI) }

  # IR's two state rows are rebuilt HERE and nowhere else: their tie lives in physical space, so it
  # needs this call's crop AND the panel box. The box is the figure minus the axis furniture, and it
  # is only an ESTIMATE (measured off a 7.0 x 6.05 in render with six columns: panels came out about
  # 0.9 x 0.76 in). Nothing breaks if it is off. It sets the proportion of the tie's bow and nothing
  # else: the arrow anchors are exact whatever the estimate, since the midpoint is computed through
  # the same map that draws the arc, and the disc is drawn in absolute units.
  FB <- .FB; FD <- .FD
  if ("IR" %in% cols) {
    gIR <- ir_gmap(xlim, yP, (wdt - 1.6) / length(cols), (hgt - 3.0) / 4)
    FB[["IR"]] <- ir_prior_panel(gIR); FD[["IR"]] <- ir_post_panel(gIR)
  }

  # windowed cumulative logL, re-anchored to 0 at the window start
  lw <- function(df) {
    d <- dplyr::arrange(win(df), step_end)
    tibble::tibble(t = c(min(d$step_start), d$step_end), cL = c(0, cumsum(d$logL)))
  }
  fE <- function(df) ggplot(lw(df), aes(t, cL)) +
    geom_line(aes(colour = "logLikelihood")) + geom_point(aes(colour = "logLikelihood")) +
    sem_scale + common_theme + guides(colour = "none")
  naive <- intersect(cols, .NAIVE)
  recg  <- setdiff(cols, .NAIVE)
  if (unify || !length(naive)) {          # one logL axis: zoomed, or no non-recursive column present
    yLn <- yLr <- range(unlist(lapply(.DAT[cols], function(x) lw(x)$cL)), na.rm = TRUE)
  } else {                                # full recording WITH the naive pair: two blocks
    yLn <- range(unlist(lapply(.DAT[naive], function(x) lw(x)$cL)), na.rm = TRUE)
    yLr <- range(unlist(lapply(.DAT[recg],  function(x) lw(x)$cL)), na.rm = TRUE)
  }

  noUpd <- function() ggplot() +          # empty open-loop panel, text centred in THIS window
    annotate("text", x = mean(xlim), y = mean(yP), label = "no update\n(open loop)", size = 2.5, colour = "grey60") +
    common_theme

  # per-cell placement. ONE coord per cell crops BOTH axes (x = window; y = the row's shared auto
  # range; coord CLIPS instead of dropping, so arrows at the edge survive).
  cell <- function(p, ylim, col, bottom = FALSE, tag = NULL, yl = NULL, title = NULL) {
    p$layers <- c(list(vgrid), p$layers)          # PREPENDED, so the boundaries sit under the marks
    p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
    if (col == first) p <- p + ylab(yl) + theme(axis.title.y = element_text(size = 8))
    else if (length(naive) && col == recg[1] && bottom && !unify) p <- p + y_no_title  # 2nd logL axis
    else p <- p + y_only_ticks
    if (bottom) {
      if (col == xcol) p <- p + xlab("acquisition interval") else p <- p + theme(axis.title.x = element_blank())
    } else p <- p + x_only_ticks
    if (!is.null(tag)) p <- p + labs(tag = tag)
    if (!is.null(title)) p <- p + ggtitle(title) +
      theme(plot.title = element_text(hjust = 0, size = 8.5, face = "plain"))
    p
  }

  # the legend rides on the first column's B panel, whichever algorithm that is
  legendise <- function(p) p + legend_arrow + legend_line +
    guides(colour = guide_legend(override.aes = list(linetype = 1, shape = NA, linewidth = 1.0)))

  rowA <- lapply(cols, function(cl) {
    p <- FB[[cl]] + chan_P
    if (cl == first) cell(legendise(p), yP, cl, tag = "A", yl = expression(P[open]), title = "Markov & prior")
    else cell(p, yP, cl)
  })
  rowB <- lapply(cols, function(cl) {
    p <- .FC[[cl]] + chan
    if (cl == first) cell(p, yI, cl, tag = "B", yl = expression("current (pA)"), title = "prediction & innovation")
    else cell(p, yI, cl)
  })
  rowC <- lapply(cols, function(cl) {
    p <- if (cl %in% .NAIVE) noUpd() else FD[[cl]] + chan_P
    if (cl == first) cell(p, yP, cl, tag = "C", yl = expression(P[open]), title = "Bayes & posterior")
    else cell(p, yP, cl)
  })
  rowD <- lapply(cols, function(cl) {
    yl_ <- if (cl %in% naive) yLn else yLr
    if (cl == first) cell(fE(.DAT[[cl]]), yl_, cl, bottom = TRUE, tag = "D",
                          yl = expression(Sigma~"logL"), title = "logLikelihood")
    else cell(fE(.DAT[[cl]]), yl_, cl, bottom = TRUE)
  })

  g <- do.call(wrap_plots, c(
    lapply(cols, function(cl) colhead(.DISP[[cl]])), rowA, rowB, rowC, rowD,
    list(ncol = length(cols), heights = c(0.01, 1, 1, 1, 1), guides = "collect"))) &
    theme(legend.position = "bottom",
          legend.text = element_text(size = 8),
          legend.key.width = unit(1, "lines"),
          legend.key.spacing.x = unit(3, "pt"),
          plot.margin = margin(1.5, 2, 1.5, 2)) &
    scale_x_continuous(breaks = xbrk, labels = xlbl) &                 # one label per interval
    theme(axis.ticks.x = element_blank())   # the dotted boundaries carry the x structure, not ticks

  ggsave(outfile, g, width = wdt, height = hgt)
  g
}
