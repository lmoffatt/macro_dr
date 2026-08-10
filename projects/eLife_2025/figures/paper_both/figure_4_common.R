## figure_4_common.R — shared machinery for the TWO figure-4 halves.
source("noise_units.R")  # the one definition of S_tilde and of when to convert
##
## Figure 4 was split in two (Luciano, 2026-07-23): one figure for the FIRST moment (bias) and one
## for the SECOND (information distortion), BOTH carrying the ellipse block and BOTH on the merged
## paper's four-column roster LSE / NR / R / IR. With four algorithms a single figure carrying both
## moments would run to sixteen map rows; split, each half is two blocks (ellipses + its own maps).
##
## Everything below the roster is shared and lives here; the two notebooks
##   figure_4_bias.Rmd        -> ellipses + distortion-induced bias   (battery_sim_G,  theta_sim)
##   figure_4_distortion.Rmd  -> ellipses + information distortion    (battery_pool_G, theta_pool)
## differ only in which map block they assemble. Edit the machinery here, not there, or they drift.
##
## TRAPS carried over from the single-figure version, none of them re-derivable from the plots:
##  - the "noise" label is NOT the physical noise: current_noise = label/1000, and the CSV column
##    noise_in_conductance_tau stores the LABEL. The y axis is the label; the caption must say so.
##  - the y axis is RAW instrumental noise, never noise/N_ch. Gating variance grows with N_ch, so
##    dividing flattens by construction the very staircase the figure can measure.
##  - group_size 10 at few channels does not merely add noise, it REVERSES the sign of the N_ch
##    direction. ELL_GROUP = 100 is the floor. The MAPS fit no estimates and are untouched.

library(tidyverse)
library(patchwork)

# ---- roster ---------------------------------------------------------------------------------
# battery-fed notebooks, so the keys are the PREFIXED form that lands in the CSV `algorithm` column.
# Order = row order inside each parameter block: the cost ladder, cheapest first.
# OPT-IN OVERRIDE: a notebook that wants a different roster or a fourth data dir declares
# FIG4_ALGOS / FIG4_ALGO_LAB / FIG4_DATA_DIRS BEFORE sourcing this file; the defaults below are
# unchanged, so every existing notebook renders exactly as it did. Added 2026-07-31 for
# figure_4_5algo.Rmd, which is figure_4.Rmd with INR as a fifth member. The variable stays set for
# the rest of the R session, so knit one notebook per session or the four-member figure inherits it.
ALGOS    <- if (exists("FIG4_ALGOS")) FIG4_ALGOS else
  c("nonlinearsqr", "macro_NR", "macro_R", "macro_IR")
# `nonlinearsqr` is averaging_approximation = 1, i.e. least squares on the INTERVAL-AVERAGED mean,
# so on the instantaneous/interval axis its name is ILSE, and LSE is the av = 0 arm that lands on
# disk as `nonlinearsqr_g`. Relabelled here 2026-08-06, following figure_1.Rmd and figure_3.Rmd,
# which made the same correction when the un-averaged runs appeared. The default ROSTER is
# untouched: no supplement grows a column, only the word on its box changes.
ALGO_LAB <- if (exists("FIG4_ALGO_LAB")) FIG4_ALGO_LAB else
  c(nonlinearsqr = "ILSE", macro_NR = "NR", macro_R = "R", macro_IR = "IR")

NCHS      <- c(10, 100, 1000, 10000)
DATA_DIRS <- if (exists("FIG4_DATA_DIRS")) FIG4_DATA_DIRS else
  c("../data/1c2ae6f", "../data/87889e6", "../data/0ffbda7")  # search path, first hit wins

# ---- file naming: the prefix is PER ALGORITHM ------------------------------------------------
# LSE (nonlinearsqr) is produced by dispatch_figure_4.sh's LSE arm, which routes to
# figure_4_LSE.macroir and stamps the prefix `figure_3_LSE_`; macro/micro get `figure_3_G_`
# (dispatch_figure_4.sh:217-218). Everything else about the two file sets is identical, including
# the Gaussian_ component names and the _G battery suffix.
# startsWith, not identical: every family-2 member writes under the LSE prefix, including
# `nonlinearsqr_g` (the averaging=0 member added 2026-08-02).
.pre  <- function(a) if (startsWith(a, "nonlinearsqr")) "figure_3_LSE_" else "figure_3_G_"
findf <- function(nch, algo, z, kind) {
  f <- file.path(DATA_DIRS,
                 sprintf("%snch_%d_nsim_10000_%s_noise_%s_%s.csv", .pre(algo), nch, algo, z, kind))
  hit <- f[file.exists(f)]
  if (!length(hit)) stop("figure_4: no file for ", algo, " N_ch ", nch, " noise ", z, " (", kind, ")")
  hit[1]
}
have_cell <- function(nch, algo, z, kind) {           # non-throwing probe, for the tolerant block A
  any(file.exists(file.path(DATA_DIRS,
    sprintf("%snch_%d_nsim_10000_%s_noise_%s_%s.csv", .pre(algo), nch, algo, z, kind))))
}

# ---- the grid, AUTO-DETECTED from disk (both batteries required) ------------------------------
# Replaces the hand-declared per-algorithm grids, which went stale every time a run landed. At render
# time we SCAN DATA_DIRS and include every (N_ch, noise) cell that has BOTH batteries — battery_pool_G
# for the distortion half AND battery_sim_G for the bias half — at nsim 10000, so a new run appears on
# the next knit with no edit here. Two rails are kept deliberately:
#   - N_ch is fixed to NCHS. The four columns are a LAYOUT choice, not auto-grown: IR alone has runs at
#     5/20/50/200/500/2000/5000, and auto-including them would silently blow the figure to 11 columns.
#     A new N_ch column is a deliberate edit to NCHS; a new NOISE level is automatic.
#   - BOTH batteries required, so a half-landed run (sim_G present, pool_G still queued — e.g. the NR
#     diagonal tail) is skipped, not hard-stopped mid-render by findf().
# COST, stated so it is not a surprise: a run that FAILS to land simply does not appear (there is no
# "expected but missing" error any more). Mitigation: the detected grid is printed loudly below on
# every render, and the data dirs are content-addressed by git hash, so at freeze the figure is still
# reproducible from a known cell set — read which cells from the render log, not from a table here.
autodetect_cells <- function(algo) {
  pat <- sprintf("^%snch_([0-9]+)_nsim_10000_%s_noise_([0-9.eE+]+)_battery_pool_G\\.csv$",
                 .pre(algo), algo)
  files <- basename(list.files(DATA_DIRS, pattern = pat))
  if (!length(files)) return(tibble::tibble(nch = integer(), z = character()))
  m  <- regmatches(files, regexec(pat, files))
  df <- tibble::tibble(nch = as.integer(vapply(m, `[`, character(1), 2)),
                       z   = vapply(m, `[`, character(1), 3)) %>%
        dplyr::distinct() %>% dplyr::filter(nch %in% NCHS)
  # keep only cells whose sim_G ALSO exists (the bias half shares this grid and findf() hard-stops)
  df[vapply(seq_len(nrow(df)),
            function(i) have_cell(df$nch[i], algo, df$z[i], "battery_sim_G"), logical(1)), , drop = FALSE]
}
GRID_BY_ALGO <- lapply(ALGOS, autodetect_cells)
names(GRID_BY_ALGO) <- ALGOS
GRID <- distinct(bind_rows(GRID_BY_ALGO))              # union, for the shared y-breaks
noise_span <- function(a) {
  z <- GRID_BY_ALGO[[a]]$z
  if (!length(z)) 0 else diff(range(log10(as.numeric(z))))
}

PIDX  <- c(on = 0, off = 1, unitary_current = 2, Current_Noise = 3, Num_ch_mean = 5)
COMPS <- c(total  = "Probit_statistics_Likelihood_Gaussian_Information_Distortion",
           sample = "Probit_statistics_Gaussian_Sample_Distortion",
           corr   = "Probit_statistics_Likelihood_Correlation_Distortion")

# the two parameter sets. k_off inflates, N_ch deflates: OPPOSITE signs, so any scalar averaged over
# parameters cancels them and reports a faithful algorithm.
PARAM_LONG <- c("off", "Num_ch_mean")
PARAM_MATH <- c(off = "log[10]~k[off]", Num_ch_mean = "log[10]~N[ch]")

# ---- colour scales ---------------------------------------------------------------------------
RAMP <- c("#021024","#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7",
          "#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F","#3d0011")
bandcols <- function(breaks) {
  n <- length(breaks) - 1; ctr <- findInterval(0, breaks); nb <- ctr - 1; nr <- n - ctr; K <- max(nb, nr, 1)
  blue <- colorRampPalette(rev(RAMP[1:7]), space = "Lab")(K + 1)
  red  <- colorRampPalette(RAMP[7:13],     space = "Lab")(K + 1)
  cols <- rep(RAMP[7], n)
  if (nb > 0) cols[seq(ctr - 1, 1)] <- blue[2:(nb + 1)]
  if (nr > 0) cols[seq(ctr + 1, n)] <- red[2:(nr + 1)]
  cols
}
fmt_big <- function(x) formatC(x, format = "g", digits = 3)
# EXTENDED 2026-07-23 for the four-column roster. The old vector topped out at 2.5, which is fine for
# the recursive family (IR 0.97-1.49, R 0.44-1.82) but NOT for NR, whose distortion runs to 1.5e3 in
# k_off and 3.4e4 in N_ch and GROWS with N_ch (median 4.6 / 11.9 / 26.3 / 74.5 across the four
# decades). On the old scale NR saturated in every cell and, where kappa also passed KAPPA_MAX, went
# grey — the mask reading as "unidentified" what is really "off the chart". The fine structure near 1
# is UNCHANGED (it is what resolves IR against R); only a geometric tail is appended, plus 0.4 at the
# bottom for R's minimum of 0.439. Same design as the archived log-det scale in
# figure_5_master.Rmd:76 ("fine near white, coarsening outward, so IR, R and NR each resolve").
# COST, deliberate: the red ramp now spans 13 bands instead of 8, so cells near 1 lose some
# saturation — R's 1.4 reads lighter than it used to. That is the price of holding NR on the same
# scale, and it is the honest one: 1.4 IS mild next to 74. The tail is kept COARSE (2, 5, 10, 30,
# 100, 1000, 1e5 — one band per half-decade to decade) precisely to limit that dilution; a finer
# tail buys nothing, since what NR has to show is a trend across decades, not a value.
DIST_LIN    <- c(0.4, 0.5, 0.67, 0.8, 0.87, 0.91, 0.94, 0.97, 0.99, 1.01, 1.03, 1.06,
                 1.1, 1.15, 1.25, 1.5, 2, 5, 10, 30, 100, 1000, 1e5)
DIST_BREAKS <- log10(DIST_LIN)
# The two validity criteria, DECLARED ONCE and fed to BOTH the colour bar (as tick marks) and the
# maps (as iso-distortion contours), so the bar and the panels cannot drift apart. Log-symmetric
# pairs: solid = a factor 1.15 either way, dashed = a factor 2 either way. Widened from 1.5 to 2
# (2026-07-23) when the scale grew a tail to 1e5: two marks crowded between 0.67 and 1.5 left the
# whole tail unreferenced. Each value MUST be a member of DIST_LIN — colorbar() places its ticks by
# match(v, lin), so a non-member vanishes silently.
DIST_SOL    <- c(0.87, 1.15)
DIST_DSH    <- c(0.5, 2)
DIST_PAL    <- bandcols(DIST_BREAKS)
DIST_LAB    <- fmt_big(DIST_LIN)
# log-spaced bands: IR's values (0.003-0.042) and R's (0.13-0.15) differ by more than a decade, and a
# linear scale would put every IR cell in the neutral band.
# EXTENDED 2026-07-23, same reason as DIST_LIN and measured the same way. The old vector stopped at
# +-0.3 and the render's own saturation report showed 3 cells below and 21 above it. The overflow is
# ENTIRELY in the N_ch direction (NR reaches +12.4, R reaches -2.79); in k_off nobody exceeds 0.055,
# so the primary k_off render used only the three innermost bands of a scale built for a range it
# never visits. Fine structure near 0 UNCHANGED (that is what resolves IR, max |0.053|); a geometric
# tail is appended on BOTH sides, symmetric, so the centre band stays white at 0.
BIAS_LIN    <- c(-15, -7, -3, -1.5, -0.7, -0.3, -0.15, -0.07, -0.03, -0.01, -0.003,
                 0.003, 0.01, 0.03, 0.07, 0.15, 0.3, 0.7, 1.5, 3, 7, 15)
BIAS_PAL    <- bandcols(BIAS_LIN)
BIAS_LAB    <- formatC(BIAS_LIN, format = "g", digits = 2)
# Same two criteria as the distortion half, expressed in the other scale: a bias of 0.07 in log10 IS
# a factor 1.175 on the parameter (~DIST_SOL's 1.15) and 0.3 IS a factor 2 (DIST_DSH). So the two
# halves of Figure 4 declare the SAME thresholds and a reader can carry one reading to the other.
# Symmetric about 0, the bias's natural centre. Members of BIAS_LIN, as colorbar() requires.
BIAS_SOL    <- c(-0.07, 0.07)
BIAS_DSH    <- c(-0.3, 0.3)

# ---- read the declared cells -------------------------------------------------------------------
cells <- imap_dfr(GRID_BY_ALGO, function(g, a)
  if (!nrow(g)) tibble(algo = character(), nch = numeric(), z = character())
  else tibble(algo = a, nch = g$nch, z = g$z))
cat("figure_4:", nrow(cells), "cells over", length(unique(cells$algo)), "algorithms with data\n")
cat("  noise levels per algorithm x N_ch (the staircase; a ragged row here is the point):\n")
print(as.data.frame(cells %>% group_by(algo, nch) %>%
  summarise(noise = paste(sort(as.numeric(z)), collapse = " "), .groups = "drop")), row.names = FALSE)
.absent <- setdiff(ALGOS, unique(cells$algo))
if (length(.absent))
  cat("  NO CELLS DECLARED (map rows will not appear):", paste(ALGO_LAB[.absent], collapse = ", "),
      "\n  -> add rows to its grid in figure_4_common.R as they land\n")

readcell <- function(kind) pmap_dfr(cells, function(algo, nch, z)
  read.csv(findf(nch, algo, z, kind), skip = 1, stringsAsFactors = FALSE) %>%
    mutate(noise = z, algo = algo))

# ---- the source data ---------------------------------------------------------------------------
# One build for the whole set, written as CSV so it doubles as the figure's source data. See
# figure_4_data.R for what makes it stale and why it is not per notebook.
source("figure_4_data.R")
.fig4_src <- fig4_source_data()

# ---- the anchor, filtered ONCE at the source ---------------------------------------------------
# The half-B products (stat / sedat / scal / spec) now carry BOTH anchors, tagged by `anchor`. Four
# supplements read `.fig4_src$scal` and `.fig4_src$sedat` DIRECTLY, so filtering here rather than at
# each consumer is what keeps them from silently doubling their rows. With ANCHOR = "pool", the
# default, every notebook sees exactly the frame it saw before the column existed.
# A notebook that wants the sim-anchored plane declares FIG4_ANCHOR <- "sim" before sourcing; one
# that wants BOTH at once (the pool-vs-sim ratio) reads `.fig4_all`, which is left untouched.
# The bias product has no anchor column: it is sim-only by construction (see figure_4_data.R).
ANCHOR <- if (exists("FIG4_ANCHOR")) FIG4_ANCHOR else "pool"
stopifnot(ANCHOR %in% c("sim", "pool"))
.fig4_all <- .fig4_src
for (.k in c("stat", "sedat", "scal", "spec"))
  .fig4_src[[.k]] <- dplyr::filter(.fig4_src[[.k]], anchor == ANCHOR)
cat("figure_4 data: anchor = ", ANCHOR, " (theta_",
    if (ANCHOR == "sim") "sim, the simulation truth" else "pool, the joint fit", ")\n", sep = "")

.withnchf <- function(d) dplyr::mutate(
  d, nchf = factor(paste0("N[ch]==", Num_ch), levels = paste0("N[ch]==", NCHS)))

# BIAS is anchored at theta_sim, DISTORTION at theta_pool; both are shaped in figure_4_data.R and
# only the roster filter happens here, so two notebooks on different rosters share one build.
bias <- .withnchf(dplyr::filter(.fig4_src$bias, algo %in% ALGOS))
stat <- .withnchf(dplyr::filter(.fig4_src$stat, algo %in% ALGOS))
# UNIDENTIFIED cells: greyed, not clamped. Off-scale AND ill-conditioned (kappa > 3e4) together
# isolate the R corner (N_ch 10000, Delta = 1), a displacement along a near-null direction whose
# magnitude is not a meaningful bias/distortion. kappa alone grows with N_ch for everyone.
kap <- .fig4_src$scal %>% dplyr::filter(param == "kappa") %>%
  dplyr::transmute(algo, Num_ch, noise, interval_in_tau, kappa = v)
KAPPA_MAX <- 3e4
bias <- bias %>% left_join(kap, by = c("algo", "Num_ch", "noise", "interval_in_tau")) %>%
  mutate(unident = (m < min(BIAS_LIN) | m > max(BIAS_LIN)) & is.finite(kappa) & kappa > KAPPA_MAX)
stat <- stat %>% left_join(kap, by = c("algo", "Num_ch", "noise", "interval_in_tau")) %>%
  mutate(unident = (m < min(DIST_LIN) | m > max(DIST_LIN)) & is.finite(kappa) & kappa > KAPPA_MAX)
cat("greyed (unidentified): bias", sum(bias$unident), "cells, distortion",
    sum(stat$unident[stat$comp == "total"]), "cells\n")

# NO SILENT CAPS. Both maps clamp to their scale, so report what the clamp swallowed.
sat <- function(v, lo, hi, what) {
  n <- sum(is.finite(v)); b <- sum(v < lo, na.rm = TRUE); t <- sum(v > hi, na.rm = TRUE)
  cat(sprintf("  %-28s %d/%d below %.3g, %d/%d above %.3g\n", what, b, n, lo, t, n, hi))
}
cat("scale saturation (0 is the only acceptable answer for a final render):\n")
sat(log10(stat$Dconf[stat$comp == "total"]), min(DIST_BREAKS), max(DIST_BREAKS), "distortion")
sat(bias$Bconf, min(BIAS_LIN), max(BIAS_LIN), "bias")

# ---- block A: the ellipses --------------------------------------------------------------------
# THREE design points, and the middle one is a control. At the failing point two things act at once:
# the algorithm's distortion and the estimator's own finite-sample departure from normality. Point 2
# holds group size, N_ch and noise FIXED and moves only the interval, so the difference between 1 and
# 2 is the algorithm alone. Point 3 is a calibrated design point, for scale.
ELL_POINTS <- tibble::tribble(
  ~lab, ~nch,   ~noise, ~dlt,  ~what,
  "1",  10,     "0.1",  0.01,  "distorted",
  "2",  10,     "0.1",  1,     "biased",
  "3",  10000,  "0.1",  0.1,   "faithful")
ELL_GROUP <- 100
PAIR <- c(off = "off", nch = "Num_ch_mean")   # the two opposite-signed directions

ell <- function(S, cx, cy, n = 200) {         # 95% ellipse, 2 dof
  if (any(!is.finite(S))) return(NULL)
  e <- eigen(S, symmetric = TRUE); if (any(e$values <= 0)) return(NULL)
  th <- seq(0, 2 * pi, length.out = n)
  P <- e$vectors %*% diag(sqrt(e$values * qchisq(0.95, 2))) %*% rbind(cos(th), sin(th))
  tibble(x = cx + P[1, ], y = cy + P[2, ])
}

panel <- function(i, algo) {
  pt   <- ELL_POINTS[i, ]
  fcl  <- findf(pt$nch, algo, pt$noise, "mle_cloud_runs")
  fbat <- findf(pt$nch, algo, pt$noise, "battery_pool_G")
  fsim <- findf(pt$nch, algo, pt$noise, "battery_sim_G")

  pts <- read.csv(fcl, skip = 1) %>%
    filter(component_path == "Model_Parameters_Hat", statistic == "value",
           group_size == ELL_GROUP, param_name %in% PAIR,
           abs(interval_in_tau - pt$dlt) < 1e-9) %>%
    select(sample_index, param_name, value) %>%
    pivot_wider(names_from = param_name, values_from = value)
  bat <- read.csv(fbat, skip = 1) %>%
    filter(statistic == "value", probit == "mean", abs(interval_in_tau - pt$dlt) < 1e-9,
           param_index %in% PIDX[PAIR], param_col %in% PIDX[PAIR])
  cov2 <- function(cp) {
    d <- bat %>% filter(component_path == cp); m <- matrix(NA_real_, 2, 2); id <- unname(PIDX[PAIR])
    for (a in 1:2) for (b in 1:2) { v <- d$value[d$param_index == id[a] & d$param_col == id[b]]
                                    if (length(v) == 1) m[a, b] <- v }
    m / ELL_GROUP
  }
  cx <- mean(pts$off, na.rm = TRUE); cy <- mean(pts$Num_ch_mean, na.rm = TRUE)
  E  <- cov(cbind(pts$off, pts$Num_ch_mean), use = "complete.obs")
  src <- bind_rows(
    mutate(ell(E, cx, cy), source = "empirical"),
    mutate(ell(cov2("Probit_statistics_Gaussian_Fisher_Covariance"), cx, cy), source = "Fisher"),
    mutate(ell(cov2("Probit_statistics_Gaussian_Distortion_Corrected_Covariance"), cx, cy),
           source = "corrected")) %>%
    mutate(source = factor(source, levels = c("empirical", "Fisher", "corrected")))
  A <- function(S) { e <- eigen(S, symmetric = TRUE)$values
                     if (any(!is.finite(S)) || any(e <= 0)) NA_real_ else sqrt(prod(e)) }
  Fm <- cov2("Probit_statistics_Gaussian_Fisher_Covariance")
  rF <- A(E) / A(Fm)
  rC <- A(E) / A(cov2("Probit_statistics_Gaussian_Distortion_Corrected_Covariance"))
  # PER-AXIS ratios, not one area. An area multiplies an inflated direction by a deflated one and
  # hides both, which is the very cancellation the maps exist to expose.
  rx <- E[1, 1] / Fm[1, 1]      # k_off, the x axis
  ry <- E[2, 2] / Fm[2, 2]      # N_ch,  the y axis
  nf <- nrow(pts)
  shrink <- function(r) { lo <- r * (nf - 1) / qchisq(.975, nf - 1)
                          hi <- r * (nf - 1) / qchisq(.025, nf - 1)
                          if (r > 1) max(1, lo) else min(1, hi) }
  bandcol <- function(r) { sr <- shrink(r)
    if (isTRUE(all.equal(sr, 1))) "grey45"
    else DIST_PAL[max(1, min(length(DIST_PAL),
           findInterval(log10(sr), DIST_BREAKS, rightmost.closed = TRUE)))] }

  # THE ELLIPSE AND THE MAP ARE DIFFERENT OBJECTS, and the panel says so: the ellipse compares
  # Cov(theta_hat) against the sandwich MARGINALISED over the other parameters; the maps compare
  # score against Fisher unmarginalised. The marginal reading is systematically milder.
  mval <- function(p) { v <- stat$Dconf[stat$comp == "total" & stat$algo == algo &
                                        stat$Num_ch == pt$nch & stat$noise == pt$noise &
                                        abs(stat$interval_in_tau - pt$dlt) < 1e-9 & stat$param == p]
                        if (length(v) == 1) v else NA_real_ }
  mx <- mval("off"); my <- mval("Num_ch_mean")

  dib <- read.csv(fsim, skip = 1) %>%
    filter(component_path == "Probit_statistics_Gaussian_Distortion_Induced_Bias",
           statistic == "value", probit == "mean", abs(interval_in_tau - pt$dlt) < 1e-9,
           param_index %in% PIDX[PAIR]) %>%
    select(param_index, value)
  gx <- 2; gy <- log10(pt$nch)     # truth: log10(k_off) = log10(100) = 2
  gv <- function(nm) { v <- dib$value[dib$param_index == PIDX[[PAIR[[nm]]]]]
                       if (length(v) == 1) v else NA_real_ }
  mk <- tibble(x = c(gx, cx, gx + gv("off")), y = c(gy, cy, gy + gv("nch")),
               kind = factor(c("truth", "empirical mean", "predicted (truth + bias)"),
                             levels = c("truth", "empirical mean", "predicted (truth + bias)")))

  hx <- 1.10 * max(abs(c(src$x, pts$off, mk$x) - gx), na.rm = TRUE)
  hy <- 1.34 * max(abs(c(src$y, pts$Num_ch_mean, mk$y) - gy), na.rm = TRUE)

  cat(sprintf(paste0("fig4A %-3s point %s (%s): N_ch %d, noise %s, Delta %.2f, n_fits %d | ",
                     "ellipse emp/Fisher  k_off %.2f  N_ch %.2f | map  k_off %.2f  N_ch %.2f | ",
                     "area emp/Fisher %.2f, emp/corrected %.2f\n"),
              ALGO_LAB[[algo]], pt$lab, pt$what, pt$nch, pt$noise, pt$dlt, nf, rx, ry,
              mx, my, rF, rC))
  list(g = ggplot() +
    geom_point(data = pts, aes(off, Num_ch_mean), colour = "grey70", size = .2, alpha = .4) +
    geom_path(data = src, aes(x, y, colour = source, linetype = source, linewidth = source)) +
    geom_point(data = mk, aes(x, y, shape = kind), size = 1.7, stroke = .55, fill = "white") +
    scale_shape_manual(values = c(truth = 3, `empirical mean` = 16, `predicted (truth + bias)` = 21)) +
    annotate("label", x = gx + hx, y = gy + hy, hjust = 1.03, vjust = 1.15, size = 2.2,
             colour = bandcol(rx), fill = "white", label.size = 0,
             label.padding = unit(0.6, "pt"), parse = TRUE,
             label = sprintf("k[off]*\" %.2f×\"", rx)) +
    annotate("label", x = gx + hx, y = gy + hy, hjust = 1.03, vjust = 2.5, size = 2.2,
             colour = bandcol(ry), fill = "white", label.size = 0,
             label.padding = unit(0.6, "pt"), parse = TRUE,
             label = sprintf("N[ch]*\" %.2f×\"", ry)) +
    scale_colour_manual(values = c(empirical = "grey35", Fisher = "#D55E00", corrected = "#0072B2")) +
    scale_linetype_manual(values = c(empirical = "solid", Fisher = "dashed", corrected = "dotted")) +
    scale_linewidth_manual(values = c(empirical = .5, Fisher = .7, corrected = .85), guide = "none") +
    labs(x = expression(log[10]~k[off]), y = expression(log[10]~N[ch]),
         title = paste0(ALGO_LAB[[algo]], "   ·   ", pt$lab, "  ", pt$what),
         # the subtitle showed `sigma == <label>`: wrong symbol (sigma is the predictive s.d.
         # in the paper) and the raw label, which is 10*S_tilde.  Both fixed 2026-08-06.
         subtitle = bquote(N[ch] == .(pt$nch) ~ "," ~
                           widetilde(S) == .(dimensionless_noise(as.numeric(pt$noise))) ~ "," ~
                           widetilde(Delta) == .(pt$dlt))) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(size = 7.5, hjust = .5, face = "bold", margin = margin(b = 1)),
          plot.subtitle = element_text(size = 6, hjust = .5, margin = margin(b = 2)),
          axis.text = element_text(size = 6), axis.title = element_text(size = 7.5),
          aspect.ratio = 1,
          legend.position = "bottom", legend.title = element_blank(), legend.box = "horizontal",
          legend.text = element_text(size = 6.5), legend.key.size = unit(.6, "lines"),
          legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"),
          plot.margin = margin(2, 4, 2, 4)),
    gx = gx, gy = gy, hx = hx, hy = hy)
}

# TOLERANT roster for block A. The maps drop an algorithm with no declared cells silently (the facet
# level is empty), but block A calls findf() per design point and would HARD STOP. So an algorithm
# whose design-point files are not all on disk is skipped HERE, named out loud — never silently.
ELL_ALGOS <- Filter(function(a)
  all(pmap_lgl(ELL_POINTS[, c("nch", "noise")], function(nch, noise)
        all(vapply(c("mle_cloud_runs", "battery_pool_G", "battery_sim_G"),
                   function(k) have_cell(nch, a, noise, k), logical(1))))), ALGOS)
if (length(ELL_ALGOS) < length(ALGOS))
  cat("block A: SKIPPING (no complete design-point cells): ",
      paste(ALGO_LAB[setdiff(ALGOS, ELL_ALGOS)], collapse = ", "),
      "\n  -> needs mle_cloud_runs + battery_pool_G + battery_sim_G at N_ch 10 and 10000, noise 0.1\n",
      sep = "")

# build all points x algorithms, then give the panels of a column the SAME window (the larger),
# so the clouds of a design point can be compared by eye rather than by reading the axis numbers.
build_gA <- function() {
  if (!length(ELL_ALGOS)) return(patchwork::plot_spacer())
  cs    <- expand_grid(algo = factor(ELL_ALGOS, levels = ELL_ALGOS), i = seq_len(nrow(ELL_POINTS)))
  built <- pmap(cs, function(algo, i) panel(i, as.character(algo)))
  gl <- lapply(seq_len(nrow(cs)), function(k) {
    mate <- which(cs$i == cs$i[k])
    hx   <- max(vapply(built[mate], function(b) b$hx, numeric(1)))
    hy   <- max(vapply(built[mate], function(b) b$hy, numeric(1)))
    b    <- built[[k]]
    b$g + coord_cartesian(xlim = c(b$gx - hx, b$gx + hx), ylim = c(b$gy - hy, b$gy + hy))
  })
  wrap_plots(gl, nrow = length(ELL_ALGOS), byrow = TRUE, guides = "collect") &
    theme(legend.position = "bottom")
}

# ---- the map block ------------------------------------------------------------------------------
X_BRK <- log10(c(.01, .1, 1))
# y breaks come from the GRID, not hardcoded: the noise axis grows when diagonal cells land, and a
# fixed set would silently stop labelling the new ones. Powers of ten, thinned to at most five.
Y_DEC <- sort(unique(as.numeric(GRID$z)))
Y_BRK <- log10(Y_DEC[seq(1, length(Y_DEC), by = max(1, ceiling(length(Y_DEC) / 5)))])
# Ticks and title in DIMENSIONLESS units (2026-08-06). The stored value is the sweep label,
# which is ten times S_tilde; joins and file matching keep using it raw and only the display
# converts, which is the rule in noise_units.R.
Y_LAB <- formatC(dimensionless_noise(10^Y_BRK), format = "g", digits = 3)

# ONE row factor, not nested facets: levels ordered parameter-major, algorithm-minor, so the four
# algorithms of a parameter are always adjacent and the ladder reads top-to-bottom in each block.
# (algorithm, parameter) pairs that must NOT be shown, because the comparison would not be
# like-for-like (Luciano, 2026-07-23). least squares x N_ch: it fits with unitary_current FIXED, so
# the mean alone pins N_ch and the N*i ridge that limits the other three simply does not exist for
# it. Side by side that reads as LSE winning the count, when what happened is that it assumed the
# hard part away. LSE therefore appears in the k_off rows only. Applies to every map built through
# rowfac(), i.e. both figure-4 halves and the corrected-SE supplement.
#
# LSE x unitary_current and LSE x Current_Noise are excluded for a SECOND and harder reason, found
# 2026-08-01: nonlinearsqr indexes its own four parameters, and the indices do NOT mean what PIDX
# says. Read off `param_name` in the files:
#     LSE    0 = on   1 = off   2 = Current_Baseline   3 = Num_ch_mean
#     macro  0 = on   1 = off   2 = unitary_current    3 = Current_Noise
#                     4 = Current_Baseline             5 = Num_ch_mean
# so the shared mapping `names(PIDX)[match(param_index, PIDX)]` labels LSE's baseline as the
# unitary current and its CHANNEL NUMBER as the instrumental noise. Nothing built so far is wrong,
# because every figure to date uses only index 1 (k_off, which agrees) and index 5 (which LSE does
# not have); but the first supplement to put index 2 or 3 on the plane would have shown LSE's N_ch
# in the sigma_noise column without a word. Excluded here rather than remapped, because the
# like-for-like argument reaches the same verdict: LSE has no unitary current at all (it is held
# fixed, which is why the Num_ch_mean exclusion above exists) and no free noise scale.
#
# BOTH least-squares arms are excluded, and for the same reasons: `nonlinearsqr_g` (LSE, av = 0)
# and `nonlinearsqr` (ILSE, av = 1) differ only in the mean model, not in what they parameterise,
# so they share the four-parameter indexing above. Verified 2026-08-06 off `param_name` in
# ../data/a202e03: LSE writes 0 = on, 1 = off, 2 = Current_Baseline, 3 = Num_ch_mean, exactly as
# ILSE does.
EXCLUDE_ROWS <- tibble::tribble(
  ~algo,            ~param,
  "nonlinearsqr",   "Num_ch_mean",
  "nonlinearsqr",   "unitary_current",
  "nonlinearsqr",   "Current_Noise",
  "nonlinearsqr_g", "Num_ch_mean",
  "nonlinearsqr_g", "unitary_current",
  "nonlinearsqr_g", "Current_Noise")

.rowlab <- function(p, a) sprintf("atop(%s, bold(\"%s\"))", PARAM_MATH[p], ALGO_LAB[a])
rowfac <- function(d, params) {
  # parameter outer, algorithm inner, so the cost ladder reads top-to-bottom inside each parameter
  keep <- tidyr::expand_grid(param = params, algo = ALGOS) %>%
    dplyr::anti_join(EXCLUDE_ROWS, by = c("algo", "param"))
  lv <- .rowlab(keep$param, keep$algo)
  d %>% filter(param %in% params) %>%
    dplyr::anti_join(EXCLUDE_ROWS, by = c("algo", "param")) %>%
    mutate(prow = factor(.rowlab(param, algo), levels = lv))
}

marks <- function(rows) ELL_POINTS %>%
  mutate(lx = log10(dlt), ly = log10(as.numeric(noise)),
         nchf = factor(paste0("N[ch]==", nch), levels = paste0("N[ch]==", NCHS))) %>%
  tidyr::crossing(prow = factor(rows, levels = rows))

# iso_sol / iso_dsh: contour levels drawn ON the panels, in the SAME space as brk (the caller
# transforms, exactly as it already does for brk). Default empty, so the bias map is unchanged.
mapblock <- function(d, zvar, pal, brk, params, top = FALSE, bottom = FALSE,
                     iso_sol = numeric(0), iso_dsh = numeric(0)) {
  if ("comp" %in% names(d)) d <- filter(d, comp == "total")   # the maps show the total, not its parts
  d  <- rowfac(d, params) %>% mutate(zz = pmin(pmax(zvar(.data), min(brk)), max(brk)))
  mk <- marks(levels(droplevels(d$prow)))
  uni <- if ("unident" %in% names(d)) filter(d, unident) else d[0, ]
  ggplot(d, aes(lx, ly, z = zz)) +
    geom_contour_filled(breaks = brk) +
    scale_fill_manual(values = pal, drop = FALSE, guide = "none") +
    # iso-distortion contours, ON TOP of the bands: the same two criteria the colour bar ticks, so a
    # reader can trace where a panel crosses them without counting bands. Drawn after the fill and
    # before the sample points. An empty level vector adds nothing.
    (if (length(iso_dsh)) geom_contour(breaks = iso_dsh, colour = "grey25", linewidth = .22,
                                       linetype = "dashed") else NULL) +
    (if (length(iso_sol)) geom_contour(breaks = iso_sol, colour = "black", linewidth = .3) else NULL) +
    geom_point(size = .08, colour = "grey35", alpha = .35) +
    geom_tile(data = uni, aes(lx, ly), inherit.aes = FALSE, fill = "grey60",
              colour = "grey30", linewidth = .25, width = 0.30, height = 1.0) +
    geom_point(data = mk, aes(lx, ly), inherit.aes = FALSE, shape = 21, fill = "white",
               colour = "black", size = 1.5, stroke = .45) +
    geom_text(data = mk, aes(lx, ly, label = lab), inherit.aes = FALSE, size = 1.9,
              fontface = "bold", vjust = -0.95) +
    # free + proportional y: an algorithm run a noise decade higher gets a proportionally taller row,
    # so a decade of noise is the SAME vertical distance everywhere; a taller row is more noise, not
    # a stretched axis. Where an algorithm has no cell the panel stays grey.
    facet_grid(prow ~ nchf, labeller = label_parsed, scales = "free_y", space = "free_y") +
    scale_x_continuous(breaks = X_BRK, labels = c(".01", ".1", "1")) +
    scale_y_continuous(breaks = Y_BRK, labels = Y_LAB) +
    labs(y = "dimensionless instrumental noise") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    # GREY PANEL, not white: grey is the absence of data, white inside the filled region is
    # calibration. Against a white panel the two would be the same picture.
    theme(panel.background = element_rect(fill = "grey85", colour = NA),
          panel.grid.minor = element_blank(), panel.spacing = unit(0.12, "cm"),
          legend.position = "none", axis.text = element_text(size = 6),
          axis.title.y = element_text(size = 8), strip.text.y = element_text(size = 6.5, angle = 0),
          plot.margin = margin(2, 2, 2, 2)) +
    (if (top) theme(strip.text.x = element_text(size = 7.5)) else theme(strip.text.x = element_blank())) +
    (if (bottom) labs(x = expression(widetilde(Delta) == Delta %.% k[off])) else
       theme(axis.title.x = element_blank(), axis.text.x = element_blank()))
}

colorbar <- function(pal, edge_lab, title, lin, sol = numeric(0), dsh = numeric(0)) {
  n <- length(pal); pos <- function(v) match(v, lin) - 0.5
  mkx <- c(
    if (length(dsh)) list(annotate("segment", x = pos(dsh), xend = pos(dsh), y = 0.5, yend = 2.0,
                                   colour = "grey25", linewidth = .35, linetype = "dashed")),
    if (length(sol)) list(annotate("segment", x = pos(sol), xend = pos(sol), y = 0.5, yend = 2.0,
                                   colour = "black", linewidth = .5)))
  ggplot(data.frame(x = seq_len(n)), aes(x, 1, fill = factor(x))) +
    geom_tile(colour = "grey92", linewidth = .15) + mkx +
    scale_fill_manual(values = pal, guide = "none") +
    scale_x_continuous(breaks = seq(0.5, n + 0.5, 1), labels = edge_lab, expand = c(0, 0)) +
    coord_cartesian(ylim = c(0.5, 2.0), clip = "off") +
    labs(y = title) + theme_void(base_family = "Helvetica") +
    theme(axis.title.y = element_text(angle = 0, size = 7, hjust = 1, vjust = .5, margin = margin(r = 4)),
          axis.text.x = element_text(size = 5.5), axis.ticks.x = element_line(colour = "grey55", linewidth = .3),
          axis.ticks.length.x = unit(2, "pt"), plot.margin = margin(3, 26, 2, 34))
}

btitle <- function(txt) ggplot() + theme_void(base_family = "Helvetica") +
  annotate("text", x = 0, y = 0, label = txt, hjust = 0, size = 2.9, fontface = "bold") +
  coord_cartesian(xlim = c(0, 1), clip = "off") + theme(plot.margin = margin(3, 2, 0, 26))

# height per noise decade, not per row: with facet space = "free_y" a row that spans more noise is
# proportionally taller, so the block height must track the SUMMED noise span of its rows or
# something gets squeezed. Only algorithms that actually have cells contribute.
map_height <- function(nparams, per_row) {
  decades <- sum(vapply(intersect(ALGOS, unique(cells$algo)), noise_span, numeric(1)))
  per_row * nparams * decades
}

# ---- TRANSPOSED map block: columns = ALGORITHM, rows = N_ch --------------------------------------
# A more space-efficient arrangement of the same plane (Luciano, 2026-07-25). One PARAMETER per grid,
# the algorithm cost ladder (LSE -> IR) across the columns and N_ch down the rows. With
# space = "free_y" the row height tracks each N_ch's noise reach: the N_ch 10 row (everyone cured by
# noise 100) is SHORT and the N_ch 10000 row (LSE/NR climb to 1e5 on the diagonal) is TALL, each
# proportional to the decades it actually spans. The algorithm-row layout, by contrast, stretched
# every low-N_ch column onto the tall shared axis of its algorithm's highest-N_ch reach. Rows read
# N_ch 10 (top) to 10000 (bottom); within a row all algorithm columns share that N_ch's noise range,
# so an R/IR cell that cured at noise 100 shows white below and grey above, which is honest.
nch_noise_span <- function(n) {
  allc <- dplyr::bind_rows(GRID_BY_ALGO)
  z <- as.numeric(unique(allc$z[allc$nch == n]))
  if (!length(z)) 0 else diff(range(log10(z)))
}
map_height_byNch <- function(nparams, per_row)
  per_row * nparams * sum(vapply(NCHS, nch_noise_span, numeric(1)))

mapblock_byNch <- function(d, zvar, pal, brk, param, top = FALSE, bottom = FALSE,
                           iso_sol = numeric(0), iso_dsh = numeric(0)) {
  if ("comp" %in% names(d)) d <- dplyr::filter(d, comp == "total")
  psel <- param                                        # avoid the arg/column name clash in filter()
  keep <- setdiff(ALGOS, EXCLUDE_ROWS$algo[EXCLUDE_ROWS$param == psel])   # e.g. drop LSE for N_ch
  d <- d %>% dplyr::filter(param == psel, algo %in% keep)
  keep <- intersect(keep, unique(d$algo))              # drop algo columns with no data for this param
  d <- d %>% dplyr::mutate(zz = pmin(pmax(zvar(.data), min(brk)), max(brk)),
                           algof = factor(ALGO_LAB[algo], levels = ALGO_LAB[keep]))
  uni <- if ("unident" %in% names(d)) dplyr::filter(d, unident) else d[0, ]
  ggplot(d, aes(lx, ly, z = zz)) +
    geom_contour_filled(breaks = brk) +
    scale_fill_manual(values = pal, drop = FALSE, guide = "none") +
    (if (length(iso_dsh)) geom_contour(breaks = iso_dsh, colour = "grey25", linewidth = .22,
                                       linetype = "dashed") else NULL) +
    (if (length(iso_sol)) geom_contour(breaks = iso_sol, colour = "black", linewidth = .3) else NULL) +
    geom_point(size = .08, colour = "grey35", alpha = .35) +
    geom_tile(data = uni, aes(lx, ly), inherit.aes = FALSE, fill = "grey60",
              colour = "grey30", linewidth = .25, width = 0.30, height = 1.0) +
    # rows = N_ch (free + proportional height), columns = algorithm (cost ladder)
    facet_grid(nchf ~ algof, labeller = labeller(nchf = label_parsed, .default = label_value),
               scales = "free_y", space = "free_y") +
    scale_x_continuous(breaks = X_BRK, labels = c(".01", ".1", "1")) +
    scale_y_continuous(breaks = Y_BRK, labels = Y_LAB) +
    labs(y = "dimensionless instrumental noise") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(panel.background = element_rect(fill = "grey85", colour = NA),
          panel.grid.minor = element_blank(), panel.spacing = unit(0.12, "cm"),
          legend.position = "none", axis.text = element_text(size = 6),
          axis.title.y = element_text(size = 8), strip.text.y = element_text(size = 6.5, angle = 0),
          plot.margin = margin(2, 2, 2, 2)) +
    (if (top) theme(strip.text.x = element_text(size = 7.5)) else theme(strip.text.x = element_blank())) +
    (if (bottom) labs(x = expression(widetilde(Delta) == Delta %.% k[off])) else
       theme(axis.title.x = element_blank(), axis.text.x = element_blank()))
}
