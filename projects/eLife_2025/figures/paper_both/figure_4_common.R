## figure_4_common.R — shared machinery for the TWO figure-4 halves.
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
ALGOS    <- c("nonlinearsqr", "macro_NR", "macro_R", "macro_IR")
ALGO_LAB <- c(nonlinearsqr = "LSE", macro_NR = "NR", macro_R = "R", macro_IR = "IR")

NCHS      <- c(10, 100, 1000, 10000)
DATA_DIRS <- c("../data/1c2ae6f", "../data/87889e6", "../data/0ffbda7")  # search path, first hit wins

# ---- file naming: the prefix is PER ALGORITHM ------------------------------------------------
# LSE (nonlinearsqr) is produced by dispatch_figure_4.sh's LSE arm, which routes to
# figure_4_LSE.macroir and stamps the prefix `figure_3_LSE_`; macro/micro get `figure_3_G_`
# (dispatch_figure_4.sh:217-218). Everything else about the two file sets is identical, including
# the Gaussian_ component names and the _G battery suffix.
.pre  <- function(a) if (identical(a, "nonlinearsqr")) "figure_3_LSE_" else "figure_3_G_"
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

# ---- the grid, DECLARED per algorithm (ragged), never a product ------------------------------
# A higher N_ch needs a higher absolute noise to reach the same regime, so the swept levels differ
# per column and per algorithm. Declaring the pairs keeps the missing-cell hard stop working on a
# ragged grid: an absent file is an error, never a silently shortened axis.
GRID_BASE <- tibble::tribble(
  ~nch,  ~z,
  10,    "0.1",  10,    "1",  10,    "10",
  100,   "0.1",  100,   "1",  100,   "10",
  1000,  "0.1",  1000,  "1",  1000,  "10",
  10000, "0.1",  10000, "1",  10000, "10")
# macro_R ONLY: the noise-100 rungs of the constant-r diagonal (r = noise/N_ch = 1 at N_ch 100,
# r = 0.1 at N_ch 1000). IR has no nsim-10000 cell there, so R's rows run one noise decade taller.
GRID_R_EXTRA <- tibble::tribble(
  ~nch,  ~z,
  100,   "100",
  1000,  "100")
# LSE: verified 2026-07-23, nonlinearsqr has ZERO complete (pool_G + sim_G) cells at nsim 10000 —
# its runs are in flight. ADD A ROW THE MOMENT A CELL LANDS; until then LSE contributes no cells,
# so nothing is looked up for it and its map rows simply do not appear (facet_grid drops the empty
# level). Target = GRID_BASE, i.e. N_ch {10,100,1000,10000} x noise {0.1,1,10} at nsim 10000.
GRID_LSE <- tibble::tribble(
  ~nch,  ~z)
GRID_BY_ALGO <- list(nonlinearsqr = GRID_LSE,
                     macro_NR     = GRID_BASE,
                     macro_R      = bind_rows(GRID_BASE, GRID_R_EXTRA),
                     macro_IR     = GRID_BASE)
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
DIST_LIN    <- c(0.5, 0.67, 0.8, 0.87, 0.91, 0.94, 0.97, 0.99, 1.01, 1.03, 1.06, 1.1, 1.15, 1.25, 1.5, 1.7, 2.5)
DIST_BREAKS <- log10(DIST_LIN)
DIST_PAL    <- bandcols(DIST_BREAKS)
DIST_LAB    <- fmt_big(DIST_LIN)
# log-spaced bands: IR's values (0.003-0.042) and R's (0.13-0.15) differ by more than a decade, and a
# linear scale would put every IR cell in the neutral band.
BIAS_LIN    <- c(-0.3, -0.15, -0.07, -0.03, -0.01, -0.003, 0.003, 0.01, 0.03, 0.07, 0.15, 0.3)
BIAS_PAL    <- bandcols(BIAS_LIN)
BIAS_LAB    <- formatC(BIAS_LIN, format = "g", digits = 2)

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

# BIAS comes from battery_sim_G, anchored at theta_sim. It cannot come from the pool files: a bias
# evaluated at the optimum is zero by construction because the score vanishes there.
bias <- readcell("battery_sim_G") %>%
  filter(component_path == "Probit_statistics_Gaussian_Distortion_Induced_Bias",
         param_index %in% PIDX[PARAM_LONG],
         statistic == "value", probit %in% c("mean", "quantile")) %>%
  mutate(key = case_when(probit == "mean" ~ "m", quantile_level == 0.025 ~ "lo",
                         quantile_level == 0.975 ~ "hi", TRUE ~ NA_character_)) %>%
  filter(!is.na(key)) %>%
  group_by(algo, noise, Num_ch, interval_in_tau, param_index, key) %>%
  summarise(v = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = key, values_from = v) %>%
  # CI-aware toward ZERO: a cell whose interval covers 0 renders white
  mutate(Bconf = ifelse(m > 0, pmax(0, lo), pmin(0, hi)),
         lx = log10(interval_in_tau), ly = log10(as.numeric(noise)),
         param = names(PIDX)[match(param_index, PIDX)],
         nchf  = factor(paste0("N[ch]==", Num_ch), levels = paste0("N[ch]==", NCHS)))

stat <- readcell("battery_pool_G") %>%
  filter(component_path %in% COMPS, param_index == param_col, param_index %in% PIDX[PARAM_LONG],
         statistic == "value", probit %in% c("mean", "quantile")) %>%
  mutate(comp = names(COMPS)[match(component_path, COMPS)],
         key  = case_when(probit == "mean" ~ "m", quantile_level == 0.025 ~ "lo",
                          quantile_level == 0.975 ~ "hi", TRUE ~ NA_character_)) %>%
  filter(!is.na(key)) %>%
  group_by(algo, comp, noise, Num_ch, interval_in_tau, param_index, key) %>%
  summarise(v = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = key, values_from = v) %>%
  # CI-aware: a cell whose 95% interval brackets 1 collapses to exactly 1 and renders white
  mutate(Dconf = ifelse(m > 1, pmax(1, lo), pmin(1, hi)),
         lx = log10(interval_in_tau), ly = log10(as.numeric(noise)),
         param = names(PIDX)[match(param_index, PIDX)],
         nchf  = factor(paste0("N[ch]==", Num_ch), levels = paste0("N[ch]==", NCHS)))

# UNIDENTIFIED cells: greyed, not clamped. Off-scale AND ill-conditioned (kappa > 3e4) together
# isolate the R corner (N_ch 10000, Delta = 1), a displacement along a near-null direction whose
# magnitude is not a meaningful bias/distortion. kappa alone grows with N_ch for everyone.
kap <- readcell("battery_pool_G") %>%
  filter(component_path == "Probit_statistics_Spectrum_Condition_Number_Gaussian_Fisher_Covariance",
         statistic == "value", probit == "mean") %>%
  group_by(algo, Num_ch, noise, interval_in_tau) %>% summarise(kappa = mean(value), .groups = "drop")
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
         subtitle = bquote(N[ch] == .(pt$nch) ~ "," ~ sigma == .(pt$noise) ~ "," ~
                           Delta %.% k[off] == .(pt$dlt))) +
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
Y_LAB <- formatC(10^Y_BRK, format = "g", digits = 3)

# ONE row factor, not nested facets: levels ordered parameter-major, algorithm-minor, so the four
# algorithms of a parameter are always adjacent and the ladder reads top-to-bottom in each block.
rowfac <- function(d, params) {
  lv <- as.vector(t(outer(params, ALGOS,
          function(p, a) sprintf("atop(%s, bold(\"%s\"))", PARAM_MATH[p], ALGO_LAB[a]))))
  d %>% filter(param %in% params) %>%
    mutate(prow = factor(sprintf("atop(%s, bold(\"%s\"))", PARAM_MATH[param], ALGO_LAB[algo]),
                         levels = lv))
}

marks <- function(rows) ELL_POINTS %>%
  mutate(lx = log10(dlt), ly = log10(as.numeric(noise)),
         nchf = factor(paste0("N[ch]==", nch), levels = paste0("N[ch]==", NCHS))) %>%
  tidyr::crossing(prow = factor(rows, levels = rows))

mapblock <- function(d, zvar, pal, brk, params, top = FALSE, bottom = FALSE) {
  if ("comp" %in% names(d)) d <- filter(d, comp == "total")   # the maps show the total, not its parts
  d  <- rowfac(d, params) %>% mutate(zz = pmin(pmax(zvar(.data), min(brk)), max(brk)))
  mk <- marks(levels(droplevels(d$prow)))
  uni <- if ("unident" %in% names(d)) filter(d, unident) else d[0, ]
  ggplot(d, aes(lx, ly, z = zz)) +
    geom_contour_filled(breaks = brk) +
    scale_fill_manual(values = pal, drop = FALSE, guide = "none") +
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
    labs(y = "noise") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    # GREY PANEL, not white: grey is the absence of data, white inside the filled region is
    # calibration. Against a white panel the two would be the same picture.
    theme(panel.background = element_rect(fill = "grey85", colour = NA),
          panel.grid.minor = element_blank(), panel.spacing = unit(0.12, "cm"),
          legend.position = "none", axis.text = element_text(size = 6),
          axis.title.y = element_text(size = 8), strip.text.y = element_text(size = 6.5, angle = 0),
          plot.margin = margin(2, 2, 2, 2)) +
    (if (top) theme(strip.text.x = element_text(size = 7.5)) else theme(strip.text.x = element_blank())) +
    (if (bottom) labs(x = expression(Delta %.% k[off])) else
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
