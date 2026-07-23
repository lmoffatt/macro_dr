# Shared machinery for Figure 4 supplements 1, 2 and 4 — the all-five-parameters R-vs-IR maps.
# One block per file (bias / distortion / correlation-fraction), so each supplement is ~9.7 in tall
# rather than one 19 in stack. Sourced, not duplicated: the palette, the readers and the row grammar
# must stay identical across the three or a reader comparing them is comparing rendering, not data.

suppressMessages({ library(tidyverse); library(patchwork) })

ALGOS    <- c("macro_R", "macro_IR")
ALGO_LAB <- c(macro_R = "R", macro_IR = "IR")
NCHS   <- c(10, 100, 1000, 10000)
NOISES <- c("0.1", "1", "10")
DATA_DIRS <- c("../data/1c2ae6f", "../data/87889e6")   # search path: R's noise 1/10 live apart

# all five, kinetic then amplitude, sigma_noise last as the shared-null control
PIDX  <- c(on = 0, off = 1, unitary_current = 2, Num_ch_mean = 5, Current_Noise = 3)
PMATH <- c(on = "log[10]~k[on]", off = "log[10]~k[off]", unitary_current = "log[10]~i",
           Num_ch_mean = "log[10]~N[ch]", Current_Noise = "log[10]~sigma[noise]")
PARAM_ORD <- names(PIDX)

DIST <- "Probit_statistics_Likelihood_Gaussian_Information_Distortion"
BIAS <- "Probit_statistics_Gaussian_Distortion_Induced_Bias"
SAMP <- "Probit_statistics_Gaussian_Sample_Distortion"
CORR <- "Probit_statistics_Likelihood_Correlation_Distortion"

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
DIST_BREAKS <- log10(DIST_LIN); DIST_PAL <- bandcols(DIST_BREAKS); DIST_LAB <- fmt_big(DIST_LIN)
# same +-0.3 bias scale figure_4 settled on, with the same three-cell clamp (R N_ch 10000, Delta=1)
BIAS_LIN <- c(-0.3, -0.15, -0.07, -0.03, -0.01, -0.003, 0.003, 0.01, 0.03, 0.07, 0.15, 0.3)
BIAS_PAL <- bandcols(BIAS_LIN); BIAS_LAB <- formatC(BIAS_LIN, format = "g", digits = 2)
# correlation FRACTION = log(corr)/log(total): 0 = all sample, 1 = all correlation. Its own 0..1 ramp.
FRAC_LIN <- seq(0, 1, by = 0.1)
FRAC_PAL <- colorRampPalette(c("#F7F7F7", "#4393C3", "#053061"), space = "Lab")(length(FRAC_LIN) - 1)
FRAC_LAB <- formatC(FRAC_LIN, format = "g", digits = 2)

findf <- function(nch, algo, z, kind) {
  f <- file.path(DATA_DIRS, sprintf("figure_3_G_nch_%d_nsim_10000_%s_noise_%s_%s.csv", nch, algo, z, kind))
  hit <- f[file.exists(f)]
  if (!length(hit)) stop("supp: no file for ", algo, " N_ch ", nch, " noise ", z, " (", kind, ")")
  hit[1]
}
CELLS <- expand_grid(algo = ALGOS, nch = NCHS, z = NOISES)
readcell <- function(kind) pmap_dfr(CELLS, function(algo, nch, z)
  read.csv(findf(nch, algo, z, kind), skip = 1, stringsAsFactors = FALSE) %>%
    mutate(noise = z, algo = algo))

# collapse the mean and the two quantiles to one CI-shrunk value; ci_toward is the null (0 or 1)
shape <- function(d, ci_toward) d %>%
  mutate(key = case_when(probit == "mean" ~ "m", quantile_level == 0.025 ~ "lo",
                         quantile_level == 0.975 ~ "hi", TRUE ~ NA_character_)) %>%
  filter(!is.na(key)) %>%
  group_by(algo, noise, Num_ch, interval_in_tau, param_index, key) %>%
  summarise(v = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = key, values_from = v) %>%
  mutate(conf = if (ci_toward == 1) ifelse(m > 1, pmax(1, lo), pmin(1, hi))
                else                 ifelse(m > 0, pmax(0, lo), pmin(0, hi)),
         lx = log10(interval_in_tau), ly = log10(as.numeric(noise)),
         param = names(PIDX)[match(param_index, PIDX)])

# UNIDENTIFIED cells: greyed, not clamped. A cell is greyed where its value falls outside the plotted
# scale AND the Gaussian Fisher covariance is ill-conditioned there (condition number > 3e4). In this
# data these coincide exactly and are only R at N_ch 10000, Delta*k_off = 1 (kappa ~ 66000): a
# displacement along a near-null direction, whose magnitude is not a meaningful bias/distortion.
# kappa alone is NOT the criterion (it grows with N_ch for everyone and would grey the whole 10000
# column); the clamp AND the high kappa together isolate the artifact. See figures_build_plan / the
# figure legends.
KAPPA_CP  <- "Probit_statistics_Spectrum_Condition_Number_Gaussian_Fisher_Covariance"
KAPPA_MAX <- 3e4
read_kappa <- function() readcell("battery_pool_G") %>%
  filter(component_path == KAPPA_CP, statistic == "value", probit == "mean") %>%
  group_by(algo, Num_ch, noise, interval_in_tau) %>%
  summarise(kappa = mean(value), .groups = "drop")
# smin/smax, NOT lo/hi: shape() already puts `lo` and `hi` quantile columns in d, and an argument
# named lo/hi would be shadowed by those columns inside mutate (the value would silently never clamp).
add_unident <- function(d, smin, smax, kap) d %>%
  left_join(kap, by = c("algo", "Num_ch", "noise", "interval_in_tau")) %>%
  mutate(unident = (m < smin | m > smax) & is.finite(kappa) & kappa > KAPPA_MAX)

# one map block: rows = parameter (outer) then algorithm (inner), columns = N_ch.
X_BRK <- log10(c(.01, .1, 1)); Y_BRK <- log10(c(.1, 1, 10))
rowfac <- function(d) {
  lv <- as.vector(t(outer(PARAM_ORD, ALGOS,
          function(p, a) sprintf("atop(%s, bold(\"%s\"))", PMATH[p], ALGO_LAB[a]))))
  d %>% mutate(prow = factor(sprintf("atop(%s, bold(\"%s\"))", PMATH[param], ALGO_LAB[algo]), levels = lv),
               nchf = factor(paste0("N[ch]==", Num_ch), levels = paste0("N[ch]==", NCHS)))
}
mapblock <- function(d, zexpr, pal, brk) {
  d <- rowfac(d) %>% mutate(zz = pmin(pmax(zexpr(.), min(brk)), max(brk)))
  uni <- if ("unident" %in% names(d)) filter(d, unident) else d[0, ]
  ggplot(d, aes(lx, ly, z = zz)) +
    geom_contour_filled(breaks = brk) +
    scale_fill_manual(values = pal, drop = FALSE, guide = "none") +
    geom_point(size = .08, colour = "grey35", alpha = .35) +
    # grey out the unidentified cells, drawn over the fill: value not meaningful there
    geom_tile(data = uni, aes(lx, ly), inherit.aes = FALSE, fill = "grey78",
              colour = "grey45", linewidth = .2, width = 0.30, height = 1.0) +
    facet_grid(prow ~ nchf, labeller = label_parsed) +
    scale_x_continuous(breaks = X_BRK, labels = c(".01", ".1", "1")) +
    scale_y_continuous(breaks = Y_BRK, labels = c(".1", "1", "10")) +
    labs(x = expression(Delta %.% k[off]), y = "noise") +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(panel.grid.minor = element_blank(), panel.spacing = unit(0.1, "cm"),
          legend.position = "none", axis.text = element_text(size = 6),
          axis.title = element_text(size = 8), strip.text.y = element_text(size = 6, angle = 0),
          strip.text.x = element_text(size = 7.5), plot.margin = margin(2, 2, 2, 2))
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
sat <- function(v, lo, hi, what) cat(sprintf("  %-12s %d/%d below %.3g, %d/%d above %.3g\n",
  what, sum(v < lo, na.rm = TRUE), sum(is.finite(v)), lo, sum(v > hi, na.rm = TRUE), sum(is.finite(v)), hi))
