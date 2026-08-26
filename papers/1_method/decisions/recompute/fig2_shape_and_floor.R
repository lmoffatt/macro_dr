# fig2_shape_and_floor.R — the second-moment SHAPE numbers of the Figure 2 subsection.
#
# Run from the repository root:  Rscript papers/1_method/decisions/recompute/fig2_shape_and_floor.R
#
# figure_2.Rmd prints two orange numbers per panel: the area ratio sqrt(det S_emp / det S_F) and the
# shape ratio sqrt(lambda_max / lambda_min) of S_F^{-1} S_emp. They are the two eigenvalues of the
# distortion written as a geometric mean and a spread, so lambda_max = area * shape and
# lambda_min = area / shape. This script recovers both, and the EIGENVECTORS the figure cannot draw,
# so that the direction each eigenvalue belongs to can be named in the text.
#
# It also simulates the finite-sample floor of the shape statistic: with a perfectly calibrated
# covariance and a finite cloud, sqrt(lambda_max/lambda_min) is still above one. The floor printed
# here is a LOWER bound on the figure's, because it treats S_F as exact whereas the figure's S_F is
# itself estimated at the pooled fit.
# Written 2026-08-06 for the anisotropy paragraph of 04_results.tex.

# UPDATED 2026-08-06 (second pass): eight members, both least-squares arms, and the panel's BLUE
# number (empirical over sandwich) is printed alongside the two orange ones.

suppressPackageStartupMessages({library(dplyr); library(tidyr); library(purrr); library(stringr)})

FIG      <- "projects/eLife_2025/figures"
DIRS     <- file.path(FIG, "data", c("a202e03", "1f7138b", "1c2ae6f", "0ffbda7", "87889e6"))
FILE_PRE <- "figure_3_G_"; BAT_SUF <- "_G"; COMP_PRE <- "Gaussian_"
REP_NCH  <- 100; REP_NSIM <- 10000; REP_NOISE <- "0.1"; REP_GROUP <- 10; REP_INT <- 0.1

ROSTER <- c("nonlinearsqr_g", "nonlinearsqr", "macro_NR", "macro_INR", "macro_R", "macro_MR", "macro_VR", "macro_IR")
LAB    <- c(nonlinearsqr_g = "LSE", nonlinearsqr = "ILSE", macro_NR = "NR", macro_INR = "INR", macro_R = "R",
            macro_MR = "MR", macro_VR = "VR", macro_IR = "IR")
IDX    <- c(on = 0, off = 1, unitary_current = 2, Current_Noise = 3, Num_ch_mean = 5)

rd   <- function(p) read.csv(p, skip = 1, stringsAsFactors = FALSE)
pre  <- function(a) if (startsWith(a, "nonlinearsqr")) "figure_3_LSE_" else FILE_PRE
stem <- function(d, a) file.path(d, sprintf("%snch_%d_nsim_%d_%s_noise_%s",
                                            pre(a), REP_NCH, REP_NSIM, a, REP_NOISE))
need <- function(s) paste0(s, c("_mle_cloud_runs.csv",
                                paste0("_battery_pool", BAT_SUF, ".csv"),
                                paste0("_battery_sim",  BAT_SUF, ".csv")))
base_for <- function(a) { for (d in DIRS) { s <- stem(d, a); if (all(file.exists(need(s)))) return(s) }
                          NA_character_ }
b <- vapply(ROSTER, base_for, character(1)); ALGOS <- ROSTER[!is.na(b)]; b <- b[!is.na(b)]

pts <- map_dfr(paste0(b, "_mle_cloud_runs.csv"), rd) %>%
  filter(component_path == "Model_Parameters_Hat", statistic == "value",
         Num_ch == REP_NCH, group_size == REP_GROUP, abs(interval_in_tau - REP_INT) < 1e-9,
         param_name %in% names(IDX)) %>%
  select(algorithm, sample_index, param_name, value) %>%
  pivot_wider(names_from = param_name, values_from = value)

bat <- map_dfr(paste0(b, "_battery_pool", BAT_SUF, ".csv"), rd) %>%
  mutate(component = str_replace(component_path, "^Probit_statistics_", ""))

# the same 2x2 the figure draws: reported covariance at the pooled fit, scaled by 1/group_size
cov_pair <- function(comp, algo, ix, jx) {
  d <- bat %>% filter(component == comp, statistic == "value", probit == "mean", algorithm == algo,
                      Num_ch == REP_NCH, abs(interval_in_tau - REP_INT) < 1e-9,
                      param_index %in% c(ix, jx), param_col %in% c(ix, jx))
  if (nrow(d) < 4) return(NULL)
  k <- c(ix, jx); M <- matrix(NA_real_, 2, 2)
  for (r in seq_len(nrow(d))) M[match(d$param_index[r], k), match(d$param_col[r], k)] <- d$value[r]
  if (any(is.na(M))) return(NULL)
  M / REP_GROUP
}

PAIRS <- list(kinetic = c("on", "off"),
              amplitude = c("Num_ch_mean", "unitary_current"),
              noise = c("Current_Noise", "Num_ch_mean"))

# Rows accumulated for Figure 2--source data 2 (added 2026-08-25: the digit-migration pass moves
# the per-pair numbers out of the Results prose; their reviewer-visible home is this CSV).
sd_rows <- list()
for (pn in names(PAIRS)) {
  px <- PAIRS[[pn]][1]; py <- PAIRS[[pn]][2]
  cat("\n===== ", pn, " pair (", px, ", ", py, ") =====\n", sep = "")
  for (a in ALGOS) {
    cl <- pts %>% filter(algorithm == a) %>% select(all_of(c(px, py))) %>% filter(complete.cases(.))
    SF <- cov_pair(paste0(COMP_PRE, "Fisher_Covariance"), a, IDX[px], IDX[py])
    if (nrow(cl) < 3 || is.null(SF)) next
    E   <- eigen(solve(SF, cov(as.matrix(cl))))
    lam <- Re(E$values); V <- Re(E$vectors)
    o <- order(lam, decreasing = TRUE); lam <- lam[o]; V <- V[, o, drop = FALSE]
    deg <- function(v) { v <- v / sqrt(sum(v^2)); 180 / pi * atan2(abs(v[2]), abs(v[1])) }
    ang <- function(v) sprintf("%5.1f deg off the %s axis", deg(v), px)
    SC  <- cov_pair(paste0(COMP_PRE, "Distortion_Corrected_Covariance"), a, IDX[px], IDX[py])
    aC  <- if (is.null(SC)) NA_real_ else
             sqrt(det(cov(as.matrix(cl)))) / sqrt(det(SC))   # the panel's BLUE number
    cat(sprintf("%-4s area %5.2f  shape %5.2f  corrected %5.2f | worst var-ratio %6.2f (%s) | best %6.2f (%s)\n",
                LAB[a], sqrt(prod(lam)), sqrt(max(lam) / min(lam)), aC,
                lam[1], ang(V[, 1]), lam[2], ang(V[, 2])))
    sd_rows[[length(sd_rows) + 1]] <- data.frame(
      pair = pn, axis_1 = px, axis_2 = py, member = unname(LAB[a]),
      magnitude = sqrt(prod(lam)), anisotropy = sqrt(max(lam) / min(lam)), corrected = aC,
      worst_var_ratio = lam[1], worst_deg_off_axis_1 = deg(V[, 1]),
      best_var_ratio = lam[2], best_deg_off_axis_1 = deg(V[, 2]))
  }
}
sd_dir <- "projects/eLife_2025/figures/figure_2_source_data"
dir.create(sd_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(dplyr::bind_rows(sd_rows) %>% mutate(across(where(is.numeric), ~round(.x, 4))),
          file.path(sd_dir, "figure_2_source_data_distortion.csv"), row.names = FALSE)
cat("\nwrote ", file.path(sd_dir, "figure_2_source_data_distortion.csv"), "\n", sep = "")

cat("\n===== finite-sample floor of the shape statistic (S_F exact, cloud of n in p dims) =====\n")
set.seed(1)
shape <- function(S) { ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values; sqrt(max(ev) / min(ev)) }
area  <- function(S) sqrt(prod(eigen(S, symmetric = TRUE, only.values = TRUE)$values))
for (p in c(2, 4)) {
  n <- 1000
  s <- replicate(20000, { X <- matrix(rnorm(n * p), n, p); c(shape(cov(X)), area(cov(X))) })
  cat(sprintf("p=%d n=%d   shape: median %.3f  q95 %.3f  q99 %.3f  |  area: median %.3f  q95 %.3f\n",
              p, n, median(s[1, ]), quantile(s[1, ], .95), quantile(s[1, ], .99),
              median(s[2, ]), quantile(s[2, ], .95)))
}
