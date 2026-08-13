# figure_4_cloud_probe.R — what the MLE cloud says about the interval a single experiment quotes.
#
# A PROBE, not a figure. It exists because the Mahalanobis supplement answers a milder question
# than the one a reader has. That figure centres on the cloud MEAN and uses one aggregate
# covariance, so it tests the second moment alone and reports IR covering 0.914 under its own
# Fisher against 0.947 under the sandwich. But every fit in the cloud carries its OWN 6x6 Fisher
# (108000 rows = 3000 fits x 36 in the high-count cell), and an experimenter with one dataset
# quotes THAT, around their own estimate, against the truth. So the statistic to report is
#
#     d2_k = (theta_hat_k - theta_true)' F_k (theta_hat_k - theta_true)   ~ chi2_6 if calibrated
#
# which carries the bias and the variance together, as a confidence interval claims to.
#
# MEASURED 2026-08-01, noise label 0.1, Delta*k_off = 1, group size 100, converged fits only,
# 100 fits per cell (so a coverage carries about +-2.2 per cent):
#
#     coverage of the nominal 95 per cent region     N_ch 10   100   1000   10000
#       R                                              0.000 0.000  0.000  0.000
#       IR                                             0.680 0.970  0.910  0.920
#
# and the reason, as the ratio of the cloud's own spread to the standard error the Fisher declares,
# and the distance of the cloud centre from the truth in units of that same standard error:
#
#                       sd/SE                          bias in SE
#              N_ch     i    k_off  k_on        N_ch      i     k_off   k_on
#       R  10   11.3   4.7    2.4    7.4        +62.1  -50.4   -2.0   -19.3
#       R 1e3   34.8  32.9    2.1    7.7       +410.9 -282.1   -8.7  -155.0
#       IR 10    8.7   4.5    2.0    4.7         +9.2   -5.1   +2.0    -2.9
#       IR 1e3  33.4  32.0    1.9    6.2          0.0   +0.1   -0.1    +0.1
#
# THREE THINGS THIS SAYS, and the third is why it is a probe and not yet a figure.
#
#   1. The recursive instantaneous filter's pooled estimate sits HUNDREDS of standard errors from
#      the truth, and N_ch and i move in OPPOSITE directions by nearly the same amount (+62/-50,
#      +411/-282). That is the N*i ridge, measured: the estimates slide along it while the product
#      stays right, and the Fisher's marginal error bar knows nothing about the slide.
#   2. Two parameters are perfectly calibrated for everyone: the baseline and the noise scale sit
#      at sd/SE = 1.00 and bias 0. So the machinery is not painting everything, and what fails is
#      exactly the gating block.
#   3. GROUP SIZE 100 MEANS 100 RECORDINGS POOLED, so the standard error is ten times smaller than
#      a single recording's and any bias counts ten times more in these units. The coverage of a
#      SINGLE recording is a different and gentler number. Sweeping group size is the missing axis
#      and it is what turns this into a figure: coverage against how much you pool.
#
# LSE IS READ OUT WITH THE WRONG PARAMETER NAMES HERE. nonlinearsqr indexes its own four (2 is the
# baseline, 3 is N_ch, where the macro members have i and the noise scale), so its per-parameter
# rows below are mislabelled and only its coverage is meaningful. Fix the mapping before using it.
#
#   Rscript tmp/cloud_sweep.R     (from the repo root; writes tmp/cloud_sweep.csv)

suppressMessages({ library(data.table) })
setwd("projects/eLife_2025/figures/data")
P <- 6; CRIT <- qchisq(0.95, P); DLT <- 1; GS <- 100; Z <- "0.1"
PN <- c("k_on", "k_off", "i", "sigma", "base", "N_ch")

cells <- rbind(
  data.table(algo = "nonlinearsqr", nch = c(10, 100, 1000)),
  data.table(algo = "macro_NR",     nch = c(10, 100, 1000, 10000)),
  data.table(algo = "macro_R",      nch = c(10, 100, 1000, 10000)),
  data.table(algo = "macro_IR",     nch = c(10, 100, 1000, 10000)))

pre <- function(a) if (a == "nonlinearsqr") "figure_3_LSE_" else "figure_3_G_"
out <- list()
for (i in seq_len(nrow(cells))) {
  a <- cells$algo[i]; n <- cells$nch[i]
  f <- Sys.glob(file.path("*", sprintf("%snch_%d_nsim_10000_%s_noise_%s_mle_cloud_runs.csv",
                                       pre(a), n, a, Z)))
  if (!length(f)) { cat(sprintf("  sin nube: %s N_ch %d\n", a, n)); next }
  d <- fread(f[1], skip = 1, showProgress = FALSE)
  if (!DLT %in% d$interval_in_tau || !GS %in% d$group_size) {
    cat(sprintf("  sin dt/gs: %s N_ch %d\n", a, n)); next }

  cv <- d[variable == "Convergence_Status_Code" & interval_in_tau == DLT & group_size == GS]
  ok <- if (nrow(cv)) cv$sample_index[cv$value == 0] else
        unique(d[variable == "Model_Parameters_Hat"]$sample_index)
  W  <- dcast(d[variable == "Model_Parameters_Hat" & interval_in_tau == DLT &
                group_size == GS & sample_index %in% ok],
              sample_index ~ param_index, value.var = "value")
  ids <- W$sample_index; T <- as.matrix(W[, -1]); rownames(T) <- as.character(ids)
  if (nrow(T) < 20) next
  pidx <- as.integer(colnames(T))

  th0 <- c(1, 2, 0, log10(as.numeric(Z) / 1000), 0, log10(n))[pidx + 1]
  dev <- sweep(T, 2, th0, "-")

  fi <- d[variable == "Gaussian_Fisher_Information" & interval_in_tau == DLT &
          group_size == GS & sample_index %in% ids]
  se <- fi[param_index == param_col, .(se = mean(1 / sqrt(value))), by = param_index]
  own <- vapply(split(seq_len(nrow(fi)), fi$sample_index), function(ix) {
    g <- fi[ix]; F <- matrix(0, P, P)
    F[cbind(g$param_index + 1, g$param_col + 1)] <- g$value
    v <- dev[as.character(g$sample_index[1]), ]
    as.numeric(t(v) %*% F[pidx + 1, pidx + 1] %*% v)
  }, numeric(1))

  out[[length(out) + 1]] <- data.table(
    algo = a, nch = n, nfit = nrow(T), conv = nrow(T) / max(1, nrow(cv)),
    cover = mean(own <= CRIT, na.rm = TRUE),
    param = PN[pidx + 1],
    ratio = apply(T, 2, sd) / se$se[match(pidx, se$param_index)],
    bias  = (colMeans(T) - th0) / se$se[match(pidx, se$param_index)])
  cat(sprintf("  %-14s N_ch %-6d n=%4d conv=%.2f cobertura=%.3f\n", a, n, nrow(T),
              nrow(T) / max(1, nrow(cv)), mean(own <= CRIT, na.rm = TRUE)))
}
R <- rbindlist(out)
fwrite(R, "../../../../tmp/cloud_sweep.csv")
cat("\n=== sd empirica / SE que declara la Fisher ===\n")
print(dcast(R, algo + nch ~ param, value.var = "ratio")[, lapply(.SD, function(x)
  if (is.numeric(x)) round(x, 2) else x)])
cat("\n=== sesgo del centro de la nube, en unidades de ese SE ===\n")
print(dcast(R, algo + nch ~ param, value.var = "bias")[, lapply(.SD, function(x)
  if (is.numeric(x)) round(x, 1) else x)])
