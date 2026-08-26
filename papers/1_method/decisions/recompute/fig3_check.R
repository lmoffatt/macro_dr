# fig3_check.R -- two measurements the Figure 3 notebook did not print when this was written: the
# accumulated standardized score bias in its DIAGONAL form, E[sum s_t]/sqrt(sum F_t), and the exact
# factorisation J_T/F_T = Ebar * inflation. Cited by figure_3_guion.md.
#
# NOTE the diagonal reading is superseded for the bias: what the body reports is the vector form
# b = F^-1 E[sum s] with the WHOLE information matrix, in fig3_score_bias.R, which at this cell
# differs in sign on the channel number. This script is kept for the factorisation.
#
# Run from the repository root:  Rscript papers/1_method/decisions/recompute/fig3_check.R
# PROMOTED 2026-08-26 from tmp/, which is gitignored; the paths changed and the saved object now
# lands beside the script instead of in tmp/.
# fig3_check.R  -- two measurements the Figure 3 notebook does not currently print:
#   (1) the ACCUMULATED standardized score bias, E[sum s_t]/sqrt(sum F_t), which is what determines
#       the estimator bias; the notebook prints only the per-interval version (row D).
#   (2) the exact factorisation of the accumulated ratio, J_T/F_T = Ebar_w * inflation, where
#       Ebar_w = sum Var(s_t) / sum F_t is the per-interval part (row E, weighted) and
#       inflation  = Var(sum s_t) / sum Var(s_t) is the correlation part (row G's content).
# Reproduces the notebook's own data prep, mask and bootstrap so the numbers are comparable.
suppressPackageStartupMessages({library(data.table)})

PSEL <- c(k_off = 1L, N_ch = 5L)
KEYS <- c("LSE", "ILSE", "NR", "INR", "R", "MR", "VR", "IR")
DTOK <- c(LSE = "LSE_av0", ILSE = "LSE", NR = "NR", INR = "INR", R = "R", MR = "MR",
          VR = "VR", IR = "IR")
B <- 400; FMIN <- 1e-6
set.seed(20260731)

colVar <- function(M) { n <- nrow(M); m <- colMeans(M); (colMeans(M * M) - m^2) * n / (n - 1) }
rowCumsum <- function(M) t(apply(M, 1, cumsum))
wide <- function(d, val) as.matrix(dcast(d, simulation_index ~ sample_index, value.var = val)[, -1])
acf_one <- function(x) as.numeric(acf(x, lag.max = 25, plot = FALSE)$acf)

out <- list()
for (k in KEYS) {
  dg <- readRDS(sprintf("projects/eLife_2025/figures/data/digest/figure_3_digest_%s.rds", DTOK[[k]]))
  P  <- as.data.table(dg$param)[param_index %in% PSEL]
  for (pn in names(PSEL)) {
    if (pn == "N_ch" && k %in% c("LSE", "ILSE")) next          # excluded, as in the notebook
    d <- P[param_index == PSEL[[pn]]]
    Sm <- wide(d, "s"); Im <- wide(d, "I")
    n  <- nrow(Sm); smp <- as.integer(colnames(Sm))
    Fb <- colMeans(Im)
    infm <- which(is.finite(Fb) & Fb > FMIN)                    # informative steps
    Tstar <- max(infm)                                          # last informative index, as dFa does
    idxT  <- seq_len(Tstar)

    Ssum <- rowSums(Sm[, idxT, drop = FALSE])                   # sum of the score per recording
    Fsum <- sum(Fb[idxT])
    zT   <- mean(Ssum) / sqrt(Fsum)                             # (1) accumulated standardized bias
    Vsum <- sum(colVar(Sm[, idxT, drop = FALSE]))
    Ebar <- Vsum / Fsum                                         # (2) per-interval part, F-weighted
    infl <- var(Ssum) / Vsum                                    # correlation part
    rC   <- var(Ssum) / Fsum                                    # the notebook's accumulated ratio

    # bootstrap over recordings, the notebook's estimator: resample rows, recompute in full
    rep_z <- rep_E <- rep_i <- numeric(B)
    for (b in seq_len(B)) {
      j  <- sample.int(n, n, replace = TRUE)
      Sb <- Sm[j, idxT, drop = FALSE]; Ib <- Im[j, idxT, drop = FALSE]
      fs <- sum(colMeans(Ib)); ss <- rowSums(Sb); vs <- sum(colVar(Sb))
      rep_z[b] <- mean(ss) / sqrt(fs); rep_E[b] <- vs / fs; rep_i[b] <- var(ss) / vs
    }
    q <- function(x) quantile(x, c(.025, .975), na.rm = TRUE)

    # lag-one ACF of the score, the notebook's estimator (per recording, then averaged)
    A <- do.call(rbind, d[, .(a = list(acf_one(s))), by = simulation_index]$a)
    rho1 <- mean(A[, 2], na.rm = TRUE)
    ar1  <- 1 + 2 * rho1 / (1 - rho1)                           # AR(1) prediction of the inflation

    out[[length(out) + 1]] <- data.table(
      mem = k, param = pn, n_inf = length(infm), Tstar = Tstar,
      zT = zT, zT_lo = q(rep_z)[1], zT_hi = q(rep_z)[2],
      Ebar = Ebar, Ebar_lo = q(rep_E)[1], Ebar_hi = q(rep_E)[2],
      infl = infl, infl_lo = q(rep_i)[1], infl_hi = q(rep_i)[2],
      rho1 = rho1, ar1 = ar1, rC = rC, prod_check = Ebar * infl)
  }
  message("done ", k)
}
R <- rbindlist(out)
setnames(R, "mem", "key"); R[, key := factor(key, levels = KEYS)]; setorder(R, key, param)
print(as.data.frame(R[, .(key, param, n_inf, zT, zT_lo, zT_hi, rC)]), digits = 4, row.names = FALSE)
cat("\n--- factorisation: rC = Ebar * inflation, and the AR(1) prediction of the inflation ---\n")
print(as.data.frame(R[, .(key, param, Ebar, Ebar_lo, Ebar_hi, infl, infl_lo, infl_hi,
                          rho1, ar1, rC, prod_check)]), digits = 4, row.names = FALSE)
saveRDS(R, "papers/1_method/decisions/recompute/fig3_check.rds")
