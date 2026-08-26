# fig3_tc.R -- an independent estimate of the total-correlation floor that fig3_kl.R infers as a
# residual, from the marginal entropies of the recorded current (Gaussian plug-in and Vasicek
# m-spacing). Feeds the same commented derivation in 04_results.tex and figure_3_guion.md, and the
# 71.9-nat floor the three open-loop members agree on.
#
# Run from the repository root:  Rscript papers/1_method/decisions/recompute/fig3_tc.R
# PROMOTED 2026-08-26 from tmp/, which is gitignored; only the paths changed.
# Independent check of  KL(p||prod q_t) = sum_t KL(p_t||q_t) + [ sum_t H(p_t) - H(p) ].
#
# The bracket was inferred as a RESIDUAL (observed gap minus the Stein term) and the only evidence
# it is the total correlation was that three open-loop members agree on it. Here it is estimated
# from the data instead.
#
#  * sum_t H(p_t): the marginal entropy of the observed current at each interval, over the 1,000
#    recordings. Estimated two ways, a Gaussian plug-in and Vasicek's m-spacing estimator, which
#    differ exactly when the marginal is not Gaussian.
#  * H(p): not estimable in 100 dimensions from 1,000 samples, but for ANY density q,
#    -E_p[log q] = H(p) + KL(p||q) >= H(p), so every member's mean total logL gives an UPPER bound
#    and the tightest available is IR's, 143.45. So TC >= sum_t H(p_t) - 143.45.
#
# Two predictions to check:
#  (1) if INR's marginals ARE the true marginals, then -logL_INR = sum_t H(p_t) exactly, and any
#      shortfall is sum_t KL(p_t||q_t^INR), which the Stein term put at 0.03 nats. This is a sharper
#      test than the Stein term because it is sensitive to the SHAPE and not only to the variance.
#  (2) sum_t H(p_t) - 143.45 should reproduce the residual 71.90.
suppressPackageStartupMessages({library(data.table)})

S <- as.data.table(readRDS("projects/eLife_2025/figures/data/digest/figure_3_digest_INR.rds")$scalar)
# the recordings are common to every member; check it rather than assume it
S2 <- as.data.table(readRDS("projects/eLife_2025/figures/data/digest/figure_3_digest_R.rds")$scalar)
setkey(S, simulation_index, sample_index); setkey(S2, simulation_index, sample_index)
cat("max |patch_current difference| between two members' digests:",
    max(abs(S$patch_current - S2$patch_current)), "\n\n")

# Vasicek m-spacing entropy, in nats
H_vasicek <- function(x, m = NULL) {
  n <- length(x); if (is.null(m)) m <- max(1, floor(sqrt(n)))
  x <- sort(x)
  i <- seq_len(n)
  hi <- pmin(i + m, n); lo <- pmax(i - m, 1)
  mean(log(n / (hi - lo) * (x[hi] - x[lo])))
}
H_gauss <- function(x) 0.5 * log(2 * pi * exp(1) * var(x))

ent <- S[, .(Hg = H_gauss(patch_current), Hv = H_vasicek(patch_current),
             sd = sd(patch_current)), by = sample_index][order(sample_index)]
sumHg <- sum(ent$Hg); sumHv <- sum(ent$Hv)

logL <- function(k) {
  s <- as.data.table(readRDS(sprintf("projects/eLife_2025/figures/data/digest/figure_3_digest_%s.rds", k))$scalar)
  mean(s[, .(t = sum(ll)), by = simulation_index]$t)
}
L_INR <- logL("INR"); L_IR <- logL("IR"); L_NR <- logL("NR"); L_LSE <- logL("LSE_av0")

cat(sprintf("sum_t H(p_t), Gaussian plug-in : %8.2f nats\n", sumHg))
cat(sprintf("sum_t H(p_t), Vasicek spacing  : %8.2f nats  (difference %+.2f = non-Gaussianity of the marginals)\n",
            sumHv, sumHv - sumHg))
cat(sprintf("-logL of INR                   : %8.2f nats\n", -L_INR))
cat(sprintf("  => sum_t KL(p_t || q_t^INR)  : %8.2f  (Gaussian)   %8.2f  (Vasicek)\n",
            -L_INR - sumHg, -L_INR - sumHv))
cat(sprintf("\nH(p) <= -logL of IR            : %8.2f nats\n", -L_IR))
cat(sprintf("TC = sum_t H(p_t) - H(p)       : >= %6.2f (Gaussian)  >= %6.2f (Vasicek)\n",
            sumHg + L_IR, sumHv + L_IR))
cat(sprintf("residual measured for INR      :    %6.2f\n", (-L_INR) - (-L_IR)))
cat(sprintf("residual measured for NR       :    %6.2f\n", 71.89))
cat("\nper-segment share of the total correlation (Gaussian estimate minus the IR reference is not\n",
    "separable per interval, so this is just where the marginal entropy sits):\n")
print(ent[, .(seg = fifelse(sample_index < 20, "pre", fifelse(sample_index < 60, "pulse", "wash"))
              )][, .N, by = seg])
print(ent[sample_index %in% c(1, 19, 20, 25, 40, 60, 70, 100)], digits = 4)
