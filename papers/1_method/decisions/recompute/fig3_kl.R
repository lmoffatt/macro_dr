# fig3_kl.R -- does the likelihood ladder of Figure 3 decompose into the two defects the figure
# measures? Feeds the commented derivation kept under the log-likelihood paragraph of
# 04_results.tex (the two-term split of the divergence, cut from the body 2026-08-12) and the
# figure_3_guion.md section that keeps it in full.
#
# Run from the repository root:  Rscript papers/1_method/decisions/recompute/fig3_kl.R
# PROMOTED 2026-08-26 from tmp/, which is gitignored; only the paths changed.
# Does the likelihood ladder decompose into the two defects the figure measures?
#
# (i) E_p[log q_m] = -H(p) - KL(p || q_m).  Averaged over recordings from the EXACT simulator, the
#     mean total logL estimates -H(p) - KL up to a constant common to every member, so DIFFERENCES
#     in the ladder are differences in KL to the true process. The ladder is a truth criterion after
#     all, up to that one shared unknown.
#
# (ii) For a member that treats the intervals as independent, q = prod_t q_t, exactly
#          KL(p || prod q_t) = sum_t KL(p_t || q_t)  +  [ sum_t H(p_t) - H(p) ],
#      the second term being the total correlation of the true observation sequence: a floor no
#      amount of marginal calibration can pay off.
#
# (iii) For a Gaussian marginal whose variance is wrong by a factor r_t = E[r_std^2]_t (mean right),
#          KL(p_t || q_t) = 0.5 * (r_t - 1 - log r_t),
#      the Stein / Itakura-Saito loss, which is exactly what row B measures interval by interval.
#
# So: predicted marginal cost = sum_t 0.5 (r_t - 1 - log r_t) from row B alone, and whatever is left
# of the observed gap is the dependence term that row G is about.
suppressPackageStartupMessages({library(data.table)})
DTOK <- c(LSE="LSE_av0", ILSE="LSE", NR="NR", INR="INR", R="R", MR="MR", VR="VR", IR="IR")

out <- rbindlist(lapply(names(DTOK), function(k) {
  s <- as.data.table(readRDS(sprintf("projects/eLife_2025/figures/data/digest/figure_3_digest_%s.rds", DTOK[[k]]))$scalar)
  s <- s[is.finite(r_std)]
  rt <- s[, .(r = mean(r_std^2)), by = sample_index]
  tot <- s[, .(tot = sum(ll)), by = simulation_index]
  data.table(mem = k, logL = mean(tot$tot), se = sd(tot$tot)/sqrt(nrow(tot)),
             stein = sum(0.5 * (rt$r - 1 - log(rt$r))),
             stein_worst = max(0.5 * (rt$r - 1 - log(rt$r))),
             worst_at = rt$sample_index[which.max(0.5 * (rt$r - 1 - log(rt$r)))])
}))
setnames(out, "mem", "key")
ref <- out[key == "IR", logL]
out[, `:=`(gap = ref - logL)]
out[, `:=`(dependence = gap - stein)]
print(as.data.frame(out[, .(key, logL, gap, stein, dependence, stein_worst, worst_at)]),
      digits = 4, row.names = FALSE)

cat("\nTwo checks the decomposition has to pass:\n")
cat(sprintf("  NR - INR observed = %.2f nats; row B alone predicts %.2f (both assume independence,\n",
    out[key=="NR", logL] * -1 + out[key=="INR", logL], out[key=="NR", stein] - out[key=="INR", stein]))
cat("    so the dependence term cancels in the difference and row B has to carry all of it)\n")
cat(sprintf("  INR - IR observed = %.2f nats, of which row B predicts %.2f: a gap that is pure\n",
    out[key=="INR", logL] - ref, out[key=="INR", stein] - out[key=="IR", stein]))
cat("    dependence, invisible in the residual variance.\n")
