# The accumulated score bias of every member of Figure 3, with a bootstrap over recordings.
#
# CITED BY 04_results.tex, the paragraph beginning "The mean of the score is what says whether the
# estimate will be displaced": the split it reports ("the four members centred there sit within a
# quarter of a standard error of zero, and the displaced four miss it on the unitary current and
# the channel number by half a standard error to two") is this table's z column.
#
# Run from the repository root:  Rscript papers/1_method/decisions/recompute/fig3_score_bias.R
# PROMOTED 2026-08-26 from tmp/fig3_bias_boot.R, which is not committed. Same computation as the
# accumulated block of figures/paper_both/figure_3.Rmd, run standalone and saved, so the numbers
# quoted in the body survive without re-rendering the figure or holding the 286 MB dumps.
#
# Per recording i:  g_i = sum_t s_t                 (vector over parameters)
#                   F_i = sum_t ( dm dm' / v_t + 0.5 dv dv' / v_t^2 )
# then  b = mean(F_i)^-1 mean(g_i)  is the first-order bias and  b_a / sqrt((F^-1)_aa)  is that
# displacement in units of the MARGINAL standard error, the one the member actually reports.
# Compared against the diagonal-only reading  mean(g_a)/sqrt(F_aa), which treats every other
# parameter as known and is what summing row D of the figure gives; the gap is not small, on N_ch
# it turns -0.34 into +0.82 for R, sign included, the channel number and the unitary current being
# strongly correlated.
suppressPackageStartupMessages({library(data.table)})
DIG  <- "projects/eLife_2025/figures/data/digest"
OUT  <- "papers/1_method/decisions/recompute/fig3_score_bias.csv"
DTOK <- c(LSE="LSE_av0", ILSE="LSE", NR="NR", INR="INR", R="R", MR="MR", VR="VR", IR="IR")
PN <- c("0"="k_on", "1"="k_off", "2"="i", "5"="N_ch")
B <- 400
set.seed(20260812)
res <- list()

for (k in names(DTOK)) {
  dg <- readRDS(file.path(DIG, sprintf("figure_3_digest_%s.rds", DTOK[[k]])))
  S  <- as.data.table(dg$scalar)[, .(simulation_index, sample_index, y_var)]
  P  <- as.data.table(dg$param)[param_index %in% c(0L,1L,2L,5L)]
  if (k %in% c("LSE","ILSE")) P <- P[param_index != 5L]   # N_ch degenerate with i there
  P  <- merge(P, S, by = c("simulation_index","sample_index"))[is.finite(y_var) & y_var > 0]
  ids <- sort(unique(P$param_index)); np <- length(ids)

  W  <- dcast(P, simulation_index + sample_index + y_var ~ param_index,
              value.var = c("s","dm","dv"))
  sims <- sort(unique(W$simulation_index)); n <- length(sims)
  DM <- as.matrix(W[, paste0("dm_", ids), with = FALSE])
  DV <- as.matrix(W[, paste0("dv_", ids), with = FALSE])
  SS <- as.matrix(W[, paste0("s_",  ids), with = FALSE])
  vv <- W$y_var; grp <- match(W$simulation_index, sims)

  # per recording: accumulated score vector and accumulated information matrix
  gi <- rowsum(SS, grp)                                   # n x np
  Fi <- array(0, c(n, np, np))
  for (a in seq_len(np)) for (b2 in seq_len(np))
    Fi[, a, b2] <- rowsum(DM[, a] * DM[, b2] / vv + 0.5 * DV[, a] * DV[, b2] / vv^2, grp)

  stat <- function(j) {
    g <- colMeans(gi[j, , drop = FALSE])
    Fm <- apply(Fi[j, , , drop = FALSE], c(2, 3), mean)
    Finv <- solve(Fm)
    bvec <- as.vector(Finv %*% g)
    list(z = bvec / sqrt(diag(Finv)), b = bvec, zd = g / sqrt(diag(Fm)))
  }
  o <- stat(seq_len(n))
  rep <- vapply(seq_len(B), function(b3) stat(sample.int(n, n, TRUE))$z, numeric(np))
  lo <- apply(rep, 1, quantile, .025); hi <- apply(rep, 1, quantile, .975)
  res[[k]] <- data.table(mem = k, param = PN[as.character(ids)],
                         z_diag = o$zd, z = o$z, lo = lo, hi = hi, bias_log10 = o$b)
  message("done ", k)
}
R <- rbindlist(res); setnames(R, "mem", "key")
R[, covers0 := lo < 0 & hi > 0]
print(as.data.frame(R), digits = 3, row.names = FALSE)
fwrite(R, OUT)
