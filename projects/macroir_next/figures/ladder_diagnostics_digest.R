# ladder_diagnostics_digest.R — per-rung time series for the ladder
# diagnostics figure. Run from the repo base:
#   Rscript projects/macroir_next/figures/ladder_diagnostics_digest.R
#
# POSITIONAL READ, as in figure_2_ladder_digest.R: files written before
# 2026-09-06 label every Moment_statistics triple mean/var/count while the
# data is count/mean/variance. Position 30 is therefore the windowed
# VARIANCE of the per-tramo evidence contribution (labeled
# count_plog_Evidence in v1 files).

suppressMessages(library(data.table))

RUN <- "projects/macroir_next/runs/e743fd5"
OUT <- "projects/macroir_next/figures/data/e743fd5"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

EX <- list.files(RUN, pattern = "^fig1_CCO_episodic_rep2_s910121_fit_CCO_.*__i_iter\\.csv$",
                 full.names = TRUE)[1]

COLS  <- c(1, 3, 5, 30, 34, 37, 39, 49, 50, 51)
NAMES <- c("iter", "i_beta", "beta", "var_plog_w", "deltaEvidence_variance",
           "acc_within", "acc_between", "ess_up", "ess_dn", "ss_count")

d <- fread(EX, select = COLS)
setnames(d, NAMES)
d <- d[is.finite(beta)]

## thin for plotting: every 8th saved event (events are every 4 iterations)
ev <- sort(unique(d$iter))
d <- d[iter %in% ev[seq(1, length(ev), by = 8)]]

fwrite(d, file.path(OUT, "ladder_diag.csv"))
cat("rows:", nrow(d), " rungs:", uniqueN(d$i_beta),
    " iter range:", range(d$iter), "\n")
cat("ess_up range:", range(d$ess_up[is.finite(d$ess_up)]), "\n")
cat("ss_count range:", range(d$ss_count), "\n")
