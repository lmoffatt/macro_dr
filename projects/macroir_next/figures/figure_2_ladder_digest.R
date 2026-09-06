# figure_2_ladder_digest.R — digests for the thermodynamic-ladder figure.
# Run from the repo base:
#   Rscript projects/macroir_next/figures/figure_2_ladder_digest.R
#
# READS BY COLUMN POSITION, not by name, on purpose: in files written before
# 2026-09-06 every Moment_statistics triple is labeled mean/var/count while
# the data is count/mean/variance (moment_statistics.h:1090; see the audit
# note). The fix renamed the titles, which leaves POSITIONS unchanged, so a
# positional read is correct for both old and new campaigns. Position 26 is
# the windowed mean of logL (labeled var_logL in v1 files, and clashing with
# the varLik column 12 under that name).

suppressMessages(library(data.table))

RUN <- "projects/macroir_next/runs/e743fd5"
OUT <- "projects/macroir_next/figures/data/e743fd5"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## The exemplary run: it both adapts the rung POSITIONS and, at iteration
## 3332, INSERTS a rung (16 -> 17), so one panel shows the two mechanisms
## the dts sampler applies to the ladder.
EX <- list.files(RUN, pattern = "^fig1_CCO_episodic_rep4_s910141_fit_CCO_.*__i_iter\\.csv$",
                 full.names = TRUE)[1]

COLS  <- c(1, 3, 4, 5, 19, 26, 37, 39, 40, 41)
NAMES <- c("iter", "i_beta", "num_beta", "beta", "plog_Evidence", "mean_logL",
           "acc_emcee", "acc_jump", "dBeta_dlogL", "plog_Evidence_ss")

d <- fread(EX, select = COLS)
setnames(d, NAMES)
d <- d[is.finite(beta)]

## (a) ladder positions over iterations -----------------------------------
## thinned for plotting; every rung keeps its own line
ev <- sort(unique(d$iter))
keep <- ev[seq(1, length(ev), by = 8)]
fwrite(d[iter %in% keep, .(iter, i_beta, num_beta, beta)],
       file.path(OUT, "ladder_trace.csv"))

## (b,c,d) the settled ladder: thermodynamic path, per-rung evidence
## contribution, mixing. Averaged over the tail so the numbers are the
## equilibrium ones, not one noisy event.
tail_stats <- d[iter >= 24000,
                .(beta = median(beta), mean_logL = median(mean_logL),
                  plog = median(plog_Evidence), plog_ss = median(plog_Evidence_ss),
                  acc_emcee = mean(acc_emcee), acc_jump = mean(acc_jump),
                  dBeta_dlogL = median(dBeta_dlogL), n = .N),
                by = i_beta][order(i_beta)]
fwrite(tail_stats, file.path(OUT, "ladder_settled.csv"))

## ladder size across the whole campaign (how often the sampler resized)
sizes <- rbindlist(lapply(list.files(RUN, pattern = "^fig1_.*__i_iter\\.csv$",
                                     full.names = TRUE), function(p) {
    x <- fread(p, select = c(1, 4))
    setnames(x, c("iter", "num_beta"))
    # nrow(x), not .N: outside a data.table call .N evaluates to 0 and the
    # column comes back empty, silently yielding a zero-row result
    data.table(run = sub("_scheme.*", "", basename(p)),
               first = x$num_beta[1], last = x$num_beta[nrow(x)],
               n_sizes = uniqueN(x$num_beta))
}))
fwrite(sizes, file.path(OUT, "ladder_sizes.csv"))
cat("runs whose ladder changed size:", sizes[n_sizes > 1, .N], "of", nrow(sizes), "\n")
cat("digests written to", OUT, "\n")
