# evidence_contributions_digest.R — per-tramo evidence contributions over
# iterations, for the three estimators. Run from the repo base:
#   Rscript projects/macroir_next/figures/evidence_contributions_digest.R
#
# POSITIONAL READ (see figure_2_ladder_digest.R): files written before
# 2026-09-06 label the Moment_statistics triples mean/var/count while the
# data is count/mean/variance, so position 29 is the WINDOWED MEAN of the
# trapezoid's per-tramo contribution (labeled var_plog_Evidence in v1).
# The telescopic columns were added later as plain scalars and are correctly
# labeled; they are taken positionally here only for uniformity.
#
# Row k carries the tramo between rungs k-1 and k; row 0 is empty.

suppressMessages(library(data.table))

RUN <- "projects/macroir_next/runs/e743fd5"
OUT <- "projects/macroir_next/figures/data/e743fd5"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

EX <- list.files(RUN, pattern = "^fig1_CCO_episodic_rep2_s910121_fit_CCO_.*__i_iter\\.csv$",
                 full.names = TRUE)[1]

COLS  <- c(1, 3, 5, 19, 29, 41, 43, 45, 47)
NAMES <- c("iter", "i_beta", "beta",
           "trap_inst", "trap_win",      # trapezoid: instantaneous, windowed
           "up_inst", "dn_inst",          # stepping stone, instantaneous
           "up_win", "dn_win")            # stepping stone, windowed

d <- fread(EX, select = COLS)
setnames(d, NAMES)
d <- d[is.finite(beta) & i_beta > 0]

ev <- sort(unique(d$iter))
d <- d[iter %in% ev[seq(1, length(ev), by = 8)]]

fwrite(d, file.path(OUT, "evidence_contrib.csv"))
cat("rows:", nrow(d), " tramos:", uniqueN(d$i_beta), "\n")
print(d[iter >= 24000, .(trap = round(median(trap_win), 3),
                         up = round(median(up_win), 3),
                         dn = round(median(dn_win), 3)), by = i_beta][order(i_beta)])
