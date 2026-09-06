# figure_1_digest.R — reduce the figure_1 campaign's raw csvs (~6 GB of
# __i_iter plus score/fim) to the small digests figure_1.Rmd consumes.
# Run from the repo base:  Rscript projects/macroir_next/figures/figure_1_digest.R
#
# Estimator per fit (settled 2026-09-06): MEDIAN over the tail (last 20% of
# iterations) of the per-window bracket midpoints (mean_log_Evidence_ss +
# _ss_dn)/2 at beta = 1. Median, not mean: rare windows are poisoned by
# heavy-tailed weights at the cold rungs (ss_ess_dn collapsing to ~10) and
# the mean swallows them (seed 910181 taught this).

suppressMessages(library(data.table))

RUN <- "projects/macroir_next/runs/e743fd5"
OUT <- "projects/macroir_next/figures/data/e743fd5"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

TAIL_FROM <- 24000L   # last 20% of 30000
MIN_LAST  <- 29000L   # a fit counts as finished past this iteration

iter_files <- list.files(RUN, pattern = "^fig1_.*__i_iter\\.csv$", full.names = TRUE)
meta <- data.table(path = iter_files)
meta[, c("truth", "prot", "rep", "seed", "fit") :=
    as.data.table(t(vapply(basename(path), function(nm) {
        m <- regmatches(nm, regexec(
            "^fig1_(CCO|COC)_(stationary|episodic)_rep(\\d+)_s(\\d+)_fit_(CCO|COC)_", nm))[[1]]
        m[2:6]
    }, character(5))))]

# ---- per-fit evidence summary ------------------------------------------------
summ <- rbindlist(lapply(seq_len(nrow(meta)), function(i) {
    dt <- fread(meta$path[i], select = c("iter", "beta",
                                         "mean_log_Evidence_ss", "mean_log_Evidence_ss_dn",
                                         "log_Evidence", "logL"))
    last <- max(dt$iter)
    tl <- dt[iter >= TAIL_FROM & beta == 1]
    if (nrow(tl) == 0L || last < MIN_LAST) return(NULL)
    data.table(meta[i, .(truth, prot, rep, seed, fit)],
               last_iter = last,
               lnZ    = median((tl$mean_log_Evidence_ss + tl$mean_log_Evidence_ss_dn) / 2),
               lnZ_up = median(tl$mean_log_Evidence_ss),
               lnZ_dn = median(tl$mean_log_Evidence_ss_dn),
               gap    = median(tl$mean_log_Evidence_ss_dn - tl$mean_log_Evidence_ss),
               trap   = median(tl$log_Evidence),
               logL1  = median(tl$logL))
}))
fwrite(summ, file.path(OUT, "evidence_by_fit.csv"))

conf <- dcast(summ, truth + prot + seed ~ fit, value.var = c("lnZ", "gap"))
conf <- conf[!is.na(lnZ_CCO) & !is.na(lnZ_COC)]
conf[, dlnZ := lnZ_CCO - lnZ_COC]
fwrite(conf, file.path(OUT, "confusion.csv"))
cat("confusion.csv:", nrow(conf), "complete pairs\n")

# ---- convergence trace of one exemplary run ---------------------------------
EX <- meta[truth == "CCO" & prot == "episodic" & seed == "910121" & fit == "CCO", path]
tr <- fread(EX[1], select = c("iter", "beta", "logL", "log_Evidence",
                              "mean_log_Evidence_ss", "mean_log_Evidence_ss_dn",
                              "log_Evidence_ss", "log_Evidence_ss_dn"))
tr <- tr[beta == 1][seq(1, .N, by = 4)]   # rows saved every 4 iters -> plot every 16
fwrite(tr, file.path(OUT, "trace_example.csv"))

# ---- example recordings (one per protocol, same truth CCO) -------------------
rec <- rbindlist(lapply(c(episodic = "fig1_CCO_episodic_rep2_s910121",
                          stationary = "fig1_CCO_stationary_rep1_s910011"), function(lb) {
    f <- list.files(RUN, pattern = paste0("^", lb, "_scheme_CCO_.*_simulation\\.csv$"),
                    full.names = TRUE)[1]
    d <- fread(f)
    setnames(d, c("i_step", "current"))
    d[, label := lb]
    d
}), idcol = "protocol")
# the evidence stage masks the stationary equilibration (steps n1..n1+n2-1);
# reproduce it here so the figure shows what the likelihood saw
rec[protocol == "stationary" & i_step >= 10 & i_step < 210, current := NA]
rec[, time_tau := i_step * 0.1]           # interval = 0.1 tau
fwrite(rec, file.path(OUT, "recordings.csv"))

# ---- score/FIM per temperature (Bartlett-style), one box fit -----------------
sc_f  <- list.files(RUN, pattern = "^fig1_CCO_episodic_rep2_s910121_fit_CCO_.*_score\\.csv$",
                    full.names = TRUE)[1]
fim_f <- list.files(RUN, pattern = "^fig1_CCO_episodic_rep2_s910121_fit_CCO_.*_fim\\.csv$",
                    full.names = TRUE)[1]
# At low beta the walkers roam the prior, where dlogL and the FIM are
# astronomically large and occasionally emitted as inf/nan (which also turns
# the csv column into character): coerce to numeric and keep finite rows only,
# so each (rung, parameter) statistic is over its finite samples.
sc <- fread(sc_f)[iter >= TAIL_FROM]
fm <- fread(fim_f)[iter >= TAIL_FROM & i_par == j_par]
sc[, dlogL := suppressWarnings(as.numeric(dlogL))]
fm[, gfi := suppressWarnings(as.numeric(gfi))]
sc <- sc[is.finite(dlogL)]
fm <- fm[is.finite(gfi)]
bart <- merge(
    sc[, .(beta = last(beta), var_score = var(dlogL), n = .N), by = .(i_beta, i_par)],
    fm[, .(mean_gfi = mean(gfi), n_gfi = .N), by = .(i_beta, i_par)],
    by = c("i_beta", "i_par"))
bart <- bart[is.finite(var_score) & is.finite(mean_gfi) & mean_gfi > 0]
fwrite(bart, file.path(OUT, "bartlett.csv"))
cat("digests written to", OUT, "\n")
