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
# The per-temperature identity is the information equality of the TEMPERED
# TARGET p_beta ~ prior*L^beta (E_q[grad log q grad log q^T] = E_q[-hess log q],
# valid for any normalized density):
# its score is s_beta = beta*dlogL + dlogprior and the
# information equality gives Var_beta(s_beta) = E_beta[beta*GFI + FIM_prior]
# at EVERY rung once that rung equilibrated (up to the Gauss-Newton gap).
# Comparing the bare likelihood score to the bare GFI conflates
# prior-dominance with failure. The prior is Gaussian and diagonal in the
# transformed space, so dlogprior_j = -(theta_j - mu_j)/var_j and
# FIM_prior_jj = 1/var_j.
# theta comes from save_Parameter, whose cadence (every ~220 iters) differs
# from save_Score's (~180): only their COMMON events can be joined, ~8 events
# x 32 walkers per rung. A future recompile should emit theta (or dlogprior)
# in the score csv itself.
# Low-beta magnitudes reach ~1e180 and are occasionally emitted as inf/nan
# (which turns csv columns into character): coerce and keep finite rows.
sc <- fread(sc_f)
fm <- fread(fim_f)[i_par == j_par]
pp_f <- list.files(RUN, pattern = "^fig1_CCO_episodic_rep2_s910121_fit_CCO_.*__i_beta__i_walker__i_par\\.csv$",
                   full.names = TRUE)[1]
pp  <- fread(pp_f)
pri <- fread(file.path(RUN, "data", "scheme_CCO_prior_N1000.csv"))
sc[, dlogL := suppressWarnings(as.numeric(dlogL))]
fm[, gfi := suppressWarnings(as.numeric(gfi))]
common <- intersect(unique(sc$iter), unique(pp$iter))
common <- common[common >= 15000]   # post-equilibrium (logL flat well before)
j <- merge(merge(sc[iter %in% common, .(iter, i_beta, beta, i_walker, i_par, dlogL)],
                 fm[iter %in% common, .(iter, i_beta, i_walker, i_par, gfi)],
                 by = c("iter", "i_beta", "i_walker", "i_par")),
           pp[iter %in% common, .(iter, i_beta, i_walker, i_par, par_value)],
           by = c("iter", "i_beta", "i_walker", "i_par"))
j <- merge(j, pri[, .(i_par, mu = transformed_mean, v = transformed_variance)], by = "i_par")
j[, s_full := beta * dlogL - (par_value - mu) / v]
j[, info   := beta * gfi + 1 / v]
j <- j[is.finite(s_full) & is.finite(info)]
bart <- j[, .(beta = last(beta), n = .N, var_s = var(s_full), mean_info = mean(info)),
          by = .(i_beta, i_par)]
bart[, ratio := var_s / mean_info]
fwrite(bart, file.path(OUT, "bartlett_tempered.csv"))
cat("digests written to", OUT, "\n")
