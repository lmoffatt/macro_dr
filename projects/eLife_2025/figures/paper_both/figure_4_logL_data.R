## figure_4_logL_data.R — the log-likelihood product for the design plane.
##
## WHAT IT IS. `Probit_statistics_Moment_statistics_Sum_logL` at the THETA_SIM anchor: the mean over
## the cell's 10,000 recordings of the log-likelihood each member assigns to a recording, summed
## over the recording's intervals. Every recording is 10 tau long in every cell of this grid
## (verified against Sum_r2_std, whose null is the interval count 10/Delta_tilde and which reads
## 999.5 / 499.8 / 199.9 / ... against 1000 / 500 / 200 / ...), so the SUM is already a quantity per
## unit time and needs no normalization. Dividing by the interval count would multiply by
## Delta_tilde/10, i.e. inject a 100-fold ramp along the horizontal axis that is nothing but the
## sample count. Measured both ways 2026-08-12: along the interval axis the median contrast varies
## 15x to 14519x per interval and 2.7x to 145x raw, so the raw form is flatter in all nine
## contrasts.
##
## WHY THETA_SIM AND NOT THETA_POOL, which is the opposite choice from the distortion half. The
## comparison is only meaningful if every member is evaluated at the SAME parameter: then
##     E[log p_A] - E[log p_B] = KL(true || p_B) - KL(true || p_A),
## a difference of Kullback-Leibler divergence rates, positive when A is the closer density. At each
## member's own optimum the difference would mix density quality with parameter fitting and would
## penalise the four-parameter least-squares arms for holding two parameters fixed.
##
## THE FILTER, and it is the one that is easy to get wrong. Each (member, cell, interval) emits SIX
## rows with statistic == "mean": one with probit == "mean", the point estimate, and five with
## probit == "quantile" at levels .025 / .16 / .5 / .84 / .975, the percentile bootstrap. Filtering
## on `statistic == "mean"` alone returns all six and they differ by about 0.05% of the sum, which
## on a sum of -7000 nats is +-4 nats, larger than several of the contrasts this file is built to
## measure. BOTH conditions are required: probit == "mean" AND statistic == "mean".
##
## THE ERROR BAR, and why it is not optional here. Each member is scored on its OWN ensemble of
## 10,000 simulated recordings (the members are separate runs), so a contrast is a difference of two
## independent means and its standard error is sqrt(v_A/n_A + v_B/n_B) with v the emitted variance
## over recordings. The recursion contrasts are of order 2 nats and the averaging contrasts of order
## 0.03, so whether "averaging buys nothing" is a measurement or a null result depends entirely on
## this number. It is carried through to the drawing, which collapses unresolved cells to zero.
##
## DATA_DIRS is figure_4.Rmd's five-directory search path, first hit wins, NOT a glob over
## figures/data: 82b956f also carries logL for IR, MR and R and is not on the figure's path.
##
## OUTPUT. ../figure_4_source_data/figure_4_source_data_logL.csv, same format as its five siblings:
## line 1 the provenance stamp, line 2 the header. Not yet folded into figure_4_data.R's single
## sweep, deliberately: doing so changes that file's mtime and forces a rebuild of all five existing
## products. Fold it in at the next rebuild those need anyway.

suppressPackageStartupMessages({library(dplyr)})

LOGL_DIRS <- c("../data/1c2ae6f", "../data/87889e6", "../data/0ffbda7", "../data/1f7138b",
               "../data/a202e03")
LOGL_MEMBERS <- c("nonlinearsqr_g", "nonlinearsqr", "macro_NR", "macro_INR", "macro_R",
                  "macro_MR", "macro_VR", "macro_IR")
LOGL_PATH <- "Probit_statistics_Moment_statistics_Sum_logL"
LOGL_OUT  <- "../figure_4_source_data/figure_4_source_data_logL.csv"

.logl_pre  <- function(a) if (startsWith(a, "nonlinearsqr")) "figure_3_LSE_" else "figure_3_G_"
.logl_read <- if (requireNamespace("data.table", quietly = TRUE))
  function(f) as.data.frame(data.table::fread(f, skip = 1, showProgress = FALSE)) else
  function(f) utils::read.csv(f, skip = 1, stringsAsFactors = FALSE)

# every battery_sim_G file on the figure's path, first directory wins per (member, N_ch, noise)
logl_files <- function() {
  out <- list()
  for (a in LOGL_MEMBERS) {
    pat <- sprintf("^%snch_([0-9]+)_nsim_10000_%s_noise_([0-9.eE+]+)_battery_sim_G\\.csv$",
                   .logl_pre(a), a)
    seen <- character(0)
    for (d in LOGL_DIRS) {
      for (b in list.files(d, pattern = pat)) {
        m <- regmatches(b, regexec(pat, b))[[1]]
        key <- paste(m[2], m[3])
        if (key %in% seen) next                       # first hit wins, same rule as findf()
        seen <- c(seen, key)
        out[[length(out) + 1]] <- data.frame(algo = a, nch = as.integer(m[2]), z = m[3],
                                             path = file.path(d, b), dir = basename(d),
                                             stringsAsFactors = FALSE)
      }
    }
  }
  bind_rows(out)
}

fig4_logL_build <- function(force = FALSE) {
  fl <- logl_files()
  stamp <- sprintf("# figure_4 logL source data | inputs=%d bytes=%.0f newest=%.0f | code=%.0f",
                   nrow(fl), sum(file.info(fl$path)$size), max(as.numeric(file.info(fl$path)$mtime)),
                   as.numeric(file.info("figure_4_logL_data.R")$mtime))
  if (!force && file.exists(LOGL_OUT) && readLines(LOGL_OUT, n = 1) == stamp) {
    cat("figure_4 logL: source data fresh, reading\n")
    return(read.csv(LOGL_OUT, skip = 1, stringsAsFactors = FALSE))
  }
  cat("figure_4 logL: rebuilding from", nrow(fl), "battery files\n")
  d <- bind_rows(lapply(seq_len(nrow(fl)), function(i) {
    x <- .logl_read(fl$path[i])
    x <- x[x$component_path == LOGL_PATH & x$probit == "mean" &
             x$statistic %in% c("mean", "variance", "count"), ]
    if (!nrow(x)) return(NULL)
    x <- x[, c("statistic", "interval_in_tau", "value")]
    x <- tidyr::pivot_wider(x, names_from = statistic, values_from = value)
    transform(x, algo = fl$algo[i], Num_ch = fl$nch[i], noise = fl$z[i], dir = fl$dir[i])
  }))
  d <- d %>% rename(logL = mean, var_logL = variance, n_sim = count) %>%
    select(algo, Num_ch, noise, interval_in_tau, logL, var_logL, n_sim, dir) %>%
    arrange(algo, Num_ch, as.numeric(noise), interval_in_tau)
  dir.create(dirname(LOGL_OUT), showWarnings = FALSE, recursive = TRUE)
  writeLines(stamp, LOGL_OUT)
  suppressWarnings(write.table(d, LOGL_OUT, sep = ",", row.names = FALSE, qmethod = "double",
                               append = TRUE))
  cat("figure_4 logL: wrote", nrow(d), "rows to", LOGL_OUT, "\n")
  d
}
