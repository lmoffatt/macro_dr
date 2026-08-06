# figure_3_digest.R — turn each ~1 GB per-step dump into a ~15 MB digest, ONCE.
#
# WHY. Every figure_3 notebook fread()s the whole dump (25 columns, ~1 GB of text, several GB once
# parsed) and then throws away all but a handful of rows. On a 16 GB laptop that is the difference
# between "run it where figures/data/ resolves" and "run it here": four dumps at once do not fit,
# and even one is uncomfortable next to a browser. The dumps are immutable, so the extraction is
# pure and can be paid once.
#
# The row filter is done in awk, streaming, so peak memory is the FILTERED table and never the
# dump. Two facts make it exact rather than approximate, both verified on the data (2026-07-31):
#   - the writer emits every evolution row TWICE, and the two copies carry byte-identical `value`;
#   - the copies differ in their METADATA: the segment_index == 0 copy carries step_middle,
#     agonist and patch_current, the blank-segment copy leaves them empty.
# So `segment_index == 0` is a lossless dedup that also keeps the columns the blank copy drops.
# The notebooks' mean-over-copies collapse gets the same numbers; this is cheaper and keeps time.
#
# OUTPUT: figures/data/digest/figure_3_digest_<ALGO>.rds, a list of two data.tables
#   $scalar  one row per (simulation_index, sample_index): t (ms), agonist, patch_current,
#            r_std, ll (per-interval logL), y_mean, y_var
#   $param   one row per (simulation_index, sample_index, param_index): s = d logL_t/d theta,
#            dm = d y_mean/d theta, dv = d y_var/d theta, I = the per-step Gaussian Fisher diagonal
#            dm^2/y_var + 0.5 dv^2/y_var^2
# param_index: 0 = k_on, 1 = k_off, 2 = unitary_current, 3 = Current_Noise, 4 = Baseline,
# 5 = Num_ch_mean. ALL SIX are kept as of 2026-07-31, where the first version kept only the four
# identified ones (0, 1, 2, 5). The extra pair costs about half again in file size and saves a
# second pass over the ~7 GB of dumps the first time any notebook wants them, which supplement 3
# does: the instrumental-noise direction is the control direction of the score-bias check, the one
# Figure 4 reports as departing in no member. Existing consumers filter on param_index and are
# unaffected by the extra rows.
#
# RUN:  Rscript figure_3_digest.R [ALGO ...]      (default: all seven)
# from projects/eLife_2025/figures/paper_both/.

suppressPackageStartupMessages({library(data.table)})

# LSE_av0 added 2026-08-05: the un-averaged least-squares arm of figure_3_time.macroir. NMR is the
# pre-2026-08-01 name of INR and is kept only so an old dump still digests.
ALL   <- c("LSE", "LSE_av0", "NR", "NMR", "INR", "R", "MR", "VR", "IR")
args  <- commandArgs(trailingOnly = TRUE)
ALGOS <- if (length(args)) args else ALL
PIDX  <- c(0L, 1L, 2L, 3L, 4L, 5L)

IN  <- "../data"
OUT <- "../data/digest"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# Column positions in the dump header (fixed by the writer; asserted below against the real header):
#   2 simulation_index  3 sample_index  4 segment_index  9 step_middle  10 agonist
#  11 patch_current    12 component_path 18 calculus     21 param_index 25 value
AWK <- paste0(
  "awk -F, 'NR>2 && $3!=\"\" && $4==\"0\" && ",
  "($12==\"r_std\"||$12==\"logL\"||$12==\"y_mean\"||$12==\"y_var\") ",
  "{print $2\",\"$3\",\"$9\",\"$10\",\"$11\",\"$12\",\"$18\",\"$21\",\"$25}'")

check_header <- function(f) {
  h <- strsplit(readLines(f, n = 2)[2], ",", fixed = TRUE)[[1]]
  want <- c(simulation_index = 2, sample_index = 3, segment_index = 4, step_middle = 9,
            agonist = 10, patch_current = 11, component_path = 12, calculus = 18,
            param_index = 21, value = 25)
  bad <- names(want)[h[want] != names(want)]
  if (length(bad)) stop("column layout moved in ", f, ": ", paste(bad, collapse = ", "))
}

digest_one <- function(a) {
  f <- file.path(IN, sprintf("figure_3_time_dlik_%s.csv", a))
  if (!file.exists(f)) { message("skip ", a, ": no dump"); return(invisible(NULL)) }
  check_header(f)
  t0 <- Sys.time()
  d <- fread(cmd = paste(AWK, shQuote(f)), header = FALSE, showProgress = FALSE,
             col.names = c("simulation_index", "sample_index", "step_middle", "agonist",
                           "patch_current", "component_path", "calculus", "param_index", "value"))

  # scalars: the per-interval primitives, one column each, plus the time map and the observation
  prim <- d[calculus == "primitive"]
  scal <- dcast(prim, simulation_index + sample_index ~ component_path, value.var = "value")
  meta <- unique(d[, .(simulation_index, sample_index, t = step_middle * 1000,
                       agonist, patch_current)])
  scal <- merge(meta, scal, by = c("simulation_index", "sample_index"))
  setnames(scal, "logL", "ll", skip_absent = TRUE)

  # per-parameter derivatives, and the Gaussian Fisher diagonal built from them. y_var comes from
  # the scalar side, so I is defined by exactly one y_var per interval and cannot drift between
  # the notebooks that rebuild it.
  der <- d[calculus == "derivative" & param_index %in% PIDX]
  parw <- dcast(der, simulation_index + sample_index + param_index ~ component_path,
                value.var = "value")
  setnames(parw, c("logL", "y_mean", "y_var"), c("s", "dm", "dv"), skip_absent = TRUE)
  parw <- merge(parw, scal[, .(simulation_index, sample_index, y_var)],
                by = c("simulation_index", "sample_index"))
  parw[, I := dm^2 / y_var + 0.5 * dv^2 / y_var^2][, y_var := NULL]

  out <- file.path(OUT, sprintf("figure_3_digest_%s.rds", a))
  saveRDS(list(scalar = scal, param = parw), out, compress = "xz")
  message(sprintf("%-4s scalar %s x %d | param %s x %d | %s | %.0f s", a,
                  format(nrow(scal), big.mark = ","), ncol(scal),
                  format(nrow(parw), big.mark = ","), ncol(parw),
                  format(structure(file.size(out), class = "object_size"), units = "MB"),
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  rm(d, prim, der, scal, parw, meta); invisible(gc(verbose = FALSE))
}

for (a in ALGOS) digest_one(a)
