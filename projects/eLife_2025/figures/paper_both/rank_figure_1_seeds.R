#!/usr/bin/env Rscript
# Rank the seeds swept by ops/local/sweep_figure_1_seeds.sh, so that the trace and the window
# for Figure 1 are chosen on a stated statistic instead of by eye over twenty PDFs.
#
# Two criteria, both per (seed, window), because a seed that is vivid at the concentration
# jump can be flat four milliseconds later:
#
#   update   how far the Bayesian step moves the end-of-window open probability,
#            |p_open(post) - p_open(prior)|, for IR. This is the row-C content of the figure.
#   contrast |update(IR) - update(R)|. Without this the trace can illustrate a filter while
#            failing to illustrate the ladder, which is the point of the figure.
#
# The quantities are read from the per-interval diagnostic dumps, same-frame prior/posterior
# pairs (the _y0 / _y1 suffix is before/after the observation):
#   IR  Algo_State_Dynamic.P_mean_0t_y0 / _y1, the 2x2 boundary table (row = start, col = end,
#       per papers/_program/notation_map.md), marginalised over the start index.
#   R   Algo_State_Dynamic.P_mean_t15_y0 / _y1, already an end-frame vector.
# The open state is identified from the data, as the one whose occupancy tracks patch_current,
# rather than assumed from the scheme's state order.
#
# Usage: Rscript rank_figure_1_seeds.R [outroot]      (default ../../ops/local/figure_1_seeds)

args <- commandArgs(trailingOnly = TRUE)
outroot <- if (length(args) >= 1) args[1] else "../../ops/local/figure_1_seeds"
stopifnot(dir.exists(outroot))

# Sample indices are the acquisition windows; a figure panel shows a consecutive pair.
WINDOWS <- list(c(0, 1), c(1, 2), c(2, 3), c(3, 4), c(4, 5))

read_dump <- function(dir, algo) {
  f <- file.path(dir, "figures", "data",
                 paste0("figure_1_likelihood_diagnostic_", algo, ".csv"))
  if (!file.exists(f)) return(NULL)
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  d <- d[d$scope == "evolution", ]
  # Every evolution row is emitted twice, once with a blank segment_index and once with
  # segment_index = 0. Summing without this filter doubles every quantity below.
  d[!is.na(d$segment_index) & d$segment_index == 0, ]
}

# Occupancy vector per sample_index for one component, marginalised over the start index when
# the component is a boundary table. Returns a data.frame(sample_index, state, p).
occupancy <- function(d, component) {
  x <- d[d$component_path == component, ]
  if (!nrow(x)) return(NULL)
  # End-state index: the column index for a boundary table (row = start), and whichever of the
  # two index columns actually varies for a plain vector.
  idx <- if (length(unique(x$value_row)) > 1 && length(unique(x$value_col)) > 1) {
    x$value_col
  } else if (length(unique(x$value_col)) > 1) {
    x$value_col
  } else {
    x$value_row
  }
  agg <- stats::aggregate(list(p = x$value),
                          by = list(sample_index = x$sample_index, state = idx),
                          FUN = sum)
  agg[order(agg$sample_index, agg$state), ]
}

# The open state is the one whose occupancy correlates with the recorded current.
open_state <- function(d, occ) {
  cur <- stats::aggregate(list(patch_current = d$patch_current),
                          by = list(sample_index = d$sample_index),
                          FUN = function(v) mean(v, na.rm = TRUE))
  states <- sort(unique(occ$state))
  r <- vapply(states, function(s) {
    z <- merge(occ[occ$state == s, ], cur, by = "sample_index")
    if (nrow(z) < 3 || stats::sd(z$p) == 0 || stats::sd(z$patch_current) == 0) return(NA_real_)
    stats::cor(z$p, z$patch_current)
  }, numeric(1))
  if (all(is.na(r))) return(states[length(states)])
  states[which.max(r)]
}

# |posterior - prior| of the open probability, per sample_index, for one algorithm.
update_size <- function(d, comp_prior, comp_post) {
  if (is.null(d)) return(NULL)
  pri <- occupancy(d, comp_prior)
  pos <- occupancy(d, comp_post)
  if (is.null(pri) || is.null(pos)) return(NULL)
  st <- open_state(d, pos)
  a <- pri[pri$state == st, c("sample_index", "p")]
  b <- pos[pos$state == st, c("sample_index", "p")]
  z <- merge(a, b, by = "sample_index", suffixes = c("_prior", "_post"))
  # Signed: the sign is itself a criterion. A pair of windows whose updates go in
  # opposite directions shows that the filter follows the sign of the innovation,
  # which two same-signed steps cannot demonstrate.
  data.frame(sample_index = z$sample_index, upd = z$p_post - z$p_prior)
}

# Per-interval log-likelihood. figure_1_panels.R plots row D as cumsum(logL), so the dumped
# quantity is per interval and summing it over a window is the window's total.
logL_by_sample <- function(d) {
  if (is.null(d)) return(NULL)
  x <- d[d$component_path == "logL", ]
  if (!nrow(x)) return(NULL)
  stats::aggregate(list(logL = x$value), by = list(sample_index = x$sample_index), FUN = sum)
}

seed_dirs <- list.files(outroot, pattern = "^seed_[0-9]+$", full.names = TRUE)
if (!length(seed_dirs)) stop("no seed_* directories under ", outroot)

# Optional second argument: a directory holding figures/data for the trace currently in the
# paper, so the incumbent appears in the table as seed 0 and the sweep has something to beat.
if (length(args) >= 2 && dir.exists(file.path(args[2], "figures", "data"))) {
  seed_dirs <- c(seed_dirs, args[2])
}

rows <- list()
for (dir in seed_dirs) {
  seed <- suppressWarnings(as.integer(sub("^seed_", "", basename(dir))))
  if (is.na(seed)) seed <- 0L        # the incumbent trace passed as the reference dir
  d_ir <- read_dump(dir, "IR")
  d_r  <- read_dump(dir, "R")
  u_ir <- update_size(d_ir, "Algo_State_Dynamic.P_mean_0t_y0",  "Algo_State_Dynamic.P_mean_0t_y1")
  u_r  <- update_size(d_r,  "Algo_State_Dynamic.P_mean_t15_y0", "Algo_State_Dynamic.P_mean_t15_y1")
  if (is.null(u_ir) || is.null(u_r)) {
    message("seed ", seed, ": incomplete dumps, skipped")
    next
  }
  u <- merge(u_ir, u_r, by = "sample_index", suffixes = c("_IR", "_R"))

  # Row D: the figure must not show least squares beating the boundary-conditioned filter on
  # its own trace. That reading is available to any reader and the paper does not claim it.
  L <- lapply(c(LSE = "LSE", NR = "NR", R = "R", IR = "IR"), function(a) {
    z <- logL_by_sample(if (a == "IR") d_ir else if (a == "R") d_r else read_dump(dir, a))
    if (is.null(z)) NULL else stats::setNames(z$logL, z$sample_index)
  })

  for (w in WINDOWS) {
    sel <- u[u$sample_index %in% w, ]
    if (nrow(sel) < length(w)) next
    Lw <- vapply(L, function(v) if (is.null(v)) NA_real_ else sum(v[as.character(w)]), numeric(1))
    # balance is the primary criterion: the weakest of the four update arrows the
    # figure actually draws in row C (R and IR, in each of the two windows). A sum
    # would let a seed with one large step and one dead step win, which is exactly
    # the failure mode: a panel where the Bayes update is invisible teaches nothing.
    sel <- sel[order(sel$sample_index), ]
    rows[[length(rows) + 1]] <- data.frame(
      seed     = seed,
      window   = paste(w, collapse = "-"),
      balance  = min(abs(c(sel$upd_IR, sel$upd_R))),
      update   = sum(abs(sel$upd_IR)),
      contrast = sum(abs(sel$upd_IR - sel$upd_R)),
      flip_IR  = sign(sel$upd_IR[1]) != sign(sel$upd_IR[2]),
      flip_R   = sign(sel$upd_R[1])  != sign(sel$upd_R[2]),
      IR_1 = sel$upd_IR[1], IR_2 = sel$upd_IR[2],
      R_1  = sel$upd_R[1],  R_2  = sel$upd_R[2],
      IR_top = !is.na(Lw[["IR"]]) && Lw[["IR"]] >= max(Lw, na.rm = TRUE) - 1e-12,
      dLogL  = Lw[["IR"]] - max(Lw[["LSE"]], Lw[["NR"]], na.rm = TRUE),
      L_LSE = Lw[["LSE"]], L_NR = Lw[["NR"]], L_R = Lw[["R"]], L_IR = Lw[["IR"]])
  }
}

res <- do.call(rbind, rows)
if (is.null(res)) stop("nothing to rank")

# Median-first: for each window, the seed nearest the median update is the defensible pick,
# and the maximum is reported only so the spread is visible. A caption that says "the median
# of N realisations" is a much stronger claim than one that says "chosen for clarity".
res$dist_to_median <- NA_real_
for (w in unique(res$window)) {
  k <- res$window == w
  res$dist_to_median[k] <- abs(res$balance[k] - stats::median(res$balance[k]))
}

out <- file.path(outroot, "figure_1_seed_ranking.csv")
utils::write.csv(res[order(res$window, -res$balance), ], out, row.names = FALSE)

fmt <- function(k) {
  k <- k[, c("window", "seed", "balance", "update", "contrast")]
  k$balance <- round(k$balance, 4); k$update <- round(k$update, 4)
  k$contrast <- round(k$contrast, 4)
  k
}
cat("\n== per window: best balance (no dead cycle) ==\n")
for (w in unique(res$window)) {
  k <- res[res$window == w, ]
  print(fmt(k[order(-k$balance), ][1:min(4, nrow(k)), ]), row.names = FALSE)
}
if (0 %in% res$seed) {
  cat("\n== the incumbent trace (seed 0 = the one currently in the paper) ==\n")
  print(fmt(res[res$seed == 0, ][order(res$window[res$seed == 0]), ]), row.names = FALSE)
}
cat("\n== sign-opposed pairs: the update reverses between the two windows ==\n")
flip <- res[res$flip_IR & res$flip_R, ]
if (!nrow(flip)) {
  cat("  none: in every (seed, window pair) both members update the same way twice\n")
} else {
  f <- flip[order(-flip$balance), ]
  f$IR <- sprintf("%+.3f/%+.3f", f$IR_1, f$IR_2)
  f$R  <- sprintf("%+.3f/%+.3f", f$R_1,  f$R_2)
  f$balance <- round(f$balance, 4)
  print(f[1:min(10, nrow(f)), c("window", "seed", "balance", "IR", "R")],
        row.names = FALSE)
}

cat("\n== all three at once: balanced updates, IR on top in row D ==\n")
ok <- res[res$IR_top, ]
if (!nrow(ok)) {
  cat("  none: no (seed, window pair) has IR best in row D\n")
} else {
  o <- ok[order(-ok$balance), ]
  o$IR <- sprintf("%+.3f/%+.3f", o$IR_1, o$IR_2)
  o$balance <- round(o$balance, 4); o$dLogL <- round(o$dLogL, 2)
  print(o[1:min(8, nrow(o)), c("window", "seed", "balance", "flip_IR", "dLogL", "IR")],
        row.names = FALSE)
}

cat("\n== per window: spread of balance ==\n")
print(do.call(rbind, lapply(split(res, res$window), function(k) data.frame(
  window = k$window[1], n = nrow(k),
  min = round(min(k$balance), 4), median = round(stats::median(k$balance), 4),
  max = round(max(k$balance), 4), best_seed = k$seed[which.max(k$balance)]))),
  row.names = FALSE)
cat("\nfull table: ", out, "\n", sep = "")
