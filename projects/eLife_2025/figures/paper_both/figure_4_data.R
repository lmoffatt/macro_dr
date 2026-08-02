## figure_4_data.R — builds the source data for the whole Figure 4 set, once, and writes it as CSV.
##
## THE PROBLEM IT SOLVES. Each notebook on this plane used to parse about a gigabyte of battery CSV
## to draw a few thousand numbers, and the six of them parsed the same gigabyte. Measured
## 2026-08-01 on the five-member roster: 307 s per notebook, 3 x 76 s of reading, to produce 4900 +
## 14700 + 4895 rows, which is 2.7 MB. Most of that was not I/O at all: `pmap_dfr` bound 140 files
## of 31675 rows into a 4.4-million-row frame and only then filtered it.
##
## THREE THINGS FIX IT, and all three are here:
##   ONE PASS. The three products come out of the same sweep, filtered PER FILE before anything is
##   bound, so nothing large is ever held.
##   THE UNION ROSTER. The sweep covers every member on disk, not the calling notebook's roster, so
##   the six notebooks share one build instead of one each. Each notebook filters afterwards, which
##   is free.
##   fread. data.table reads a 7.1 MB battery file in 0.03 s against read.csv's 0.25 s. Falls back
##   to read.csv if data.table is not installed.
##
## THE OUTPUT IS CSV, NOT AN OPAQUE CACHE, and that is deliberate: these files ARE the figure's
## source data. eLife asks for one per figure and this is it, already in the shape the figure
## draws, with a provenance line on top. They can ship with the paper unchanged.
##
## FORMAT. Line 1 is the stamp, line 2 the header, the rest data — the same layout the battery
## files themselves use, so `read.csv(f, skip = 1)` reads them and `readLines(f, n = 1)` reads the
## provenance. The stamp is what decides staleness: number of input files, their total size, the
## newest mtime among them, and the mtimes of the two code files that build the products. A new run
## landing changes the count or the newest mtime; editing the filters changes the code mtime.
##
## FLOW. Rendering any figure sources figure_4_common.R, which calls fig4_source_data() here. If
## nothing has moved it reads three CSVs and returns in about a second. If anything has, it
## rebuilds all three and writes them, and every other figure then finds them fresh. There is no
## separate step to remember and no ordering between notebooks.

SRC_DIR <- "../figure_4_source_data"
SRC <- c(bias  = "figure_4_source_data_bias.csv",
         stat  = "figure_4_source_data_distortion.csv",
         sedat = "figure_4_source_data_standard_error.csv",
         scal  = "figure_4_source_data_cell_scalars.csv")

# Every member with runs on disk. NOT the calling notebook's roster: the point is that one build
# serves all six. A member added here costs one sweep the first time and nothing after.
FIG4_ALL_MEMBERS <- c("nonlinearsqr", "macro_NR", "macro_INR", "macro_R",
                      "macro_MR", "macro_VR", "macro_IR")

.fast_read <- local({
  if (requireNamespace("data.table", quietly = TRUE)) {
    function(f) as.data.frame(data.table::fread(f, skip = 1, showProgress = FALSE))
  } else {
    function(f) utils::read.csv(f, skip = 1, stringsAsFactors = FALSE)
  }
})

fig4_stamp <- function(paths) {
  code <- c("figure_4_common.R", "figure_4_layout.R", "figure_4_data.R")
  code <- code[file.exists(code)]
  i <- file.info(paths); j <- file.info(code)
  sprintf("# figure_4 source data | inputs=%d bytes=%.0f newest=%.0f | code=%s",
          length(paths), sum(i$size), max(as.numeric(i$mtime)),
          paste(sprintf("%s:%.0f", basename(code), as.numeric(j$mtime)), collapse = ","))
}

# the union cell set, detected the same way figure_4_common.R detects the roster's own
fig4_all_cells <- function() {
  do.call(rbind, lapply(FIG4_ALL_MEMBERS, function(a) {
    g <- autodetect_cells(a)
    if (!nrow(g)) NULL else data.frame(algo = a, nch = g$nch, z = g$z,
                                       stringsAsFactors = FALSE)
  }))
}

KEEP_COLS <- c("component_path", "variable", "param_index", "param_col", "statistic", "probit",
               "quantile_level", "interval_in_tau", "Num_ch", "value")

# The per-cell scalars: no parameter index, one number per (member, cell, interval). They were each
# costing a FOURTH and FIFTH full sweep -- the condition number for the greying in
# figure_4_common.R, the residual lag for the lag/kappa supplement -- so they ride this pass too.
KAPPA_PATH <- "Probit_statistics_Spectrum_Condition_Number_Gaussian_Fisher_Covariance"

# BOTH ANCHORS for the half-B products (2026-08-02). The distortion, the corrected standard error
# and the per-cell scalars are extracted from battery_sim_G AND battery_pool_G and tagged with
# `anchor`, because GIDM at theta_sim and GIDM at theta_pool are DIFFERENT quantities: the
# numerator J and the denominator G_bar both move with the anchor, and whether that displacement
# matters is an open question (figure_4_sim.Rmd draws the sim-anchored plane; the ratio
# G_pool^-1/2 G_sim G_pool^-1/2 is the next step). Costs nothing in I/O: the sim file was already
# being read for the bias, so this is three more subsets of a frame already in memory.
#
# The BIAS is deliberately NOT extracted at both anchors. GDIB = G_bar^-1 * mean score, and at
# theta_pool the pooled score vanishes by the stationarity condition of the joint fit, so the
# pool-anchored bias is zero by construction for every member and carries no information.
.extract_half_b <- function(d, a, z, anchor) list(
  stat = transform(
    d[d$component_path %in% COMPS & d$param_index == d$param_col &
      d$param_index %in% PIDX & d$statistic == "value" &
      d$probit %in% c("mean", "quantile"), , drop = FALSE],
    algo = a, noise = z, anchor = anchor),
  se = transform(
    d[d$component_path == "Probit_statistics_Gaussian_Distortion_Corrected_Covariance" &
      d$param_index == d$param_col & d$param_index %in% PIDX &
      d$statistic == "value" & d$probit == "mean", , drop = FALSE],
    algo = a, noise = z, anchor = anchor),
  scal = transform(
    d[(d$component_path == KAPPA_PATH & d$statistic == "value" & d$probit == "mean") |
      (d$variable == "Report_integral_r_std" &
         d$statistic == "integral_correlation_lag" & d$probit == "mean") |
      (d$variable == "r2_std" & d$statistic == "mean" & d$probit == "mean"), , drop = FALSE],
    algo = a, noise = z, anchor = anchor))

fig4_build <- function(cl) {
  t0 <- Sys.time()
  bias_l <- stat_l <- se_l <- scal_l <- vector("list", nrow(cl))
  for (i in seq_len(nrow(cl))) {
    a <- cl$algo[i]; n <- cl$nch[i]; z <- cl$z[i]

    ds <- .fast_read(findf(n, a, z, "battery_sim_G"))[, KEEP_COLS]
    bias_l[[i]] <- transform(
      ds[ds$component_path == "Probit_statistics_Gaussian_Distortion_Induced_Bias" &
         ds$param_index %in% PIDX & ds$statistic == "value" &
         ds$probit %in% c("mean", "quantile"), , drop = FALSE],
      algo = a, noise = z)
    bs <- .extract_half_b(ds, a, z, "sim")

    dp <- .fast_read(findf(n, a, z, "battery_pool_G"))[, KEEP_COLS]
    bp <- .extract_half_b(dp, a, z, "pool")

    stat_l[[i]] <- rbind(bs$stat, bp$stat)
    se_l[[i]]   <- rbind(bs$se,   bp$se)
    scal_l[[i]] <- rbind(bs$scal, bp$scal)

    if (i %% 25 == 0) cat(sprintf("  figure_4 data: %d/%d cells\n", i, nrow(cl)))
  }
  cat(sprintf("  figure_4 data: swept %d cells in %.0f s\n", nrow(cl),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  list(bias = fig4_shape_bias(dplyr::bind_rows(bias_l)),
       stat = fig4_shape_stat(dplyr::bind_rows(stat_l)),
       sedat = fig4_shape_sedat(dplyr::bind_rows(se_l)),
       scal  = fig4_shape_scal(dplyr::bind_rows(scal_l)))
}

# ---- the shaping, kept here so the CSVs are already what the blocks draw --------------------
# BIAS is anchored at theta_sim: it cannot come from the pool files, where a bias evaluated at the
# optimum is zero by construction because the score vanishes there.
fig4_shape_bias <- function(d) d %>%
  dplyr::mutate(key = dplyr::case_when(probit == "mean" ~ "m", quantile_level == 0.025 ~ "lo",
                                       quantile_level == 0.975 ~ "hi", TRUE ~ NA_character_)) %>%
  dplyr::filter(!is.na(key)) %>%
  dplyr::group_by(algo, noise, Num_ch, interval_in_tau, param_index, key) %>%
  dplyr::summarise(v = mean(value), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = key, values_from = v) %>%
  # CI-aware toward ZERO: a cell whose interval covers 0 renders white
  dplyr::mutate(Bconf = ifelse(m > 0, pmax(0, lo), pmin(0, hi)),
                lx = log10(interval_in_tau), ly = log10(as.numeric(noise)),
                param = names(PIDX)[match(param_index, PIDX)])

fig4_shape_stat <- function(d) d %>%
  dplyr::mutate(comp = names(COMPS)[match(component_path, COMPS)],
                key  = dplyr::case_when(probit == "mean" ~ "m", quantile_level == 0.025 ~ "lo",
                                        quantile_level == 0.975 ~ "hi", TRUE ~ NA_character_)) %>%
  dplyr::filter(!is.na(key)) %>%
  dplyr::group_by(algo, anchor, comp, noise, Num_ch, interval_in_tau, param_index, key) %>%
  dplyr::summarise(v = mean(value), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = key, values_from = v) %>%
  # CI-aware: a cell whose 95% interval brackets 1 collapses to exactly 1 and renders white
  dplyr::mutate(Dconf = ifelse(m > 1, pmax(1, lo), pmin(1, hi)),
                lx = log10(interval_in_tau), ly = log10(as.numeric(noise)),
                param = names(PIDX)[match(param_index, PIDX)])

fig4_shape_sedat <- function(d) d %>%
  dplyr::group_by(algo, anchor, noise, Num_ch, interval_in_tau, param_index) %>%
  dplyr::summarise(v = mean(value), .groups = "drop") %>%
  dplyr::filter(is.finite(v), v > 0) %>%
  dplyr::mutate(sev = sqrt(v), lx = log10(interval_in_tau), ly = log10(as.numeric(noise)),
                param = names(PIDX)[match(param_index, PIDX)])

# r2_std is emitted as the SUM over the recording's intervals, so its null is the NUMBER of them,
# 10 / interval_in_tau, and not 1; the emitted `count` is 10000, the number of simulations, and is
# NOT the divisor. Normalized here so the CSV carries a quantity whose null is 1 like the others.
fig4_shape_scal <- function(d) d %>%
  dplyr::mutate(param = dplyr::case_when(component_path == KAPPA_PATH ~ "kappa",
                                         variable == "r2_std" ~ "r2", TRUE ~ "lag"),
                value = ifelse(param == "r2", value / (10 / interval_in_tau), value)) %>%
  # every evolution row is emitted TWICE (blank-segment and segment == 0 copies): collapse first
  dplyr::group_by(algo, anchor, param, noise, Num_ch, interval_in_tau) %>%
  dplyr::summarise(v = mean(value), .groups = "drop") %>%
  dplyr::mutate(lx = log10(interval_in_tau), ly = log10(as.numeric(noise)))

# ---- the entry point ---------------------------------------------------------------------------
fig4_source_data <- function() {
  cl <- fig4_all_cells()
  paths <- c(vapply(seq_len(nrow(cl)),
                    function(i) findf(cl$nch[i], cl$algo[i], cl$z[i], "battery_sim_G"), ""),
             vapply(seq_len(nrow(cl)),
                    function(i) findf(cl$nch[i], cl$algo[i], cl$z[i], "battery_pool_G"), ""))
  want <- fig4_stamp(paths)
  fs   <- file.path(SRC_DIR, SRC)

  if (all(file.exists(fs)) &&
      all(vapply(fs, function(f) identical(readLines(f, n = 1), want), TRUE))) {
    cat("figure_4 data: source CSVs are current, reading them\n")
    return(lapply(setNames(fs, names(SRC)),
                  function(f) tibble::as_tibble(utils::read.csv(f, skip = 1,
                                                                stringsAsFactors = FALSE))))
  }

  cat("figure_4 data: inputs changed, rebuilding the source CSVs\n")
  out <- fig4_build(cl)
  if (!dir.exists(SRC_DIR)) dir.create(SRC_DIR, recursive = TRUE)
  for (k in names(SRC)) {
    f <- file.path(SRC_DIR, SRC[[k]])
    # write.csv has no `append`; write.table does, and the stamp has to be line 1
    writeLines(want, f)
    suppressWarnings(utils::write.table(out[[k]], f, sep = ",", row.names = FALSE,
                                        col.names = TRUE, append = TRUE, qmethod = "double"))
  }
  out
}
