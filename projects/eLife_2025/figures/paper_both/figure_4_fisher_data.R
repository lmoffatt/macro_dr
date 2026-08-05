## figure_4_fisher_data.R — the numerical-against-analytic Fisher comparison, built on its own.
##
## WHY THIS IS A SEPARATE BUILDER and not three more lines in figure_4_data.R. Everything the
## Figure 4 set draws comes from the `_G` battery files, which are produced in GAUSSIAN-ONLY mode:
## there the numerical Fisher F_b is never computed (likelihood.cpp:3396-3402), so
## Likelihood_Gaussian_Fisher_Distortion is emitted as a slot and is `nan` in every one of the 240
## cells on disk. The quantity only exists in the NON-`_G` runs, which are a different file family
## with a different member list, a different noise reach and one extra column. Sharing a builder
## would mean one product whose staleness stamp mixes two file sets and whose rows are half empty;
## this file keeps the two apart.
##
## WHAT IT MEASURES. Likelihood_Gaussian_Fisher_Distortion = inv(sqrt(G_b)) * F_b * inv(sqrt(G_b)),
## with G_b the cheap analytic Gaussian-formula Fisher and F_b the finite-difference one, in the
## same frame the information distortion uses, so the identity is the null. It answers a question
## about the DEVICE and not about calibration: does the closed-form derivative reproduce what
## differencing the likelihood actually gives? That is a different question from GIDM, which asks
## whether either of them matches the score covariance.
##
## The two numbers are the same pair the spectrum product carries, and for the same reason:
##   mag    exp(mean(log lambda)), the typical factor by which the analytic formula misses
##   aniso  exp(sd(log lambda)),   whether it misses by the same factor in every direction
## Natural logs in, variance over p, and affine^2 = p*(mag_log^2 + aniso_log^2) holds here too.
##
## COVERAGE, measured 2026-08-02, and it is the reason this figure is provisional. Members with a
## non-`_G` run at nsim 10000: NR, R, MR, NMR, IR. Noise levels: 0.1, 1 and 10, three against the
## nine to twelve the `_G` sweep reaches. NOT present at all: nonlinearsqr (least squares) and
## macro_INR. Those two are exactly the members the question is about, least squares because that is
## where the analytic Fisher is reported to fail and INR because it is the high-noise comparison, so
## what this builder can draw today is the control group and not the case.
##
## READ NMR WITH CARE. The non-`_G` runs predate the fix at HEAD ("restore the interval conductance
## variance in the non-recursive path; NMR -> INR"), so the NMR column here is the version missing
## the N*ms term and is NOT the current INR. It is drawn because it is the largest departure in the
## file and worth seeing, and it must be labelled as superseded wherever it appears.
##
## THE EXTRA COLUMN. Non-`_G` files carry `axis_h_fim` (the finite-difference step, 1e-5) between
## `simulation_algorithm` and `value`, so `value` is column 30 where the `_G` files have it at 29.
## Selecting by NAME, as below, makes that a non-issue; positional readers get the step size instead
## of the number and it looks like a plausible constant, which is how it was nearly reported as a
## result.

FISHER_SRC <- "../figure_4_source_data/figure_4_source_data_gaussian_fisher.csv"

GFD_PATHS <- c(
  mag   = "Probit_statistics_Mean_Log_Eigenvalue_Likelihood_Gaussian_Fisher_Distortion",
  aniso = "Probit_statistics_Log_Eigenvalue_Variance_Likelihood_Gaussian_Fisher_Distortion")

# every member that could have a non-`_G` run, not the notebook's roster: one build serves any
# roster, and a member with no file simply contributes no rows
FISHER_MEMBERS <- c("nonlinearsqr", "macro_NR", "macro_INR", "macro_R", "macro_MR",
                    "macro_VR", "macro_NMR", "macro_IR")

.fisher_pat <- function(a, kind)
  sprintf("^%snch_([0-9]+)_nsim_10000_%s_noise_([0-9.eE+]+)_%s\\.csv$",
          if (startsWith(a, "nonlinearsqr")) "figure_3_LSE_" else "figure_3_", a, kind)

# the cells with BOTH anchors present, so a member never appears at one anchor only
fisher_cells <- function(algo) {
  pat   <- .fisher_pat(algo, "battery_sim")
  files <- basename(list.files(DATA_DIRS, pattern = pat))
  if (!length(files)) return(tibble::tibble(nch = integer(), z = character()))
  m  <- regmatches(files, regexec(pat, files))
  tibble::tibble(nch = as.integer(vapply(m, `[`, character(1), 2)),
                 z   = vapply(m, `[`, character(1), 3)) %>%
    dplyr::distinct() %>% dplyr::filter(nch %in% NCHS)
}

fisher_findf <- function(nch, algo, z, kind) {
  f <- file.path(DATA_DIRS,
                 sprintf("%snch_%d_nsim_10000_%s_noise_%s_%s.csv",
                         if (startsWith(algo, "nonlinearsqr")) "figure_3_LSE_" else "figure_3_",
                         nch, algo, z, kind))
  hit <- f[file.exists(f)]
  if (!length(hit)) NA_character_ else hit[1]
}

FISHER_KEEP <- c("component_path", "statistic", "probit", "quantile_level",
                 "interval_in_tau", "Num_ch", "value")

fisher_stamp <- function(paths) {
  i <- file.info(paths)
  j <- file.info("figure_4_fisher_data.R")
  sprintf("# figure_4 gaussian-vs-numerical Fisher | inputs=%d bytes=%.0f newest=%.0f | code=%.0f",
          length(paths), sum(i$size), max(as.numeric(i$mtime)), as.numeric(j$mtime))
}

# Same shaping as the spectrum product, and deliberately so: the two figures then read the same
# columns and a colour means the same thing on both. CI-aware toward 1 for the same reason, since
# a differenced Fisher estimated over a bootstrap spreads even where the two matrices agree.
fisher_shape <- function(d) d %>%
  dplyr::mutate(param = names(GFD_PATHS)[match(component_path, GFD_PATHS)],
                key = dplyr::case_when(probit == "mean" ~ "m", quantile_level == 0.025 ~ "lo",
                                       quantile_level == 0.975 ~ "hi", TRUE ~ NA_character_)) %>%
  dplyr::filter(!is.na(key), !is.na(param)) %>%
  dplyr::group_by(algo, anchor, param, noise, Num_ch, interval_in_tau, key) %>%
  dplyr::summarise(v = mean(value), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = key, values_from = v) %>%
  dplyr::mutate(dplyr::across(dplyr::any_of(c("m", "lo", "hi")),
                              ~ exp(ifelse(param == "aniso", sqrt(pmax(0, .x)), .x)))) %>%
  dplyr::mutate(Sconf = if ("lo" %in% names(.)) ifelse(m > 1, pmax(1, lo), pmin(1, hi)) else m,
                lx = log10(interval_in_tau), ly = log10(as.numeric(noise)))

fisher_source_data <- function() {
  cl <- do.call(rbind, lapply(FISHER_MEMBERS, function(a) {
    g <- fisher_cells(a)
    if (!nrow(g)) NULL else data.frame(algo = a, nch = g$nch, z = g$z, stringsAsFactors = FALSE)
  }))
  if (is.null(cl) || !nrow(cl)) stop("figure_4 fisher: no non-_G battery files under ",
                                     paste(DATA_DIRS, collapse = ", "))
  paths <- na.omit(unlist(lapply(c("battery_sim", "battery_pool"), function(k)
    vapply(seq_len(nrow(cl)), function(i) fisher_findf(cl$nch[i], cl$algo[i], cl$z[i], k), ""))))
  paths <- paths[nzchar(paths)]
  want  <- fisher_stamp(paths)

  if (file.exists(FISHER_SRC) && identical(readLines(FISHER_SRC, n = 1), want)) {
    cat("figure_4 fisher: source CSV is current, reading it\n")
    return(tibble::as_tibble(utils::read.csv(FISHER_SRC, skip = 1, stringsAsFactors = FALSE)))
  }

  cat("figure_4 fisher: inputs changed, rebuilding\n")
  out <- list()
  for (i in seq_len(nrow(cl))) {
    for (k in c("battery_sim", "battery_pool")) {
      f <- fisher_findf(cl$nch[i], cl$algo[i], cl$z[i], k)
      if (is.na(f)) next
      d <- .fast_read(f)
      d <- d[d$component_path %in% GFD_PATHS & d$statistic == "value" &
             d$probit %in% c("mean", "quantile"), intersect(FISHER_KEEP, names(d)), drop = FALSE]
      if (!nrow(d)) next
      out[[length(out) + 1]] <- transform(d, algo = cl$algo[i], noise = cl$z[i],
                                          anchor = if (k == "battery_sim") "sim" else "pool")
    }
  }
  res <- fisher_shape(dplyr::bind_rows(out))
  if (!dir.exists(dirname(FISHER_SRC))) dir.create(dirname(FISHER_SRC), recursive = TRUE)
  writeLines(want, FISHER_SRC)
  suppressWarnings(utils::write.table(res, FISHER_SRC, sep = ",", row.names = FALSE,
                                      col.names = TRUE, append = TRUE, qmethod = "double"))
  res
}

# ---- re-point the plane ------------------------------------------------------------------------
# figure_4_layout.R takes XLIM and YLIM from `bias` and `stat`, and the row heights from
# GRID_BY_ALGO through nch_noise_span. Those all describe the `_G` grid, which reaches noise 1e7;
# this figure has three noise levels and would be three slivers in an empty plane. So the four
# objects the layout reads for GEOMETRY are re-pointed at the Fisher grid here, between sourcing
# figure_4_common.R and sourcing figure_4_layout.R. Nothing else about common.R is disturbed: the
# scales, the criteria, the roster and the block builder are used exactly as they are.
fisher_repoint <- function(fis) {
  g <- dplyr::distinct(dplyr::transmute(fis, algo, nch = Num_ch, z = as.character(noise)))
  GRID_BY_ALGO <<- lapply(setNames(ALGOS, ALGOS),
                          function(a) dplyr::select(dplyr::filter(g, algo == a), nch, z))
  GRID   <<- dplyr::distinct(dplyr::select(g, nch, z))
  # Y_BRK is built in figure_4_common.R from GRID, which at that point is still the _G one reaching
  # 1e7. Left alone, the only break landing inside this plane's two decades is the 0.5 the IR sweep
  # contributes, so every row gets exactly ONE tick and it is a level this figure never ran. Rebuilt
  # from the Fisher grid by the same rule, powers of ten thinned to at most five.
  Y_DEC  <<- sort(unique(as.numeric(GRID$z)))
  Y_BRK  <<- log10(Y_DEC[seq(1, length(Y_DEC), by = max(1, ceiling(length(Y_DEC) / 5)))])
  Y_LAB  <<- formatC(10^Y_BRK, format = "g", digits = 3)
  # only Num_ch / lx / ly are read off these two; `comp` is there because the layout filters on it
  stat  <<- dplyr::mutate(dplyr::filter(fis, algo %in% ALGOS), comp = "total")
  bias  <<- dplyr::select(stat, Num_ch, ly)
  invisible(NULL)
}
