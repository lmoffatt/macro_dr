# fig2_first_moment.R — the first-moment numbers of the Figure 2 subsection of 04_results.tex.
#
# Run from the repository root:  Rscript papers/1_method/decisions/recompute/fig2_first_moment.R
#
# Reads exactly the files figures/paper_both/figure_2.Rmd reads (same search path, same cell,
# same filters) and reports, per member and per parameter:
#   - the measured bias of the MLE cloud, mean(theta_hat) - theta_sim, with a 2000-resample
#     percentile bootstrap interval (the cloud's own filled marker against its cross);
#   - the bias of the amplitude PRODUCT log10(N_ch * i), which is what the mean current fixes;
#   - the predicted distortion-induced bias from battery_sim (the figure's open marker), against
#     the measured one.
# Written 2026-08-06 because the Figure 2 subsection asserted a kinetic-pair unbiasedness that the
# markers contradict. The numbers it prints are the ones quoted in 04_results.tex.

# UPDATED 2026-08-06 (second pass): the roster is EIGHT members. The un-averaged least-squares arm
# (nonlinearsqr_g, averaging=0, labelled LSE) landed from dirac under figures/data/a202e03; the arm
# that was previously labelled "LSE" is nonlinearsqr, averaging=1, and is relabelled ILSE
# (projects/eLife_2025/ops/slurm/dispatch_figure_3_LSE.sh:144-146). a202e03's engine change is av=0
# only, so no other column moves.

suppressPackageStartupMessages({library(dplyr); library(tidyr); library(purrr); library(stringr)})

FIG      <- "projects/eLife_2025/figures"
DIRS     <- file.path(FIG, "data", c("a202e03", "1f7138b", "1c2ae6f", "0ffbda7", "87889e6"))
FILE_PRE <- "figure_3_G_"; BAT_SUF <- "_G"; COMP_PRE <- "Gaussian_"
REP_NCH  <- 100; REP_NSIM <- 10000; REP_NOISE <- "0.1"; REP_GROUP <- 10; REP_INT <- 0.1

ROSTER <- c("nonlinearsqr_g", "nonlinearsqr", "macro_NR", "macro_INR", "macro_R", "macro_MR", "macro_VR", "macro_IR")
LAB    <- c(nonlinearsqr_g = "LSE", nonlinearsqr = "ILSE", macro_NR = "NR", macro_INR = "INR", macro_R = "R",
            macro_MR = "MR", macro_VR = "VR", macro_IR = "IR")
IDX    <- c(on = 0, off = 1, unitary_current = 2, Current_Noise = 3, Num_ch_mean = 5)
TRUTH  <- c(on = log10(10), off = log10(100), unitary_current = log10(1),
            Current_Noise = log10(as.numeric(REP_NOISE) / 1000), Num_ch_mean = log10(REP_NCH))

rd   <- function(p) read.csv(p, skip = 1, stringsAsFactors = FALSE)
pre  <- function(a) if (startsWith(a, "nonlinearsqr")) "figure_3_LSE_" else FILE_PRE
stem <- function(d, a) file.path(d, sprintf("%snch_%d_nsim_%d_%s_noise_%s",
                                            pre(a), REP_NCH, REP_NSIM, a, REP_NOISE))
need <- function(s) paste0(s, c("_mle_cloud_runs.csv",
                                paste0("_battery_pool", BAT_SUF, ".csv"),
                                paste0("_battery_sim",  BAT_SUF, ".csv")))
base_for <- function(a) { for (d in DIRS) { s <- stem(d, a); if (all(file.exists(need(s)))) return(s) }
                          NA_character_ }
b <- vapply(ROSTER, base_for, character(1))
if (any(is.na(b))) cat("MISSING, dropped:", paste(ROSTER[is.na(b)], collapse = ", "), "\n")
ALGOS <- ROSTER[!is.na(b)]; b <- b[!is.na(b)]

pts <- map_dfr(paste0(b, "_mle_cloud_runs.csv"), rd) %>%
  filter(component_path == "Model_Parameters_Hat", statistic == "value",
         Num_ch == REP_NCH, group_size == REP_GROUP, abs(interval_in_tau - REP_INT) < 1e-9,
         param_name %in% names(TRUTH)) %>%
  select(algorithm, sample_index, param_name, value)

set.seed(7)
bootci <- function(v, B = 2000) {
  m <- replicate(B, mean(sample(v, length(v), TRUE)))
  unname(quantile(m, c(.025, .975), type = 6))     # type=6 matches the program's Probit_statistics
}

res <- pts %>% group_by(algorithm, param_name) %>%
  summarise(n = n(), bias = mean(value) - TRUTH[first(param_name)],
            lo = bootci(value)[1] - TRUTH[first(param_name)],
            hi = bootci(value)[2] - TRUTH[first(param_name)], .groups = "drop") %>%
  mutate(algo = LAB[algorithm], resolved = ifelse(lo > 0 | hi < 0, "*", "")) %>%
  select(algo, param_name, n, bias, lo, hi, resolved)
res$algo <- factor(res$algo, levels = unname(LAB[ALGOS]))

for (p in names(TRUTH)) {
  cat("\n== ", p, "  (truth log10 = ", TRUTH[p], ") ==\n", sep = "")
  print(res %>% filter(param_name == p) %>% arrange(algo) %>%
          mutate(across(c(bias, lo, hi), ~round(.x, 4))) %>% as.data.frame(), row.names = FALSE)
}

cat("\n\n===== the amplitude split: the product is recovered, the factors are not =====\n")
w <- pts %>% filter(param_name %in% c("unitary_current", "Num_ch_mean")) %>%
  pivot_wider(names_from = param_name, values_from = value) %>%
  filter(is.finite(unitary_current), is.finite(Num_ch_mean)) %>%
  mutate(prod = unitary_current + Num_ch_mean)
t_prod <- TRUTH["unitary_current"] + TRUTH["Num_ch_mean"]
print(w %>% group_by(algorithm) %>%
        summarise(n = n(), bias_i = mean(unitary_current) - TRUTH["unitary_current"],
                  bias_N = mean(Num_ch_mean) - TRUTH["Num_ch_mean"],
                  bias_prod = mean(prod) - t_prod,
                  lo = bootci(prod)[1] - t_prod, hi = bootci(prod)[2] - t_prod, .groups = "drop") %>%
        mutate(algo = LAB[algorithm]) %>% select(algo, n, bias_i, bias_N, bias_prod, lo, hi) %>%
        mutate(across(where(is.numeric), ~round(.x, 4))) %>% as.data.frame(), row.names = FALSE)

cat("\n\n===== predicted bias (open marker) against measured bias (filled marker) =====\n")
dib <- map_dfr(paste0(b, "_battery_sim", BAT_SUF, ".csv"), rd) %>%
  mutate(component = str_replace(component_path, "^Probit_statistics_", "")) %>%
  filter(component == paste0(COMP_PRE, "Distortion_Induced_Bias"), statistic == "value",
         probit == "mean", Num_ch == REP_NCH, abs(interval_in_tau - REP_INT) < 1e-9) %>%
  mutate(param_name = names(IDX)[match(param_index, IDX)], algo = LAB[algorithm]) %>%
  filter(!is.na(param_name)) %>% select(algo, param_name, predicted = value)

cmp <- res %>% select(algo, param_name, bias, lo, hi) %>%
  left_join(dib, by = c("algo", "param_name")) %>%
  mutate(ratio   = predicted / bias,
         covered = ifelse(is.na(predicted), "", ifelse(predicted >= lo & predicted <= hi, "yes", "NO")))
print(cmp %>% arrange(param_name, algo) %>%
        mutate(across(where(is.numeric), ~round(.x, 4))) %>% as.data.frame(), row.names = FALSE)

# Figure 2--source data 1 (added 2026-08-25: the digit-migration pass moves the first-moment
# numbers out of the Results prose; their reviewer-visible home is this CSV).
sd_dir <- "projects/eLife_2025/figures/figure_2_source_data"
dir.create(sd_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(cmp %>% arrange(param_name, algo) %>%
            mutate(across(where(is.numeric), ~round(.x, 4))),
          file.path(sd_dir, "figure_2_source_data_bias.csv"), row.names = FALSE)
cat("\nwrote ", file.path(sd_dir, "figure_2_source_data_bias.csv"), "\n", sep = "")
