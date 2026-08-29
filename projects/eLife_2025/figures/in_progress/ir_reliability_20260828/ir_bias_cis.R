setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(purrr);library(readr)})
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/data")
DIRS <- c("1c2ae6f","0ffbda7","87889e6","82b956f")
files <- unlist(lapply(DIRS, function(d) list.files(d, pattern="macro_IR.*battery_sim_G\\.csv$", full.names=TRUE)))
meta <- tibble(f=files) |>
  mutate(nch = as.numeric(sub(".*nch_([0-9]+)_.*","\\1",f)),
         z   = as.numeric(sub(".*noise_([0-9.]+)_battery.*","\\1",f))) |>
  group_by(nch,z) |> slice(1) |> ungroup()
GDIB <- "Probit_statistics_Gaussian_Distortion_Induced_Bias"
GFC  <- "Probit_statistics_Gaussian_Fisher_Covariance"
one <- function(f, nch, z){
  d <- suppressWarnings(read_csv(f, skip=1, show_col_types=FALSE, progress=FALSE)) |>
       filter(statistic=="value")
  b <- d |> filter(component_path==GDIB,
                   probit=="mean" | (probit=="quantile" & quantile_level %in% c(0.025,0.975))) |>
       mutate(key = ifelse(probit=="mean","m", ifelse(quantile_level<0.5,"lo","hi"))) |>
       select(interval_in_tau, param_index, key, value) |>
       pivot_wider(names_from=key, values_from=value)
  v <- d |> filter(component_path==GFC, probit=="mean", value_row==value_col) |>
       transmute(interval_in_tau, param_index=value_row, se=sqrt(value))
  b |> inner_join(v, by=c("interval_in_tau","param_index")) |> mutate(nch=nch, z=z)
}
res <- pmap_dfr(meta, function(f,nch,z) one(f,nch,z))
write_csv(res, "/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828/digest/ir_bias_cis.csv")
cat("rows:", nrow(res), " cells:", nrow(distinct(res, nch, z, interval_in_tau)),
    " finite CI:", sum(is.finite(res$lo) & is.finite(res$hi)), "\n")
