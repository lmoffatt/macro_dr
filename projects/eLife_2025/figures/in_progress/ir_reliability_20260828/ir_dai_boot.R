# Bootstrap quantiles of the EMITTED affine distance of the GIDM (d_AI of each resample's matrix), per cell,
# to put error bars on the matrix-computed D_AI. Also compares emitted mean/median to the matrix value.
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(purrr);library(readr)})
DATA <- "/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/data"
DIRS <- c("1c2ae6f","0ffbda7","87889e6","82b956f")
files <- unlist(lapply(DIRS, function(d) list.files(file.path(DATA,d), pattern="macro_IR.*battery_sim_G\\.csv$", full.names=TRUE)))
meta <- tibble(f=files) |> mutate(nch=as.numeric(sub(".*nch_([0-9]+)_.*","\\1",f)), z=as.numeric(sub(".*noise_([0-9.]+)_battery.*","\\1",f))) |>
  group_by(nch,z) |> slice(1) |> ungroup()
TAG <- "Probit_statistics_Affine_Invariant_Distance_Likelihood_Gaussian_Information_Distortion"
one <- function(f,nch,z){
  d <- suppressWarnings(read_csv(f, skip=1, show_col_types=FALSE, progress=FALSE)) |> filter(component_path==TAG, statistic=="value")
  d |> mutate(key=ifelse(probit=="mean","mean", paste0("q", sub("^0\\.","",as.character(quantile_level))))) |>
    select(interval=interval_in_tau, key, value) |> pivot_wider(names_from=key, values_from=value) |> mutate(nch=nch, z=z)
}
res <- pmap_dfr(meta, function(f,nch,z) one(f,nch,z))
write_csv(res, "digest/ir_dai_boot.csv")
m <- read_csv("digest/ir_collapse_cells.csv", show_col_types=FALSE) |> filter(anchor=="sim") |> select(nch,z,interval,dAI)
cmp <- inner_join(res, m, by=c("nch","z","interval"))
cat("cells:", nrow(cmp), " columns:", paste(names(res), collapse=" "), "\n")
cat("emitted mean / matrix dAI, quantiles of the ratio:"); print(round(quantile(cmp$mean/cmp$dAI, c(.05,.5,.95)),2))
cat("bracket rate of matrix dAI by [q025,q975]:", round(mean(cmp$dAI >= cmp$q025 & cmp$dAI <= cmp$q975),2), "  by [q16,q84]:", round(mean(cmp$dAI >= cmp$q16 & cmp$dAI <= cmp$q84),2), "\n")
cat("half-width (q975-q025)/2 relative to matrix dAI, by dAI band:\n")
print(cmp |> mutate(hw=(q975-q025)/2, band=cut(dAI, c(0,0.1,0.2,0.4,1))) |> group_by(band) |> summarise(n=n(), med_hw=round(median(hw),3), med_rel=round(median(hw/dAI),2)) |> as.data.frame(), row.names=FALSE)
