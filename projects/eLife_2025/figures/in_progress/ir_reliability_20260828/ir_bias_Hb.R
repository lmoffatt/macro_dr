# H*b per cell (H = inverse of the Gaussian Fisher covariance), for a delta-method SE of b'Hb:
#   var(b'Hb) ~ 4 * sum_k (Hb)_k^2 * var(b_k), with var(b_k) from the per-parameter bootstrap CI of b (ir_bias_cis).
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(purrr);library(readr)})
DATA <- "/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/data"
DIRS <- c("1c2ae6f","0ffbda7","87889e6","82b956f")
files <- unlist(lapply(DIRS, function(d) list.files(file.path(DATA,d), pattern="macro_IR.*battery_sim_G\\.csv$", full.names=TRUE)))
meta <- tibble(f=files) |> mutate(nch=as.numeric(sub(".*nch_([0-9]+)_.*","\\1",f)), z=as.numeric(sub(".*noise_([0-9.]+)_battery.*","\\1",f))) |>
  group_by(nch,z) |> slice(1) |> ungroup()
GDIB <- "Probit_statistics_Gaussian_Distortion_Induced_Bias"; GFC <- "Probit_statistics_Gaussian_Fisher_Covariance"
one <- function(f,nch,z){
  d <- suppressWarnings(read_csv(f, skip=1, show_col_types=FALSE, progress=FALSE)) |> filter(statistic=="value", probit=="mean")
  b <- d |> filter(component_path==GDIB) |> select(interval_in_tau, param_index, value)
  C <- d |> filter(component_path==GFC) |> select(interval_in_tau, value_row, value_col, value)
  map_dfr(sort(unique(b$interval_in_tau)), function(iv){
    bb <- b |> filter(interval_in_tau==iv) |> arrange(param_index) |> pull(value)
    cc <- C |> filter(interval_in_tau==iv); M <- matrix(NA_real_,6,6)
    for(k in seq_len(nrow(cc))) M[cc$value_row[k]+1, cc$value_col[k]+1] <- cc$value[k]
    if(length(bb)!=6 || any(is.na(M))) return(NULL)
    M <- (M+t(M))/2; e <- eigen(M, symmetric=TRUE); if(any(e$values<=0)) return(NULL)
    Hb <- as.vector(e$vectors %*% ((t(e$vectors) %*% bb)/e$values))
    tibble(nch=nch, z=z, interval=iv, param_index=0:5, Hb=Hb)
  })
}
res <- pmap_dfr(meta, function(f,nch,z) one(f,nch,z))
ci <- read_csv("digest/ir_bias_cis.csv", show_col_types=FALSE) |> mutate(se_b=(hi-lo)/3.92) |> select(nch,z,interval=interval_in_tau,param_index,m,se_b)
out <- inner_join(res, ci, by=c("nch","z","interval","param_index")) |> group_by(nch,z,interval) |>
  summarise(q=sum(m*Hb), se_q=2*sqrt(sum(Hb^2*se_b^2)), .groups="drop")
write_csv(out, "digest/ir_bias_Hb_se.csv")
mh <- read_csv("digest/ir_bias_mahal.csv", show_col_types=FALSE)
chk <- inner_join(out, mh, by=c("nch","z","interval"))
cat("cells:", nrow(out), " max |q - mahal^2| =", signif(max(abs(chk$q - chk$mahal^2)),3), "\n")
cat("relative SE of b'Hb (se_q/q), by q band:\n")
print(out |> mutate(band=cut(q, c(0,0.003,0.01,0.03,1))) |> group_by(band) |> summarise(n=n(), med_rel=round(median(se_q/q),2), q90_rel=round(quantile(se_q/q,.9),2)) |> as.data.frame(), row.names=FALSE)
