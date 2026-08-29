setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
# Full Gaussian-Fisher covariance per cell -> Mahalanobis bias norm sqrt(b' C^-1 b) (coordinate-free)
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(purrr);library(readr)})
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/data")
DIRS <- c("1c2ae6f","0ffbda7","87889e6","82b956f")
files <- unlist(lapply(DIRS, function(d) list.files(d, pattern="macro_IR.*battery_sim_G\\.csv$", full.names=TRUE)))
meta <- tibble(f=files) |> mutate(nch=as.numeric(sub(".*nch_([0-9]+)_.*","\\1",f)),
                                  z=as.numeric(sub(".*noise_([0-9.]+)_battery.*","\\1",f))) |>
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
    M <- (M+t(M))/2; e <- eigen(M, symmetric=TRUE)
    if(any(e$values<=0)) return(tibble(nch=nch,z=z,interval=iv,mahal=NA_real_,note="C not PD"))
    Ginv_b <- e$vectors %*% ((t(e$vectors) %*% bb)/e$values)   # C^-1 b
    tibble(nch=nch, z=z, interval=iv, mahal=sqrt(sum(bb*Ginv_b)), note="")
  })
}
res <- pmap_dfr(meta, function(f,nch,z) one(f,nch,z))
write_csv(res, "/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828/digest/ir_bias_mahal.csv")
cat("cells:", nrow(res), " not PD:", sum(res$note!=""), "\n")
print(round(quantile(res$mahal, c(0,.5,.9,.99,1), na.rm=TRUE),4))
