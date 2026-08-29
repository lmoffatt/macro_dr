setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
# Bias magnitude + anisotropy: GDIB vector against the full Gaussian Fisher covariance eigenbasis.
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
GID  <- "Probit_statistics_Likelihood_Gaussian_Information_Distortion"
one <- function(f, nch, z){
  d <- suppressWarnings(read_csv(f, skip=1, show_col_types=FALSE, progress=FALSE)) |>
       filter(statistic=="value", probit=="mean")
  map_dfr(sort(unique(d$interval_in_tau)), function(iv){
    di <- d |> filter(interval_in_tau==iv)
    bv <- di |> filter(component_path==GDIB) |> arrange(param_index)
    if(nrow(bv)!=6) return(NULL)
    b <- bv$value
    Cm <- matrix(NA_real_,6,6); s <- di |> filter(component_path==GFC)
    for(k in seq_len(nrow(s))) Cm[s$value_row[k]+1, s$value_col[k]+1] <- s$value[k]
    Gm <- matrix(NA_real_,6,6); g <- di |> filter(component_path==GID)
    for(k in seq_len(nrow(g))) Gm[g$value_row[k]+1, g$value_col[k]+1] <- g$value[k]
    if(any(is.na(Cm))) return(NULL)
    Cm <- (Cm+t(Cm))/2; e <- eigen(Cm, symmetric=TRUE)   # covariance: e$values desc = softest first
    if(any(e$values<=0)) return(NULL)
    ck <- as.vector(t(e$vectors) %*% b) / sqrt(e$values) # bias in sigmas along each principal axis
    mag <- sqrt(sum(ck^2))
    tr_gid <- if(!any(is.na(Gm))) sum(diag((Gm+t(Gm))/2)) else 6
    tibble(nch=nch, z=z, interval=iv, mag=mag,
           share_soft = ck[1]^2/sum(ck^2),               # softest (largest variance) direction
           share_top2 = sum(ck[1:2]^2)/sum(ck^2),
           share_stiff= ck[6]^2/sum(ck^2),
           c_soft = ck[1],
           floor = sqrt(tr_gid/10000))                   # ensemble-noise Mahalanobis floor
  })
}
res <- pmap_dfr(meta, function(f,nch,z) one(f,nch,z))
write_csv(res, "/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828/digest/ir_bias_direction.csv")
cat("cells:", nrow(res), "\n")
