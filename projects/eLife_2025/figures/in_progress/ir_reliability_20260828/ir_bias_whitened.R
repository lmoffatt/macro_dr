setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
# Whitened bias v = G^{1/2} b (symmetric sqrt of the Gaussian Fisher, G = C^-1), 6 components per cell.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(purrr);library(readr)})
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/data")
DIRS <- c("1c2ae6f","0ffbda7","87889e6","82b956f")
files <- unlist(lapply(DIRS, function(d) list.files(d, pattern="macro_IR.*battery_sim_G\\.csv$", full.names=TRUE)))
meta <- tibble(f=files) |> mutate(nch=as.numeric(sub(".*nch_([0-9]+)_.*","\\1",f)),
                                  z=as.numeric(sub(".*noise_([0-9.]+)_battery.*","\\1",f))) |>
  group_by(nch,z) |> slice(1) |> ungroup()
GDIB <- "Probit_statistics_Gaussian_Distortion_Induced_Bias"; GFC <- "Probit_statistics_Gaussian_Fisher_Covariance"
PN <- c("k_on","k_off","i","S","baseline","N_ch")
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
    Ghalf <- e$vectors %*% diag(1/sqrt(e$values)) %*% t(e$vectors)     # G^{1/2} = C^{-1/2}
    v <- as.vector(Ghalf %*% bb)
    tibble(nch=nch, z=z, interval=iv, param=PN, v=v, b=bb, se=sqrt(diag(M)))
  })
}
res <- pmap_dfr(meta, function(f,nch,z) one(f,nch,z))
write_csv(res, "/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828/digest/ir_bias_whitened.csv")
chk <- res |> group_by(nch,z,interval) |> summarise(norm=sqrt(sum(v^2)), .groups="drop")
mh <- read_csv("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828/digest/ir_bias_mahal.csv", show_col_types=FALSE)
cat("cells:", nrow(chk), "  max |norm - mahal| =", signif(max(abs(chk$norm - inner_join(chk, mh, by=c("nch","z","interval"))$mahal)),3), "\n")
