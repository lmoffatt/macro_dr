setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(purrr);library(readr)})
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/data")
DIRS <- c("1c2ae6f","0ffbda7","87889e6","82b956f")
files <- unlist(lapply(DIRS, function(d) list.files(d, pattern="macro_IR.*battery_sim_G\\.csv$", full.names=TRUE)))
meta <- tibble(f=files) |> mutate(nch=as.numeric(sub(".*nch_([0-9]+)_.*","\\1",f)),
         z=as.numeric(sub(".*noise_([0-9.]+)_battery.*","\\1",f))) |> group_by(nch,z) |> slice(1) |> ungroup()
GDIB<-"Probit_statistics_Gaussian_Distortion_Induced_Bias"; GFC<-"Probit_statistics_Gaussian_Fisher_Covariance"
one <- function(f, nch, z){
  d <- suppressWarnings(read_csv(f, skip=1, show_col_types=FALSE, progress=FALSE)) |>
       filter(statistic=="value", probit=="mean")
  map_dfr(sort(unique(d$interval_in_tau)), function(iv){
    di <- d |> filter(interval_in_tau==iv)
    bv <- di |> filter(component_path==GDIB) |> arrange(param_index); if(nrow(bv)!=6) return(NULL)
    Cm <- matrix(NA_real_,6,6); s <- di |> filter(component_path==GFC)
    for(k in seq_len(nrow(s))) Cm[s$value_row[k]+1,s$value_col[k]+1] <- s$value[k]
    if(any(is.na(Cm))) return(NULL); Cm <- (Cm+t(Cm))/2
    e <- eigen(Cm, symmetric=TRUE); if(any(e$values<=0)) return(NULL)
    ck <- as.vector(t(e$vectors)%*%bv$value)/sqrt(e$values)
    tibble(nch=nch,z=z,interval=iv,rank=1:6,share=ck^2/sum(ck^2),c=ck,
           mag=sqrt(sum(ck^2)), sigma=sqrt(e$values),
           # dominant parameter of this eigenvector (for interpretation)
           dompar=apply(abs(e$vectors),2,which.max))
  })
}
res <- pmap_dfr(meta, function(f,nch,z) one(f,nch,z))
write_csv(res, "/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828/digest/ir_bias_rank.csv")
d2 <- read.csv("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828/digest/ir_bias_direction.csv") |>
      transmute(nch,z,interval,floor)
r <- res |> inner_join(d2, by=c("nch","z","interval")) |> mutate(real = mag > 2*floor)
cat("median share by eigen-rank (1=softest), REAL-bias cells vs null cells (isotropy 0.167):\n")
print(r |> group_by(real, rank) |> summarise(med=round(median(share),3), .groups="drop") |>
      pivot_wider(names_from=rank, values_from=med))
PN <- c("on","off","unitary","C_Noise","Baseline","Nch_mean")
cat("\nwhich rank dominates each real cell (rank of max share), counts:\n")
print(r |> filter(real) |> group_by(nch,z,interval) |> slice_max(share,n=1) |> ungroup() |> count(rank))
cat("\nfor real cells, dominant parameter of the WINNING eigenvector:\n")
print(r |> filter(real) |> group_by(nch,z,interval) |> slice_max(share,n=1) |> ungroup() |> count(dompar) |>
      mutate(param=PN[dompar]))
cat("\nsign consistency of winning projection c:\n")
print(r |> filter(real) |> group_by(nch,z,interval) |> slice_max(share,n=1) |> ungroup() |>
      mutate(sgn=sign(c)) |> count(dompar,sgn) |> mutate(param=PN[dompar]))
