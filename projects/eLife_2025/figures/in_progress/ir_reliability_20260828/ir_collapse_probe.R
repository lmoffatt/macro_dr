setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(purrr); library(readr)})
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/data")
DIRS <- c("1c2ae6f","0ffbda7","87889e6","82b956f")
files <- unlist(lapply(DIRS, function(d) list.files(d, pattern="macro_IR.*battery_(sim|pool)_G\\.csv$", full.names=TRUE)))
meta <- tibble(f=files) |>
  mutate(nch = as.numeric(sub(".*nch_([0-9]+)_.*","\\1",f)),
         z   = as.numeric(sub(".*noise_([0-9.]+)_battery.*","\\1",f)),
         anchor = ifelse(grepl("battery_sim",f),"sim","pool"),
         dir = dirname(f)) |>
  group_by(nch,z,anchor) |> slice(1) |> ungroup()

GID <- "Probit_statistics_Likelihood_Gaussian_Information_Distortion"
GDIB<- "Probit_statistics_Gaussian_Distortion_Induced_Bias"
GFC <- "Probit_statistics_Gaussian_Fisher_Covariance"

one <- function(f, nch, z, anchor, dir){
  d <- suppressWarnings(read_csv(f, skip=1, show_col_types=FALSE, progress=FALSE))
  d <- d |> filter(probit=="mean", statistic=="value")
  gi <- d |> filter(component_path==GID)
  bi <- d |> filter(component_path==GDIB)
  fc <- d |> filter(component_path==GFC, value_row==value_col)
  ivs <- sort(unique(gi$interval_in_tau))
  map_dfr(ivs, function(iv){
    M <- matrix(NA_real_,6,6); s <- gi |> filter(interval_in_tau==iv)
    if(nrow(s)==0) return(NULL)
    for(k in seq_len(nrow(s))) M[s$value_row[k]+1, s$value_col[k]+1] <- s$value[k]
    if(any(is.na(M))) return(NULL)
    M <- (M+t(M))/2; e <- eigen(M, symmetric=TRUE)$values
    if(any(e<=0)) return(NULL)
    L <- log(e); ebar <- mean(L); vlog <- mean((L-ebar)^2)
    b <- bi |> filter(interval_in_tau==iv) |> arrange(param_index) |> pull(value)
    v <- fc |> filter(interval_in_tau==iv) |> arrange(value_row) |> pull(value)
    zmax <- if(length(b)==6 && length(v)==6) max(abs(b)/sqrt(v)) else NA_real_
    znrm <- if(length(b)==6 && length(v)==6) sqrt(sum((b/sqrt(v))^2)) else NA_real_
    tibble(nch=nch, z=z, anchor=anchor, dir=basename(dir), interval=iv,
           dAI=sqrt(sum(L^2)), ebar=ebar, aniso=sqrt(vlog), zmax=zmax, znrm=znrm)
  })
}
res <- pmap_dfr(meta |> select(f,nch,z,anchor,dir), one)
write_csv(res, "/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828/digest/ir_collapse_cells.csv")
cat("cells:", nrow(res), " anchors:", paste(unique(res$anchor),collapse=" "), "\n")
print(res |> count(anchor, nch) |> pivot_wider(names_from=anchor, values_from=n), n=30)
