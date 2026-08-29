# The paper's pair (m = e^ebar, a = e^s) and the worst direction max|log lambda_i| from the AGGREGATE GIDM, plus the
# EMITTED bootstrap quantiles (over 100 resamples) of Mean_Log_Eigenvalue, Log_Eigenvalue_Variance, Max/Min_Eigenvalue
# for error bars. Sim anchor, all 560 IR cells.
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(purrr);library(readr)})
DATA <- "/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/data"
DIRS <- c("1c2ae6f","0ffbda7","87889e6","82b956f")
files <- unlist(lapply(DIRS, function(d) list.files(file.path(DATA,d), pattern="macro_IR.*battery_sim_G\\.csv$", full.names=TRUE)))
meta <- tibble(f=files) |> mutate(nch=as.numeric(sub(".*nch_([0-9]+)_.*","\\1",f)), z=as.numeric(sub(".*noise_([0-9.]+)_battery.*","\\1",f))) |>
  group_by(nch,z) |> slice(1) |> ungroup()
GID <- "Probit_statistics_Likelihood_Gaussian_Information_Distortion"
SC <- c(ebar="Probit_statistics_Mean_Log_Eigenvalue_Likelihood_Gaussian_Information_Distortion",
        s2="Probit_statistics_Log_Eigenvalue_Variance_Likelihood_Gaussian_Information_Distortion",
        lmax="Probit_statistics_Max_Eigenvalue_Likelihood_Gaussian_Information_Distortion",
        lmin="Probit_statistics_Min_Eigenvalue_Likelihood_Gaussian_Information_Distortion")
one <- function(f,nch,z){
  d <- suppressWarnings(read_csv(f, skip=1, show_col_types=FALSE, progress=FALSE)) |> filter(statistic=="value")
  g <- d |> filter(component_path==GID, probit=="mean")
  mat <- map_dfr(sort(unique(g$interval_in_tau)), function(iv){ s <- g |> filter(interval_in_tau==iv); M <- matrix(NA_real_,6,6)
    for(k in seq_len(nrow(s))) M[s$value_row[k]+1, s$value_col[k]+1] <- s$value[k]
    if(any(is.na(M))) return(NULL); M <- (M+t(M))/2; e <- eigen(M, symmetric=TRUE)$values; if(any(e<=0)) return(NULL)
    L <- log(e); tibble(interval=iv, ebar=mean(L), s=sqrt(mean((L-mean(L))^2)), wmax=max(abs(L)), lmax=max(e), lmin=min(e)) })
  sc <- d |> filter(component_path %in% SC, probit %in% c("mean","quantile")) |> filter(probit=="mean" | quantile_level %in% c(0.025,0.975)) |>
    mutate(what=names(SC)[match(component_path, SC)], key=ifelse(probit=="mean","m", ifelse(quantile_level<0.5,"lo","hi"))) |>
    transmute(interval=interval_in_tau, col=paste(what,key,sep="_"), value) |> pivot_wider(names_from=col, values_from=value)
  inner_join(mat, sc, by="interval") |> mutate(nch=nch, z=z)
}
res <- pmap_dfr(meta, function(f,nch,z) one(f,nch,z))
write_csv(res, "digest/ir_pair_boot.csv")
cat("cells:", nrow(res), " cols:", paste(names(res), collapse=" "), "\n")
res <- res |> mutate(St=0.1*z, x=nch*St*sqrt(interval))
cat("worst direction max|log lambda| vs s: ratio quantiles"); print(round(quantile(res$wmax/res$s, c(.05,.5,.95)),2))
cat("subordination |ebar| <= s:", sum(abs(res$ebar) <= res$s), "of", nrow(res), "\n")
fl <- function(v) max(v[res$x > 1e4])
for(c0 in c(0.4,0.45,0.5,0.55,0.6)) cat(sprintf("  wmax <= max(%.2f x^-1/5, floor %.3f): violations %d\n", c0, fl(res$wmax), sum(res$wmax > pmax(c0*res$x^-0.2, fl(res$wmax)))))
cat(sprintf("  s <= max(0.21 x^-1/5, floor %.3f): violations %d ; |ebar| <= same: %d\n", fl(res$s), sum(res$s > pmax(0.21*res$x^-0.2, fl(res$s))), sum(abs(res$ebar) > pmax(0.21*res$x^-0.2, fl(res$s)))))
cat("emitted vs matrix: ebar diff q50/q95 |", round(quantile(abs(res$ebar_m-res$ebar), c(.5,.95)),4), "| s (sqrt of emitted var) ratio", round(quantile(sqrt(res$s2_m)/res$s, c(.05,.5,.95)),2), "\n")
cat("half-widths (95% band): ebar", round(median((res$ebar_hi-res$ebar_lo)/2),3), " s2->s approx", round(median((sqrt(res$s2_hi)-sqrt(res$s2_lo))/2),3), " log lmax", round(median((log(res$lmax_hi)-log(res$lmax_lo))/2),3), "\n")
