setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
# Test of the insight: whitened bias  w = G^{1/2} b  projected on the GIDM eigenbasis W,
# beta_k = W_k' w, against log lambda_k of the GIDM. Orientation-free tests use beta_k^2.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(purrr);library(readr)})
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/data")
DIRS <- c("1c2ae6f","0ffbda7","87889e6","82b956f")
files <- unlist(lapply(DIRS, function(d) list.files(d, pattern="macro_IR.*battery_sim_G\\.csv$", full.names=TRUE)))
meta <- tibble(f=files) |> mutate(nch=as.numeric(sub(".*nch_([0-9]+)_.*","\\1",f)),
  z=as.numeric(sub(".*noise_([0-9.]+)_battery.*","\\1",f))) |> group_by(nch,z) |> slice(1) |> ungroup()
GDIB<-"Probit_statistics_Gaussian_Distortion_Induced_Bias"; GFC<-"Probit_statistics_Gaussian_Fisher_Covariance"
GID<-"Probit_statistics_Likelihood_Gaussian_Information_Distortion"
mat <- function(s){ M<-matrix(NA_real_,6,6); for(k in seq_len(nrow(s))) M[s$value_row[k]+1,s$value_col[k]+1]<-s$value[k]; if(any(is.na(M))) return(NULL); (M+t(M))/2 }
one <- function(f,nch,z){
  d <- suppressWarnings(read_csv(f, skip=1, show_col_types=FALSE, progress=FALSE)) |> filter(statistic=="value", probit=="mean")
  map_dfr(sort(unique(d$interval_in_tau)), function(iv){
    di <- d |> filter(interval_in_tau==iv)
    bv <- di |> filter(component_path==GDIB) |> arrange(param_index); if(nrow(bv)!=6) return(NULL)
    C <- mat(di |> filter(component_path==GFC)); M <- mat(di |> filter(component_path==GID))
    if(is.null(C)||is.null(M)) return(NULL)
    ec <- eigen(C, symmetric=TRUE); if(any(ec$values<=0)) return(NULL)
    Ghalf <- ec$vectors %*% diag(1/sqrt(ec$values)) %*% t(ec$vectors)   # G^{1/2} = C^{-1/2}
    w <- as.vector(Ghalf %*% bv$value)
    em <- eigen(M, symmetric=TRUE); if(any(em$values<=0)) return(NULL)
    beta <- as.vector(t(em$vectors) %*% w)
    tibble(nch=nch,z=z,interval=iv,k=1:6, loglam=log(em$values), beta=beta,
           share=beta^2/sum(beta^2), mag=sqrt(sum(beta^2)))
  })
}
res <- pmap_dfr(meta, function(f,nch,z) one(f,nch,z))
write_csv(res, "/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828/digest/ir_bias_vs_logeig.csv")
fl <- read.csv("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828/digest/ir_bias_direction.csv") |> select(nch,z,interval,floor)
r <- res |> inner_join(fl, by=c("nch","z","interval")) |> mutate(real = mag > 2*floor)
cat("cells:", nrow(distinct(r,nch,z,interval)), " real:", nrow(distinct(r |> filter(real),nch,z,interval)), "\n")
# within-cell Spearman: share of bias^2 vs |log lambda|
wc <- r |> group_by(nch,z,interval,real) |> summarise(rho_abs = cor(share, abs(loglam), method="spearman"),
        rho_sgn = cor(share, loglam, method="spearman"), .groups="drop")
cat("\nwithin-cell Spearman(share_k, |log lambda_k|): median real / null:",
    round(median(wc$rho_abs[wc$real]),2), "/", round(median(wc$rho_abs[!wc$real]),2), "\n")
cat("within-cell Spearman(share_k, log lambda_k signed): median real / null:",
    round(median(wc$rho_sgn[wc$real]),2), "/", round(median(wc$rho_sgn[!wc$real]),2), "\n")
cat("fraction of real cells with rho_abs>0:", round(mean(wc$rho_abs[wc$real]>0),2), "\n")
# where does the bias sit: in lambda>1 (variance understated) or lambda<1 directions?
cat("\nreal cells: share of bias^2 in directions with lambda>1 (isotropy would give the fraction of such directions):\n")
print(r |> filter(real) |> group_by(nch,z,interval) |> summarise(sh_gt1 = sum(share[loglam>0]), frac_gt1 = mean(loglam>0), .groups="drop") |>
      summarise(med_share_gt1=round(median(sh_gt1),2), med_frac_dirs_gt1=round(median(frac_gt1),2)))
# pooled: share by rank of |log lambda| (1 = most distorted direction)
cat("\nmedian share by rank of |log lambda| (1 = most distorted), real cells:\n")
print(r |> filter(real) |> group_by(nch,z,interval) |> mutate(rk = rank(-abs(loglam))) |> ungroup() |>
      group_by(rk) |> summarise(med=round(median(share),3)) |> pivot_wider(names_from=rk, values_from=med))
cat("same, null cells:\n")
print(r |> filter(!real) |> group_by(nch,z,interval) |> mutate(rk = rank(-abs(loglam))) |> ungroup() |>
      group_by(rk) |> summarise(med=round(median(share),3)) |> pivot_wider(names_from=rk, values_from=med))
# magnitude relation across cells: mag_bias vs d_AI
cc <- r |> group_by(nch,z,interval,real,mag) |> summarise(dAI=sqrt(sum(loglam^2)), .groups="drop")
m <- lm(log10(mag)~log10(dAI), data=cc |> filter(real))
cat(sprintf("\nreal cells: log|bias|_M vs log d_AI slope %.2f R2 %.2f ; median ratio |bias|_M/d_AI = %.2f\n",
    coef(m)[2], summary(m)$r.squared, median(cc$mag[cc$real]/cc$dAI[cc$real])))
