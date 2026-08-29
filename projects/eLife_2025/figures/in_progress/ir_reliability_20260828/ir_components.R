# ebar / s for the three distortion matrices (total, sample, correlation), from the aggregate (probit-mean) matrices.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(purrr);library(readr)})
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/data")
DIRS <- c("1c2ae6f","0ffbda7","87889e6","82b956f")
files <- unlist(lapply(DIRS, function(d) list.files(d, pattern="macro_IR.*battery_sim_G\\.csv$", full.names=TRUE)))
meta <- tibble(f=files) |> mutate(nch=as.numeric(sub(".*nch_([0-9]+)_.*","\\1",f)), z=as.numeric(sub(".*noise_([0-9.]+)_battery.*","\\1",f))) |>
  group_by(nch,z) |> slice(1) |> ungroup()
TAGS <- c(total="Probit_statistics_Likelihood_Gaussian_Information_Distortion",
          sample="Probit_statistics_Gaussian_Sample_Distortion",
          corr="Probit_statistics_Likelihood_Correlation_Distortion")
eigstats <- function(sub){ M <- matrix(NA_real_,6,6); for(k in seq_len(nrow(sub))) M[sub$value_row[k]+1, sub$value_col[k]+1] <- sub$value[k]
  if(any(is.na(M))) return(c(NA,NA,NA)); M <- (M+t(M))/2; e <- eigen(M, symmetric=TRUE)$values; if(any(e<=0)) return(c(NA,NA,NA))
  L <- log(e); c(mean(L), sqrt(mean((L-mean(L))^2)), sqrt(sum(L^2))) }
one <- function(f,nch,z){
  d <- suppressWarnings(read_csv(f, skip=1, show_col_types=FALSE, progress=FALSE)) |> filter(probit=="mean", statistic=="value", component_path %in% TAGS)
  map_dfr(names(TAGS), function(cn) d |> filter(component_path==TAGS[[cn]]) |> group_by(interval_in_tau) |>
    group_modify(~{ s <- eigstats(.x); tibble(ebar=s[1], s=s[2], dAI=s[3]) }) |> ungroup() |> mutate(comp=cn)) |>
  rename(interval=interval_in_tau) |> mutate(nch=nch, z=z)
}
res <- pmap_dfr(meta, function(f,nch,z) one(f,nch,z))
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828"); write_csv(res, "digest/ir_components_cells.csv")
cat("rows:", nrow(res), " cells per comp:"); print(table(res$comp)); cat("NA:", sum(is.na(res$s)), "\n")
wide <- res |> select(nch,z,interval,comp,ebar,s) |> pivot_wider(names_from=comp, values_from=c(ebar,s))
cat("identity check ebar_total = ebar_sample + ebar_corr : max |diff| =", signif(max(abs(wide$ebar_total - wide$ebar_sample - wide$ebar_corr), na.rm=TRUE),3), "\n")
cat("s quantiles by component:\n"); print(res |> group_by(comp) |> summarise(q05=round(quantile(s,.05),3), med=round(median(s),3), q95=round(quantile(s,.95),3), max=round(max(s),3)) |> as.data.frame(), row.names=FALSE)
cat("ebar quantiles by component:\n"); print(res |> group_by(comp) |> summarise(q05=round(quantile(ebar,.05),3), med=round(median(ebar),3), q95=round(quantile(ebar,.95),3)) |> as.data.frame(), row.names=FALSE)
# which variable does each component want? free exponents on s above its floor
cc <- res |> mutate(St=0.1*z)
for(cn in names(TAGS)){ s2 <- cc |> filter(comp==cn); fl <- quantile(s2$s[s2$nch>=5000 | s2$z>=100], .9); s2 <- s2 |> filter(s > 1.3*fl)
  m <- lm(log10(s) ~ log10(nch)+log10(St)+log10(interval), data=s2); co <- coef(m)
  cat(sprintf("  s[%-6s] floor~%.3f n=%3d  N %+.3f  S %+.3f  D %+.3f  R2 %.3f   -> N * S^%.2f * D^%.2f\n", cn, fl, nrow(s2), co[2],co[3],co[4], summary(m)$r.squared, co[3]/co[2], co[4]/co[2])) }
