# Shape witnesses per cell: per-parameter J_pp/G_pp, residual autocorrelation tau_int, residual excess kurtosis.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(purrr);library(readr)})
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/data")
DIRS <- c("1c2ae6f","0ffbda7","87889e6","82b956f")
files <- unlist(lapply(DIRS, function(d) list.files(d, pattern="macro_IR.*battery_sim_G\\.csv$", full.names=TRUE)))
meta <- tibble(f=files) |> mutate(nch=as.numeric(sub(".*nch_([0-9]+)_.*","\\1",f)), z=as.numeric(sub(".*noise_([0-9.]+)_battery.*","\\1",f))) |>
  group_by(nch,z) |> slice(1) |> ungroup()
PN <- c("k_on","k_off","i","S","baseline","N_ch")
one <- function(f,nch,z){
  d <- suppressWarnings(read_csv(f, skip=1, show_col_types=FALSE, progress=FALSE)) |> filter(probit=="mean")
  sc <- d |> filter(component_path=="Probit_statistics_Moment_statistics_Sum_dlogL", statistic=="covariance", value_row==value_col) |>
    transmute(interval=interval_in_tau, p=value_row, J=value)
  sm <- d |> filter(component_path=="Probit_statistics_Moment_statistics_Sum_dlogL", statistic=="mean") |>
    transmute(interval=interval_in_tau, p=value_row, smean=value)
  st <- d |> filter(component_path=="Probit_statistics_Sum_Moment_statistics_dlogL_true", statistic=="covariance", value_row==value_col) |>
    transmute(interval=interval_in_tau, p=value_row, Jstep=value)
  sc <- sc |> inner_join(sm, by=c("interval","p")) |> inner_join(st, by=c("interval","p"))
  G  <- d |> filter(component_path=="Probit_statistics_Moment_statistics_Sum_Gaussian_Fisher_Information", statistic=="mean", value_row==value_col) |>
    transmute(interval=interval_in_tau, p=value_row, G=value)
  pp <- inner_join(sc, G, by=c("interval","p")) |> mutate(param=PN[p+1], nch=nch, z=z)
  ti <- d |> filter(component_path=="Probit_statistics_Report_integral_r_std", statistic=="integral_correlation_lag") |> transmute(interval=interval_in_tau, tau_int=value)
  r2 <- d |> filter(component_path=="Probit_statistics_Moment_statistics_Sum_r2_std", statistic %in% c("mean","variance")) |>
    transmute(interval=interval_in_tau, statistic, value) |> pivot_wider(names_from=statistic, values_from=value) |> rename(r2m=mean, r2v=variance)
  r1 <- d |> filter(component_path=="Probit_statistics_Moment_statistics_Sum_r_std", statistic %in% c("mean","variance")) |>
    transmute(interval=interval_in_tau, statistic, value) |> pivot_wider(names_from=statistic, values_from=value) |> rename(r1m=mean, r1v=variance)
  r1s <- d |> filter(component_path=="Probit_statistics_Sum_Moment_statistics_r_std_false", statistic=="variance") |> transmute(interval=interval_in_tau, r1v_step=value)
  cw <- ti |> inner_join(r2, by="interval") |> inner_join(r1, by="interval") |> inner_join(r1s, by="interval") |> mutate(nch=nch, z=z)
  list(pp=pp, cw=cw)
}
res <- pmap(meta, function(f,nch,z) one(f,nch,z))
pp <- bind_rows(map(res,"pp")); cw <- bind_rows(map(res,"cw")) |> mutate(n=10/interval, kurt=r2v/n-2, acf1=r1v/r1v_step-1, r2bar=r2m/n-1)
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828"); write_csv(pp, "digest/ir_score_JG_perparam.csv"); write_csv(cw, "digest/ir_residual_shape.csv")
d <- read.csv("digest/ir_collapse_cells.csv") |> filter(anchor=="sim") |> select(nch,z,interval,ebar,aniso)
cell <- cw |> inner_join(d, by=c("nch","z","interval")) |> mutate(St=0.1*z, x=nch*St*sqrt(interval),
  Dcl=cut(interval, c(0,0.05,0.2,1.01), labels=c("short","mid","long")))
cat("cells:", nrow(cell), "\n")
cat("\nresidual witnesses in the corner x<1, medians by interval class (null: tau_int 1, kurt 0, acf-term 0, r2bar 0)\n")
print(cell |> filter(x<1) |> group_by(Dcl) |> summarise(n=n(), ebar=round(median(ebar),3), tau_int=round(median(tau_int),3), kurt=round(median(kurt),3),
  acf_term=round(median(acf1),3), r2bar=round(median(r2bar),4), .groups="drop") |> as.data.frame(), row.names=FALSE)
big <- cell |> filter(x<1)
cat("\ncorrelations with ebar in the corner: tau_int", round(cor(big$ebar, big$tau_int),2), " kurt", round(cor(big$ebar, big$kurt),2),
    " acf_term", round(cor(big$ebar, big$acf1),2), " r2bar", round(cor(big$ebar, big$r2bar),2), "\n")
sh <- big |> filter(Dcl=="short")
cat("short-interval corner only (n=", nrow(sh), "): cor(ebar, tau_int)", round(cor(sh$ebar, sh$tau_int),2), " cor(ebar, kurt)", round(cor(sh$ebar, sh$kurt),2), " cor(ebar, acf_term)", round(cor(sh$ebar, sh$acf1),2), "\n")
# per-parameter J/G
pq <- pp |> mutate(lr=log(J/G)) |> inner_join(cell |> select(nch,z,interval,x,Dcl,ebar), by=c("nch","z","interval"))
cat("\nper-parameter log(J_pp/G_pp), corner x<1, median by interval class (0 = calibrated in that coordinate)\n")
print(pq |> filter(x<1) |> group_by(param, Dcl) |> summarise(lr=round(median(lr),3), .groups="drop") |> pivot_wider(names_from=Dcl, values_from=lr) |> as.data.frame(), row.names=FALSE)
cat("\nDECOMPOSITION log(J/G) = log(Jstep/G) [per-step shape] + log(J/Jstep) [temporal correlation], corner x<1, SHORT intervals, median by param\n")
print(pq |> filter(x<1, Dcl=="short") |> mutate(shape=log(Jstep/G), tcorr=log(J/Jstep)) |> group_by(param) |>
  summarise(total=round(median(lr),3), shape=round(median(shape),3), tcorr=round(median(tcorr),3), .groups="drop") |> as.data.frame(), row.names=FALSE)
cat("same for LONG intervals\n")
print(pq |> filter(x<1, Dcl=="long") |> mutate(shape=log(Jstep/G), tcorr=log(J/Jstep)) |> group_by(param) |>
  summarise(total=round(median(lr),3), shape=round(median(shape),3), tcorr=round(median(tcorr),3), .groups="drop") |> as.data.frame(), row.names=FALSE)
cat("\nsame, all cells with x>100 (should be ~0):\n")
print(pq |> filter(x>100) |> group_by(param) |> summarise(lr=round(median(lr),3), q05=round(quantile(lr,.05),3), q95=round(quantile(lr,.95),3)) |> as.data.frame(), row.names=FALSE)
