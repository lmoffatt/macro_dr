setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
suppressPackageStartupMessages({library(dplyr);library(tidyr)})
d <- read.csv("digest/ir_collapse_cells.csv") |> filter(anchor=="sim") |> mutate(S=0.1*z)
cat("=== signed amplitude ebar = mean(log lambda) : sign and size ===\n")
print(as.data.frame(d |> filter(interval==0.05) |> select(nch,z,ebar) |> mutate(ebar=round(ebar,3)) |>
  pivot_wider(names_from=z, values_from=ebar) |> arrange(nch)), row.names=FALSE)
cat("\n=== anisotropy s = sd(log lambda) ===\n")
print(as.data.frame(d |> filter(interval==0.05) |> select(nch,z,aniso) |> mutate(aniso=round(aniso,3)) |>
  pivot_wider(names_from=z, values_from=aniso) |> arrange(nch)), row.names=FALSE)
cat("\n=== z*sqrt(N_ch)  (bias/SE times sqrt N) at low noise (label<=1) ===\n")
print(as.data.frame(d |> filter(z<=1) |> mutate(k=round(zmax*sqrt(nch),2)) |> select(nch,interval,k,z) |>
  group_by(nch) |> summarise(median_k=median(k), q25=quantile(k,.25), q75=quantile(k,.75), n=n())), row.names=FALSE)
cat("\n=== dAI*sqrt(N_ch) ===\n")
print(as.data.frame(d |> filter(z<=1) |> group_by(nch) |> summarise(median=median(round(dAI*sqrt(nch),2)))), row.names=FALSE)
