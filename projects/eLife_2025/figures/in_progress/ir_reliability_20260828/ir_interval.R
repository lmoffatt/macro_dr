setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
suppressPackageStartupMessages({library(dplyr);library(tidyr)})
d <- read.csv("digest/ir_collapse_cells.csv") |> filter(anchor=="sim")
cat("dAI vs interval (rows N_ch, cols interval) at noise label 0.1\n")
print(as.data.frame(d |> filter(z==0.1) |> select(nch,interval,dAI) |> mutate(dAI=round(dAI,3)) |>
  pivot_wider(names_from=interval, values_from=dAI) |> arrange(nch)), row.names=FALSE)
cat("\ndAI vs interval at noise label 1\n")
print(as.data.frame(d |> filter(z==1) |> select(nch,interval,dAI) |> mutate(dAI=round(dAI,3)) |>
  pivot_wider(names_from=interval, values_from=dAI) |> arrange(nch)), row.names=FALSE)
cat("\nzmax vs interval at noise label 0.1\n")
print(as.data.frame(d |> filter(z==0.1) |> select(nch,interval,zmax) |> mutate(zmax=round(zmax,3)) |>
  pivot_wider(names_from=interval, values_from=zmax) |> arrange(nch)), row.names=FALSE)
