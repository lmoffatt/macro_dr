setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
suppressPackageStartupMessages({library(dplyr);library(tidyr)})
d <- read.csv("digest/ir_collapse_cells.csv") |> filter(anchor=="sim")
for(iv in c(0.1,1)){
 for(y in c("dAI","zmax")){
  cat("\n### ",y," at interval",iv," (rows N_ch, cols noise label)\n")
  s <- d |> filter(interval==iv) |> select(nch, z, v=all_of(y)) |>
       mutate(v=round(v,3)) |> pivot_wider(names_from=z, values_from=v) |> arrange(nch)
  print(as.data.frame(s), row.names=FALSE)
 }
}
