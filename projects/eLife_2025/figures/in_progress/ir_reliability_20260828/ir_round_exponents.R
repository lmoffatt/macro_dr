suppressPackageStartupMessages({library(dplyr);library(tidyr)})
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
comp <- read.csv("digest/ir_components_cells.csv")
bias <- read.csv("digest/ir_bias_whitened.csv") |> group_by(nch,z,interval) |> summarise(ebar=mean(v), s=sqrt(mean((v-mean(v))^2)), .groups="drop") |> mutate(comp="bias")
all <- bind_rows(bias, comp |> select(nch,z,interval,comp,ebar,s)) |> mutate(St=0.1*z, D=interval, xref=nch*St*sqrt(D))
fl <- all |> filter(xref > 1e4) |> group_by(comp) |> summarise(floor=max(s), .groups="drop"); all <- all |> inner_join(fl, by="comp")
grid <- expand_grid(eS=c(0.5,0.75,1,1.5,2), eD=c(0,0.5,1))
for(cn in c("bias","total","sample","corr")){
  d <- all |> filter(comp==cn, s > 1.5*floor)
  res <- grid |> rowwise() |> mutate(fit=list({ x <- d$nch*d$St^eS*d$D^eD; m <- lm(log10(d$s) ~ log10(x)); c(slope=coef(m)[[2]], R2=summary(m)$r.squared, spread=median(10^abs(resid(m)))) })) |>
    ungroup() |> mutate(slope=sapply(fit,`[[`,"slope"), R2=sapply(fit,`[[`,"R2"), spread=sapply(fit,`[[`,"spread")) |> select(-fit) |> arrange(-R2)
  cat(sprintf("\n=== %s (n=%d) top 5 of %d ===\n", cn, nrow(d), nrow(grid))); print(as.data.frame(head(res,5) |> mutate(across(c(slope,R2,spread), ~round(.,3)))), row.names=FALSE)
}
