setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
suppressPackageStartupMessages({library(dplyr)})
d <- read.csv("digest/ir_collapse_cells.csv") |> filter(anchor=="sim") |>
  mutate(S=0.1*z, x2=S*nch*interval, x1=S*nch)
# above-floor cells for tightness comparison
s <- d |> filter(aniso > 0.055)
for(v in c("x1","x2")){
  m <- lm(log10(aniso) ~ log10(s[[v]]), data=s)
  r <- 10^abs(resid(m))
  cat(sprintf("%s: slope %+.3f  R2 %.3f  spread median x%.2f  q90 x%.2f\n",
      v, coef(m)[2], summary(m)$r.squared, median(r), quantile(r,.9)))
}
# envelope for x2: smallest a with 0 violations at exponents -0.2 and -0.25, floor 0.045
for(b in c(0.2, 0.25)){
  a <- max((d$aniso/pmax(d$x2^-b, 0.045/1))[d$aniso > 0.045] * 1)  # crude; refine below
  amin <- max(ifelse(d$aniso > 0.045, d$aniso/d$x2^-b, 0))
  # proper: a must satisfy aniso <= max(a*x2^-b, 0.045) -> only cells above floor constrain
  cells <- d |> filter(aniso > 0.045)
  amin <- max(cells$aniso * cells$x2^b)
  cat(sprintf("exponent -%.2f: min envelope constant a = %.3f\n", b, amin))
}
# same for ebar envelope vs x2
cells <- d |> filter(abs(ebar) > 0.02)
cat(sprintf("ebar: exponent -0.2 min a = %.3f\n", max(abs(cells$ebar)*cells$x2^0.2)))
