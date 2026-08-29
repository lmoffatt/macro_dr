setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
suppressPackageStartupMessages({library(dplyr);library(tidyr)})
d <- read.csv("digest/ir_collapse_cells.csv") |> filter(anchor=="sim") |>
  mutate(S = 0.1*z,                 # dimensionless noise  S~ = 0.1 * label
         vin = S/interval,          # instrumental variance per sample, units i^2
         vgate = nch/4,             # gating variance per sample, units i^2
         V = vgate + vin)
FLOOR_D <- 0.07; FLOOR_Z <- 0.015
d <- d |> mutate(exD = sqrt(pmax(0, dAI^2 - FLOOR_D^2)),
                 exZ = sqrt(pmax(0, zmax^2 - FLOOR_Z^2)))
cat("floor check: quantiles of dAI in the flat corner (nch>=5000 | z>=100):\n")
print(round(quantile(d |> filter(nch>=5000 | z>=100) |> pull(dAI), c(0,.25,.5,.75,1)),3))
cat("quantiles of zmax there:\n")
print(round(quantile(d |> filter(nch>=5000 | z>=100) |> pull(zmax), c(0,.25,.5,.75,1)),4))

s <- d |> filter(exD > 0.05)          # cells clearly above the floor
cat("\ncells above floor:", nrow(s), "of", nrow(d), "\n")
y <- log10(s$exD)
mods <- list(
 "A  logN + logS + logDelta"      = lm(y ~ log10(s$nch) + log10(s$S) + log10(s$interval)),
 "B  logN + log(S/Delta)=log vin" = lm(y ~ log10(s$nch) + log10(s$vin)),
 "C  log(V)= log(N/4 + S/Delta)"  = lm(y ~ log10(s$V)),
 "C' log(V) + logN"               = lm(y ~ log10(s$V) + log10(s$nch)),
 "D  log(S*N)"                    = lm(y ~ log10(s$S*s$nch)),
 "D' log(vin*N)"                  = lm(y ~ log10(s$vin*s$nch)),
 "E  log(S/N)  [ratio]"           = lm(y ~ log10(s$S/s$nch)),
 "E' log(vin/N) [ratio, per-sample]" = lm(y ~ log10(s$vin/s$nch)),
 "F  logN only"                   = lm(y ~ log10(s$nch)),
 "G  log(N/V^1.5) [skewness]"     = lm(y ~ log10(s$nch/s$V^1.5)),
 "H  log(N/V^2)   [kurtosis]"     = lm(y ~ log10(s$nch/s$V^2))
)
for(nm in names(mods)){ m <- mods[[nm]]
  cat(sprintf("  %-34s R2=%.3f  rmse(dex in log10)=%.3f  coef= %s\n", nm, summary(m)$r.squared,
      sqrt(mean(resid(m)^2)), paste(sprintf("%+.3f", coef(m)[-1]), collapse=" ")))}
