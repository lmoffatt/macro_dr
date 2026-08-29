setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
suppressPackageStartupMessages({library(dplyr);library(tidyr)})
b <- read.csv("digest/ir_bias_cis.csv") |>
  mutate(se_boot = (hi-lo)/(2*1.96),          # bootstrap SE of the bias estimate itself
         zb = m/se_boot,                       # how resolved the bias is (measurement z)
         zr = m/se,                            # bias in units of the REPORTED parameter SE
         S = 0.1*z)
PN <- c("on","off","unitary_current","Current_Noise","Current_Baseline","Num_ch_mean")
b$param <- PN[b$param_index+1]

# 1) Bonferroni joint test per cell: any |zb| > 2.64 (alpha .05 over 6 params)
cell <- b |> group_by(nch, z, interval_in_tau) |>
  summarise(any_sig = any(abs(zb) > 2.64), zmax = max(abs(zr)), .groups="drop")
cat("fraction of cells with Bonferroni-significant bias, by N_ch (null expectation ~0.05):\n")
print(cell |> group_by(nch) |> summarise(frac = round(mean(any_sig),2), n=n()), n=12)

# 2) sign consistency across independent noise files at high N: per param, mean zr across cells
cat("\nhigh-N (>=2000) per-param bias in reported-SE units: mean zr, t-stat across cells:\n")
hi <- b |> filter(nch >= 2000)
print(hi |> group_by(param) |> summarise(mean_zr = round(mean(zr),4), sd = round(sd(zr),4),
   t = round(mean(zr)/ (sd(zr)/sqrt(n())),1), n = n()) |> arrange(-abs(t)), n=8)
cat("\nsame at N_ch = 10000 only:\n")
print(b |> filter(nch==10000) |> group_by(param) |> summarise(mean_zr=round(mean(zr),4),
   t=round(mean(zr)/(sd(zr)/sqrt(n())),1), n=n()) |> arrange(-abs(t)), n=8)
cat("\nand mid-N (100-1000):\n")
print(b |> filter(nch>=100, nch<=1000) |> group_by(param) |> summarise(mean_zr=round(mean(zr),4),
   t=round(mean(zr)/(sd(zr)/sqrt(n())),1), n=n()) |> arrange(-abs(t)), n=8)

# 3) measurement floor: with nsim=10000, zr noise per param ~ 1/sqrt(10000)=0.01
cat("\nmedian |zr| at high N per param (measurement noise alone would give ~0.008):\n")
print(hi |> group_by(param) |> summarise(med = round(median(abs(zr)),4)), n=8)
