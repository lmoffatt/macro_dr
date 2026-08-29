setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
suppressPackageStartupMessages({library(dplyr);library(tidyr)})
PN <- c("k_on","k_off","i","S","baseline","N_ch")
ci <- read.csv("digest/ir_bias_cis.csv") |> mutate(S=0.1*z, D=interval_in_tau, x=nch*S*sqrt(D),
        sig=(lo>0|hi<0), zp=abs(m)/se, pname=PN[param_index+1])
for(p in PN){ s <- ci |> filter(pname==p)
  cat(sprintf("%-9s |bias| q50 %.4f q90 %.4f max %.3f   |bias|/SE q50 %.3f max %.3f   SE q50 %.3f max %.2f   sig %d/%d\n",
   p, median(abs(s$m)), quantile(abs(s$m),.9), max(abs(s$m)), median(s$zp), max(s$zp), median(s$se), max(s$se), sum(s$sig), nrow(s)))}
phys <- ci |> filter(pname %in% c("k_on","k_off","i","N_ch")) |> group_by(nch,z,D,x,S) |>
  summarise(babs=max(abs(m)), p=pname[which.max(abs(m))], sig=sig[which.max(abs(m))], .groups="drop") |>
  mutate(relerr = 10^babs - 1)
cat("\nPHYSICAL params only: max|bias| quantiles (log10):\n"); print(round(quantile(phys$babs, c(.5,.9,.99,1)),4))
cat("relative error (%) quantiles:\n"); print(round(100*quantile(phys$relerr, c(.5,.9,.99,1)),2))
cat("carrier:\n"); print(table(phys$p))
s <- phys |> filter(sig); cat("significant:", nrow(s), "\n")
fit <- function(y, xs, lab){ m <- lm(log10(y) ~ log10(xs)); r <- 10^abs(resid(m))
  cat(sprintf("  %-16s slope %+.3f  R2 %.3f  spread x%.2f\n", lab, coef(m)[2], summary(m)$r.squared, median(r))) }
fit(s$babs, s$x, "x=N S sqrtD"); fit(s$babs, s$nch, "N"); fit(s$babs, s$nch*s$S, "N S")
m3 <- lm(log10(babs) ~ log10(nch)+log10(S)+log10(D), data=s)
cat(sprintf("  free: N %+.3f S %+.3f D %+.3f R2 %.3f\n", coef(m3)[2],coef(m3)[3],coef(m3)[4],summary(m3)$r.squared))
cat("\nworst physical cells:\n"); print(phys |> arrange(-babs) |> head(8) |> mutate(across(c(babs,relerr), ~round(.,3))) |> as.data.frame(), row.names=FALSE)
