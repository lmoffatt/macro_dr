suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(patchwork)})
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
d <- read.csv("digest/ir_components_cells.csv") |> filter(comp=="sample") |> mutate(q=6*(ebar^2+s^2), St=0.1*z, D=interval, xref=nch*St*sqrt(D))
floor <- max(d$q[d$xref>1e4]); cat("sample floor (max q at xref>1e4):", signif(floor,3), "\n")
# 1. where is the peak? per (N,S) cell with all 7 intervals, argmax over D of q; and the profile normalised to the cell max
prof <- d |> group_by(nch,z) |> filter(n()==7, max(q) > 3*floor) |> mutate(qn=q/max(q), Dmax=D[which.max(q)]) |> ungroup()
cat("cells with a resolved peak:", n_distinct(prof$nch, prof$z), "\n"); cat("argmax D distribution:\n"); print(table(prof |> distinct(nch,z,Dmax) |> pull(Dmax)))
cat("median normalised profile q/qmax by D:\n"); print(prof |> group_by(D) |> summarise(med=round(median(qn),2), q25=round(quantile(qn,.25),2), q75=round(quantile(qn,.75),2)) |> as.data.frame(), row.names=FALSE)
# 2. fits: log10 q ~ log10 N + log10 S + g(D), for two g families and a grid of D0
s <- d |> filter(q > 1.5*floor); cat("fit cells:", nrow(s), "\n")
fitg <- function(g, lab){ m <- lm(log10(q) ~ log10(nch) + log10(St) + g, data=s); co <- coef(m)
  m2 <- lm(log10(q) ~ log10(nch) + log10(St) + g + factor(D), data=s)
  c(R2=summary(m)$r.squared, gainD=summary(m2)$r.squared-summary(m)$r.squared, bN=co[2], bS=co[3], bG=co[4]) }
cat("\nreference (power law in D):\n"); r0 <- fitg(log10(s$D), "pow"); cat(sprintf("  R2 %.3f gainD %.3f  N %+.3f S %+.3f D %+.3f\n", r0[1], r0[2], r0[3], r0[4], r0[5]))
cat("\nGAUSSIAN in log D:  g = (log10(D/D0))^2\n")
for(D0 in c(0.03,0.05,0.07,0.1,0.14,0.2)){ r <- fitg((log10(s$D/D0))^2); cat(sprintf("  D0=%.2f  R2 %.3f gainD %.3f  N %+.3f S %+.3f  coef(g) %+.3f  -> x = N * S^%.2f * 10^(%.2f*(log10(D/D0))^2)\n", D0, r[1], r[2], r[3], r[4], r[5], r[4]/r[3], r[5]/r[3])) }
cat("\nSUM OF POWERS:  g = log10(D0/D + D/D0)\n")
for(D0 in c(0.03,0.05,0.07,0.1,0.14,0.2)){ r <- fitg(log10(D0/s$D + s$D/D0)); cat(sprintf("  D0=%.2f  R2 %.3f gainD %.3f  N %+.3f S %+.3f  coef(g) %+.3f  -> x = N * S^%.2f * (D0/D + D/D0)^%.2f\n", D0, r[1], r[2], r[3], r[4], r[5], r[4]/r[3], r[5]/r[3])) }
