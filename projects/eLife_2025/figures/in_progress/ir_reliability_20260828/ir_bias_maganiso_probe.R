setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
suppressPackageStartupMessages({library(dplyr);library(tidyr)})
d <- read.csv("digest/ir_bias_direction.csv") |> mutate(S=0.1*z, x2=S*nch*interval, x1=S*nch)
cat("floor (Mahalanobis noise, sqrt(tr GIDM/nsim)) quantiles:\n"); print(round(quantile(d$floor, c(0,.5,.9,1)),4))
cat("\nmag quantiles by nch (null ~ floor ~0.025):\n")
print(d |> group_by(nch) |> summarise(med=round(median(mag),3), q90=round(quantile(mag,.9),3),
      med_over_floor=round(median(mag/floor),1), n=n()), n=12)
# scaling of the excess magnitude
s <- d |> filter(mag > 2*floor)
cat("\ncells with mag > 2*floor:", nrow(s), "\n")
m <- lm(log10(mag) ~ log10(nch)+log10(S)+log10(interval), data=s)
cat(sprintf("exponents: N %+.3f  S %+.3f  Delta %+.3f  R2=%.3f\n",
    coef(m)[2],coef(m)[3],coef(m)[4], summary(m)$r.squared))
for(v in c("x1","x2","nch")){
  mm <- lm(log10(s$mag) ~ log10(s[[v]]))
  cat(sprintf("  x=%-4s slope %+.3f R2=%.3f\n", v, coef(mm)[2], summary(mm)$r.squared))}
# envelope vs nch
cells <- s
cat(sprintf("\nenvelope mag <= a/sqrt(N): min a = %.2f ; vs a*x2^-b fit b=%.2f\n",
    max(cells$mag*sqrt(cells$nch)), -coef(lm(log10(mag)~log10(x2), data=s))[2]))
# anisotropy: share in softest direction (isotropy = 1/6 = 0.167)
cat("\nshare of |b|^2 in SOFTEST covariance direction, by nch (real-bias cells mag>2*floor):\n")
print(s |> group_by(nch) |> summarise(med_soft=round(median(share_soft),2),
      med_top2=round(median(share_top2),2), med_stiff=round(median(share_stiff),3), n=n()), n=12)
cat("\nnull-region cells (mag<2*floor) for comparison:\n")
print(d |> filter(mag<=2*floor) |> summarise(med_soft=round(median(share_soft),2),
      med_top2=round(median(share_top2),2), med_stiff=round(median(share_stiff),3), n=n()))
cat("\nsign consistency of the softest-direction projection c_soft at low N (real cells):\n")
print(s |> filter(nch<=100) |> group_by(nch) |> summarise(pos=sum(c_soft>0), n=n()), n=8)
