suppressPackageStartupMessages({library(dplyr)})
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
d <- read.csv("digest/ir_components_cells.csv") |> filter(comp=="sample") |> mutate(q=6*(ebar^2+s^2), St=0.1*z, D=interval, xref=nch*St*sqrt(D))
cf <- 0.00139
# additive-floor model fitted on log scale, all 560 cells:  q = A * N^a * S^b * (D0/D + D/D0)^(-k) + cf*D
f <- function(p, dd) with(dd, exp(p[1]) * nch^p[2] * St^p[3] * (exp(p[4])/D + D/exp(p[4]))^(-p[5]) + cf*D)
obj <- function(p, dd) sum((log(dd$q) - log(f(p, dd)))^2)
p0 <- c(log(0.05), -0.27, -0.7, log(0.1), 2.5)
o <- optim(p0, obj, dd=d, method="BFGS", control=list(maxit=2000))
p <- o$par; R2 <- 1 - o$value/sum((log(d$q)-mean(log(d$q)))^2)
cat(sprintf("ALL 560 cells, model q = A N^a S^b (D0/D+D/D0)^-k + floor(D):\n  A=%.3g  a=%+.3f  b=%+.3f  D0=%.3f  k=%.2f   R2(log) %.3f  spread x%.2f\n",
  exp(p[1]), p[2], p[3], exp(p[4]), p[5], R2, exp(median(abs(log(d$q)-log(f(p,d)))))))
cat(sprintf("  normalised: x = N * S^%.2f * (D0/D + D/D0)^%.2f , slope %+.3f\n", p[3]/p[2], p[5]/p[2]*(-1)*(-1), p[2]))
# same model with the floor FIXED to zero (to show what the floor buys)
f0 <- function(p, dd) with(dd, exp(p[1]) * nch^p[2] * St^p[3] * (exp(p[4])/D + D/exp(p[4]))^(-p[5]))
o0 <- optim(p0, function(p,dd) sum((log(dd$q)-log(f0(p,dd)))^2), dd=d, method="BFGS", control=list(maxit=2000))
cat(sprintf("  without floor term: R2(log) %.3f\n", 1 - o0$value/sum((log(d$q)-mean(log(d$q)))^2)))
# signal-dominated cells: model signal > 3*floor -> collapse quality of (q - floor) vs x there
d$sig <- f(p,d) - cf*d$D; s <- d |> filter(sig > 3*cf*D) |> mutate(qc = q - cf*D, x = nch * St^(p[3]/p[2]) * (exp(p[4])/D + D/exp(p[4]))^(-p[5]/p[2]))
s <- s |> filter(qc > 0); m1 <- lm(log10(qc) ~ log10(x), data=s); m2 <- lm(log10(qc) ~ log10(x) + factor(D), data=s)
cat(sprintf("signal-dominated cells n=%d: 1-D slope %+.3f R2 %.3f gainD %.3f spread x%.2f\n", nrow(s), coef(m1)[2], summary(m1)$r.squared, summary(m2)$r.squared-summary(m1)$r.squared, median(10^abs(resid(m1)))))
cat("residual by N:\n"); s$res <- resid(m1); print(s |> group_by(nch) |> summarise(med=round(median(res),2), n=n()) |> as.data.frame(), row.names=FALSE)
cat("residual by S:\n"); print(s |> group_by(St) |> summarise(med=round(median(res),2), n=n()) |> as.data.frame(), row.names=FALSE)
saveRDS(list(p=p, cf=cf), "digest/ir_sample_nls.rds")
