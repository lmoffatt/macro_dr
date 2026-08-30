# TEST of the resolvability-window hypothesis (Luciano, 2026-08-29) against the empirical power law.
#  H: distortion is a function of the obstacle to resolving one opening, O = S~/Dhat + c*N*Dhat
#     (per-sample noise power over the unitary step, plus transitions per sample), i.e. of min over the
#     two edges of the window S~ < Dhat < 1/(cN). Predictions: (1) one variable O collapses the distortion
#     at least as well as x = N S~ sqrt(Dhat); (2) the worst interval moves as Dhat* = sqrt(S~/(cN)).
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
suppressPackageStartupMessages({library(dplyr);library(tidyr)})
r <- read.csv("digest/ir_pair_boot.csv") |> mutate(St=0.1*z, D=interval, x=nch*St*sqrt(D))
cs <- read.csv("digest/ir_components_cells.csv") |> mutate(q=6*(ebar^2+s^2)) |> select(nch,z,interval,comp,q) |>
      pivot_wider(names_from=comp, values_from=q, names_prefix="q_") |> rename(D=interval)
r <- inner_join(r, cs, by=c("nch","z","D"))
fl <- max(r$s[r$x>1e4]); ab <- r |> filter(s > 1.5*fl)
R2 <- function(m) summary(m)$r.squared
gain <- function(f1, f2, d) R2(lm(f2, d)) - R2(lm(f1, d))
cat("=== PREDICTION 1: one-variable collapse of the total anisotropy s (cells above 1.5x floor, n=", nrow(ab), ") ===\n")
m_x  <- lm(log10(s) ~ log10(x), ab);  cat(sprintf("  x = N S sqrtD (1 param)      R2 %.3f   residual gain with D as factor %.3f\n", R2(m_x), gain(log10(s)~log10(x), log10(s)~log10(x)+factor(D), ab)))
m_f  <- lm(log10(s) ~ log10(nch)+log10(St)+log10(D), ab); cat(sprintf("  free exponents (3 params)     R2 %.3f\n", R2(m_f)))
best <- NULL
for(cc in 10^seq(-3, 2, 0.1)){ ab$O <- ab$St/ab$D + cc*ab$nch*ab$D; m <- lm(log10(s) ~ log10(O), ab); if(is.null(best) || R2(m) > best$R2) best <- list(c=cc, R2=R2(m), slope=coef(m)[2]) }
ab$O <- ab$St/ab$D + best$c*ab$nch*ab$D
cat(sprintf("  O = S/D + c N D, best c = %.3g (1 fitted exponent + c)  R2 %.3f  slope %+.3f   residual gain with D factor %.3f\n", best$c, best$R2, best$slope, gain(log10(s)~log10(O), log10(s)~log10(O)+factor(D), ab)))
# generalised obstacle with free relative exponent: O2 = (S/D)^p + c (N D)^q  -> equivalent freedom to the 3-param law; report for fairness
best2 <- NULL
for(p in c(0.5,0.75,1)) for(qq in c(0.5,0.75,1)) for(cc in 10^seq(-3,2,0.25)){ O2 <- (ab$St/ab$D)^p + cc*(ab$nch*ab$D)^qq; m <- lm(log10(ab$s) ~ log10(O2)); if(is.null(best2) || R2(m) > best2$R2) best2 <- list(p=p,q=qq,c=cc,R2=R2(m)) }
cat(sprintf("  O2 = (S/D)^p + c (N D)^q, best p=%.2f q=%.2f c=%.3g  R2 %.3f\n", best2$p, best2$q, best2$c, best2$R2))
# does O carry information beyond x?  (partial)
cat(sprintf("  x + log O together: R2 %.3f ; x + free D: %.3f\n", R2(lm(log10(s) ~ log10(x)+log10(O), ab)), R2(lm(log10(s) ~ log10(x)+log10(D), ab))))

cat("\n=== PREDICTION 2: does the worst interval move with sqrt(S/N)?  (sample component, argmax over the 7 intervals per (N,S) cell) ===\n")
flq <- max(r$q_sample[r$x>1e4]); pk <- r |> group_by(nch, St) |> filter(n()==7, max(q_sample) > 3*flq) |>
  summarise(Dstar = D[which.max(q_sample)], ratio = first(St)/first(nch), qmax=max(q_sample), .groups="drop")
cat("  cells with a resolved peak:", nrow(pk), "\n")
m2 <- lm(log10(Dstar) ~ log10(ratio), pk); cat(sprintf("  log10 Dstar ~ log10(S/N): slope %+.3f +- %.3f  (window predicts +0.5; fixed timescale predicts 0)   R2 %.2f\n", coef(m2)[2], summary(m2)$coef[2,2], R2(m2)))
print(pk |> mutate(lr=round(log10(ratio))) |> group_by(lr) |> summarise(n=n(), med_Dstar=median(Dstar), pred_window=signif(sqrt(10^first(lr)/best$c),2)) |> as.data.frame(), row.names=FALSE)
# same test on the TOTAL anisotropy: is there a peak at all, and where?
pk2 <- r |> group_by(nch, St) |> filter(n()==7, max(s) > 3*fl) |> summarise(Dstar=D[which.max(s)], ratio=first(St)/first(nch), .groups="drop")
cat(sprintf("  total anisotropy: %d cells; argmax at the shortest interval in %d of them; slope of log Dstar on log(S/N) %+.2f\n", nrow(pk2), sum(pk2$Dstar==0.01), coef(lm(log10(Dstar)~log10(ratio), pk2))[2]))
cat("\n=== the noise edge: at fixed (N,S), does s DECREASE toward short intervals where Dhat < S~ ?  (cells with S~ >= 0.05, compare Dhat=0.01 vs 0.05) ===\n")
ne <- r |> filter(St >= 0.05, D %in% c(0.01,0.05)) |> select(nch,St,D,s) |> pivot_wider(names_from=D, values_from=s, names_prefix="s_") |> mutate(drop = s_0.01/s_0.05)
print(ne |> group_by(St) |> summarise(n=n(), median_ratio_s001_over_s005=round(median(drop),2), share_below_1=round(mean(drop<1),2)) |> as.data.frame(), row.names=FALSE)
