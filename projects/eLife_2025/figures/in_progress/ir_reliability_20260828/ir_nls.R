setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
suppressPackageStartupMessages({library(dplyr)})
d <- read.csv("digest/ir_collapse_cells.csv") |> filter(anchor=="sim") |> mutate(S=0.1*z, D=interval)
fitnls <- function(y, lab, start){
  dd <- d |> filter(is.finite(.data[[y]])); dd$Y <- dd[[y]]^2
  m <- nls(Y ~ f2 + exp(2*lA)*nch^(2*a)*S^(2*b)*D^(2*c), data=dd,
           start=c(f2=start$f2, lA=start$lA, a=-0.4, b=-0.2, c=-0.2),
           control=nls.control(maxiter=500, warnOnly=TRUE))
  co <- coef(m); f <- sqrt(max(co["f2"],0))
  cat(sprintf("\n%s   n=%d\n", lab, nrow(dd)))
  cat(sprintf("  floor  = %.4f   (fitted, free)\n", f))
  cat(sprintf("  signal = %.3f * N^%.3f * S^%.3f * Delta^%.3f\n", exp(co["lA"]), co["a"], co["b"], co["c"]))
  cat(sprintf("  normalised: x = N * S^%.2f * Delta^%.2f ,  signal = %.3f * x^%.3f\n",
              co["b"]/co["a"], co["c"]/co["a"], exp(co["lA"]), co["a"]))
  pr <- sqrt(pmax(exp(2*co["lA"])*dd$nch^(2*co["a"])*dd$S^(2*co["b"])*dd$D^(2*co["c"]),0))
  ob <- sqrt(pmax(dd$Y - co["f2"],0))
  ok <- ob > f
  cat(sprintf("  R2 on sqrt scale (cells above floor, n=%d): %.3f ; median |ratio| = %.2f\n",
      sum(ok), 1-var(log10(ob[ok]/pr[ok]))/var(log10(ob[ok])), 10^median(abs(log10(ob[ok]/pr[ok])))))
  invisible(co)
}
cD <- fitnls("dAI","DISTORTION (affine of GIDM)", list(f2=0.005,lA=0))
cZ <- fitnls("zmax","BIAS max|b|/SE", list(f2=2e-4,lA=-1))
# ---- one-variable versions, floor free, exponents tied to a single composite
one <- function(y, xf, lab, start=list(f2=0.005,lA=0)){
  dd <- d |> filter(is.finite(.data[[y]])); dd$Y <- dd[[y]]^2; dd$X <- xf(dd)
  m <- nls(Y ~ f2 + exp(2*lA)*X^(2*a), data=dd, start=c(f2=start$f2,lA=start$lA,a=-0.2),
           control=nls.control(maxiter=500, warnOnly=TRUE)); co <- coef(m)
  f <- sqrt(max(co["f2"],0)); pr <- exp(co["lA"])*dd$X^co["a"]; ob <- sqrt(pmax(dd$Y-co["f2"],0)); ok <- ob>f
  cat(sprintf("  %-22s floor %.3f  slope %+.3f  R2(above floor) %.3f  spread x%.2f\n", lab, f, co["a"],
      1-var(log10(ob[ok]/pr[ok]))/var(log10(ob[ok])), 10^median(abs(log10(ob[ok]/pr[ok])))))
}
cat("\n--- single-variable collapses, DISTORTION ---\n")
one("dAI", function(z) z$nch*z$S, "x = N*S")
one("dAI", function(z) z$nch*z$S*z$D, "x = N*S*Delta")
one("dAI", function(z) z$S/z$nch, "x = S/N (ratio)")
one("dAI", function(z) z$nch, "x = N")
one("dAI", function(z) z$nch*z$S^0.6*z$D^0.73, "x = N*S^.6*D^.73")
cat("\n--- single-variable collapses, BIAS ---\n")
one("zmax", function(z) z$nch, "x = N", list(f2=2e-4,lA=-1))
one("zmax", function(z) z$nch*z$S, "x = N*S", list(f2=2e-4,lA=-1))
one("zmax", function(z) z$S/z$nch, "x = S/N (ratio)", list(f2=2e-4,lA=-1))
one("zmax", function(z) z$nch*z$S^0.2, "x = N*S^0.2", list(f2=2e-4,lA=-1))
