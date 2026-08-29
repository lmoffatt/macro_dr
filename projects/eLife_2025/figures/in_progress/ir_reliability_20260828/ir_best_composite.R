setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
suppressPackageStartupMessages({library(dplyr);library(tidyr)})
d <- read.csv("digest/ir_collapse_cells.csv") |> filter(anchor=="sim") |>
  mutate(S=0.1*z, D=interval)
fit_one <- function(y, floor, lab){
  s <- d |> mutate(ex = sqrt(pmax(0, .data[[y]]^2 - floor^2))) |> filter(ex > 0.6*floor)
  m <- lm(log10(ex) ~ log10(nch)+log10(S)+log10(D), data=s); co <- coef(m)
  a <- co[2]; b <- co[3]; g <- co[4]
  cat(sprintf("\n%s  (floor %.3f, n=%d above)\n", lab, floor, nrow(s)))
  cat(sprintf("  exponents: N %+.3f   S %+.3f   Delta %+.3f    R2=%.3f  rmse=%.3f dex\n",
              a,b,g,summary(m)$r.squared, sqrt(mean(resid(m)^2))))
  cat(sprintf("  normalised (N=1): S^%.2f  Delta^%.2f   -> x = N * S^%.2f * D^%.2f\n", b/a, g/a, b/a, g/a))
  cands <- list("N*S (product)"=log10(s$nch*s$S), "S/N (ratio)"=log10(s$S/s$nch),
                "N alone"=log10(s$nch), "best composite"=log10(s$nch*s$S^(b/a)*s$D^(g/a)),
                "N*S*D"=log10(s$nch*s$S*s$D), "N*sqrt(S)"=log10(s$nch*sqrt(s$S)))
  for(nm in names(cands)){
    mm <- lm(log10(s$ex) ~ cands[[nm]])
    cat(sprintf("     x=%-16s R2=%.3f  rmse=%.3f dex (factor %.2f)\n", nm, summary(mm)$r.squared,
        sqrt(mean(resid(mm)^2)), 10^sqrt(mean(resid(mm)^2))))
  }
  invisible(list(a=a,b=b,g=g,s=s))
}
D1 <- fit_one("dAI", 0.070, "DISTORTION  (affine distance of GIDM to I, 6 params)")
D2 <- fit_one("dAI", 0.050, "DISTORTION  (lower floor, sensitivity check)")
B1 <- fit_one("zmax", 0.015, "BIAS  (max |bias| / SE, per recording)")
